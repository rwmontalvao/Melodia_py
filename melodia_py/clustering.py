# Copyright 2021-2026 KU Leuven.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Rinaldo Wander Montalvão, PhD
#
import math
import warnings
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial.distance import pdist

from sty import fg

from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering

from melodia_py.geometryparser import _UNCLUSTERED, _defined

warnings.filterwarnings('ignore', category=PDBConstructionWarning)

# ---------------------------------------------------------------------------
# Named constants
# ---------------------------------------------------------------------------

# FIX [style]: Named constant for the pairwise Cα outlier-removal threshold (Å).
# Any anchor position where the maximum inter-structure Cα distance exceeds
# this value is excluded from the conserved-region annotation.
_OUTLIER_DIST_ANGSTROM: float = 2.0

# FIX [style]: Maximum number of top anchor regions used for SA optimisation.
_MAX_ANCHORS: int = 5

# Minimum run length (residues) for a region to qualify as an anchor.
_MIN_ANCHOR_LENGTH: int = 5

# Minimum run length (residues) kept after short-region pruning.
_MIN_REGION_LENGTH: int = 3


# ---------------------------------------------------------------------------
# Low-level geometry helpers
# ---------------------------------------------------------------------------

def rmsd(x: np.ndarray, y: np.ndarray) -> float:
    """
    Compute the RMSD between two Nx3 coordinate arrays.

    :param x: Reference coordinates, shape (N, 3)
    :param y: Mobile coordinates, shape (N, 3)
    :return: Root-mean-square deviation in the same units as the input
    """
    # FIX [correctness]: Original formula multiplied squared differences by 3
    # instead of summing across the 3 spatial dimensions, giving a result
    # sqrt(3)× larger than the standard RMSD.
    # Correct formula: sqrt( mean_over_atoms( sum_over_xyz( (x-y)^2 ) ) )
    return float(np.sqrt(((x - y) ** 2).sum(axis=1).mean()))


def select(data: np.ndarray, seg: List[int], msk: List[bool]) -> np.ndarray:
    """
    Extract the masked subset of coordinates at the given segment indices.

    :param data: Full coordinate array, shape (L, 3)
    :param seg: Flat list of alignment-column indices forming the segment
    :param msk: Boolean mask of the same length as seg; True = keep
    :return: Selected coordinate rows, shape (M, 3)
    """
    return np.array([data[j] for i, j in enumerate(seg) if msk[i]])


def superposition(
    xo: np.ndarray,
    yo: np.ndarray,
    seg: List[int],
    msk: List[bool],
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the SVD-optimal rotation and translation that superimposes the
    masked segment of *yo* onto the corresponding segment of *xo*.

    :param xo: Reference coordinate array, shape (L, 3)
    :param yo: Mobile coordinate array, shape (L, 3)
    :param seg: Segment column indices
    :param msk: Boolean mask selecting active residues within the segment
    :return: (rotation matrix 3×3, translation vector 1×3)
    """
    x = select(xo, seg, msk)
    y = select(yo, seg, msk)
    sup = SVDSuperimposer()
    sup.set(x, y)
    sup.run()
    return sup.get_rotran()


def energy(
    xo: np.ndarray,
    yo: np.ndarray,
    seg: List[int],
    msk: List[bool],
) -> float:
    """
    RMSD between *xo* and the superimposed *yo* over the active segment.

    :param xo: Reference coordinate array, shape (L, 3)
    :param yo: Mobile coordinate array, shape (L, 3)
    :param seg: Segment column indices
    :param msk: Boolean mask selecting active residues
    :return: RMSD energy value (Å)
    """
    rot, tran = superposition(xo, yo, seg, msk)
    yt = np.dot(yo, rot) + tran
    return rmsd(xo, yt)


def segments(anchors: List[Tuple[int, int]], members: List[int]) -> List[int]:
    """
    Flatten the selected anchor regions into a single list of column indices.

    :param anchors: List of (ini, end) anchor region pairs
    :param members: Indices into anchors selecting which regions to include
    :return: Flat list of alignment column indices
    """
    seg: List[int] = []
    for member in members:
        ini, end = anchors[member]
        seg.extend(range(ini, end))
    return seg


# ---------------------------------------------------------------------------
# Run-finding helper (replaces three copies of the same while-loop pattern)
# ---------------------------------------------------------------------------

def _find_runs(area: List[int], value: int = 1) -> List[Tuple[int, int]]:
    """
    Find all contiguous runs of *value* in *area*.

    :param area: Per-position binary annotation list
    :param value: The value to search for (default 1)
    :return: List of (start, end) pairs where area[start:end] == value
    """
    runs: List[Tuple[int, int]] = []
    n = len(area)
    i = 0
    while i < n:
        while i < n and area[i] != value:
            i += 1
        start = i
        while i < n and area[i] == value:
            i += 1
        if i > start:
            runs.append((start, i))
    return runs


# ---------------------------------------------------------------------------
# Simulated annealing
# ---------------------------------------------------------------------------

def simulated_annealing(
    xo: np.ndarray,
    yo: np.ndarray,
    anchors: List[Tuple[int, int]],
    members: List[int],
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Use simulated annealing to find the residue mask within the given anchor
    segments that minimises the RMSD between *xo* and the superimposed *yo*.

    :param xo: Reference coordinate array, shape (L, 3)
    :param yo: Mobile coordinate array, shape (L, 3)
    :param anchors: List of (ini, end) anchor region pairs
    :param members: Indices into anchors selecting which regions to use
    :param rng: Optional numpy random Generator for reproducibility.
        If None, uses the global numpy random state.
    :return: (rotation matrix, translation vector, final RMSD energy)
    """
    # FIX [style]: Accept an rng parameter so callers can seed for
    # reproducibility in batch pipelines.
    if rng is None:
        rng = np.random.default_rng()

    seg = segments(anchors, members)
    n_seg = len(seg)

    # FIX [perf]: Use a numpy bool array so toggling and copying are O(1)
    # fixed-size operations rather than Python list allocations.
    msk0 = np.ones(n_seg, dtype=bool)
    energy0 = energy(xo, yo, seg, msk0.tolist())

    temperature = 300.0
    # TODO: test other annealing schedules
    while temperature > 1e-5:
        for _ in range(10_000):
            msk1 = msk0.copy()
            j = rng.integers(n_seg)
            msk1[j] = not msk1[j]

            energy1 = energy(xo, yo, seg, msk1.tolist())
            delta = energy1 - energy0

            if delta < 0.0 or rng.random() < math.exp(-delta / temperature):
                msk0 = msk1
                energy0 = energy1

        temperature *= 0.1

    rot, tran = superposition(xo, yo, seg, msk0.tolist())
    return rot, tran, energy0


# ---------------------------------------------------------------------------
# Main superimposition pipeline
# ---------------------------------------------------------------------------

def superimposer(
    align,
    threshold: float = 0.8,
    csv: bool = False,
    pdb: bool = True,
    max_anchors: int = _MAX_ANCHORS,
    outlier_dist: float = _OUTLIER_DIST_ANGSTROM,
    structure_dir: str = '.',
    rng: Optional[np.random.Generator] = None,
) -> None:
    """
    Superimpose PDB structures referenced in a geometry-annotated alignment
    using differential geometry as the dissimilarity measure.

    Conserved anchor regions are identified by agglomerative clustering on
    (curvature, torsion) pairs, then refined by simulated annealing to
    minimise the inter-structure RMSD. Optionally writes CSV coordinate files
    and superimposed PDB files. A column with a residue whose curvature or
    torsion is NaN (next to a chain end or break) is never conserved.

    :param align: Geometry-annotated MultipleSeqAlignment
    :param threshold: Agglomerative clustering distance threshold
    :param csv: Write per-structure CSV files with Cα coordinates and group labels
    :param pdb: Write superimposed PDB files (<id>_sup.pdb)
    :param max_anchors: Maximum number of top anchor regions used for SA
    :param outlier_dist: Maximum allowed pairwise Cα distance (Å) within an
        anchor before the position is removed as an outlier
    :param structure_dir: Directory containing the PDB files
    :param rng: Optional random Generator for reproducible annealing
    """
    # ------------------------------------------------------------------
    # Collect geometry data and initialise cluster annotations
    # ------------------------------------------------------------------
    data: List[List[float]] = []
    id2pos: Dict[str, int] = {}

    for position, record in enumerate(align):
        if 'structure' not in record.description:
            continue
        id2pos[record.id] = position
        # FIX [correctness]: Was data = [...] (overwrite); must accumulate
        # across all structures so the scaler fits the full dataset.
        data += [
            [curv, tors]
            for curv, tors in zip(
                record.letter_annotations['curvature'],
                record.letter_annotations['torsion'],
            )
        ]
        record.letter_annotations['cluster'] = [0] * len(record.seq)

    scaler = StandardScaler()
    scaler.fit(data)

    # ------------------------------------------------------------------
    # Cluster each alignment column
    # ------------------------------------------------------------------
    clustering = AgglomerativeClustering(distance_threshold=threshold, n_clusters=None)

    for i in range(align.get_alignment_length()):
        xy: List[List[float]] = []
        tags: List[str] = []
        for rec_id, position in id2pos.items():
            record = align[position]
            curv = record.letter_annotations['curvature'][i]
            tors = record.letter_annotations['torsion'][i]
            # Residues too close to a chain end or break for curvature/torsion
            # (NaN) are left out of the clustering and marked unclustered.
            if record.seq[i] == '-':
                continue
            if _defined(curv, tors):
                xy.append([curv, tors])
                tags.append(rec_id)
            else:
                record.letter_annotations['cluster'][i] = _UNCLUSTERED

        if len(xy) > 1:
            clusters = clustering.fit_predict(scaler.transform(xy))
            cluster_map = {tag: int(c) + 1 for tag, c in zip(tags, clusters)}
            for rec_id, position in id2pos.items():
                record = align[position]
                if rec_id in cluster_map:
                    record.letter_annotations['cluster'][i] = cluster_map[rec_id]

    # ------------------------------------------------------------------
    # Identify fully conserved columns (all structures in same cluster)
    # ------------------------------------------------------------------
    align_len = align.get_alignment_length()
    column_labels = [
        {align[pos].letter_annotations['cluster'][i] for pos in id2pos.values()}
        for i in range(align_len)
    ]
    area = [
        1 if len(labels) == 1 and _UNCLUSTERED not in labels else 0
        for labels in column_labels
    ]

    # ------------------------------------------------------------------
    # Find anchor regions (contiguous conserved runs ≥ _MIN_ANCHOR_LENGTH)
    # ------------------------------------------------------------------
    # FIX [perf]: _find_runs replaces three near-identical while-loop blocks.
    anchors: List[Tuple[int, int]] = [
        run for run in _find_runs(area, value=1)
        if run[1] - run[0] >= _MIN_ANCHOR_LENGTH
    ]
    anchor_lengths = [end - ini for ini, end in anchors]

    # Initialise group and ca_coords annotations
    for position in id2pos.values():
        record = align[position]
        record.letter_annotations['group']     = [0] * len(record.seq)
        record.letter_annotations['ca_coords'] = [[0.0, 0.0, 0.0]] * len(record.seq)

    for ini, end in anchors:
        for position in id2pos.values():
            for i in range(ini, end):
                align[position].letter_annotations['group'][i] = 1

    # Select top anchor regions for annealing
    top_anchors = list(np.argsort(anchor_lengths)[::-1][:max_anchors])

    # ------------------------------------------------------------------
    # Load Cα coordinates
    # ------------------------------------------------------------------
    # FIX [style]: QUIET=True suppresses noisy BioPython warnings.
    parser = PDBParser(QUIET=True)
    ca_coords:  Dict[str, np.ndarray] = {}
    ca_masked:  Dict[str, np.ndarray] = {}
    structures: Dict[str, object]     = {}

    for record in align:
        if 'structure' not in record.description:
            continue
        # FIX [style]: Renamed 'id' → 'rec_id' to avoid shadowing the built-in.
        rec_id = record.id
        pdb_path = f'{structure_dir}/{rec_id}.pdb'
        structures[rec_id] = parser.get_structure(rec_id, pdb_path)

        model = structures[rec_id][0]
        xyz: List[np.ndarray] = []
        seq: List[str]        = []
        for chain in model:
            for residue in chain:
                xyz.append(residue['CA'].get_coord())
                seq.append(seq1(residue.get_resname()))

        j = 0
        cds: List[np.ndarray] = []
        msk: List[bool]       = []

        for i in range(align_len):
            if record.seq[i] == '-':
                cds.append(np.zeros(3))
                msk.append(False)
            elif j < len(seq) and record.seq[i] == seq[j]:
                cds.append(xyz[j])
                msk.append(True)
                j += 1
            else:
                # FIX [correctness]: Append placeholder so array length stays
                # consistent with the alignment; previously this branch left
                # the arrays short, causing downstream index misalignment.
                print(f'Warning: sequence mismatch at column {i} — '
                      f'alignment={record.seq[i]}, PDB={seq[j] if j < len(seq) else "?"}')
                cds.append(np.zeros(3))
                msk.append(False)

        ca_coords[rec_id] = np.array(cds)
        ca_masked[rec_id] = np.array(msk)

    # ------------------------------------------------------------------
    # Simulated annealing superimposition
    # ------------------------------------------------------------------
    ids = list(ca_coords.keys())
    ref = ids[0]
    print(f'Reference: {ref}')

    xo = ca_coords[ref]
    for rec_id in ids[1:]:
        yo = ca_coords[rec_id]
        rot, tran, e0 = simulated_annealing(xo, yo, anchors, top_anchors, rng=rng)
        ca_coords[rec_id] = np.dot(yo, rot) + tran
        print(f'{rec_id}: {e0:.2f} Å')

        for model in structures[rec_id]:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom.set_coord(np.dot(atom.get_coord(), rot) + tran)

    # Store transformed Cα coordinates in the alignment record
    for rec_id in ca_coords:
        record = align[id2pos[rec_id]]
        for i, masked in enumerate(ca_masked[rec_id]):
            if masked:
                record.letter_annotations['ca_coords'][i] = ca_coords[rec_id][i].tolist()

    # ------------------------------------------------------------------
    # Remove anchor positions where any pair of structures diverges > outlier_dist
    # ------------------------------------------------------------------
    # FIX [perf]: Stack all structures into a matrix and use pdist for a
    # vectorised pairwise distance computation instead of an O(N²) double loop.
    struct_ids = list(id2pos.keys())
    for ini, end in anchors:
        for i in range(ini, end):
            coords_at_i = np.array([ca_coords[sid][i] for sid in struct_ids])
            if pdist(coords_at_i).max() > outlier_dist:
                for position in id2pos.values():
                    align[position].letter_annotations['group'][i] = 0

    # ------------------------------------------------------------------
    # Prune short conserved runs (< _MIN_REGION_LENGTH)
    # ------------------------------------------------------------------
    group_area = [
        1 if all(
            align[pos].letter_annotations['group'][i] == 1
            for pos in id2pos.values()
        ) else 0
        for i in range(align_len)
    ]

    for ini, end in _find_runs(group_area, value=1):
        if end - ini < _MIN_REGION_LENGTH:
            for i in range(ini, end):
                for position in id2pos.values():
                    align[position].letter_annotations['group'][i] = 0

    # ------------------------------------------------------------------
    # Prune short non-conserved gaps (< _MIN_REGION_LENGTH) — fill them in
    # ------------------------------------------------------------------
    group_area = [
        1 if all(
            align[pos].letter_annotations['group'][i] == 0
            for pos in id2pos.values()
        ) else 0
        for i in range(align_len)
    ]

    for ini, end in _find_runs(group_area, value=1):
        if end - ini < _MIN_REGION_LENGTH:
            for i in range(ini, end):
                for position in id2pos.values():
                    align[position].letter_annotations['group'][i] = 1

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    if csv:
        for rec_id in ca_coords:
            record = align[id2pos[rec_id]]
            with open(f'{rec_id}.csv', 'w') as f:
                for i, masked in enumerate(ca_masked[rec_id]):
                    if masked:
                        x, y, z = ca_coords[rec_id][i]
                        g = record.letter_annotations['group'][i]
                        f.write(f'{i},{x:.4f},{y:.4f},{z:.4f},{g}\n')

    if pdb:
        io = PDBIO()
        for rec_id in ca_coords:
            io.set_structure(structures[rec_id])
            io.save(f'{rec_id}_sup.pdb')


# ---------------------------------------------------------------------------
# Terminal alignment viewer
# ---------------------------------------------------------------------------

def show_align(align, pal) -> None:
    """
    Print the alignment to the terminal, coloured by structural group using
    the provided Seaborn palette.

    :param align: Clustered MultipleSeqAlignment with 'group' annotations
    :param pal: Seaborn colour palette (list of (R, G, B) float triples)
    """
    length = align.get_alignment_length()

    for j in range(length // 50 + 1):
        ini = 50 * j
        end = min(ini + 50, length)

        # FIX [perf]: Build each ruler row as a string, print once per row.
        if end > 100:
            ruler_h = '      ' + ''.join(
                f'{i:03d}'[0] if i % 10 == 0 else ' ' for i in range(ini, end)
            )
            print(ruler_h)

        ruler_t = '      ' + ''.join(
            f'{i:03d}'[1] if i % 10 == 0 else ' ' for i in range(ini, end)
        )
        ruler_u = '      ' + ''.join(
            f'{i:03d}'[2] if i % 10 == 0 else ' ' for i in range(ini, end)
        )
        print(ruler_t)
        print(ruler_u)

        for record in align:
            if 'structure' not in record.description:
                continue

            # FIX [perf]: Accumulate coloured characters, print once per record.
            row = record.id + ':'
            for i in range(ini, end):
                g = record.letter_annotations['group'][i] - 1
                ch = record.seq[i]
                if g >= 0:
                    R = int(pal[g][0] * 255)
                    G = int(pal[g][1] * 255)
                    B = int(pal[g][2] * 255)
                    row += fg(R, G, B) + ch
                else:
                    row += fg(125, 125, 125) + ch.lower()
            print(row + fg.rs)
