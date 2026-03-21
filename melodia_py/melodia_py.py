# Copyright 2021-2024 KU Leuven.
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
import os
from typing import Dict, List, Optional, Tuple, Any

import Bio.Align
import pandas as pd
import nglview as nv
import seaborn as sns

from ipywidgets import Box
from Bio import AlignIO
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure
from joblib import Parallel, delayed

from melodia_py.geometryparser import GeometryParser

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering
from importlib import resources as importlib_resources


# ---------------------------------------------------------------------------
# Ramachandran bin definitions for PropensityTable.get_score
# Bins derived from: Kleywegt & Jones (1996) Ramachandran revisited.
# Each entry is (phi_min, phi_max, psi_min, psi_max).
# ---------------------------------------------------------------------------
_RAMA_BINS: List[Tuple[float, float, float, float]] = [
    (-180.0,   0.0,  -90.0,   45.0),
    (-110.0,   0.0,  100.0,  180.0),
    (-110.0,   0.0, -180.0,  -90.0),
    (-180.0, -110.0,  100.0,  180.0),
    (-180.0, -110.0, -180.0,  -90.0),
    (-180.0,   0.0,   45.0,  100.0),
    (  20.0, 140.0,  -40.0,   80.0),
    (   0.0, 180.0, -180.0,  -40.0),
    (   0.0, 180.0,   80.0,  180.0),
]
_RAMA_TAB_MAP: List[int] = [0, 1, 1, 2, 2, 3, 4, 5, 5]


# ---------------------------------------------------------------------------
# Structure → DataFrame helpers
# ---------------------------------------------------------------------------

def geometry_from_structure_file(
    file_name: str,
    n_jobs: int = 1,
    rna_atom: str = "C4'",
    deg: bool = True,
) -> pd.DataFrame:
    """
    Parse a PDB file and compute per-residue geometric properties.

    For multi-model files (NMR ensembles, MD trajectories), set ``n_jobs=-1``
    to parallelise across all models × chains using all available CPUs.

    :param file_name: Path to the PDB file
    :param n_jobs: Number of parallel worker processes.
        ``1`` = sequential (default). ``-1`` = all CPUs.
    :param rna_atom: Backbone atom for RNA chains (default ``"C4'"``)
    :param deg: Return phi/psi in degrees (True) or radians (False)
    :return: DataFrame with curvature, torsion, arc_length, writhing, phi, psi
    """
    name, _ = os.path.splitext(os.path.basename(file_name))
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(name, file_name)
    return proc_chains(structure, n_jobs=n_jobs, rna_atom=rna_atom, deg=deg)


def geometry_from_structure(
    structure: Structure,
    n_jobs: int = 1,
    rna_atom: str = "C4'",
    deg: bool = True,
) -> pd.DataFrame:
    """
    Compute per-residue geometric properties from a BioPython Structure object.

    :param structure: BioPython PDB Structure
    :param n_jobs: Number of parallel worker processes.
        ``1`` = sequential (default). ``-1`` = all CPUs.
    :param rna_atom: Backbone atom for RNA chains (default ``"C4'"``)
    :param deg: Return phi/psi in degrees (True) or radians (False)
    :return: DataFrame with curvature, torsion, arc_length, writhing, phi, psi
    """
    return proc_chains(structure, n_jobs=n_jobs, rna_atom=rna_atom, deg=deg)


# ---------------------------------------------------------------------------
# Model coordinate extraction (runs in parent process, BioPython-safe)
# ---------------------------------------------------------------------------

def _extract_model_data(
    model,
    pdb_code: str,
    rna_atom: str,
) -> Optional[dict]:
    """
    Extract all coordinate data for every chain in one model into plain
    Python/numpy objects that are cheap to pickle.

    Runs in the **parent** process. Workers receive only the returned dict —
    no BioPython objects cross the process boundary.

    :return: Dict with keys: model_id, pdb_code, chains (list of chain dicts)
             Each chain dict contains: chain_id, rna, coords, residue_meta,
             backbone (list of (N, CA, C) triples or empty for RNA).
    """
    from melodia_py.geometryparser import _RNA_RESIDUE_NAMES

    model_id = int(model.id)
    chain_data = []

    for chain in model:
        all_res = [r for r in chain.get_residues() if r.id[0] == ' ']
        if not all_res:
            continue

        rna  = all_res[0].get_resname() in _RNA_RESIDUE_NAMES
        atom = rna_atom if rna else 'CA'

        needs_fallback = rna and atom in ("P", "C5'")
        last_pos = all_res[-1].id[1] if needs_fallback else None

        coords:       List[List[float]] = []
        residue_meta: List[Tuple[str, int]] = []
        backbone:     List[Optional[Tuple]] = []

        for res in all_res:
            pos = res.id[1]

            if atom in res:
                coord = res[atom].get_coord().tolist()
            elif needs_fallback and pos == last_pos:
                coord = next(res.get_atoms()).get_coord().tolist()
            else:
                print(f'Warning: missing {atom} at {res.get_resname()} {pos} '
                      f'chain {chain.id} model {model_id} — skipped')
                continue

            coords.append(coord)
            residue_meta.append((res.get_resname(), pos))

            if not rna:
                try:
                    n  = res['N'].get_coord().tolist()
                    ca = res['CA'].get_coord().tolist()
                    c  = res['C'].get_coord().tolist()
                    backbone.append((n, ca, c))
                except KeyError:
                    backbone.append(None)

        if len(coords) < 2:
            continue

        chain_data.append({
            'chain_id':     chain.id,
            'rna':          rna,
            'coords':       coords,
            'residue_meta': residue_meta,
            'backbone':     backbone,
        })

    if not chain_data:
        return None

    return {
        'model_id':  model_id,
        'pdb_code':  pdb_code,
        'chains':    chain_data,
    }


# ---------------------------------------------------------------------------
# Worker function — processes all chains of one model (no BioPython objects)
# ---------------------------------------------------------------------------

def _compute_model_geometry(task: dict, deg: bool) -> List[dict]:
    """
    Compute geometry for every chain in one model.

    Receives only plain dicts/lists — no BioPython objects.
    Safe to run in a separate loky worker process.

    :param task: Output of _extract_model_data()
    :param deg: Return dihedral angles in degrees
    :return: Flat list of per-residue record dicts for all chains
    """
    import numpy as np
    from scipy.interpolate import CubicSpline
    from melodia_py.geometryparser import (
        GeometryParser, _NORM_EPS, _calc_writhing_jit,
    )

    model_id = task['model_id']
    pdb_code = task['pdb_code']
    rows: List[dict] = []

    for chain in task['chains']:
        chain_id     = chain['chain_id']
        rna          = chain['rna']
        coords       = chain['coords']
        residue_meta = chain['residue_meta']
        backbone     = chain['backbone']

        n_res = len(coords)
        t = list(range(n_res))
        x = [c[0] for c in coords]
        y = [c[1] for c in coords]
        z = [c[2] for c in coords]

        xt = CubicSpline(t, x, bc_type='natural')
        yt = CubicSpline(t, y, bc_type='natural')
        zt = CubicSpline(t, z, bc_type='natural')

        xa = np.ascontiguousarray(x, dtype=np.float64)
        ya = np.ascontiguousarray(y, dtype=np.float64)
        za = np.ascontiguousarray(z, dtype=np.float64)

        ini, end = 0, n_res - 1

        for i, (resname, pos) in enumerate(residue_meta):
            p_curv = float(t[1]) if i == ini else (float(t[-2]) if i == end else float(i))

            curvature, torsion = GeometryParser.calc_curvature_torsion(
                p=p_curv, t=t, xt=xt, yt=yt, zt=zt
            )
            arc_len  = GeometryParser.calc_arc_length(p=float(i), xt=xt, yt=yt, zt=zt)
            writhing = (
                float(_calc_writhing_jit(i, xa, ya, za, _NORM_EPS))
                if not rna else 0.0
            )

            row: dict = {
                'id':         i,
                'model':      model_id,
                'code':       pdb_code,
                'chain':      chain_id,
                'order':      pos,
                'name':       resname,
                'curvature':  curvature,
                'torsion':    torsion,
                'arc_length': arc_len,
                'writhing':   writhing,
            }

            if not rna:
                phi: Optional[float] = None
                psi: Optional[float] = None

                bb = backbone[i]
                if bb is not None:
                    n_coord  = np.array(bb[0])
                    ca_coord = np.array(bb[1])
                    c_coord  = np.array(bb[2])

                    if i > 0 and backbone[i - 1] is not None:
                        prev_c = np.array(backbone[i - 1][2])
                        phi = GeometryParser.calc_dihedral_torsion(
                            p1=prev_c, p2=n_coord, p3=ca_coord, p4=c_coord, deg=deg
                        )
                    if i < n_res - 1 and backbone[i + 1] is not None:
                        next_n = np.array(backbone[i + 1][0])
                        psi = GeometryParser.calc_dihedral_torsion(
                            p1=n_coord, p2=ca_coord, p3=c_coord, p4=next_n, deg=deg
                        )
                row['phi'] = phi
                row['psi'] = psi

            rows.append(row)

    return rows


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def proc_chains(
    structure: Structure,
    n_jobs: int = 1,
    rna_atom: str = "C4'",
    deg: bool = True,
) -> pd.DataFrame:
    """
    Compute per-residue geometry for every model and chain in *structure*,
    returning a tidy DataFrame.

    Parallelism is **across models** — each worker process handles all chains
    of one model. This is the natural unit of work for multi-model files
    (NMR ensembles, MD trajectories) where models are independent conformations.

    BioPython Chain objects are never sent to workers. Coordinates are
    extracted as plain lists in the parent process and dispatched as small
    dicts, keeping pickle overhead minimal.

    :param structure: BioPython PDB Structure
    :param n_jobs: Number of parallel worker processes.
        ``1`` = sequential (default). ``-1`` = all CPUs.
        For a 300-model ensemble with 8 cores, expect ~8× speedup.
    :param rna_atom: Backbone atom for RNA chains (default ``"C4'"``)
    :param deg: Return phi/psi in degrees (True) or radians (False)
    :return: DataFrame with one row per residue
    """
    pdb_code = structure.id.upper()

    # Extract coordinate data for every model in the parent process.
    # BioPython objects never leave this process.
    tasks = []
    for model in structure:
        task = _extract_model_data(model, pdb_code, rna_atom)
        if task is not None:
            tasks.append(task)

    if not tasks:
        return pd.DataFrame()

    # Dispatch one task per model to the loky worker pool.
    # Each worker processes all chains of its model sequentially.
    nested: List[List[dict]] = Parallel(n_jobs=n_jobs)(
        delayed(_compute_model_geometry)(task, deg)
        for task in tasks
    )

    records = [row for model_rows in nested for row in model_rows]
    return pd.DataFrame(records)



def _compute_model_gp(task: dict, deg: bool) -> List[Tuple[str, dict]]:
    """
    Worker: compute geometry for all chains of one model and return
    (key, residues_dict) pairs that can be used to reconstruct results.
    Receives only plain dicts — no BioPython objects.
    """
    import numpy as np
    from scipy.interpolate import CubicSpline
    from melodia_py.geometryparser import (
        GeometryParser, ResidueGeometry, _NORM_EPS, _calc_writhing_jit,
    )

    model_id = task['model_id']
    results  = []

    for chain in task['chains']:
        chain_id     = chain['chain_id']
        rna          = chain['rna']
        coords       = chain['coords']
        residue_meta = chain['residue_meta']
        backbone     = chain['backbone']

        n_res = len(coords)
        t = list(range(n_res))
        x = [c[0] for c in coords]
        y = [c[1] for c in coords]
        z = [c[2] for c in coords]

        xt = CubicSpline(t, x, bc_type='natural')
        yt = CubicSpline(t, y, bc_type='natural')
        zt = CubicSpline(t, z, bc_type='natural')

        xa = np.ascontiguousarray(x, dtype=np.float64)
        ya = np.ascontiguousarray(y, dtype=np.float64)
        za = np.ascontiguousarray(z, dtype=np.float64)

        ini, end = 0, n_res - 1
        residues:     dict = {}
        residues_map: dict = {}

        for i, (resname, pos) in enumerate(residue_meta):
            p_curv = float(t[1]) if i == ini else (float(t[-2]) if i == end else float(i))

            curvature, torsion = GeometryParser.calc_curvature_torsion(
                p=p_curv, t=t, xt=xt, yt=yt, zt=zt
            )
            arc_len  = GeometryParser.calc_arc_length(p=float(i), xt=xt, yt=yt, zt=zt)
            writhing = (
                float(_calc_writhing_jit(i, xa, ya, za, _NORM_EPS))
                if not rna else 0.0
            )

            rg = ResidueGeometry(
                name=resname,
                chain=chain_id,
                res_num=i,
                res_order=pos,
                curvature=curvature,
                torsion=torsion,
                arc_len=arc_len,
                writhing=writhing,
            )

            if not rna:
                bb = backbone[i]
                if bb is not None:
                    n_coord  = np.array(bb[0])
                    ca_coord = np.array(bb[1])
                    c_coord  = np.array(bb[2])
                    if i > 0 and backbone[i - 1] is not None:
                        prev_c = np.array(backbone[i - 1][2])
                        rg.phi = GeometryParser.calc_dihedral_torsion(
                            p1=prev_c, p2=n_coord, p3=ca_coord, p4=c_coord, deg=deg
                        )
                    if i < n_res - 1 and backbone[i + 1] is not None:
                        next_n = np.array(backbone[i + 1][0])
                        rg.psi = GeometryParser.calc_dihedral_torsion(
                            p1=n_coord, p2=ca_coord, p3=c_coord, p4=next_n, deg=deg
                        )

            residues[i]   = rg
            residues_map[pos] = i

        key = f'{model_id}:{chain_id}'
        results.append((key, residues, residues_map, rna))

    return results


def geometry_dict_from_structure(
    structure: Structure,
    n_jobs: int = 1,
    rna_atom: str = "C4'",
    deg: bool = True,
) -> Dict[str, GeometryParser]:
    """
    Build a mapping of 'model_id:chain_id' → GeometryParser for all chains
    that contain at least one standard residue.

    For multi-model structures (NMR ensembles, MD trajectories) set
    ``n_jobs=-1`` to parallelise across models.

    :param structure: BioPython PDB Structure
    :param n_jobs: Number of parallel worker processes (default 1)
    :param rna_atom: Backbone atom for RNA chains (default ``"C4'"``)
    :param deg: Return phi/psi in degrees (True) or radians (False)
    :return: Dict of GeometryParser objects keyed by 'model_id:chain_id'
    """
    pdb_code = structure.id.upper()

    # Extract coordinate data in the parent process
    tasks = []
    for model in structure:
        task = _extract_model_data(model, pdb_code, rna_atom)
        if task is not None:
            tasks.append(task)

    if not tasks:
        return {}

    # Dispatch to workers — one task per model
    nested: List[List[tuple]] = Parallel(n_jobs=n_jobs)(
        delayed(_compute_model_gp)(task, deg)
        for task in tasks
    )

    # Reconstruct GeometryParser objects from the worker results.
    # We bypass __init__ and set internals directly since we already
    # have all computed data — no need to reparse the BioPython chain.
    chains: Dict[str, GeometryParser] = {}
    for model_results in nested:
        for key, residues, residues_map, rna in model_results:
            gp = object.__new__(GeometryParser)
            object.__setattr__(gp, '_GeometryParser__residues',     residues)
            object.__setattr__(gp, '_GeometryParser__residues_map', residues_map)
            object.__setattr__(gp, '_GeometryParser__degrees',      deg)
            object.__setattr__(gp, '_GeometryParser__gap_list',     [])
            object.__setattr__(gp, '_GeometryParser__anomaly_list', [])
            object.__setattr__(gp, '_GeometryParser__rna_atom',     rna_atom)
            object.__setattr__(gp, 'RNA',                           rna)
            chains[key] = gp

    return chains



# ---------------------------------------------------------------------------
# B-factor mapping
# ---------------------------------------------------------------------------

def bfactor_from_geo(
    structure: Structure,
    attribute: str,
    geo: Optional[Dict[str, GeometryParser]] = None,
    n_jobs: int = 1,
    rna_atom: str = "C4'",
) -> None:
    """
    Set every atom's B-factor to a geometric property value.

    :param structure: BioPython PDB Structure (mutated in place)
    :param attribute: One of 'curvature', 'torsion', 'custom'
    :param geo: Pre-computed geometry dict; computed from structure if None
    :param n_jobs: Number of parallel workers used if geo must be computed
    :param rna_atom: Backbone atom for RNA chains if geo must be computed
    """
    if geo is None:
        geo = geometry_dict_from_structure(structure, n_jobs=n_jobs, rna_atom=rna_atom)

    def _get_val(gp: GeometryParser, res_idx: int) -> float:
        res = gp.residues[res_idx]
        if attribute == 'curvature':
            return res.curvature
        if attribute == 'torsion':
            return res.torsion
        if attribute == 'custom':
            return res.custom
        return 0.0

    # Determine global minimum for default fill
    min_value = min(
        _get_val(gp, res_idx)
        for gp in geo.values()
        for res_idx in gp.residues
    ) if attribute in ('curvature', 'torsion', 'custom') else 0.0

    for atom in structure.get_atoms():
        if atom.is_disordered():
            for da in atom.disordered_get_list():
                da.set_bfactor(min_value)
        else:
            atom.set_bfactor(min_value)

    for model in structure:
        for chain in model:
            for atom in chain.get_atoms():
                het_flag, sequence_id, _ = atom.get_parent().id
                if het_flag[0] == ' ':
                    key = f'{model.id}:{chain.id}'
                    res_idx = geo[key].residues_map[sequence_id]
                    atom.set_bfactor(_get_val(geo[key], res_idx))


# ---------------------------------------------------------------------------
# NGL viewers
# ---------------------------------------------------------------------------

def _make_view(structure: Structure, representation: dict, width: int, height: int) -> Box:
    view = nv.show_biopython(structure)
    view.representations = [representation]
    view.layout.width = '100%'
    view.layout.height = '100%'
    box = Box([view])
    box.layout.width = f'{width}px'
    box.layout.height = f'{height}px'
    return box


def view_putty(structure: Structure, radius_scale: float = 1.0, width: int = 1200, height: int = 600) -> Box:
    """Display PDB structure as a putty (tube-radius-by-bfactor) model."""
    return _make_view(structure, {
        'type': 'tube',
        'params': {
            'sele': 'protein',
            'radius': 'bfactor',
            'radiusScale': radius_scale,
            'color': 'bfactor',
            'colorScale': 'RdYlBu',
        },
    }, width, height)


def view_cartoon(structure: Structure, width: int = 1200, height: int = 600) -> Box:
    """Display PDB structure as a cartoon coloured by B-factor."""
    return _make_view(structure, {
        'type': 'cartoon',
        'params': {'sele': 'protein', 'color': 'bfactor', 'colorScale': 'RdYlBu'},
    }, width, height)


def view_tube(structure: Structure, width: int = 1200, height: int = 600) -> Box:
    """Display PDB structure as a tube coloured by B-factor."""
    return _make_view(structure, {
        'type': 'tube',
        'params': {'sele': 'protein', 'color': 'bfactor', 'colorScale': 'RdYlBu'},
    }, width, height)


# ---------------------------------------------------------------------------
# PIR alignment parser
# ---------------------------------------------------------------------------

# Standard 3-letter to 1-letter amino acid code mapping
_C321: Dict[str, str] = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
}


def parser_pir_file(
    pir_file: str,
    structure_dir: str = '.',
) -> Bio.Align.MultipleSeqAlignment:
    """
    Parse a PIR alignment file and annotate each structure record with
    per-position geometric properties.

    :param pir_file: Path to the PIR alignment file
    :param structure_dir: Directory containing the PDB files referenced in the
        alignment. Defaults to the current working directory.
    :return: Annotated MultipleSeqAlignment
    """
    # FIX [style]: Accept structure_dir so callers are not forced to run from
    # the directory that contains the PDB files.
    align = AlignIO.read(pir_file, 'pir')
    parser = PDBParser(QUIET=True)

    for record in align:
        if not record.description.startswith('structure'):
            continue

        pdb_path = os.path.join(structure_dir, f'{record.id}.pdb')
        structure = parser.get_structure(record.id, pdb_path)

        # Renumber residues starting at 1 (two-pass to avoid id collisions)
        residue_number = 1
        for model in structure:
            for chain in model:
                for residue in chain:
                    residue.id = (residue.id[0], residue_number, 'Z')
                    residue_number += 1
                for residue in chain:
                    residue.id = (residue.id[0], residue.id[1], ' ')

        geo = geometry_dict_from_structure(structure)

        # Build a flat ordered index of (chain_key, res_idx) pairs
        idx: List[Tuple[str, int]] = [
            (key, res)
            for key in geo
            for res in geo[key].residues
        ]

        # FIX [correctness]: Use None for gap positions so consumers can
        # distinguish "no angle computed" from a genuine angle of 0.0.
        curvature: List[Optional[float]] = []
        torsion:   List[Optional[float]] = []
        arc_length: List[Optional[float]] = []
        writhing:  List[Optional[float]] = []
        phi:       List[Optional[float]] = []
        psi:       List[Optional[float]] = []

        j = 0
        for letter in record.seq:
            if letter not in ('-', '/'):
                curr_chain, curr_residue = idx[j]
                res = geo[curr_chain].residues[curr_residue]
                try:
                    res_code = _C321[res.name]
                except KeyError:
                    raise NameError(
                        f'Alignment error: {record.id} {curr_chain} '
                        f'residue {curr_residue + 1} — unknown residue {res.name}'
                    )
                if letter != res_code:
                    raise NameError(
                        f'Alignment error: {record.id} seq={letter} pdb={res_code} '
                        f'at chain={curr_chain} res={curr_residue + 1}'
                    )
                curvature.append(res.curvature)
                torsion.append(res.torsion)
                arc_length.append(res.arc_len)
                writhing.append(res.writhing)
                phi.append(res.phi)
                psi.append(res.psi)
                j += 1
            else:
                curvature.append(None)
                torsion.append(None)
                arc_length.append(None)
                writhing.append(None)
                phi.append(None)
                psi.append(None)

        record.letter_annotations['curvature']  = curvature
        record.letter_annotations['torsion']    = torsion
        record.letter_annotations['arc_length'] = arc_length
        record.letter_annotations['writhing']   = writhing
        record.letter_annotations['phi']        = phi
        record.letter_annotations['psi']        = psi

    return align


def dataframe_from_alignment(
    align: Bio.Align.MultipleSeqAlignment,
    keys: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Build a DataFrame from an annotated alignment, with one column per
    (annotation × structure) combination.

    :param align: Annotated MultipleSeqAlignment (from parser_pir_file)
    :param keys: Annotation keys to include; defaults to all available keys
    :return: DataFrame
    """
    data: Dict[str, list] = {}
    for record in align:
        if not record.description.startswith('structure'):
            continue
        data[f'seq_{record.id}'] = list(record.seq)
        items = keys if keys is not None else list(record.letter_annotations.keys())
        for key in items:
            data[f'{key}_{record.id}'] = record.letter_annotations[key]
    return pd.DataFrame.from_dict(data)


# ---------------------------------------------------------------------------
# Propensity table
# ---------------------------------------------------------------------------

class PropensityTable:
    """
    Residue-pair propensity scores binned by Ramachandran (phi, psi) region.
    Data is loaded from the bundled luthier.dat resource file.
    """
    __slots__ = ('__data',)

    def __init__(self) -> None:
        """Load and parse the propensity table from the package data file."""
        ref = importlib_resources.files('melodia_py').joinpath('data/luthier.dat')
        with ref.open('rb') as fp:
            raw = fp.readlines()

        lines = [line.decode('utf-8').rstrip('\n') for line in raw]
        # FIX [style]: Parsing extracted to a private static method so __init__
        # has no inner functions closing over mutable state.
        self.__data = PropensityTable._parse_blocks(lines)

    @staticmethod
    def _parse_blocks(lines: List[str]) -> Dict[int, Dict[str, int]]:
        """Parse all > blocks from the luthier.dat line list."""
        data: Dict[int, Dict[str, int]] = {}
        i = 0
        n = len(lines)

        while i < n:
            # Seek next block header
            while i < n and (not lines[i] or lines[i][0] != '>'):
                i += 1
            if i >= n:
                break

            tag = lines[i].split()
            block_key = int(tag[1])
            i += 1

            rows: List[str] = []
            while i < n and lines[i] and lines[i][0] != 'U':
                rows.append(lines[i])
                i += 1
            if i < n:
                rows.append(lines[i])  # include the 'U...' terminator row
                i += 1

            head = rows[0].split()
            dct: Dict[str, int] = {}
            for row in rows[1:]:
                parts = row.split()
                row_key = parts[0]
                for j, col_key in enumerate(head[1:], start=1):
                    dct[f'{row_key},{col_key}'] = int(parts[j])
            data[block_key] = dct

        return data

    def get_score(
        self,
        target: str,
        residue: str,
        phi: float,
        psi: float,
    ) -> Optional[int]:
        """
        Return the propensity score for a residue pair in a given
        Ramachandran bin, or None if the (phi, psi) point falls outside
        all defined bins.

        :param target: Target residue 1-letter code
        :param residue: Query residue 1-letter code
        :param phi: Phi angle in degrees
        :param psi: Psi angle in degrees
        :return: Integer score, or None if no bin matches
        """
        # FIX [correctness]: Return None instead of 0 when no bin matches, so
        # callers can distinguish a genuine score of 0 from a missed bin.
        key = f'{target.upper()},{residue.upper()}'
        for i, (phi_min, phi_max, psi_min, psi_max) in enumerate(_RAMA_BINS):
            if phi_min <= phi < phi_max and psi_min <= psi < psi_max:
                return self.__data[3][key]
        return None


# ---------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------

def cluster_alignment(
    align: Bio.Align.MultipleSeqAlignment,
    threshold: float = 0.7,
    long: bool = False,
) -> None:
    """
    Cluster alignment positions by structural similarity in (curvature, torsion)
    space and annotate each structure record with a 'cluster' letter_annotation.

    :param align: Annotated MultipleSeqAlignment (mutated in place)
    :param threshold: Agglomerative clustering distance threshold
    :param long: If True, remove clusters shorter than 3 positions
    """
    # FIX [perf]: Pre-index structure records once instead of scanning align
    # on every column iteration.
    structure_records = {
        record.id: record
        for record in align
        if 'structure' in record.description
    }

    all_pairs: List[List[float]] = []
    for record in structure_records.values():
        record.letter_annotations['cluster'] = [0] * len(record.seq)
        for curv, tors in zip(
            record.letter_annotations['curvature'],
            record.letter_annotations['torsion'],
        ):
            all_pairs.append([curv, tors])

    scaler = StandardScaler()
    scaler.fit(all_pairs)

    clustering = AgglomerativeClustering(distance_threshold=threshold, n_clusters=None)

    for i in range(align.get_alignment_length()):
        xy: List[List[float]] = []
        tags: List[str] = []
        for rec_id, record in structure_records.items():
            if record.seq[i] != '-':
                xy.append([
                    record.letter_annotations['curvature'][i],
                    record.letter_annotations['torsion'][i],
                ])
                tags.append(rec_id)

        if len(xy) > 1:
            clusters = clustering.fit_predict(scaler.transform(xy))
            cluster_map = dict(zip(tags, clusters))
            for rec_id, record in structure_records.items():
                if rec_id in cluster_map:
                    record.letter_annotations['cluster'][i] = cluster_map[rec_id]
        else:
            for record in structure_records.values():
                record.letter_annotations['cluster'][i] = 0

    # Propagate consistent cluster labels across adjacent columns
    last_cluster = max(
        max(record.letter_annotations['cluster'])
        for record in structure_records.values()
    )

    for i in range(align.get_alignment_length() - 1):
        j = i + 1
        left: Dict[int, set] = {}
        right: Dict[int, set] = {}

        for k, record in enumerate(align):
            if 'structure' not in record.description:
                continue
            ca = record.letter_annotations['cluster'][i]
            cb = record.letter_annotations['cluster'][j]
            left.setdefault(ca, set()).add(k)
            right.setdefault(cb, set()).add(k)

        for right_key, right_members in right.items():
            found_key = next(
                (lk for lk, lm in left.items() if not lm.symmetric_difference(right_members)),
                None,
            )
            if found_key is None:
                last_cluster += 1
                for k in right_members:
                    align[k].letter_annotations['cluster'][j] = last_cluster
            else:
                for k in right_members:
                    align[k].letter_annotations['cluster'][j] = found_key

    if long:
        data, idx = get_idx(align)
        last_cluster = 0
        for j in idx:
            cluster, ini, end, size = data[j]
            new_label = -1 if size < 3 else last_cluster
            for i in range(ini, end + 1):
                for record in align:
                    if 'structure' not in record.description:
                        continue
                    if record.letter_annotations['cluster'][i] == cluster:
                        record.letter_annotations['cluster'][i] = new_label
            data[j] = (new_label, ini, end, size)
            if size >= 3:
                last_cluster += 1


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def save_pymol_script(
    align: Bio.Align.MultipleSeqAlignment,
    pml_file: str,
    palette: str = 'Dark2',
    colors: int = 7,
) -> None:
    """
    Write a PyMOL .pml script that loads, superimposes, and colours structures
    by cluster assignment.

    :param align: Clustered MultipleSeqAlignment
    :param pml_file: Output filename (without .pml extension)
    :param palette: Seaborn colour palette name
    :param colors: Number of colours to draw from the palette
    """
    data, idx = get_idx(align)
    tags = [r.id for r in align if 'structure' in r.description]
    pal = sns.color_palette(palette, colors).as_hex()

    with open(f'{pml_file}.pml', 'w') as f:
        f.write('# Script generated by Melodia\n\n')
        f.write('# load structures\n')
        for tag in tags:
            f.write(f'load {tag}.pdb\n')
        f.write('\n# superimpose structures\n')
        for tag in tags[1:]:
            f.write(f'super {tag}, {tags[0]}\n')
        f.write('\n# non-conserved cluster color\ncolor gray40\n\n# cluster colors\n')
        for i in idx:
            cluster, ini, end, _ = data[i]
            if cluster >= 0:
                for record in align:
                    if 'structure' not in record.description:
                        continue
                    if record.letter_annotations['cluster'][ini] == cluster:
                        color = f'0x{pal[cluster % colors][1:]}'
                        f.write(f'color {color}, {record.id} and resi {ini + 1}-{end + 1}\n')
                f.write('\n')
        f.write('\ncenter\n')


def get_idx(
    align: Bio.Align.MultipleSeqAlignment,
) -> Tuple[List[Tuple[Any, int, int, int]], List[int]]:
    """
    Build a sorted index of (cluster, ini, end, size) blocks from a clustered
    alignment.

    :param align: Clustered MultipleSeqAlignment
    :return: (data list, sorted index into data list)
    """
    clusters: Dict[int, List[int]] = {}
    for record in align:
        if 'structure' not in record.description:
            continue
        for j, cluster in enumerate(record.letter_annotations['cluster']):
            clusters.setdefault(cluster, []).append(j)

    data: List[Tuple[Any, int, int, int]] = []
    block_init: List[int] = []
    for key, positions in clusters.items():
        lo, hi = min(positions), max(positions)
        data.append((key, lo, hi, hi - lo + 1))
        block_init.append(lo)

    idx = sorted(range(len(block_init)), key=block_init.__getitem__)
    return data, idx


def save_align_to_ps(
    align: Bio.Align.MultipleSeqAlignment,
    ps_file: str,
    palette: str = 'Dark2',
    colors: int = 7,
) -> None:
    """
    Write a PostScript visualisation of the clustered alignment.

    :param align: Clustered MultipleSeqAlignment
    :param ps_file: Output filename (without .ps extension)
    :param palette: Seaborn colour palette name
    :param colors: Number of colours to draw from the palette
    """
    pal = sns.color_palette(palette, colors)
    rgb = [f'{c[0]:4.2f} {c[1]:4.2f} {c[2]:4.2f}' for c in pal]
    black = '0.00 0.00 0.00'
    grey  = '0.50 0.50 0.50'

    length = align.get_alignment_length()
    count  = len(align)
    total  = length // 50
    block  = count + 4
    blocks_per_page = int(76 / block) + 1

    out_file = f'{ps_file}.ps'
    # FIX [perf]: Accumulate lines into a buffer and write once per logical
    # section rather than calling f.write() for every single character.
    buf: List[str] = []

    def flush(f):
        f.write(''.join(buf))
        buf.clear()

    with open(out_file, 'w') as ps:
        buf += [
            '%%!PS-Adobe-3.0\n',
            '%%%%Pages: 1\n',
            '%%%%Creator: Melodia 1.0\n',
            '%%%%CreationDate:\n',
            '%%%%EndComments\n',
            '%%%%Page: 1 1\n',
            '/Courier-Regular findfont  16.0 scalefont  setfont\n',
            '0.00 0.00 0.83 setrgbcolor\n',
            f'72.0 735.0 moveto ({ps_file}) show\n',
        ]
        flush(ps)

        page = 1
        line = 705.0
        blocks = 0
        position = 10

        for j in range(total + 1):
            ini = 50 * j
            end = min(ini + 49, length - 1)

            column = 203.0
            buf.append('/Courier-Regular findfont  8.0 scalefont  setfont\n')
            buf.append('0.00 0.00 0.83 setrgbcolor\n')

            for k in range(1, 6):
                if position <= length:
                    x = column if position < 100 else column - 2.0
                    buf.append(f'{x:5.1f} {line:5.1f} moveto ({position}) show\n')
                    position += 10
                    column += 80.0

            for seq_i, record in enumerate(align, start=1):
                column = 52.0
                line -= 10.0
                row_chars: List[str] = [
                    '/Courier-Regular findfont  10.0 scalefont  setfont\n',
                    '0.00 0.00 0.00 setrgbcolor\n',
                    f'{column:5.1f} {line:5.1f} moveto ({seq_i:3d}) show\n',
                ]
                column += 20.0
                row_chars.append(f'{column:5.1f} {line:5.1f} moveto ({record.id}) show\n')

                column = 131.1
                bold = False
                last_color = black
                is_structure = 'structure' in record.description

                for cur_res in range(ini, end + 1):
                    ch = record.seq[cur_res]
                    if ch == '-':
                        if bold:
                            bold = False
                            row_chars.append('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                        if last_color != black:
                            last_color = black
                            row_chars.append(f'{black} setrgbcolor\n')
                        display_ch = ch
                    elif is_structure:
                        k = record.letter_annotations['cluster'][cur_res]
                        if k >= 0:
                            c = k % colors
                            color = rgb[c]
                            if not bold:
                                bold = True
                                row_chars.append('/Courier-Bold findfont  10.0 scalefont  setfont\n')
                            if last_color != color:
                                last_color = color
                                row_chars.append(f'{color} setrgbcolor\n')
                            display_ch = ch
                        else:
                            if bold:
                                bold = False
                                row_chars.append('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                            if last_color != grey:
                                last_color = grey
                                row_chars.append(f'{grey} setrgbcolor\n')
                            display_ch = ch.lower()
                    else:
                        if bold:
                            bold = False
                            row_chars.append('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                        if last_color != black:
                            last_color = black
                            row_chars.append(f'{black} setrgbcolor\n')
                        display_ch = ch.lower()

                    row_chars.append(f'{column:5.1f} {line:5.1f} moveto ({display_ch}) show\n')
                    column += 8.0

                buf.extend(row_chars)

            line -= 20.0
            blocks += 1
            flush(ps)

            if blocks == blocks_per_page and j != total:
                blocks = 0
                line = 705.0
                page += 1
                ps.write(f'showpage\n%%%%Page: {page} {page}\n')
                ps.write('/Courier-Regular findfont  16.0 scalefont  setfont\n')
                ps.write(f'0.00 0.00 0.83 setrgbcolor\n72.0 735.0 moveto ({out_file}) show\n')

        ps.write('showpage\n')
