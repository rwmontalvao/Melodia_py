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
import math
import warnings

import numpy as np

from sty import fg

from numpy.random import rand, randint

from Bio.PDB import PDBParser, PDBIO
from Bio.SeqUtils import seq1
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering

warnings.filterwarnings('ignore', category=PDBConstructionWarning)


def rmsd(x: np.ndarray, y: np.ndarray) -> float:
    """
         Compute the RMSD between the two numpy arrays.
         :param x: The first numpy array
         :type x: np.ndarray
         :param y: The second numpy array
         :type y: np.ndarray
         :return: The RMSD between the two numpy arrays
         :rtype: float
    """

    return np.sqrt((((x - y) ** 2) * 3).mean())


def select(data: np.ndarray, seg: list, msk: list) -> np.ndarray:
    """
    Select the data in the given segment from the mask.
    :param data: The data to select
    :type data: np.ndarray
    :param seg: The segment to select
    :type seg: list
    :param msk: The mask for aligned segment
    :type msk: list
    :return: The selected data
    :rtype: np.ndarray
    """
    tmp = []
    for i, j in enumerate(seg):
        if msk[i]:
            tmp.append(data[j])

    return np.array(tmp)


def superposition(xo: np.ndarray, yo: np.ndarray, seg: list, msk: list) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the superposition of the two numpy arrays.
    :param xo: The first numpy array
    :type xo: np.ndarray
    :param yo: The second numpy array
    :type yo: np.ndarray
    :param seg: The segment to superimpose
    :type seg: list
    :param msk: The mask for the aligned segment
    :type msk: list
    :return: The superposition of the two numpy arrays
    :rtype: tuple(np.ndarray, np.ndarray)
    """
    x = select(xo, seg, msk)
    y = select(yo, seg, msk)

    sup = SVDSuperimposer()

    sup.set(x, y)
    sup.run()

    return sup.get_rotran()


def energy(xo: np.ndarray, yo: np.ndarray, seg: list, msk: list) -> float:
    """
    Compute the energy between two numpy arrays.
    :param xo: The first numpy array
    :type xo: np.ndarray
    :param yo: The second numpy array
    :type yo: np.ndarray
    :param seg: The segment to superimpose
    :type seg: list
    :param msk: The mask for the aligned segment
    :type msk: list
    :return: The energy between two numpy arrays
    :rtype: float
    """
    rot, tran = superposition(xo, yo, seg, msk)

    yt = np.dot(yo, rot) + tran

    return rmsd(xo, yt)


def segments(anchors: list, members: list) -> list:
    """
    Return a list of segments for each member of the protein alignment of an anchor.
    :param anchors: The list of anchors
    :type anchors: list
    :param members: The list of members
    :type members: list
    :return: The of segments
    :rtype: list
    """
    seg = []
    for member in members:
        ini, end = anchors[member]
        seg.extend([x for x in range(ini, end)])

    return seg


def simulated_annealing(xo: np.ndarray,
                        yo: np.ndarray,
                        anchors: list,
                        members: list) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Compute the rotation and translation matrices that reduces the RMSD
    between two protein segments.
    :param xo: The first numpy array
    :type xo: np.ndarray
    :param yo: The second numpy array
    :type yo: np.ndarray
    :param anchors: The list of anchors
    :type anchors: list
    :param members: The list of members
    :type members: list
    :return: The rotation and translation matrices and the final RMSD
    :rtype: tuple(np.ndarray, np.ndarray, float)
    """
    seg = segments(anchors, members)

    msk0 = [True for _ in range(len(seg))]

    temperature = 300
    energy0 = energy(xo, yo, seg, msk0)
    while temperature > 0.00001:
        for i in range(10000):
            msk1 = msk0.copy()
            j = randint(len(msk1))
            if msk1[j]:
                msk1[j] = False
            else:
                msk1[j] = True

            energy1 = energy(xo, yo, seg, msk1)

            energy_delta = energy1 - energy0

            if energy_delta < 0.0 or rand() < math.exp(-energy_delta / temperature):
                msk0 = msk1.copy()
                energy0 = energy1

        # TODO: test other schedules
        temperature *= 0.1

    rot, tran = superposition(xo, yo, seg, msk0)
    return rot, tran, energy0


def superimposer(align, threshold=0.8, csv=False, pdb=True):
    """
    Given an alignment, superimpose the PDB structure using Differential Geometry
    as dissimilarity measure. It can output CSV files with the alpha carbon coordinates
    and cluster groups, and also save the superimposed PDB structures.
    :param align: Alignment to superimpose
    :type align: Bio.Align.MultipleSeqAlignment
    :param threshold: Threshold for cluster selection
    :type threshold: float
    :param csv: If true, output the CSV files
    :type csv: bool
    :param pdb: If true, output the superimposed PDB structures
    :type pdb: bool
    """
    # Collect Differential Geometry data for initial clustering
    data = []
    id2pos = {}

    for position, record in enumerate(align):
        if 'structure' in record.description:
            id2pos[record.id] = position
            pairs = zip(
                record.letter_annotations['curvature'],
                record.letter_annotations['torsion'],
            )
            data = [list(pair) for pair in pairs]
            record.letter_annotations['cluster'] = [0 for _ in record.seq]

    scaler = StandardScaler()
    scaler.fit(data)

    # Clustering to find maximal conserved regions
    clustering = AgglomerativeClustering(distance_threshold=threshold, n_clusters=None)
    for i in range(align.get_alignment_length()):
        xy = []
        tags = []
        for identity, position in id2pos.items():
            record = align[position]
            if record.seq[i] != '-':
                xy.append(
                    [
                        record.letter_annotations['curvature'][i],
                        record.letter_annotations['torsion'][i],
                    ]
                )
                tags.append(identity)

        if len(xy) > 1:
            xy = scaler.transform(xy)

            clusters = clustering.fit_predict(xy)

            map_of_clusters = {pair[0]: pair[1] for pair in zip(tags, clusters)}

            for identity, position in id2pos.items():
                record = align[position]
                if identity in map_of_clusters:
                    record.letter_annotations['cluster'][i] = (
                            map_of_clusters[identity] + 1
                    )
                    # print(id, record.letter_annotations['cluster'][i])
        else:
            for identity, position in id2pos.items():
                record = align[position]
                record.letter_annotations['cluster'][i] = 0

    # Find the highly conserved areas
    ini = 0
    end = align.get_alignment_length()

    area = []
    for i in range(ini, end):
        s = set()
        for identity, position in id2pos.items():
            record = align[position]
            j = record.letter_annotations['cluster'][i]
            s = s.union({j})
        # Only regions conserved across all proteins
        if len(s) == 1:
            area.append(1)
        else:
            area.append(0)

    # Find all anchor conserved regions
    anchors = []
    anchor_lengths = []
    new_area = [0 for _ in area]

    ini, end = 0, 0
    while ini < len(area):
        while ini < len(area):
            if area[ini] == 1:
                break
            else:
                ini += 1

        end = ini
        while end < len(area):
            if area[end] != 1:
                break
            else:
                end += 1

        # Minimal length for an anchor regions is 5 residues long
        if end - ini >= 5:
            anchor_lengths.append(end - ini)
            anchors.append((ini, end))
            for i in range(ini, end):
                new_area[i] = 1
        ini = end

    # Group is the annotation for the new clusters
    for identity, position in id2pos.items():
        record = align[position]
        record.letter_annotations['group'] = [0 for _ in record.seq]
        record.letter_annotations['ca_coords'] = [[0.0, 0.0, 0.0] for _ in record.seq]

    for anchor in anchors:
        for identity, position in id2pos.items():
            record = align[position]
            ini, end = anchor
            for i in range(ini, end):
                record.letter_annotations['group'][i] = 1

    # Select the top 5 anchor regions
    top_anchors = np.argsort(anchor_lengths)[::-1]
    if len(top_anchors) > 5:
        top_anchors = top_anchors[:5]

    # Read the CA coordinates for all the files
    parser = PDBParser()

    ca_coords = {}
    ca_masked = {}
    structures = {}
    for record in align:
        if 'structure' not in record.description:
            continue
        id = record.id
        structures[id] = parser.get_structure(id, f'{id}.pdb')

        model = structures[id][0]
        xyz = []
        seq = []
        for chain in model:
            for residue in chain:
                xyz.append(residue['CA'].get_coord())
                seq.append(seq1(residue.get_resname()))
        j = 0
        cds = []
        msk = []
        for i in range(align.get_alignment_length()):
            # Insure coordinate alignment
            if record.seq[i] == '-':
                cds.append(np.array([0.0, 0.0, 0.0]))
                msk.append(False)
            elif record.seq[i] == seq[j]:
                cds.append(xyz[j])
                msk.append(True)
                j += 1
            else:
                print(f'Error: {i} {record.seq[i]} - {j} {seq[j]}')
        ca_coords[id] = np.array(cds)
        ca_masked[id] = np.array(msk)

    # Optimise the anchor positions for minimal RMSD
    # TODO: Use a dendrogram for the order?
    ids = list(ca_coords.keys())
    ref = ids.pop(0)
    print(f'{ref}')

    xo = ca_coords[ref]
    for idx in ids:
        yo = ca_coords[idx]
        rot, tran, E0 = simulated_annealing(xo, yo, anchors, top_anchors)
        ca_coords[idx] = np.dot(yo, rot) + tran
        print(f'{idx}: {E0:.2f} Å')

        for model in structures[idx]:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        ya = atom.get_coord()
                        atom.set_coord(np.dot(ya, rot) + tran)

    # Store the CA coordinates in the record
    for idx in ca_coords.keys():
        for i, m in enumerate(ca_masked[idx]):
            if m:
                record = align[id2pos[idx]]
                x = ca_coords[idx][i][0]
                y = ca_coords[idx][i][1]
                z = ca_coords[idx][i][2]
                record.letter_annotations['ca_coords'][i] = [x, y, z]

    # Remove group outliers
    for anchor in anchors:
        ini, end = anchor
        for i in range(ini, end):
            max_dist = 0.0
            for j in id2pos.keys():
                cj = ca_coords[j][i]
                for k in id2pos.keys():
                    if k != j:
                        ck = ca_coords[k][i]
                        dist = np.linalg.norm(ck - cj)
                        if dist > max_dist:
                            max_dist = dist

            # Any region with at least one distance
            # over the threshold is removed
            if max_dist > 2.0:
                for j in id2pos.values():
                    record = align[j]
                    record.letter_annotations['group'][i] = 0

    # Remove conserved regions with less than 3 residues
    ini = 0
    end = align.get_alignment_length()

    area = []
    for i in range(ini, end):
        s = set()
        for position in id2pos.values():
            record = align[position]
            j = record.letter_annotations['group'][i]
            if j == 1:
                s = s.union({j})
        if len(s) == 1:
            area.append(1)
        else:
            area.append(0)

    regions = []
    regions_lengths = []
    new_area = [0 for _ in area]

    ini, end = 0, 0
    while ini < len(area):
        while ini < len(area):
            if area[ini] == 1:
                break
            else:
                ini += 1

        end = ini
        while end < len(area):
            if area[end] != 1:
                break
            else:
                end += 1
        delta = end - ini
        if 0 < delta < 3:
            regions_lengths.append(delta)
            regions.append((ini, end))
            for i in range(ini, end):
                new_area[i] = 1
        ini = end

    for region in regions:
        ini, end = region
        for i in range(ini, end):
            for position in id2pos.values():
                record = align[position]
                record.letter_annotations['group'][i] = 0

    # Remove non-conserved regions with less than 3 residues
    ini = 0
    end = align.get_alignment_length()

    area = []
    for i in range(ini, end):
        s = set()
        for position in id2pos.values():
            record = align[position]
            j = record.letter_annotations['group'][i]
            if j == 0:
                s = s.union({j})
        if len(s) == 1:
            area.append(1)
        else:
            area.append(0)

    regions = []
    regions_lengths = []
    new_area = [0 for _ in area]

    ini, end = 0, 0
    while ini < len(area):
        while ini < len(area):
            if area[ini] == 1:
                break
            else:
                ini += 1

        end = ini
        while end < len(area):
            if area[end] != 1:
                break
            else:
                end += 1
        delta = end - ini
        if 0 < delta < 3:
            regions_lengths.append(delta)
            regions.append((ini, end))
            for i in range(ini, end):
                new_area[i] = 1
        ini = end

    for region in regions:
        ini, end = region
        for i in range(ini, end):
            for position in id2pos.values():
                record = align[position]
                record.letter_annotations['group'][i] = 1

    # Output CSV coordinate files
    if csv:
        for idx in ca_coords.keys():
            record = align[id2pos[idx]]
            with open(f'{idx}.csv', 'w') as f:
                for i, m in enumerate(ca_masked[idx]):
                    if m:
                        x = ca_coords[idx][i][0]
                        y = ca_coords[idx][i][1]
                        z = ca_coords[idx][i][2]
                        g = record.letter_annotations['group'][i]
                        f.write(f'{i},{x:.4f},{y:.4f},{z:.4f},{g}\n')

    # Output superposed PDB files
    if pdb:
        io = PDBIO()
        for idx in ca_coords.keys():
            io.set_structure(structures[idx])
            io.save(f'{idx}_sup.pdb')

    return


def show_align(align, pal):
    """
    Show the alignment colored by cluster using a Seaborn pallet
    :param align: Alignment to show
    :type align: Bio.Align.MultipleSeqAlignment
    :param pal: Palette
    :type pal: seaborn.palettes._ColorPalette
    """
    length = align.get_alignment_length()
    total = int(length / 50.0)
    for j in range(total + 1):
        ini = 0 + 50 * j
        end = 49 + 50 * j + 1

        if end >= length:
            end = length - 1

        if end > 100:
            print('      ', end='')
            for i in range(ini, end):
                s = f'{i:03n}'
                if i % 10 == 0:
                    print(s[0], end='')
                else:
                    print(' ', end='')
            print()

        print('      ', end='')
        for i in range(ini, end):
            s = f'{i:03n}'
            if i % 10 == 0:
                print(s[1], end='')
            else:
                print(' ', end='')
        print()

        print('      ', end='')
        for i in range(ini, end):
            s = f'{i:03n}'
            if i % 10 == 0:
                print(s[2], end='')
            else:
                print(' ', end='')
        print()

        for record in align:
            if 'structure' not in record.description:
                continue

            print(record.id, end=':')
            for i in range(ini, end):
                j = record.letter_annotations['group'][i] - 1
                if j >= 0:
                    R = int(pal[j][0] * 255)
                    G = int(pal[j][1] * 255)
                    B = int(pal[j][2] * 255)
                    print(fg(R, G, B) + record.seq[i], end='')
                else:
                    print(fg(125, 125, 125) + record.seq[i].lower(), end='')
            print(fg.rs)
    return
