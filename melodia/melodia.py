# Copyright 2021 KU Leuven.
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
# TODO: Better error messages
#
import os
from typing_extensions import SupportsIndex

import Bio.Align
import pkg_resources

import pandas as pd
import nglview as nv
import seaborn as sns

from typing import Dict, Tuple, List, Any
from ipywidgets import Box

from Bio import AlignIO
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure

from melodia.geometryparser import GeometryParser

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering


def geometry_from_structure_file(file_name: str) -> pd.DataFrame:
    """
         Function used to compute the geometric properties around residues
         for (BioPython compatible) structures.
         It computes curvature, torsion, arc-length and writhing number
         :param file_name: Protein file name
         :type file_name: String
         :return: Pandas dataframe
         :rtype: pd.DataFrame
    """
    parser = PDBParser()
    name, ext = os.path.splitext(file_name)
    structure = parser.get_structure(name, file_name)

    base = proc_chains(structure)

    return base


def geometry_from_structure(structure: Structure) -> pd.DataFrame:
    """
         Function used to compute the geometric properties around residues
         for (BioPython compatible) structures.
         It computes curvature, torsion, arc-length and writhing number
         :param structure: BioPython PDB Structure class
         :type structure: Structure
         :return: Panda dataframe
         :rtype: pd.DataFrame
    """
    base = proc_chains(structure)

    return base


def proc_chains(structure: Structure) -> pd.DataFrame:
    """
         Function used to compute the geometric properties around residues
         for (BioPython compatible) structures.
         It computes curvature, torsion, arc-length and writhing number
         :param structure: BioPython PDB Structure class
         :type structure: Structure
         :return: Panda dataframe
         :rtype: pd.DataFrame
    """
    pdb_code = structure.id.upper()
    chains = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            key = f'{pdb_code}:{model.id}:{chain_id}'
            chains[key] = GeometryParser(chain)

    # Generate a Pandas Dataframe
    identity = []
    model = []
    chain = []
    order = []
    name = []
    curvature = []
    torsion = []
    arc_length = []
    writhing = []
    angle_phi = []
    angle_psi = []
    for key, value in chains.items():
        # print(key)
        tokens = key.split(':')

        for res_id in value.residues:
            model.append(int(tokens[1]))
            chain.append(tokens[2])

            identity.append(value.residues[res_id].res_num)
            order.append(value.residues[res_id].res_order)
            name.append(value.residues[res_id].name)
            curvature.append(value.residues[res_id].curvature)
            torsion.append(value.residues[res_id].torsion)
            arc_length.append(value.residues[res_id].arc_len)
            writhing.append(value.residues[res_id].writhing)
            angle_phi.append(value.residues[res_id].phi)
            angle_psi.append(value.residues[res_id].psi)

    base = {'id': identity,
            'model': model,
            'code': pdb_code,
            'chain': chain,
            'order': order,
            'name': name,
            'curvature': curvature,
            'torsion': torsion,
            'arc_length': arc_length,
            'writhing': writhing,
            'phi': angle_phi,
            'psi': angle_psi}

    return pd.DataFrame(base)


def geometry_dict_from_structure(structure: Structure) -> Dict[str, GeometryParser]:
    """
         Function used to compute the geometric properties around residues
         for (BioPython compatible) structures.
         It computes curvature, torsion, arc-length and writhing number
         :param structure: BioPython PDB Structure class
         :type structure: Structure
         :return: Geometry Dictionary
         :rtype: Dict[str, GeometryParser]
    """
    chains = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id
            key = f'{model.id}:{chain_id}'
            residues = [res.id[1] for res in list(chain.get_residues()) if res.id[0] == ' ']
            if residues:
                chains[key] = GeometryParser(chain)

    return chains


def bfactor_from_geo(structure: Structure, attribute: str, geo: Dict[str, GeometryParser] = None) -> None:
    """
        Set the PDB bfactor to geometric values
         :param structure: BioPython PDB Structure class
         :type structure: Structure
         :param attribute: geometric attribute
         :type attribute: str
        :param geo: Geometry Dictionary
        :type: Dict[str, GeometryParser]

    """
    if geo is None:
        geo = geometry_dict_from_structure(structure)

    for model in structure.get_list():
        for chain in model.get_list():
            for atom in chain.get_atoms():
                het_flag, sequence_id, insertion_code = atom.get_parent().id
                if het_flag[0] == " ":
                    key = f'{model.id}:{chain.id}'
                    res = geo[key].residues_map[sequence_id]

                    if attribute == 'curvature':
                        value = geo[key].residues[res].curvature
                    elif attribute == 'torsion':
                        value = geo[key].residues[res].torsion
                    elif attribute == 'custom':
                        value = geo[key].residues[res].custom
                    else:
                        value = 0.0

                    atom.set_bfactor(value)


def view_putty(structure: Structure, radius_scale=1.0, width=1200, height=600) -> Box:
    """
        Display PDB structure as a putty model
        :param structure: BioPython PDB Structure class
        :type structure: Structure
        :param radius_scale: Radius Scale
        :type width: float
        :param width: Widget width
        :type width: int
        :param height: Widget height
        :type height: int
    """
    view = nv.show_biopython(structure)
    view.representations = [
        {"type": "tube",
         "params": {
             "sele": "protein",
             "radius": "bfactor",
             "radiusScale": radius_scale,
             "color": "bfactor",
             "colorScale": "RdYlBu"
         }}
    ]
    view.layout.width = 'auto'
    view.layout.height = 'auto'

    box = Box([view])
    box.layout.width = f'{width}px'
    box.layout.height = f'{height}px'

    return box


def view_cartoon(structure: Structure, width=1200, height=600) -> Box:
    """
        Display PDB structure as a cartoon model
        :param structure: BioPython PDB Structure class
        :type structure: Structure
        :param width: Widget width
        :type width: int
        :param height: Widget height
        :type height: int
    """
    view = nv.show_biopython(structure)
    view.representations = [
        {"type": "cartoon",
         "params": {
             "sele": "protein",
             "color": "bfactor",
             "colorScale": "RdYlBu"
         }}
    ]
    view.layout.width = 'auto'
    view.layout.height = 'auto'

    box = Box([view])
    box.layout.width = f'{width}px'
    box.layout.height = f'{height}px'

    return box


def view_tube(structure: Structure, width=1200, height=600) -> Box:
    """
        Display PDB structure as a tube model
        :param structure: BioPython PDB Structure class
        :type structure: Structure
        :param width: Widget width
        :type width: int
        :param height: Widget height
        :type height: int
    """
    view = nv.show_biopython(structure)
    view.representations = [
        {"type": "tube",
         "params": {
             "sele": "protein",
             "color": "bfactor",
             "colorScale": "RdYlBu"
         }}
    ]
    view.layout.width = 'auto'
    view.layout.height = 'auto'

    box = Box([view])
    box.layout.width = f'{width}px'
    box.layout.height = f'{height}px'

    return box


def parser_pir_file(pir_file: str) -> Bio.Align.MultipleSeqAlignment:
    """
        Parser for PIR alignment file and associated PDB structures
        :param pir_file: BioPython PDB Structure class
        :type pir_file: str
        :return: Protein sequence alignment
        :rtype: Bio.Align.MultipleSeqAlignment
    """
    c321 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    align = AlignIO.read(pir_file, 'pir')

    parser = PDBParser()
    for record in align:
        # Only select structures
        if record.description[:9] != 'structure':
            continue

        file_name = f'{record.id}.pdb'

        structure = parser.get_structure(record.id, file_name)

        # Renumber the residues to start with 1
        residue_number = 1
        for model in structure:
            for chain in model:
                for residue in chain:
                    # print(residue.id)
                    residue.id = (residue.id[0], residue_number, 'Z')
                    residue_number += 1
                for residue in chain:
                    residue.id = (residue.id[0], residue.id[1], ' ')
                    # print('----', residue.id)

        geo = geometry_dict_from_structure(structure)

        curvature = []
        torsion = []
        arc_length = []
        writhing = []
        phi = []
        psi = []

        idx = []
        for key in geo:
            for res in geo[key].residues:
                idx.append((key, res))
        j = 0
        for i, letter in enumerate(record.seq):
            if letter != '-' and letter != '/':
                try:
                    curr_chain = idx[j][0]
                    curr_residue = idx[j][1]
                    res_code = c321[geo[curr_chain].residues[curr_residue].name]
                except KeyError:
                    raise NameError(f'Alignment error: {record.id} {curr_chain} residue {curr_residue + 1} {letter}')
                if letter != res_code:
                    raise NameError(f'Alignment error: {record.id} res={i + 1} seq={letter} pdb={res_code}')
                curvature.append(geo[curr_chain].residues[curr_residue].curvature)
                torsion.append(geo[curr_chain].residues[curr_residue].torsion)
                arc_length.append(geo[curr_chain].residues[curr_residue].arc_len)
                writhing.append(geo[curr_chain].residues[curr_residue].writhing)
                phi.append(geo[curr_chain].residues[curr_residue].phi)
                psi.append(geo[curr_chain].residues[curr_residue].psi)
                j += 1
            else:
                curvature.append(0.0)
                torsion.append(0.0)
                arc_length.append(0.0)
                writhing.append(0.0)
                phi.append(0.0)
                psi.append(0.0)

        record.letter_annotations['curvature'] = curvature
        record.letter_annotations['torsion'] = torsion
        record.letter_annotations['arc_length'] = arc_length
        record.letter_annotations['writhing'] = writhing
        record.letter_annotations['phi'] = phi
        record.letter_annotations['psi'] = psi

    return align


def dataframe_from_alignment(align: Bio.Align.MultipleSeqAlignment, keys=None) -> pd.DataFrame:
    """
         Create a Pandas Dataframe from an annotated alignment
         :param align: Protein alignment
         :type align: Bio.Align.MultipleSeqAlignment
         :param keys: keys for the annotations
         :type keys: list
         :return: Pandas dataframe
         :rtype: pd.DataFrame
    """
    data = {}

    for record in align:
        if record.description[:9] == 'structure':
            data[f'seq_{record.id}'] = list(record.seq)
            if keys is None:
                items = record.letter_annotations.keys()
            else:
                items = keys
            for key in items:
                data[f'{key}_{record.id}'] = record.letter_annotations[key]

    return pd.DataFrame.from_dict(data)


class PropensityTable:
    """
    Class for the Protein Propensity Table data
    """
    __slots__ = ('__ini', '__end', '__lines', '__data')

    def __init__(self):
        """
        Load the Propensity Table data
        """

        def find_block():
            found = False
            while self.__ini < self.__end:
                if self.__lines[self.__ini][0] == '>':
                    found = True
                    break
                self.__ini += 1
            return found

        def load_block():
            tag = self.__lines[self.__ini].split()
            self.__ini += 1
            rows = []
            while self.__lines[self.__ini][0] != 'U':
                rows.append(self.__lines[self.__ini])
                self.__ini += 1
            rows.append(self.__lines[self.__ini])
            return int(tag[1]), rows

        stream = pkg_resources.resource_stream(__name__, 'data/luthier.dat')
        lines = stream.readlines()

        # Remove newlines
        for i, line in enumerate(lines):
            lines[i] = lines[i].decode('utf-8')
            if lines[i][-1] == '\n':
                lines[i] = lines[i][:-1]

        # TODO: trim data after use
        self.__ini = 0
        self.__end = len(lines)
        self.__lines = lines
        self.__data = {}

        while find_block():
            key, table = load_block()

            dct = {}
            head = table[0].split()
            for i in range(1, len(table)):
                row = table[i].split()
                key1 = row[0]
                for j in range(1, len(row)):
                    key2 = head[j]
                    key12 = f'{key1},{key2}'
                    value = int(row[j])
                    dct[key12] = value

            self.__data[key] = dct

    def get_score(self, target: str, residue: str, phi: float, psi: float) -> int:
        """
        Find anomalies in the protein's chain
        :param target: Target residue 1 letter code
        :type target: str
        :param residue: Residue 1 letter code
        :type residue: str
        :param phi: phi angle (degrees)
        :type phi: float
        :param psi: psi angle (degrees)
        :type psi: float
        :return: Score
        :rtype: int
        """
        max_phi = [0.0, 0.0, 0.0, -110.0, -110.0, 0.0, 140.0, 180.0, 180.0]
        min_phi = [-180.0, -110.0, -110.0, -180.0, -180.0, -180.0, 20.0, 0.0, 0.0]
        max_psi = [45.0, 180.0, -90.0, 180.0, -90.0, 100.0, 80.0, -40.0, 180.0]
        min_psi = [-90.0, 100.0, -180.0, 100.0, -180.0, 45.0, -40.0, -180.0, 80.0]

        tab_map = [0, 1, 1, 2, 2, 3, 4, 5, 5]

        score = 0
        for i in range(len(tab_map)):
            if (min_phi[i] <= phi < max_phi[i]) and (min_psi[i] <= psi < max_psi[i]):
                key = f'{target.upper()},{residue.upper()}'
                score = self.__data[3][key]
                break

        return score


def cluster_alignment(align: Bio.Align.MultipleSeqAlignment, threshold=0.7, long=False) -> None:
    """
         Cluster the alignments by structural similarity
         :param align: Protein alignment
         :type align: Bio.Align.MultipleSeqAlignment
         :param threshold: similarity threshold
         :type threshold: float
         :param long: output only long cluster (length > 3)
         :type long: bool
    """
    data = []
    id2pos = {}
    for position, record in enumerate(align):
        if 'structure' in record.description:
            id2pos[record.id] = position
            for pair in zip(record.letter_annotations['curvature'], record.letter_annotations['torsion']):
                data.append(list(pair))
            record.letter_annotations['cluster'] = [0 for _ in record.seq]

    scaler = StandardScaler()
    scaler.fit(data)

    clustering = AgglomerativeClustering(distance_threshold=threshold, n_clusters=None)
    for i in range(align.get_alignment_length()):
        xy = []
        tags = []
        for identity, position in id2pos.items():
            record = align[position]
            if record.seq[i] != '-':
                xy.append([record.letter_annotations['curvature'][i], record.letter_annotations['torsion'][i]])
                tags.append(identity)

        xy = scaler.transform(xy)

        clusters = clustering.fit_predict(xy)

        map_of_clusters = {pair[0]: pair[1] for pair in zip(tags, clusters)}

        for identity, position in id2pos.items():
            record = align[position]
            if identity in map_of_clusters:
                record.letter_annotations['cluster'][i] = map_of_clusters[identity]
                # print(id, record.letter_annotations['cluster'][i])

    last_cluster = 0
    for record in align:
        if 'structure' not in record.description:
            continue
        max_cluster = max(record.letter_annotations['cluster'])
        if max_cluster > last_cluster:
            last_cluster = max_cluster

    for i in range(align.get_alignment_length() - 1):
        j = i + 1
        left = {}
        right = {}
        for k, record in enumerate(align):
            if 'structure' not in record.description:
                continue

            ca = record.letter_annotations['cluster'][i]
            cb = record.letter_annotations['cluster'][j]

            if ca not in left:
                left[ca] = [k]
            else:
                left[ca].append(k)

            if cb not in right:
                right[cb] = [k]
            else:
                right[cb].append(k)

        for key, value in left.items():
            left[key] = set(value)

        for key, value in right.items():
            right[key] = set(value)

        for right_key in right:
            found_key = None
            for left_key in left:
                if right[right_key].symmetric_difference(left[left_key]) == set():
                    found_key = left_key

            if found_key is None:
                last_cluster += 1
                for k in right[right_key]:
                    record = align[k]
                    record.letter_annotations['cluster'][j] = last_cluster
            else:
                for k in right[right_key]:
                    record = align[k]
                    record.letter_annotations['cluster'][j] = found_key

    data, idx = get_idx(align)

    if long:
        last_cluster = 0
        for j in idx:
            cluster, ini, end, size = data[j]
            if size < 3:
                for i in range(ini, end + 1):
                    for record in align:
                        if 'structure' not in record.description:
                            continue
                        if record.letter_annotations['cluster'][i] == cluster:
                            record.letter_annotations['cluster'][i] = -1
                data[j] = (-1, ini, end, size)
            else:
                # print(f'{ini}-{end} {cluster}->{p}')
                for i in range(ini, end + 1):
                    for record in align:
                        if 'structure' not in record.description:
                            continue
                        if record.letter_annotations['cluster'][i] == cluster:
                            record.letter_annotations['cluster'][i] = last_cluster

                data[j] = (last_cluster, ini, end, size)
                last_cluster += 1
    return


def save_pymol_script(align: Bio.Align.MultipleSeqAlignment, pml_file: str, palette='Dark2', colors=7) -> None:
    """
         Cluster the alignments by structural similarity
         :param align: Protein alignment
         :type align: Bio.Align.MultipleSeqAlignment
         :param pml_file: Pymol script file name (pml)
         :type pml_file: str
         :param palette: colour palette
         :type palette: str
         :param colors: number of colours in the palette
         :type colors: int
    """
    data, idx = get_idx(align)

    tags = []
    pal = sns.color_palette(palette, colors).as_hex()
    with open(f'{pml_file}.pml', 'w') as f:
        f.write('# Script generated by Melodia\n\n')
        f.write('# load structures\n')
        for record in align:
            if 'structure' not in record.description:
                continue
            tags.append(record.id)
            f.write(f'load {record.id}.pdb\n')
        f.write('\n')

        f.write('# superimpose structures\n')
        for i in range(1, len(tags)):
            f.write(f'super {tags[i]}, {tags[0]}\n')
        f.write('\n')

        f.write('# non-conserved cluster color\n')
        f.write(f'color gray40\n\n')
        f.write('# cluster colors\n')
        for i in idx:
            positions = data[i]
            cluster, ini, end, size = positions
            if cluster >= 0:
                for record in align:
                    if 'structure' not in record.description:
                        continue
                    if record.letter_annotations['cluster'][ini] == cluster:
                        colour = f'0x{pal[cluster % colors][1:]}'
                        f.write(f'color {colour}, {record.id} and resi {ini + 1}-{end + 1}\n')
                f.write('\n')
    return


def get_idx(align: Bio.Align.MultipleSeqAlignment) -> Tuple[List[Tuple[Any, int, int, int]], List[SupportsIndex]]:
    """
    Internal function for clustering
    :param align: Protein alignment
    :return: Indexed clusters
    """
    clusters = {}
    for i, record in enumerate(align):
        if 'structure' not in record.description:
            continue
        for j, cluster in enumerate(record.letter_annotations['cluster']):
            if cluster not in clusters:
                clusters[cluster] = [j]
            else:
                clusters[cluster].append(j)
    data = []
    block_init = []
    for key, value in clusters.items():
        data.append((key, min(value), max(value), max(value) - min(value) + 1))
        block_init.append(min(value))
    idx = sorted(range(len(block_init)), key=block_init.__getitem__)
    return data, idx


def save_align_to_ps(align: Bio.Align.MultipleSeqAlignment, ps_file: str, palette='Dark2', colors=7) -> None:
    """
         Cluster the alignments by structural similarity
         :param align: Protein alignment
         :type align: Bio.Align.MultipleSeqAlignment
         :param ps_file: post-script file name
         :type ps_file: str
         :param palette: colour palette
         :type palette: str
         :param colors: number of colours in the palette
         :type colors: int
    """

    pal = sns.color_palette(palette, colors)
    rgb = [f'{c[0]:4.2f} {c[1]:4.2f} {c[2]:4.2f}' for c in pal]

    black = '0.00 0.00 0.00'
    grey = '0.50 0.50 0.50'

    out_file = f'{ps_file}.ps'
    with open(out_file, 'w') as ps:
        ps.write('%%!PS-Adobe-3.0\n')
        ps.write('%%%%Pages: 1\n')
        ps.write('%%%%Creator: Melodia 1.0\n')
        ps.write('%%%%CreationDate:\n')
        ps.write('%%%%EndComments\n')
        ps.write('%%%%Page: 1 1\n')
        ps.write('/Courier-Regular findfont  16.0 scalefont  setfont\n')
        ps.write('0.00 0.00 0.83 setrgbcolor\n')
        ps.write(f'72.0 735.0 moveto ({ps_file}) show\n')

        length = align.get_alignment_length()
        count = len(align)
        total = int(length / 50)
        block = count + 4
        blocks_per_page = int(76 / block) + 1

        page = 1
        line = 705.0
        blocks = 0
        position = 10

        for j in range(total + 1):
            ini = 0 + 50 * j
            end = 49 + 50 * j

            if end >= length:
                end = length - 1

            column = 203.0

            ps.write('/Courier-Regular findfont  8.0 scalefont  setfont\n')
            ps.write('0.00 0.00 0.83 setrgbcolor\n')

            for k in range(1, 6):
                if position < 100:
                    ps.write('%5.1f %5.1f moveto (%d) show\n' % (column, line, position))
                else:
                    ps.write('%5.1f %5.1f moveto (%d) show\n' % (column - 2.0, line, position))

                position += 10
                column += 80.0

                if position > length:
                    break

            i = 1
            for record in align:
                column = 52.0
                line -= 10.0
                ps.write('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                ps.write('0.00 0.00 0.00 setrgbcolor\n')
                ps.write('%5.1f %5.1f moveto (%3d) show\n' % (column, line, i))

                column += 20.0

                ps.write('%5.1f %5.1f moveto (%s) show\n' % (column, line, record.id))

                column = 131.1
                bold = False
                last_color = black
                for cur_res in range(ini, end + 1):
                    if record.seq[cur_res] == '-':
                        if bold:
                            bold = False
                            ps.write('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                        if last_color != black:
                            last_color = black
                            ps.write(f'{last_color} setrgbcolor\n')
                    else:
                        if 'structure' in record.description:
                            k = record.letter_annotations['cluster'][cur_res]
                            if k >= 0:
                                c = k % colors
                                if not bold:
                                    bold = True
                                    ps.write('/Courier-Bold findfont  10.0 scalefont  setfont\n')
                                if last_color != rgb[c]:
                                    last_color = rgb[c]
                                    ps.write(f'{last_color} setrgbcolor\n')
                            else:
                                if bold:
                                    bold = False
                                    ps.write('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                                if last_color != grey:
                                    last_color = grey
                                    ps.write(f'{last_color} setrgbcolor\n')
                        else:
                            if bold:
                                bold = False
                                ps.write('/Courier-Regular findfont  10.0 scalefont  setfont\n')
                            if last_color != black:
                                last_color = black
                                ps.write(f'{last_color} setrgbcolor\n')
                    if 'structure' in record.description:
                        if record.letter_annotations['cluster'][cur_res] < 0:
                            ps.write(f'%5.1f %5.1f moveto (%c) show\n' % (column, line, record.seq[cur_res].lower()))
                        else:
                            ps.write(f'%5.1f %5.1f moveto (%c) show\n' % (column, line, record.seq[cur_res]))
                    else:
                        ps.write(f'%5.1f %5.1f moveto (%c) show\n' % (column, line, record.seq[cur_res].lower()))
                    column += 8.0
                i += 1

            line -= 20.0
            blocks += 1

            if blocks == blocks_per_page and j != total:
                blocks = 0
                line = 705.0
                page += 1

                ps.write('showpage\n')
                ps.write(f'%%%%Page: {page} {page}\n')
                ps.write('/Courier-Regular findfont  16.0 scalefont  setfont\n')
                ps.write('0.00 0.00 0.83 setrgbcolor\n')
                ps.write(f'72.0 735.0 moveto (out_file) show\n')

        ps.write('showpage\n')

    return
