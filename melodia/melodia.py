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
import os
import Bio.Align
import pkg_resources

import pandas as pd
import nglview as nv

from typing import Dict
from ipywidgets import Box

from Bio import AlignIO
from Bio.PDB import PDBParser
from Bio.PDB.Structure import Structure

from melodia.geometryparser import GeometryParser


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
    for model in structure.get_list():
        for chain in model.get_list():
            chain_id = chain.id
            key = f'{pdb_code}:{model.id}:{chain_id}'
            chains[key] = GeometryParser(chain)

    # Generate a Pandas Dataframe
    identity = []
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
    for model in structure.get_list():
        for chain in model.get_list():
            chain_id = chain.id
            key = f'{model.id}:{chain_id}'
            chains[key] = GeometryParser(chain)

    return chains


def bfactor_from_geo(structure: Structure, attribute: str) -> None:
    """
        Set the PDB bfactor to torsion values
         :param structure: BioPython PDB Structure class
         :type structure: Structure
         :param attribute: geometric attribute
         :type attribute: str
    """
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
                    else:
                        value = 0.0

                    atom.set_bfactor(value)


def view_putty(structure: Structure, radius_scale=1.0, width=1200, height=600) -> Box:
    """
        Set the PDB bfactor to torsion values
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

        chain = record.description.split(':')[3]

        file_name = f'{record.id}.pdb'

        structure = parser.get_structure(record.id, file_name)

        geo = geometry_dict_from_structure(structure)

        curvature = []
        torsion = []
        arc_length = []
        writhing = []
        phi = []
        psi = []

        j = 0
        for i, letter in enumerate(record.seq):
            if letter != '-':
                res_code = c321[geo[f'0:{chain}'].residues[j].name]
                if letter != res_code:
                    raise NameError(f'Alignment error: res={i+1} seq={letter} pdb={res_code}')
                curvature.append(geo[f'0:{chain}'].residues[j].curvature)
                torsion.append(geo[f'0:{chain}'].residues[j].torsion)
                arc_length.append(geo[f'0:{chain}'].residues[j].arc_len)
                writhing.append(geo[f'0:{chain}'].residues[j].writhing)
                phi.append(geo[f'0:{chain}'].residues[j].phi)
                psi.append(geo[f'0:{chain}'].residues[j].psi)
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
