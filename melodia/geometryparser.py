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
import math
import numpy as np

from typing import Dict, List, Tuple
from Bio.PDB import Chain
from dataclasses import dataclass
from numpy.polynomial import chebyshev
from scipy.interpolate import CubicSpline


@dataclass
class ResidueGeometry:
    # Residue information
    name: str = None
    chain: str = None
    res_num: int = None
    res_order: int = None

    # Differential geometry
    curvature: float = None
    torsion: float = None
    arc_len: float = None

    # Knot theory invariant
    writhing: float = None

    # Dihedral angles
    phi: float = None
    psi: float = None

    # Residue annotation
    res_ann: Dict[str, str] = None


class GeometryParser:
    """
    Class for parsing the geometrical properties of a protein chain
    """
    __slots__ = ('__residues',
                 '__residues_map',
                 '__degrees',
                 '__gap_list',
                 '__anomaly_list')

    def __init__(self, chain: Chain, deg: bool = True) -> None:
        """
        :param chain: Protein chain
        :type chain: Chain
        :param deg: Degree?
        :type deg: bool
        """
        residues, residues_map = GeometryParser.calc_geometry(chain=chain, deg=deg)
        self.__residues = residues
        self.__residues_map = residues_map
        self.__degrees = deg
        self.__gap_list = GeometryParser.find_gaps(chain=chain)
        self.__anomaly_list = GeometryParser.find_anomalies(chain=chain)

    @property
    def residues(self) -> Dict[int, ResidueGeometry]:
        """
        Access to residue geometry
        :return: dataclass for residue geometry
        :rtype: dict[int, ResidueGeometry]
        """
        return self.__residues

    @property
    def residues_map(self) -> Dict[int, ResidueGeometry]:
        """
        Maps residues number
        :return: Dictionary for the residue map
        :rtype: dict[int, int]
        """
        return self.__residues_map

    @property
    def deg(self) -> bool:
        """
        Are phi and psi in degrees?
        :return: True if phi and psi are in degrees, False for radians
        :rtype: bool
        """
        return self.__degrees

    @staticmethod
    def find_gaps(chain: Chain) -> List[Tuple[int, int]]:
        """
        Find gaps in the protein's chain
        :param chain: Protein chain
        :type chain: Chain
        :return: List of gaps
        :rtype: list[tuple[int, int]]
        """
        gaps = []

        i = 0
        residue = list(chain.get_residues())[i]
        het_flag, prev, insertion_code = residue.id

        for i in range(len(chain)):
            # get chain and pos
            residue = list(chain.get_residues())[i]
            het_flag, pos, insertion_code = residue.id

            if pos - prev > 1:
                gaps.append((prev, pos))

            prev = pos

        return gaps

    @property
    def gaps(self) -> List[Tuple[int, int]]:
        """
        Access to chain's gaps
        :return:  Chain's gaps
        :rtype: list[tuple[int, int]]
        """
        return self.__gap_list

    @property
    def gap(self) -> bool:
        """
        Are there  gaps in the chain?
        :return: True if gaps were found
        :rtype: bool
        """
        return len(self.__gap_list) > 0

    @staticmethod
    def find_anomalies(chain: Chain) -> List[str]:
        """
        Find anomalies in the protein's chain
        :param chain: Protein chain
        :type chain: Chain
        :return: List of anomalies
        :rtype: list[str]
        """
        # TODO: insert anomalies
        anomalies = []

        return anomalies

    @property
    def anomalies(self) -> List[str]:
        """
        Access to chain's anomalies
        :return:  Chain's anomalies
        :rtype: list[str]
        """
        return self.__anomaly_list

    @property
    def anomaly(self) -> bool:
        """
        Are there anomalies in the chain?
        :return: True if anomalies were found
        :rtype: bool
        """
        return len(self.__anomaly_list) > 0

    @staticmethod
    def calc_curvature_torsion(p: float,
                               t: List[float],
                               xt: CubicSpline,
                               yt: CubicSpline,
                               zt: CubicSpline) -> (float, float):
        """
        Function to compute Curvature and Torsion
        :param p: point of calculation
        :type p: float
        :param t: list os parameters
        :type t: list[float]
        :param xt: x(t)
        :type xt: CubicSpline
        :param yt: y(t)
        :type yt: CubicSpline
        :param zt: z(t)
        :type zt: CubicSpline
        :return: curvature, torsion
        :rtype: (float, float)
        """
        numpts = 51

        # TODO: Improve balance method
        mn = np.min(t)
        mx = np.max(t)
        for dt in range(1, 4):
            delta = float(dt)

            ini = p - delta
            end = p + delta

            if ini < mn:
                offset = mn - ini
            elif end > mx:
                offset = mx - end
            else:
                offset = 0.0

            ini += offset
            end += offset

            tp = np.linspace(ini, end, numpts)

            cxt, resxt = chebyshev.chebfit(tp, xt(tp), deg=10, full=True)
            cyt, resyt = chebyshev.chebfit(tp, yt(tp), deg=10, full=True)
            czt, reszt = chebyshev.chebfit(tp, zt(tp), deg=10, full=True)

            if resxt[0].size != 0 and resyt[0].size != 0 and reszt[0].size != 0:
                break

        cxtd1 = chebyshev.chebder(cxt, m=1)
        cytd1 = chebyshev.chebder(cyt, m=1)
        cztd1 = chebyshev.chebder(czt, m=1)

        cxtd2 = chebyshev.chebder(cxt, m=2)
        cytd2 = chebyshev.chebder(cyt, m=2)
        cztd2 = chebyshev.chebder(czt, m=2)

        cxtd3 = chebyshev.chebder(cxt, m=3)
        cytd3 = chebyshev.chebder(cyt, m=3)
        cztd3 = chebyshev.chebder(czt, m=3)

        xtd1 = chebyshev.chebval(p, cxtd1)
        ytd1 = chebyshev.chebval(p, cytd1)
        ztd1 = chebyshev.chebval(p, cztd1)

        xtd2 = chebyshev.chebval(p, cxtd2)
        ytd2 = chebyshev.chebval(p, cytd2)
        ztd2 = chebyshev.chebval(p, cztd2)

        xtd3 = chebyshev.chebval(p, cxtd3)
        ytd3 = chebyshev.chebval(p, cytd3)
        ztd3 = chebyshev.chebval(p, cztd3)

        # Compute curvature
        v1 = np.array([xtd1, ytd1, ztd1])
        v2 = np.array([xtd2, ytd2, ztd2])

        rs = np.cross(v1, v2)
        r1 = np.dot(rs, rs)
        r2 = np.dot(v1, v1)

        curvature = math.sqrt(r1) / math.sqrt(r2) ** 3

        # Compute torsion
        det = -xtd3 * ytd2 * ztd1
        det += xtd2 * ytd3 * ztd1
        det += xtd3 * ytd1 * ztd2
        det -= xtd1 * ytd3 * ztd2
        det -= xtd2 * ytd1 * ztd3
        det += xtd1 * ytd2 * ztd3

        torsion = det / r1

        return curvature, torsion

    @staticmethod
    def calc_arc_length(p: float, xt: CubicSpline, yt: CubicSpline, zt: CubicSpline) -> float:
        """
        Compute the arc length of a 3-residues long curve
        :param p: point around the curve is calculated
        :type p: float
        :param xt: x(t)
        :type xt: CubicSpline
        :param yt: y(t)
        :type yt: CubicSpline
        :param zt: z(t)
        :type zt: CubicSpline
        :return: arc length
        :rtype: float
        """
        arc_len = 0.0

        i = p - 1.0
        while i < (p + 1.0):
            dx = xt(i + 0.1) - xt(i)
            dy = yt(i + 0.1) - yt(i)
            dz = zt(i + 0.1) - zt(i)

            dist = np.array([dx, dy, dz])

            arc_len += math.sqrt(np.dot(dist, dist))

            i += 0.1

        return arc_len

    @staticmethod
    def calc_writhing(i: int, t: List[float], x: List[float], y: List[float], z: List[float]) -> float:
        """
        Compute the writhing number in a 5-residue long window
        :param i: residue postion
        :type i: int
        :param t: curve's parameters
        :type t: list[float]
        :param x: x(t)
        :type x: list[float]
        :param y: y(t}
        :type y: list[float]
        :param z: z(t)
        :type z: list[float]
        :return: writhing number
        :rtype: float
        """
        start = i - 2
        stop = i + 2

        ini = 0
        end = len(t) - 1
        if start < ini:
            offset = ini - start
        elif stop > end:
            offset = end - stop
        else:
            offset = 0

        start += offset
        stop += offset

        rij = np.zeros(3)
        ri1j = np.zeros(3)
        rij1 = np.zeros(3)
        rjj1 = np.zeros(3)
        rii1 = np.zeros(3)
        ri1j1 = np.zeros(3)

        # Return the number's sign
        def sgn(v: np.ndarray) -> float:
            return v and (1.0, -1.0)[v < 0.0]

        total = 0.0
        for ii in range(start, stop - 2):
            for jj in range(ii + 2, stop):
                rij[0] = x[jj] - x[ii]
                rij[1] = y[jj] - y[ii]
                rij[2] = z[jj] - z[ii]

                ri1j[0] = x[jj] - x[ii + 1]
                ri1j[1] = y[jj] - y[ii + 1]
                ri1j[2] = z[jj] - z[ii + 1]

                rij1[0] = x[jj + 1] - x[ii]
                rij1[1] = y[jj + 1] - y[ii]
                rij1[2] = z[jj + 1] - z[ii]

                ri1j1[0] = x[jj + 1] - x[ii + 1]
                ri1j1[1] = y[jj + 1] - y[ii + 1]
                ri1j1[2] = z[jj + 1] - z[ii + 1]

                rjj1[0] = x[jj + 1] - x[jj]
                rjj1[1] = y[jj + 1] - y[jj]
                rjj1[2] = z[jj + 1] - z[jj]

                rii1[0] = x[ii + 1] - x[ii]
                rii1[1] = y[ii + 1] - y[ii]
                rii1[2] = z[ii + 1] - z[ii]

                aij = (np.cross(rij, rij1) / np.linalg.norm(np.cross(rij, rij1)))
                bij = (np.cross(rij1, ri1j1) / np.linalg.norm(np.cross(rij1, ri1j1)))
                cij = (np.cross(ri1j1, ri1j) / np.linalg.norm(np.cross(ri1j1, ri1j)))
                dij = (np.cross(ri1j, rij) / np.linalg.norm(np.cross(ri1j, rij)))

                omegaij = (math.asin(np.dot(aij, bij)) +
                           math.asin(np.dot(bij, cij)) +
                           math.asin(np.dot(cij, dij)) +
                           math.asin(np.dot(dij, aij))) * sgn(np.dot(np.cross(rjj1, rii1), rij1))

                total += omegaij / (4.0 * math.pi)
        writhing = 2.0 * total
        return writhing

    @staticmethod
    def calc_geometry(chain: Chain, deg: bool) -> Dict[int, ResidueGeometry]:
        """
        Function used to compute the geometric properties around residues.
        It computes curvature, torsion, arc-length and writhing number
        :param chain: Protein main-chain
        :type chain: Chain
        :param deg: angle in degrees?
        :type deg: bool
        :return:  Residue dictionary
        :rtype: Dict[int, ResidueGeometry]
        """
        t = []
        x = []
        y = []
        z = []

        residues_map = {}

        residues = {}
        num = 0
        for residue in chain:
            # Skip invalid residues
            res_type, model, chain_id, res_id = residue.get_full_id()
            het_flag, pos, insertion_code = res_id
            if het_flag[0] != " ":
                continue

            try:
                coord = residue["CA"].get_coord()
            except KeyError:
                print(f'ERROR: missing CA atom at {residue.get_resname()} - {residue.get_full_id()}!')
                raise

            t.append(float(num))
            x.append(coord[0])
            y.append(coord[1])
            z.append(coord[2])

            residues[int(num)] = ResidueGeometry(name=residue.get_resname(),
                                                 chain=chain_id,
                                                 res_num=num,
                                                 res_order=pos,
                                                 curvature=0.0,
                                                 torsion=0.0,
                                                 arc_len=0.0,
                                                 writhing=0.0,
                                                 res_ann={})
            residues_map[pos] = num

            num += 1

        # Fit the alpha-carbons with a cubic-spline
        xt = CubicSpline(t, x, bc_type='natural')
        yt = CubicSpline(t, y, bc_type='natural')
        zt = CubicSpline(t, z, bc_type='natural')

        ini = 0
        end = len(t) - 1
        for i, tp in enumerate(t):
            # Compute curvature and torsion
            if ini < i < end:
                curvature, torsion = GeometryParser.calc_curvature_torsion(p=tp, t=t, xt=xt, yt=yt, zt=zt)
            elif i == ini:
                curvature, torsion = GeometryParser.calc_curvature_torsion(p=t[+1], t=t, xt=xt, yt=yt, zt=zt)
            else:
                curvature, torsion = GeometryParser.calc_curvature_torsion(p=t[-2], t=t, xt=xt, yt=yt, zt=zt)

            # Compute the arc length
            arc_len = GeometryParser.calc_arc_length(p=tp, xt=xt, yt=yt, zt=zt)

            # Compute the writhing number
            writhing = GeometryParser.calc_writhing(i=i, t=t, x=x, y=y, z=z)

            residues[int(tp)].curvature = curvature
            residues[int(tp)].torsion = torsion
            residues[int(tp)].arc_len = arc_len
            residues[int(tp)].writhing = writhing

        GeometryParser.calc_dihedral_angles(chain=chain, residues=residues, deg=deg)

        return residues, residues_map

    @staticmethod
    def calc_dihedral_torsion(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray, deg: bool) -> float:
        """
        Compute the dihedral angles between four vectors
        :param p1: vector 1
        :type p1: np.ndarray
        :param p2: vector 2
        :type p2:  np.ndarray
        :param p3: vector 3
        :type p3: np.ndarray
        :param p4: vector 4
        :type p4: np.ndarray
        :param deg: dihedral angle
        :type deg: bool
        :return: Dihedral angle (in degrees if deg=True, radians otherwise)
        :rtype: float
        """
        b1 = p2 - p1
        b2 = p2 - p3
        b3 = p4 - p3

        # Normalize a vector
        def norm_vec(v: np.ndarray) -> np.ndarray:
            return v/np.linalg.norm(v)

        n1 = norm_vec(np.cross(b1, b2))
        n2 = norm_vec(np.cross(b2, b3))

        m1 = np.cross(n1, norm_vec(b2))

        x = np.dot(n1, n2)
        y = np.dot(m1, n2)

        if deg:
            theta = math.degrees(math.atan2(y, x))
        else:
            theta = math.atan2(y, x)
        return theta

    @staticmethod
    def calc_dihedral_angles(chain: Chain, residues: Dict[int, ResidueGeometry], deg: bool) -> None:
        """
        Compute the dihedral angles phi and psi.
        :param chain: Protein chain
        :type chain: Chain
        :param residues:  Residue dictionary
        :type residues: Dict[int, ResidueGeometry]
        :param deg: angle in degrees?
        :type deg: bool
        """
        idx = [res.id[1] for res in list(chain.get_residues()) if res.id[0] == ' ']

        ini = idx[+0]
        end = idx[-2]

        num = 0
        for i, residue in enumerate(chain):
            # Skip invalid residues
            het_flag, pos, insertion_code = residue.id
            if het_flag[0] != ' ':
                continue

            # core atoms
            atom_n = residue['N'].get_coord()
            atom_ca = residue['CA'].get_coord()
            atom_c = residue['C'].get_coord()

            # phi
            if pos > ini:
                het_flag, prev_res, insertion_code = chain[idx[i - 1]].id

                if abs(prev_res - pos) <= 1:
                    p1 = chain[pos - 1]['C'].get_coord()
                    p2 = atom_n
                    p3 = atom_ca
                    p4 = atom_c
                    phi = GeometryParser.calc_dihedral_torsion(p1=p1, p2=p2, p3=p3, p4=p4, deg=deg)
                else:
                    phi = 0.0
            else:
                phi = 0.0

            # psi
            if pos < end:
                het_flag, next_res, insertion_code = chain[idx[i + 1]].id

                if abs(next_res - pos) <= 1:
                    p1 = atom_n
                    p2 = atom_ca
                    p3 = atom_c
                    p4 = chain[pos + 1]['N'].get_coord()
                    psi = GeometryParser.calc_dihedral_torsion(p1=p1, p2=p2, p3=p3, p4=p4, deg=deg)
                else:
                    psi = 0.0
            else:
                psi = 0.0

            residues[num].phi = phi
            residues[num].psi = psi

            num += 1
