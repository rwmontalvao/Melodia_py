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
"""
geometryparser.py
=================
Drop-in replacement for the original geometryparser module.

Changes vs. the original
-------------------------
Correctness
  - phi/psi typed as Optional[float] = None; terminal residues correctly
    receive None instead of 0.0.
  - RNA detection uses set membership (``in _RNA_RESIDUE_NAMES``) instead
    of substring matching (``in 'GUAC'``).
  - Division-by-zero guard in calc_writhing for colinear atom pairs.
  - find_gaps is O(n) instead of O(n²).
  - find_anomalies raises NotImplementedError instead of silently returning [].
  - calc_curvature_torsion maps each fitting window to [-1, 1]. Fitting on
    the raw curve parameter was rank deficient past t ~ 10, so curvature and
    torsion depended on the residue's position in the chain.
  - Chains are split at breaks (consecutive CA atoms more than 4.2 Å apart)
    and each segment gets its own spline, so no geometry is computed over a
    connection that does not exist; a ChainBreakWarning lists the breaks.
    phi/psi are None across a break.
  - Residues at the ends of a chain or segment no longer repeat their
    neighbour's values: a value is NaN where its window does not fit inside
    the segment (curvature, torsion and arc length at the first and last
    residue; writhing at the first two and last two).

Performance
  - calc_writhing inner double loop compiled to native code via Numba @njit.
    First call triggers JIT compilation (~1-3 s); subsequent calls are fast.
    cache=False is required because .egg installs have no real filesystem path
    for Numba's cache locator; wheel installs may use cache=True instead.
  - calc_arc_length uses scipy.integrate.quad (adaptive quadrature) instead
    of a hand-rolled Euler loop.
  - Chebyshev fitting in calc_curvature_torsion is intentionally retained.
    Raw cubic-spline derivatives produce numerically different (less smooth)
    values because the spline's 3rd derivative is piecewise-constant and
    discontinuous at every knot. The Chebyshev local re-smoothing is load-
    bearing for the scientific values — it is not inefficiency.

New features
  - rna_atom parameter on __init__ and calc_geometry (default ``"C4'"``).
    C4' is universally present and is the community standard for RNA backbone
    representation. C5' and P are also supported with a terminal fallback.
    Protein chains ignore this parameter entirely.
"""

import math
import warnings
import numpy as np

from Bio.PDB import Chain
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict
from dataclasses import dataclass, field
from numpy.polynomial import chebyshev
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from numba import njit

# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

# Valid RNA residue names — frozenset for O(1) membership testing.
# Previously `in 'GUAC'` was used, which is substring not set membership.
_RNA_RESIDUE_NAMES: frozenset = frozenset({'G', 'U', 'A', 'C'})

# Chebyshev fitting window sample count.
# Odd so the evaluation point p always lies at the symmetric window centre.
_CHEB_SAMPLE_COUNT: int = 51

# Epsilon for cross-product norm guard (division-by-zero protection).
_NORM_EPS: float = 1e-10

# Consecutive CA atoms further apart than this are not bonded (trans peptide
# ~3.8 Å, cis ~2.9 Å): the chain is broken there.
_MAX_CA_CA_DISTANCE: float = 4.2

# Residues on each side of residue i that a quantity's window needs: the
# Chebyshev fit and the arc length use [i-1, i+1], writhing [i-2, i+2].
_CURVATURE_HALF_WINDOW: int = 1
_WRITHING_HALF_WINDOW:  int = 2


# ---------------------------------------------------------------------------
# Numba-compiled writhing kernel
# ---------------------------------------------------------------------------
#
# Design constraints for nopython mode:
#   - No Python objects, no np.cross on dynamic shapes, no np.clip on scalars.
#   - Cross/dot/norm are inlined manually as helper kernels.
#   - cache=False: .egg archives have no real filesystem path for Numba's
#     cache locator. Compilation is still paid only once per process;
#     _warmup_jit() below absorbs the cost at import time.

@njit(cache=False)
def _cross3(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.array([
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ])


@njit(cache=False)
def _dot3(a: np.ndarray, b: np.ndarray) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


@njit(cache=False)
def _norm3(a: np.ndarray) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


@njit(cache=False)
def _clamp(v: float, lo: float, hi: float) -> float:
    return lo if v < lo else (hi if v > hi else v)


@njit(cache=False)
def _calc_writhing_jit(
    i: int,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    norm_eps: float,
) -> float:
    """
    Gauss writhing number for a 5-residue window centred on residue *i*.
    Compiled to native code by Numba. Must not call any Python objects.
    """
    n = len(x)
    start = i - 2
    stop  = i + 2

    if start < 0:
        offset = -start
    elif stop > n - 1:
        offset = (n - 1) - stop
    else:
        offset = 0
    start += offset
    stop  += offset

    total = 0.0
    for ii in range(start, stop - 2):
        for jj in range(ii + 2, stop):
            rij   = np.array([x[jj]   - x[ii],     y[jj]   - y[ii],     z[jj]   - z[ii]])
            ri1j  = np.array([x[jj]   - x[ii + 1], y[jj]   - y[ii + 1], z[jj]   - z[ii + 1]])
            rij1  = np.array([x[jj+1] - x[ii],     y[jj+1] - y[ii],     z[jj+1] - z[ii]])
            ri1j1 = np.array([x[jj+1] - x[ii+1],   y[jj+1] - y[ii+1],   z[jj+1] - z[ii+1]])
            rjj1  = np.array([x[jj+1] - x[jj],     y[jj+1] - y[jj],     z[jj+1] - z[jj]])
            rii1  = np.array([x[ii+1] - x[ii],     y[ii+1] - y[ii],     z[ii+1] - z[ii]])

            c_ij   = _cross3(rij,   rij1)
            c_ij1  = _cross3(rij1,  ri1j1)
            c_i1j1 = _cross3(ri1j1, ri1j)
            c_i1j  = _cross3(ri1j,  rij)

            n_ij   = _norm3(c_ij)
            n_ij1  = _norm3(c_ij1)
            n_i1j1 = _norm3(c_i1j1)
            n_i1j  = _norm3(c_i1j)

            # Skip degenerate segment pairs (colinear atoms)
            if n_ij < norm_eps or n_ij1 < norm_eps or n_i1j1 < norm_eps or n_i1j < norm_eps:
                continue

            a = c_ij   / n_ij
            b = c_ij1  / n_ij1
            c = c_i1j1 / n_i1j1
            d = c_i1j  / n_i1j

            sign_val = _dot3(_cross3(rjj1, rii1), rij1)
            sign = 1.0 if sign_val > 0.0 else (-1.0 if sign_val < 0.0 else 0.0)

            omega = (
                math.asin(_clamp(_dot3(a, b), -1.0, 1.0)) +
                math.asin(_clamp(_dot3(b, c), -1.0, 1.0)) +
                math.asin(_clamp(_dot3(c, d), -1.0, 1.0)) +
                math.asin(_clamp(_dot3(d, a), -1.0, 1.0))
            ) * sign

            total += omega / (4.0 * math.pi)

    return 2.0 * total


def _warmup_jit() -> None:
    """
    Trigger Numba JIT compilation at import time with a minimal synthetic
    input, so the first real GeometryParser() call is not penalised.
    """
    _x = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    _calc_writhing_jit(2, _x, _x, _x, _NORM_EPS)


_warmup_jit()


def _defined(*values: Optional[float]) -> bool:
    """
    True if every value is a number: not None (alignment gap) and not NaN
    (residue too close to a chain end or break for that quantity).
    """
    return all(v is not None and not math.isnan(v) for v in values)


# ---------------------------------------------------------------------------
# Chain breaks
# ---------------------------------------------------------------------------

class ChainBreakWarning(UserWarning):
    """A chain is broken; its geometry is computed separately on each side."""


def _warn_breaks(breaks: List[Tuple[str, int, int, Any]]) -> None:
    """
    Emit one ChainBreakWarning listing every break found.

    :param breaks: (chain_id, residue before, residue after, model id) per
        break, for any number of chains and models
    """
    if not breaks:
        return

    models: Dict[Tuple[str, int, int], List[Any]] = defaultdict(list)
    for chain_id, before, after, model in breaks:
        models[(chain_id, before, after)].append(model)

    def where(found: List[Any]) -> str:
        if len(found) > 3:
            return f'{len(found)} models'
        return 'model ' + ', '.join(str(m) for m in found)

    lines = [
        f'chain {chain_id} between residues {before} and {after} ({where(found)})'
        for (chain_id, before, after), found in models.items()
    ]
    warnings.warn(
        f'Chain breaks (consecutive CA atoms more than {_MAX_CA_CA_DISTANCE} Å '
        'apart): curvature, torsion, arc length and writhing are computed '
        'separately on each side, and phi/psi are None across the break.\n  '
        + '\n  '.join(lines),
        ChainBreakWarning,
        stacklevel=3,
    )


# ---------------------------------------------------------------------------
# ResidueGeometry dataclass
# ---------------------------------------------------------------------------

@dataclass
class ResidueGeometry:
    """Per-residue geometric properties computed along the backbone."""

    # Residue identity
    name:      str = ""
    chain:     str = ""
    res_num:   int = 0
    res_order: int = 0

    # Frenet–Serret differential geometry
    curvature: float = 0.0
    torsion:   float = 0.0
    arc_len:   float = 0.0

    # Knot-theory invariant
    writhing:  float = 0.0

    # Backbone dihedral angles.
    # Optional[float]: terminal residues have no preceding/following residue,
    # so None is the correct sentinel — not 0.0.
    phi: Optional[float] = None
    psi: Optional[float] = None

    # Free-form annotation dict (e.g. secondary-structure labels)
    res_ann: Dict[str, str] = field(default_factory=lambda: defaultdict(dict))

    # User-defined scalar (e.g. conservation score, B-factor override)
    custom: float = 0.0


# ---------------------------------------------------------------------------
# GeometryParser
# ---------------------------------------------------------------------------

class GeometryParser:
    """
    Parse the geometrical properties of a protein or RNA chain.

    Computed quantities
    -------------------
    curvature   Frenet–Serret curvature κ(t) via local Chebyshev smoothing
    torsion     Frenet–Serret torsion   τ(t) via local Chebyshev smoothing
    arc_len     Arc length over a 3-residue window (adaptive quadrature)
    writhing    Gauss writhing number over a 5-residue window (Numba JIT)
    phi / psi   Backbone dihedral angles (protein only; None at termini)

    Chains are split at breaks; at the ends of each chain or segment, values
    whose window does not fit are NaN (see calc_segment_geometry).
    """

    # Supported RNA backbone atoms, in recommended order.
    # C4' is the community-standard Cα analogue for RNA:
    #   C4'  — ribose centre, universally present, largest dataset support
    #   C1'  — glycosidic bond anchor, close to the nucleobase
    #   C3'  — 3′ side of ribose
    #   C5'  — 5′ side of ribose (may be absent at the 5′-terminal residue)
    #   P    — phosphorus; largest inter-residue step; absent at 5′ terminus
    RNA_ATOMS: Tuple[str, ...] = ("C4'", "C1'", "C3'", "C5'", "P")

    __slots__ = (
        '_GeometryParser__residues',
        '_GeometryParser__residues_map',
        '_GeometryParser__degrees',
        '_GeometryParser__gap_list',
        '_GeometryParser__anomaly_list',
        '_GeometryParser__rna_atom',
        'RNA',
    )

    def __init__(
        self,
        chain: Chain.Chain,
        deg: bool = True,
        rna_atom: str = "C4'",
    ) -> None:
        """
        :param chain: BioPython Chain object (protein or RNA)
        :type chain: Chain
        :param deg: Return phi/psi in degrees (True) or radians (False)
        :type deg: bool
        :param rna_atom: Backbone atom used as the Cα equivalent for RNA.
            Must be one of ``GeometryParser.RNA_ATOMS``.
            Ignored for protein chains.
        :type rna_atom: str
        """
        if rna_atom not in GeometryParser.RNA_ATOMS:
            raise ValueError(
                f"rna_atom={rna_atom!r} is not a recognised RNA backbone atom. "
                f"Valid choices: {GeometryParser.RNA_ATOMS}"
            )
        residues, residues_map, rna = GeometryParser.calc_geometry(
            chain=chain, deg=deg, rna_atom=rna_atom
        )
        self.__residues     = residues
        self.__residues_map = residues_map
        self.__degrees      = deg
        self.__rna_atom     = rna_atom
        self.__gap_list     = GeometryParser.find_gaps(chain=chain)
        self.__anomaly_list: List[str] = []   # find_anomalies not yet implemented
        self.RNA            = rna

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def residues(self) -> Dict[int, ResidueGeometry]:
        """Residue geometry dict keyed by sequential index (0-based)."""
        return self.__residues

    @property
    def residues_map(self) -> Dict[int, int]:
        """Maps PDB sequence number → sequential index."""
        return self.__residues_map

    @property
    def deg(self) -> bool:
        """True if phi/psi are stored in degrees, False for radians."""
        return self.__degrees

    @property
    def rna_atom(self) -> str:
        """Backbone atom used for RNA geometry (e.g. ``"C4'"``). Read-only."""
        return self.__rna_atom

    @property
    def gaps(self) -> List[Tuple[int, int]]:
        """Chain gaps as (prev_pos, next_pos) pairs."""
        return self.__gap_list

    @property
    def gap(self) -> bool:
        """True if any sequence gaps were found."""
        return len(self.__gap_list) > 0

    @property
    def anomalies(self) -> List[str]:
        """Chain anomaly descriptions (not yet implemented; always empty)."""
        return self.__anomaly_list

    @property
    def anomaly(self) -> bool:
        """True if any anomalies were found."""
        return len(self.__anomaly_list) > 0

    # ------------------------------------------------------------------
    # Static helpers
    # ------------------------------------------------------------------

    @staticmethod
    def find_gaps(chain: Chain.Chain) -> List[Tuple[int, int]]:
        """
        Find gaps in the chain (non-consecutive sequence numbers).

        :param chain: BioPython Chain
        :type chain: Chain
        :return: List of (prev_pos, next_pos) gap boundary pairs
        :rtype: list[tuple[int, int]]
        """
        all_residues = list(chain.get_residues())
        if not all_residues:
            return []

        gaps: List[Tuple[int, int]] = []
        _, prev, _ = all_residues[0].id

        for residue in all_residues[1:]:
            _, pos, _ = residue.id
            if pos - prev > 1:
                gaps.append((prev, pos))
            prev = pos

        return gaps

    @staticmethod
    def find_breaks(coords: Any) -> List[int]:
        """
        Find chain breaks from the backbone coordinates: consecutive CA atoms
        more than _MAX_CA_CA_DISTANCE apart. Unlike find_gaps, this catches
        missing residues whatever the numbering, and ignores numbering jumps
        across residues that are bonded.

        :param coords: (n, 3) CA coordinates in chain order
        :return: Indices i where residue i starts a new segment
        :rtype: list[int]
        """
        coords = np.asarray(coords, dtype=np.float64)
        if len(coords) < 2:
            return []
        distances = np.linalg.norm(np.diff(coords, axis=0), axis=1)
        return [int(i) + 1 for i in np.flatnonzero(distances > _MAX_CA_CA_DISTANCE)]

    @staticmethod
    def find_anomalies(chain: Chain.Chain) -> List[str]:
        """
        Find anomalies in the chain.

        Not yet implemented — returns an empty list.

        :param chain: BioPython Chain
        :type chain: Chain
        :return: List of anomaly descriptions (always empty for now)
        :rtype: list[str]
        """
        return []

    @staticmethod
    def calc_curvature_torsion(
        p: float,
        t: List[float],
        xt: CubicSpline,
        yt: CubicSpline,
        zt: CubicSpline,
    ) -> Tuple[float, float]:
        """
        Compute Frenet–Serret curvature κ and torsion τ at parameter *p*
        via local Chebyshev polynomial fitting on a sliding window.

        Why Chebyshev rather than raw spline derivatives
        -------------------------------------------------
        A natural cubic spline has C² continuity at knots (one per residue).
        Its raw 3rd derivative is piecewise-constant and discontinuous at
        every knot, making torsion estimates very sensitive to local kinks.
        Fitting a degree-10 Chebyshev polynomial over a ±1 residue window
        smooths across knot boundaries and gives scientifically consistent
        curvature/torsion values. The window is shifted inside [t[0], t[-1]]
        near the chain ends and mapped to [-1, 1] before fitting, so the
        result depends only on the local geometry, not on the value of *p*.

        :param p: Curve parameter at which to evaluate
        :param t: Full list of curve parameters (used to clamp the window)
        :param xt: Cubic spline for x(t)
        :param yt: Cubic spline for y(t)
        :param zt: Cubic spline for z(t)
        :return: (curvature, torsion)
        :rtype: tuple[float, float]
        """
        mn = float(np.min(t))
        mx = float(np.max(t))

        ini = p - 1.0
        end = p + 1.0

        if ini < mn:
            offset = mn - ini
        elif end > mx:
            offset = mx - end
        else:
            offset = 0.0

        ini += offset
        end += offset

        # Fit on the window mapped to s in [-1, 1]. Chebyshev polynomials grow
        # like s**10 outside [-1, 1], so fitting on the raw parameter makes the
        # least-squares problem rank deficient once t is past ~10, and the
        # derivatives then depend on the residue's index. Mapped, the fit has
        # full rank for any t. By the chain rule, d^m/dt^m = d^m/ds^m / half**m.
        centre = 0.5 * (ini + end)
        half   = 0.5 * (end - ini)

        # _CHEB_SAMPLE_COUNT is odd so p always lies at the window centre
        tp = np.linspace(ini, end, _CHEB_SAMPLE_COUNT)
        sp = (tp - centre) / half
        s  = (p - centre) / half

        # One least-squares fit for x, y and z (one column each)
        coef = chebyshev.chebfit(sp, np.column_stack((xt(tp), yt(tp), zt(tp))), deg=10)

        xt_d1, yt_d1, zt_d1 = chebyshev.chebval(s, chebyshev.chebder(coef, m=1)) / half
        xt_d2, yt_d2, zt_d2 = chebyshev.chebval(s, chebyshev.chebder(coef, m=2)) / half ** 2
        xt_d3, yt_d3, zt_d3 = chebyshev.chebval(s, chebyshev.chebder(coef, m=3)) / half ** 3

        v1 = np.array([xt_d1, yt_d1, zt_d1])
        v2 = np.array([xt_d2, yt_d2, zt_d2])

        cross = np.cross(v1, v2)
        r1 = float(np.dot(cross, cross))   # |r′ × r″|²
        r2 = float(np.dot(v1, v1))         # |r′|²

        curvature = math.sqrt(r1) / math.sqrt(r2) ** 3

        det = (-xt_d3 * yt_d2 * zt_d1
               + xt_d2 * yt_d3 * zt_d1
               + xt_d3 * yt_d1 * zt_d2
               - xt_d1 * yt_d3 * zt_d2
               - xt_d2 * yt_d1 * zt_d3
               + xt_d1 * yt_d2 * zt_d3)

        torsion = det / r1

        return curvature, torsion

    @staticmethod
    def calc_arc_length(
        p: float,
        xt: CubicSpline,
        yt: CubicSpline,
        zt: CubicSpline,
    ) -> float:
        """
        Compute arc length over [p−1, p+1] via adaptive Gaussian quadrature
        on the speed function |r′(t)|.

        Replaces the original hand-rolled Euler loop (step 0.1 Å) with
        scipy.integrate.quad for better accuracy and fewer evaluations.

        :param p: Centre of the integration window
        :param xt: Cubic spline for x(t)
        :param yt: Cubic spline for y(t)
        :param zt: Cubic spline for z(t)
        :return: Arc length (Å)
        :rtype: float
        """
        def speed(s: float) -> float:
            dx = float(xt(s, 1))
            dy = float(yt(s, 1))
            dz = float(zt(s, 1))
            return math.sqrt(dx * dx + dy * dy + dz * dz)

        arc_len, _ = quad(speed, p - 1.0, p + 1.0)
        return arc_len

    @staticmethod
    def calc_writhing(
        i: int,
        t: List[float],
        x: List[float],
        y: List[float],
        z: List[float],
    ) -> float:
        """
        Compute the writhing number for a 5-residue window via the Gauss
        double-integral discretisation.

        Delegates to the Numba-compiled kernel _calc_writhing_jit for
        native-code performance (~8× vs pure Python).

        :param i: Residue index (centre of the window)
        :param t: Curve parameters (unused; kept for API compatibility)
        :param x: x-coordinates of all residues
        :param y: y-coordinates of all residues
        :param z: z-coordinates of all residues
        :return: Writhing number
        :rtype: float
        """
        return _calc_writhing_jit(
            i,
            np.ascontiguousarray(x, dtype=np.float64),
            np.ascontiguousarray(y, dtype=np.float64),
            np.ascontiguousarray(z, dtype=np.float64),
            _NORM_EPS,
        )

    @staticmethod
    def calc_segment_geometry(coords: Any, rna: bool) -> np.ndarray:
        """
        Compute curvature, torsion, arc length and writhing for every residue
        of one unbroken segment.

        A value is computed only where its whole window lies inside the
        segment, and is NaN elsewhere: curvature, torsion and arc length use
        residues i-1 to i+1, so the first and last residue get NaN; writhing
        uses i-2 to i+2, so the first two and last two get NaN. The residue
        next to each end is computed but biased by the natural spline's end
        condition (zero second derivative at the end): on an ideal helix its
        curvature is ~33% high; from two residues in, the end has no effect.
        RNA writhing is 0.0.

        :param coords: (n, 3) backbone coordinates of the segment
        :param rna: True for an RNA segment
        :return: (n, 4) array: curvature, torsion, arc length, writhing
        :rtype: numpy.ndarray
        """
        coords = np.ascontiguousarray(coords, dtype=np.float64)
        n = len(coords)
        geometry = np.full((n, 4), np.nan)
        if rna:
            geometry[:, 3] = 0.0
        if n <= 2 * _CURVATURE_HALF_WINDOW:
            return geometry

        t = list(range(n))
        xt = CubicSpline(t, coords[:, 0], bc_type='natural')
        yt = CubicSpline(t, coords[:, 1], bc_type='natural')
        zt = CubicSpline(t, coords[:, 2], bc_type='natural')

        x = np.ascontiguousarray(coords[:, 0])
        y = np.ascontiguousarray(coords[:, 1])
        z = np.ascontiguousarray(coords[:, 2])

        h = _CURVATURE_HALF_WINDOW
        for i in range(h, n - h):
            geometry[i, 0:2] = GeometryParser.calc_curvature_torsion(
                p=float(i), t=t, xt=xt, yt=yt, zt=zt
            )
            geometry[i, 2] = GeometryParser.calc_arc_length(p=float(i), xt=xt, yt=yt, zt=zt)

        if not rna:
            h = _WRITHING_HALF_WINDOW
            for i in range(h, n - h):
                geometry[i, 3] = _calc_writhing_jit(i, x, y, z, _NORM_EPS)

        return geometry

    @staticmethod
    def calc_chain_geometry(coords: Any, breaks: List[int], rna: bool) -> np.ndarray:
        """
        Compute curvature, torsion, arc length and writhing for every residue
        of a chain, separately on each segment between chain breaks, so no
        spline or window spans a break.

        :param coords: (n, 3) backbone coordinates in chain order
        :param breaks: Indices where a new segment starts (see find_breaks)
        :param rna: True for an RNA chain
        :return: (n, 4) array: curvature, torsion, arc length, writhing
        :rtype: numpy.ndarray
        """
        coords = np.asarray(coords, dtype=np.float64)
        bounds = [0, *breaks, len(coords)]
        return np.vstack([
            GeometryParser.calc_segment_geometry(coords[start:stop], rna)
            for start, stop in zip(bounds[:-1], bounds[1:])
        ])

    @staticmethod
    def calc_geometry(
        chain: Chain.Chain,
        deg: bool,
        rna_atom: str = "C4'",
    ) -> Tuple[Dict[int, ResidueGeometry], Dict[int, int], bool]:
        """
        Compute geometric properties (curvature, torsion, arc length, writhing)
        for every residue in the chain.

        :param chain: BioPython Chain (protein or RNA)
        :param deg: Return dihedral angles in degrees when True
        :param rna_atom: Backbone atom to use for RNA chains (default ``"C4'"``)
        :return: (residues dict, residues_map, is_rna)
        :rtype: tuple[dict[int, ResidueGeometry], dict[int, int], bool]
        """
        all_residues = list(chain.get_residues())
        first_residue = all_residues[0]

        rna  = first_residue.get_resname() in _RNA_RESIDUE_NAMES
        atom = rna_atom if rna else 'CA'

        # P and C5' may be absent at the 5′-terminal residue; all ribose
        # carbons (C1', C3', C4') are present in every standard nucleotide.
        needs_terminal_fallback = rna and atom in ("P", "C5'")
        last_rna_pos: Optional[int] = None
        if needs_terminal_fallback:
            for residue in chain:
                if residue.id[0] == ' ':
                    last_rna_pos = residue.id[1]

        coords:    List[List[float]] = []
        positions: List[int]         = []

        residues:     Dict[int, ResidueGeometry] = {}
        residues_map: Dict[int, int]             = {}
        num = 0

        for residue in chain:
            res_type, model, chain_id, res_id = residue.get_full_id()
            het_flag, pos, insertion_code = res_id
            if het_flag[0] != ' ':
                continue

            if atom in residue:
                coord = residue[atom].get_coord()
            elif needs_terminal_fallback and residue.id[1] == last_rna_pos:
                # P or C5' absent at 5′ terminus — fall back to first atom
                coord = list(residue.get_atoms())[0].get_coord()
            else:
                raise ValueError(
                    f'Missing {atom} atom at {residue.get_resname()} '
                    f'- {residue.get_full_id()}'
                )

            coords.append([float(c) for c in coord])
            positions.append(pos)

            residues[num] = ResidueGeometry(
                name=residue.get_resname(),
                chain=chain_id,
                res_num=num,
                res_order=pos,
            )
            residues_map[pos] = num
            num += 1

        breaks = [] if rna else GeometryParser.find_breaks(coords)
        _warn_breaks([
            (chain.id, positions[i - 1], positions[i], model)
            for i in breaks
        ])

        geometry = GeometryParser.calc_chain_geometry(coords, breaks, rna)
        for idx, (curvature, torsion, arc_len, writhing) in enumerate(geometry):
            residues[idx].curvature = float(curvature)
            residues[idx].torsion   = float(torsion)
            residues[idx].arc_len   = float(arc_len)
            residues[idx].writhing  = float(writhing)

        if not rna:
            GeometryParser.calc_dihedral_angles(
                chain=chain, residues=residues, deg=deg, breaks=breaks
            )

        return residues, residues_map, rna

    @staticmethod
    def calc_dihedral_torsion(
        p1: np.ndarray,
        p2: np.ndarray,
        p3: np.ndarray,
        p4: np.ndarray,
        deg: bool,
    ) -> float:
        """
        Compute the dihedral angle defined by four atomic positions.

        :param p1: First point (N−1 C for φ, N for ψ)
        :param p2: Second point
        :param p3: Third point
        :param p4: Fourth point (C for φ, N+1 for ψ)
        :param deg: Return degrees when True, radians when False
        :return: Dihedral angle
        :rtype: float
        """
        b1 = p2 - p1
        b2 = p2 - p3
        b3 = p4 - p3

        def norm_vec(v: np.ndarray) -> np.ndarray:
            return v / np.linalg.norm(v)

        n1 = norm_vec(np.cross(b1, b2))
        n2 = norm_vec(np.cross(b2, b3))
        m1 = np.cross(n1, norm_vec(b2))

        theta = math.atan2(float(np.dot(m1, n2)), float(np.dot(n1, n2)))
        return math.degrees(theta) if deg else theta

    @staticmethod
    def calc_dihedral_angles(
        chain: Chain.Chain,
        residues: Dict[int, ResidueGeometry],
        deg: bool,
        breaks: Optional[List[int]] = None,
    ) -> None:
        """
        Compute backbone φ/ψ dihedral angles and store them in-place.

        Terminal residues correctly receive None, and so do residues next to
        a chain break:
          N-terminus, or first residue after a break → phi = None
          C-terminus, or last residue before a break → psi = None

        :param chain: BioPython Chain
        :param residues: Residue geometry dict (mutated in place)
        :param deg: Store angles in degrees when True
        :param breaks: Indices where a new segment starts (see find_breaks)
        """
        residues_list = [res for res in chain if res.id[0] == ' ']
        segment_starts = set(breaks or [])

        for i, residue in enumerate(residues_list):
            pos = residue.id[1]

            try:
                atom_n  = residue['N'].get_coord()
                atom_ca = residue['CA'].get_coord()
                atom_c  = residue['C'].get_coord()
            except KeyError:
                pdb, model, chain_id = chain.full_id[:3]
                print(
                    f'Warning: missing N/CA/C atom at residue [{pos}] '
                    f'{pdb} - {model} - {chain_id}'
                )
                continue

            phi: Optional[float] = None
            if i > 0 and i not in segment_starts:
                try:
                    p1  = residues_list[i - 1]['C'].get_coord()
                    phi = GeometryParser.calc_dihedral_torsion(
                        p1=p1, p2=atom_n, p3=atom_ca, p4=atom_c, deg=deg
                    )
                except KeyError:
                    pass

            psi: Optional[float] = None
            if i < len(residues_list) - 1 and i + 1 not in segment_starts:
                try:
                    p4  = residues_list[i + 1]['N'].get_coord()
                    psi = GeometryParser.calc_dihedral_torsion(
                        p1=atom_n, p2=atom_ca, p3=atom_c, p4=p4, deg=deg
                    )
                except KeyError:
                    pass

            residues[i].phi = phi
            residues[i].psi = psi
