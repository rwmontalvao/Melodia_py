"""Curvature and torsion depend only on the local geometry of the chain.

An ideal alpha helix has the same geometry at every residue, so its interior
residues must all get the same curvature and torsion, and neither the value of
the curve parameter nor residues added before the first one may change them.
"""

import io

import numpy as np
import pytest
from Bio.PDB import PDBParser
from scipy.interpolate import CubicSpline

import melodia_py as mel
from melodia_py.geometryparser import GeometryParser


def helix_ca(first, n):
    """CA trace of an ideal alpha helix (radius 2.3 A, 100 deg and 1.5 A per residue)."""
    i = np.arange(first, first + n)
    angle = np.radians(100.0 * i)
    return np.column_stack((2.3 * np.cos(angle), 2.3 * np.sin(angle), 1.5 * i))


def helix_geometry(first, n):
    """Melodia curvature and torsion of the helix residues first .. first + n - 1."""
    lines = [
        f"ATOM  {k + 1:5d}  CA  ALA A{k + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        for k, (x, y, z) in enumerate(helix_ca(first, n))
    ]
    pdb = "\n".join(lines + ["END"])
    structure = PDBParser(QUIET=True).get_structure("helix", io.StringIO(pdb))
    df = mel.geometry_from_structure(structure)
    return df.sort_values("order")[["curvature", "torsion"]].to_numpy()


def test_constant_along_ideal_helix():
    geometry = helix_geometry(0, 120)
    interior = geometry[10:110]
    assert np.ptp(interior[:, 0]) < 1e-3
    assert np.ptp(interior[:, 1]) < 1e-3


def test_unchanged_by_residues_added_at_start():
    # The 120 residues of the first helix are the last 120 of the second; away from the
    # first helix's start (natural-spline end effect) the values must agree.
    geometry = helix_geometry(0, 120)
    extended = helix_geometry(-30, 150)[30:]
    np.testing.assert_allclose(extended[10:115], geometry[10:115], rtol=0, atol=1e-6)


@pytest.mark.parametrize("shift", [10.0, 100.0, 1000.0])
def test_unchanged_by_shifting_curve_parameter(shift):
    ca = helix_ca(0, 40)
    t = np.arange(40, dtype=float)

    def curvature_torsion(t):
        splines = [CubicSpline(t, ca[:, k], bc_type="natural") for k in range(3)]
        return np.array(
            [GeometryParser.calc_curvature_torsion(p, list(t), *splines) for p in t[1:-1]]
        )

    np.testing.assert_allclose(
        curvature_torsion(t + shift), curvature_torsion(t), rtol=0, atol=1e-8
    )
