"""At the ends of a chain, values whose window does not fit are NaN, not copies.

Curvature, torsion and arc length use residues i-1 to i+1, so the first and
last residue have none; writhing uses i-2 to i+2, so the first two and last two
have none. Every other residue is computed from its own window.
"""

import io

import numpy as np
from Bio.PDB import PDBParser

import melodia_py as mel


def helix_geometry(n):
    """Melodia geometry of an ideal alpha helix CA trace (2.3 A, 100 deg, 1.5 A)."""
    i = np.arange(n)
    angle = np.radians(100.0 * i)
    xyz = np.column_stack((2.3 * np.cos(angle), 2.3 * np.sin(angle), 1.5 * i))
    pdb = "\n".join(
        f"ATOM  {k + 1:5d}  CA  ALA A{k + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        for k, (x, y, z) in enumerate(xyz)
    )
    structure = PDBParser(QUIET=True).get_structure("helix", io.StringIO(pdb + "\nEND\n"))
    return mel.geometry_from_structure(structure).sort_values("order").reset_index(drop=True)


def test_values_at_chain_ends_are_nan():
    g = helix_geometry(20)
    for column in ("curvature", "torsion", "arc_length"):
        assert g[column].isna().tolist() == [True] + [False] * 18 + [True]
    assert g["writhing"].isna().tolist() == [True] * 2 + [False] * 16 + [True] * 2

