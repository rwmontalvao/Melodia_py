"""Chains are split at breaks: no geometry is computed across a missing connection.

A break is two consecutive CA atoms more than 4.2 A apart. Each side of a break
must get exactly the geometry it would get as a chain of its own.
"""

import io
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.PDB import PDBParser

import melodia_py as mel
from melodia_py.geometryparser import ChainBreakWarning, GeometryParser

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "2lj5.pdb"
GEOMETRY = ["curvature", "torsion", "arc_length", "writhing"]


def first_model_lines():
    """ATOM records of the first model of 2LJ5 (76 residues, chain A, no breaks)."""
    lines = []
    with open(EXAMPLE) as f:
        for line in f:
            if line.startswith("ATOM"):
                lines.append(line.rstrip("\n"))
            elif line.startswith("ENDMDL"):
                break
    return lines


def structure(lines, keep=lambda resseq: True):
    pdb = "\n".join(line for line in lines if keep(int(line[22:26]))) + "\nEND\n"
    return PDBParser(QUIET=True).get_structure("2lj5", io.StringIO(pdb))


def geometry(s):
    """Melodia geometry indexed by residue number, ChainBreakWarnings collected."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df = mel.geometry_from_structure(s)
    breaks = [w for w in caught if issubclass(w.category, ChainBreakWarning)]
    return df.set_index("order"), breaks


@pytest.fixture(scope="module")
def lines():
    return first_model_lines()


def test_intact_chain_has_no_breaks(lines):
    _, breaks = geometry(structure(lines))
    assert breaks == []


def test_each_side_of_a_break_is_computed_on_its_own(lines):
    gapped, breaks = geometry(structure(lines, lambda r: not 36 <= r <= 40))
    before, _ = geometry(structure(lines, lambda r: r <= 35))
    after, _ = geometry(structure(lines, lambda r: r >= 41))

    np.testing.assert_allclose(gapped.loc[1:35, GEOMETRY], before[GEOMETRY], rtol=0, atol=1e-12)
    np.testing.assert_allclose(gapped.loc[41:76, GEOMETRY], after[GEOMETRY], rtol=0, atol=1e-12)

    assert len(breaks) == 1
    assert "chain A between residues 35 and 41" in str(breaks[0].message)


def test_no_dihedrals_across_a_break(lines):
    intact, _ = geometry(structure(lines))
    gapped, _ = geometry(structure(lines, lambda r: not 36 <= r <= 40))

    # the DataFrame stores None as NaN
    assert pd.notna(intact.loc[35, "psi"]) and pd.notna(intact.loc[41, "phi"])
    assert pd.isna(gapped.loc[35, "psi"])
    assert pd.isna(gapped.loc[41, "phi"])
    assert gapped.loc[35, "phi"] == pytest.approx(intact.loc[35, "phi"])
    assert gapped.loc[41, "psi"] == pytest.approx(intact.loc[41, "psi"])


def test_numbering_jump_without_missing_residues_is_not_a_break(lines):
    # Renumber residues 41-76 to 141-176: same atoms, still bonded.
    renumbered = [
        line[:22] + f"{int(line[22:26]) + 100:4d}" + line[26:] if int(line[22:26]) > 40 else line
        for line in lines
    ]
    intact, _ = geometry(structure(lines))
    jumped, breaks = geometry(structure(renumbered))
    assert breaks == []
    np.testing.assert_allclose(jumped[GEOMETRY].to_numpy(), intact[GEOMETRY].to_numpy(), rtol=0, atol=0)


def test_short_segments_get_nan(lines):
    # Residues 30-31 (2 residues) and 34-37 (4 residues) between gaps.
    gapped, _ = geometry(structure(lines, lambda r: r <= 27 or 30 <= r <= 31 or 34 <= r <= 37 or r >= 40))

    two = gapped.loc[30:31]
    assert two[["curvature", "torsion", "writhing"]].isna().all().all()
    assert two["arc_length"].notna().all()

    four = gapped.loc[34:37]
    assert four[["curvature", "torsion", "arc_length"]].notna().all().all()
    assert four["writhing"].isna().all()


def test_all_entry_points_agree(lines):
    s = structure(lines, lambda r: not 36 <= r <= 40)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ChainBreakWarning)
        df = mel.geometry_from_structure(s).set_index("order")
        from_dict = mel.geometry_dict_from_structure(s)["0:A"]
        from_class = GeometryParser(s[0]["A"])

    for gp in (from_dict, from_class):
        values = np.array([
            [r.curvature, r.torsion, r.arc_len, r.writhing]
            for _, r in sorted(gp.residues.items())
        ])
        np.testing.assert_allclose(values, df[GEOMETRY].to_numpy(), rtol=0, atol=0)
        phi = [r.phi for _, r in sorted(gp.residues.items())]
        assert phi[35] is None  # residue 41, first after the break
