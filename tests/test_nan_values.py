"""Functions that use curvature and torsion accept residues where they are NaN.

Residues in segments too short for a quantity (between chain breaks) get NaN;
the B-factor mapping and the alignment clustering must treat them like
residues with no value, not fail or propagate NaN.
"""

import io
import math
import warnings
from pathlib import Path

import pytest
from Bio.PDB import PDBParser

import melodia_py as mel
from melodia_py.geometryparser import ChainBreakWarning

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@pytest.fixture(scope="module")
def short_segment_structure():
    """2LJ5 model 1 without residues 3-5: residues 1-2 form a 2-residue segment."""
    lines = []
    with open(EXAMPLES / "2lj5.pdb") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and not 3 <= int(line[22:26]) <= 5:
                lines.append(line)
    pdb = "".join(lines) + "END\n"
    return PDBParser(QUIET=True).get_structure("2lj5", io.StringIO(pdb))


def test_bfactor_from_geo_keeps_fill_value_for_nan(short_segment_structure):
    s = short_segment_structure
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ChainBreakWarning)
        geo = mel.geometry_dict_from_structure(s)
        mel.bfactor_from_geo(s, "curvature", geo=geo)

    curvature = [r.curvature for r in geo["0:A"].residues.values()]
    assert math.isnan(curvature[0]) and math.isnan(curvature[1])
    fill = min(c for c in curvature if not math.isnan(c))

    bfactors = {
        residue.id[1]: {atom.get_bfactor() for atom in residue}
        for residue in s[0]["A"]
    }
    assert not any(math.isnan(b) for values in bfactors.values() for b in values)
    assert bfactors[1] == {fill} and bfactors[2] == {fill}
    assert bfactors[10] == {geo["0:A"].residues[6].curvature}  # residue 10 is the 7th left


def test_cluster_alignment_skips_nan_residues(monkeypatch):
    monkeypatch.chdir(EXAMPLES)
    align = mel.parser_pir_file("model.ali")
    record = next(r for r in align if r.description.startswith("structure"))
    column = next(i for i, letter in enumerate(record.seq) if letter != "-")

    # The same residue with no curvature/torsion, as for a too-short segment
    for key in ("curvature", "torsion"):
        values = list(record.letter_annotations[key])
        values[column] = float("nan")
        record.letter_annotations[key] = values
    mel.cluster_alignment(align, threshold=1.1)

    for clustered in align:
        if clustered.description.startswith("structure"):
            assert len(clustered.letter_annotations["cluster"]) == len(clustered.seq)
