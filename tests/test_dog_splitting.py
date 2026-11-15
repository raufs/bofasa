"""
Regression tests for phylogenetic splitting of coarse domain ortholog groups.

The expected results (tests/data/dog_trees.tar.gz) were produced by the bofasa
v1.3.0 implementation (before splitting was optimised and parallelised) on 129
synthetic gene-family trees; see tests/data/make_test_fixtures.py. Cluster ids
and memberships must match exactly.
"""
import json
import logging
import os

import pytest

import make_test_fixtures as fx
from bofasa import core

LOG = logging.getLogger("bofasa_test")
TREES = sorted(f[:-len(".tre")] for f in os.listdir(fx.TREE_DIR) if f.endswith(".tre"))
CONFIG_IDS = [fx.config_tag(*c) for c in fx.CONFIGS]


@pytest.fixture(scope="module")
def expected():
    with open(fx.EXPECTED_FILE) as fh:
        return json.load(fh)


def test_fixture_covers_all_trees(expected):
    assert len(TREES) == 129
    for tag in CONFIG_IDS:
        assert sorted(expected["dog_splits"][tag]) == TREES
    assert sorted(expected["protein_splits"]) == TREES


@pytest.mark.parametrize("skip,rs,fic", fx.CONFIGS, ids=CONFIG_IDS)
def test_split_dogs_serial_matches_v1_3_0(expected, tmp_path, skip, rs, fic):
    tag = fx.config_tag(skip, rs, fic)
    mismatches = []
    for og in TREES:
        out = tmp_path / f"{og}.txt"
        core.split_dogs([og, os.path.join(fx.TREE_DIR, og + ".tre"), str(out),
                         skip, rs, fic, 1, LOG])
        if fx.read_split_file(str(out)) != expected["dog_splits"][tag][og]:
            mismatches.append(og)
    assert mismatches == []


@pytest.mark.parametrize(
    "skip,rs,fic",
    [(False, 1, 0.1), (False, 5, 0.1)],
    ids=["midpoint_only", "rooting_seeds_5"],
)
def test_split_dogs_in_parallel_matches_v1_3_0(expected, tmp_path, skip, rs, fic):
    tag = fx.config_tag(skip, rs, fic)
    core.split_dogs_in_parallel(
        fx.TREE_DIR, str(tmp_path), LOG,
        skip_merge_back_flag=skip, rooting_seeds=rs,
        fixation_index_cutoff=fic, threads=2,
    )
    mismatches = [
        og for og in TREES
        if fx.read_split_file(str(tmp_path / f"{og}.txt")) != expected["dog_splits"][tag][og]
    ]
    assert mismatches == []


def test_protein_level_splitting_matches_v1_3_0(expected):
    """recursive_splitting(are_proteins=True) + disjoint splitting, as used by split_njt."""
    mismatches = []
    for og in TREES:
        t = fx.protein_style_tree(os.path.join(fx.TREE_DIR, og + ".tre"))
        samples = set(l.name.split("|")[0] for l in t.get_leaves())
        sp = core.recursive_splitting(t.copy(), samples, are_proteins=True)
        sp = core.further_split_disjoint_dog_partitions(t, sp)
        if [sorted(x) for x in sp] != expected["protein_splits"][og]:
            mismatches.append(og)
    assert mismatches == []


def test_split_dog_rooting_failure_raises_runtime_error(tmp_path):
    """Worker errors must surface as RuntimeError (a SystemExit would hang the pool)."""
    with pytest.raises(RuntimeError):
        core.split_dog_rooting(
            ["OG_missing", str(tmp_path / "missing.tre"), 0, None, False, 0.1, LOG]
        )


def test_split_dogs_in_parallel_empty_dir(tmp_path):
    core.split_dogs_in_parallel(str(tmp_path), str(tmp_path), LOG, threads=2)
    assert list(tmp_path.iterdir()) == []
