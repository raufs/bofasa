"""
Unit tests for the tree helpers used in DOG splitting, each checked against a
straightforward brute-force reference (the pre-optimisation formulation).
"""
import copy
import itertools
import math
import os
import random
import statistics
from collections import defaultdict

import pytest
from ete3 import Tree

import make_test_fixtures as fx
from bofasa import core

TREE_FILES = sorted(
    os.path.join(fx.TREE_DIR, f) for f in os.listdir(fx.TREE_DIR) if f.endswith(".tre")
)
# Keep the brute-force references fast: every 3rd tree.
SAMPLE_TREES = TREE_FILES[::3]
# The O(n^2)-per-node pairwise merge-back reference only on trees <= 200 leaves.
SMALL_TREES = [tf for tf in SAMPLE_TREES if len(Tree(tf)) <= 200]


def _midpoint_rooted(tre_file):
    t = Tree(tre_file)
    t.set_outgroup(t.get_midpoint_outgroup())
    return t


def _samples(t):
    return set(l.name.split("|")[0] for l in t.get_leaves())


def _random_partitions(leaves, rng):
    leaves = list(leaves)
    rng.shuffle(leaves)
    parts, i = [], 0
    while i < len(leaves):
        k = rng.choice([1, 1, 2, 3, 5, 8])
        parts.append(set(leaves[i:i + k]))
        i += k
    return parts


# --------------------------------------------------------------------------- #
# Brute-force references
# --------------------------------------------------------------------------- #
def _reference_score_partitions(rooted_t, partitions):
    pscore_sum, bl_sum = 0, 0.0
    for p in partitions:
        ptree = copy.deepcopy(rooted_t)
        ptree.prune(list(p), preserve_branch_length=True)
        pscore_sum += core.determine_tree_score(ptree, set(x.split("|")[0] for x in p))
        bl_sum += core.determine_tree_score_bl(ptree)
    return pscore_sum, bl_sum


def _reference_merge_back(t, split_partitions, fixation_index_cutoff):
    """Pairwise-distance formulation of merge_back_dog_partitions (bofasa <= 1.3.0)."""
    leaf_to_clade, clade_leaves, all_sps = {}, defaultdict(set), set()
    for i, sp in enumerate(split_partitions):
        for l in sp:
            clade_leaves[i].add(l)
            leaf_to_clade[l] = i
            all_sps.add(i)
    if len(clade_leaves) == 1:
        return split_partitions
    leaves = t.get_leaves()
    dist = {}
    for a, b in itertools.combinations(leaves, 2):
        dist[(a, b)] = dist[(b, a)] = t.get_distance(a, b)
    merge_sps, paired = [], set()
    for n in t.traverse("postorder"):
        if n.is_leaf():
            continue
        by_clade = defaultdict(list)
        for l in n.get_leaves():
            by_clade[leaf_to_clade[l.name]].append(l)
        if len(by_clade) == 1:
            continue
        within = [dist[(a, b)] for c in sorted(by_clade)
                  for a, b in itertools.combinations(by_clade[c], 2)]
        between = [dist[(a, b)]
                   for c1, c2 in itertools.combinations(sorted(by_clade), 2)
                   for a in by_clade[c1] for b in by_clade[c2]]
        if not within:
            continue
        pw, pb = statistics.mean(within), statistics.mean(between)
        fix = 0.0 if pb == 0 else (pb - pw) / pb
        if fix > fixation_index_cutoff:
            continue
        for c1, c2 in itertools.combinations(sorted(by_clade), 2):
            merge_sps.append([c1, c2])
            paired.update((c1, c2))
    if not merge_sps:
        return split_partitions
    merged = core.single_linkage_cluster(merge_sps, all_sps, paired)
    return [[l for sp in msp for l in clade_leaves[sp]] for msp in merged]


# --------------------------------------------------------------------------- #
# Tests
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("are_proteins", [False, True])
def test_subtree_duplication_scores_match_determine_tree_score(are_proteins):
    for tf in SAMPLE_TREES:
        t = _midpoint_rooted(tf)
        all_samples = _samples(t)
        scores, named = core._subtree_duplication_scores(t, all_samples, are_proteins)
        for n in t.traverse():
            assert scores[n] == core.determine_tree_score(n, all_samples, are_proteins)
            assert named[n] == len(n.get_leaves())


def test_score_partitions_matches_prune_reference():
    rng = random.Random(5)
    for tf in SAMPLE_TREES:
        t = _midpoint_rooted(tf)
        for _ in range(3):
            parts = _random_partitions(t.get_leaf_names(), rng)
            got_score, got_bl = core.score_partitions(t, parts)
            ref_score, ref_bl = _reference_score_partitions(t, parts)
            assert got_score == ref_score
            assert math.isclose(got_bl, ref_bl, rel_tol=1e-9, abs_tol=1e-12)


def test_score_partitions_single_leaf_uses_root_to_leaf_depth():
    # ete3's prune keeps the original root, so a 1-leaf partition scores its depth.
    t = Tree("((A|a|1|2:1,B|b|1|2:2):3,(C|c|1|2:4,D|d|1|2:5):6);")
    score, bl = core.score_partitions(t, [{"A|a|1|2"}, {"B|b|1|2", "C|c|1|2", "D|d|1|2"}])
    assert score == 0
    assert bl == pytest.approx(4.0 + (2 + 3 + 4 + 6 + 5))


@pytest.mark.parametrize("cutoff", [0.1, 0.5, 0.9])
def test_merge_back_matches_pairwise_reference(cutoff):
    rng = random.Random(11)
    assert len(SMALL_TREES) >= 20
    for tf in SMALL_TREES:
        t = _midpoint_rooted(tf)
        sp = core.recursive_splitting(t, _samples(t))
        sp = core.further_split_disjoint_dog_partitions(t, sp)
        for parts in (sp, _random_partitions(t.get_leaf_names(), rng)):
            got = core.merge_back_dog_partitions(t, parts, fixation_index_cutoff=cutoff)
            ref = _reference_merge_back(t, parts, cutoff)
            assert [sorted(x) for x in got] == [sorted(x) for x in ref]


def test_merge_back_identical_sequences_merge():
    # All-zero branch lengths: partitions are undifferentiated and merged back.
    t = Tree("((A|a|1|2:0,B|b|1|2:0):0,(A|c|1|2:0,B|d|1|2:0):0);")
    parts = [{"A|a|1|2", "B|b|1|2"}, {"A|c|1|2", "B|d|1|2"}]
    got = core.merge_back_dog_partitions(t, parts, fixation_index_cutoff=0.1)
    assert [sorted(x) for x in got] == [sorted(parts[0] | parts[1])]


def test_exact_depths_match_float_depths():
    for tf in SAMPLE_TREES:
        t = _midpoint_rooted(tf)
        depth, scale = core._exact_depths(t)
        for n in t.traverse():
            assert depth[n] / (1 << scale) == pytest.approx(
                t.get_distance(n) if not n.is_root() else 0.0, rel=1e-12, abs=1e-15
            )


def test_further_split_disjoint_monophyletic_and_cover():
    t = Tree("((A|1:1,B|2:1):1,(C|3:1,(D|4:1,E|5:1):1):1);")
    # Monophyletic partition kept; non-monophyletic one covered by maximal subtrees.
    out = core.further_split_disjoint_dog_partitions(
        t, [{"D|4", "E|5"}, {"A|1", "B|2", "C|3"}]
    )
    assert [sorted(x) for x in out] == [["D|4", "E|5"], ["A|1", "B|2"], ["C|3"]]


def test_plan_dog_rootings():
    t = core._load_named_dog_tree(TREE_FILES[-1])
    assert core.plan_dog_rootings("OG", t, 1) == [None]
    assert core.plan_dog_rootings("OG", t, 0) == [None]
    plan = core.plan_dog_rootings("OG", t, 5)
    assert len(plan) == 5 and plan[-1] is None
    names = set(n.name for n in t.traverse())
    assert all(p in names for p in plan[:-1])
    # Deterministic per OG, independent of call order.
    assert core.plan_dog_rootings("OG", core._load_named_dog_tree(TREE_FILES[-1]), 5) == plan


def test_count_newick_leaves():
    for tf in SAMPLE_TREES:
        assert core._count_newick_leaves(tf) == len(Tree(tf))
