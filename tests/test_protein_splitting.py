"""
Regression tests for protein-level splitting (split_njt).

tests/data/protein_njt.tar.gz holds a synthetic protein dataset (DOG counts per
protein, coarse clusters), the FastME tree bofasa v1.3.0 built for every
cluster, and the split v1.3.0 produced. A stand-in ``fastme`` on PATH returns
the stored tree, so the tests run without FastME and independently of its version.
"""
import io
import json
import logging
import os
import random
import shutil
import subprocess
import sys
import zlib

import pytest
from ete3 import Tree
from scipy.spatial import distance

import make_test_fixtures as fx
from bofasa import core
from bofasa.utils import create_pool

LOG = logging.getLogger("bofasa_test")
TREE_DIR = os.path.join(fx.PROTEIN_DIR, "trees")


@pytest.fixture(scope="module")
def dataset():
    with open(os.path.join(fx.PROTEIN_DIR, "dataset.json")) as fh:
        data = json.load(fh)
    with open(os.path.join(fx.PROTEIN_DIR, "expected_splits.json")) as fh:
        data["expected"] = json.load(fh)
    return data


@pytest.fixture
def fake_fastme(tmp_path, monkeypatch):
    """Put a `fastme` on PATH that copies the stored tree for the input matrix."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    helper = bin_dir / "fake_fastme.py"
    helper.write_text(
        "import os, shutil, sys\n"
        "a = sys.argv[1:]\n"
        "src = a[a.index('-i') + 1]\n"
        "og = os.path.basename(src)[:-len('.phylip')]\n"
        f"shutil.copy(os.path.join({TREE_DIR!r}, og + '.tre'), a[a.index('-o') + 1])\n"
    )
    exe = bin_dir / "fastme"
    exe.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{helper}" "$@"\n')
    exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return exe


def _reference_cosine_phylip(members, prot_dogs):
    """The v1.3.0 per-pair scipy formulation of the FastME input matrix."""
    all_dogs = sorted(set(d for p in members for d, c in prot_dogs[p].items() if c > 0))
    vec = {p: [prot_dogs[p].get(d, 0) for d in all_dogs] for p in members}
    lines = [str(len(members))]
    for i, p1 in enumerate(sorted(members)):
        row = ['{:.10f}'.format(distance.cosine(vec[p1], vec[p2])) for p2 in sorted(members)]
        lines.append(str(i) + ' ' * (10 - len(str(i))) + ' '.join(row))
    return '\n'.join(lines) + '\n'


def _reference_outlier_split(all_proteins, sps, prot_dogs):
    """The v1.3.0 per-pair weighted-Jaccard formulation of outlier splitting."""
    def jacc(a, b):
        da = {d: c for d, c in prot_dogs.get(a, {}).items() if c != 0}
        db = {d: c for d, c in prot_dogs.get(b, {}).items() if c != 0}
        union = inter = 0
        for d in set(da) | set(db):
            x, y = da.get(d, 0), db.get(d, 0)
            union += x + y - min(x, y)
            inter += min(x, y)
        return inter / float(union) if union > 0 else None

    out = []
    for sp in sps:
        singletons = set()
        for p1 in sorted(sp):
            internal = external = 0.0
            for p2 in sorted(sp):
                if p1 != p2:
                    j = jacc(p1, p2)
                    if j is not None and j > internal:
                        internal = j
            for p2 in sorted(all_proteins):
                if p2 not in sp:
                    j = jacc(p1, p2)
                    if j is not None and j > external:
                        external = j
            if internal <= external:
                singletons.add(p1)
                out.append({p1})
        rest = set(sp) - singletons
        if rest:
            out.append(rest)
    return out


def _run_all(dataset, work_dir, pass_dogs):
    mismatches = []
    for ci, members in enumerate(dataset["clusters"]):
        dogs = {p: dataset["protein_dogs"][p] for p in members} if pass_dogs else None
        inp = fx.njt_inputs(str(work_dir), ci, members, LOG, prot_dogs=dogs)
        core.split_njt(inp)
        if fx.read_njt_split(inp[4]) != dataset["expected"][f"CoarseOG_{ci}"]:
            mismatches.append(ci)
    return mismatches


def test_fixture_contents(dataset):
    assert len(dataset["clusters"]) == 86
    assert len(dataset["expected"]) == 86
    assert len(os.listdir(TREE_DIR)) == 86


def test_split_njt_global_protein_dogs_matches_v1_3_0(dataset, fake_fastme, tmp_path, monkeypatch):
    monkeypatch.setattr(core, "protein_dogs", dataset["protein_dogs"])
    assert _run_all(dataset, tmp_path, pass_dogs=False) == []


def test_split_njt_passed_protein_dogs_matches_v1_3_0(dataset, fake_fastme, tmp_path, monkeypatch):
    # Worker processes do not see the parent's protein_dogs; jobs carry their own.
    monkeypatch.setattr(core, "protein_dogs", {})
    assert _run_all(dataset, tmp_path, pass_dogs=True) == []


def test_split_njt_jobs_in_process_pool_match_v1_3_0(dataset, fake_fastme, tmp_path):
    jobs = []
    for ci, members in enumerate(dataset["clusters"]):
        dogs = {p: dataset["protein_dogs"][p] for p in members}
        jobs.append(fx.njt_inputs(str(tmp_path), ci, members, LOG, prot_dogs=dogs))
    pool = create_pool(2)
    try:
        results = dict(pool.imap_unordered(core._split_njt_job, jobs))
    finally:
        pool.close()
        pool.join()
    assert all(err is None for err in results.values())
    mismatches = [
        ci for ci in range(len(jobs))
        if fx.read_njt_split(jobs[ci][4]) != dataset["expected"][f"CoarseOG_{ci}"]
    ]
    assert mismatches == []


def test_cosine_distance_matrix_matches_scipy(dataset, tmp_path):
    clusters = dataset["clusters"]
    for ci in list(range(0, len(clusters), 7)) + [max(range(len(clusters)), key=lambda i: len(clusters[i]))]:
        members = clusters[ci]
        out = tmp_path / f"{ci}.phylip"
        naming = core.write_protein_cosine_distance_matrix(
            set(members), dataset["protein_dogs"], str(out), block_size=37)
        assert out.read_text() == _reference_cosine_phylip(members, dataset["protein_dogs"])
        assert naming == {str(i): p for i, p in enumerate(sorted(members))}


def test_outlier_split_matches_reference(dataset):
    rng = random.Random(9)
    for ci in range(0, len(dataset["clusters"]), 5):
        members = sorted(dataset["clusters"][ci])
        t = Tree(os.path.join(TREE_DIR, f"CoarseOG_{ci}.tre"))
        for leaf in t.get_leaves():
            leaf.name = members[int(leaf.name)]
        leaves = list(members)
        rng.shuffle(leaves)
        sps, i = [], 0
        while i < len(leaves):
            k = rng.choice([1, 2, 4, 9, 20])
            sps.append(set(leaves[i:i + k]))
            i += k
        got = core.further_split_outlier_artifact_groups(
            t, sps, prot_dogs=dataset["protein_dogs"], block_size=3)
        ref = _reference_outlier_split(members, sps, dataset["protein_dogs"])
        assert got == ref


def test_split_njt_job_reports_errors(tmp_path):
    og, err = core._split_njt_job(
        ["CoarseOG_x", str(tmp_path / "missing.txt"), str(tmp_path / "x.phylip"),
         str(tmp_path / "x.tre"), str(tmp_path / "x.split"), LOG, 1, {}])
    assert og == "CoarseOG_x"
    assert err is not None


# --------------------------------------------------------------------------- #
# Protein pair listing for MCL (repeated pairs removed)
# --------------------------------------------------------------------------- #
def _dog_member_lists(prot_dogs):
    members = {}
    for p, dogs in prot_dogs.items():
        for d, c in dogs.items():
            if c > 0:
                members.setdefault(d, set()).add(p)
    return [sorted(members[d]) for d in sorted(members)]


def _write_pairs(dog_lists, num_batches):
    handles = {b: io.StringIO() for b in range(num_batches)}
    counts = core.write_unique_protein_pairs(dog_lists, handles)
    return ''.join(handles[b].getvalue() for b in range(num_batches)), counts


def _reference_pair_listing(dog_lists, num_batches):
    """v1.3.0 listing: every pair once per shared DOG, round-robin over batches."""
    batches = [[] for _ in range(num_batches)]
    n = 0
    for members in dog_lists:
        for j, p1 in enumerate(members):
            for p2 in members[j + 1:]:
                batches[n % num_batches].append(p1 + '\t' + p2 + '\n')
                n += 1
    return ''.join(''.join(b) for b in batches), n


def _drop_repeats(listing):
    seen, out = set(), []
    for line in listing.splitlines(keepends=True):
        if line not in seen:
            seen.add(line)
            out.append(line)
    return ''.join(out)


@pytest.mark.parametrize("num_batches", [1, 2, 3, 7])
def test_unique_protein_pairs_are_old_listing_minus_repeats(dataset, num_batches):
    dog_lists = _dog_member_lists(dataset["protein_dogs"])
    got, (total, written) = _write_pairs(dog_lists, num_batches)
    ref, ref_total = _reference_pair_listing(dog_lists, num_batches)
    assert total == ref_total
    assert got == _drop_repeats(ref)
    assert written == got.count('\n') < total


def test_unique_protein_pairs_random_overlaps():
    rng = random.Random(4)
    prots = [f"S{i % 5}|p{i}" for i in range(40)]
    for _ in range(20):
        dog_lists = [sorted(rng.sample(prots, rng.randint(2, 12))) for _ in range(rng.randint(1, 15))]
        for nb in (1, 4, 5):
            got, _ = _write_pairs(dog_lists, nb)
            assert got == _drop_repeats(_reference_pair_listing(dog_lists, nb)[0])


@pytest.mark.skipif(shutil.which("mcl") is None, reason="mcl not installed")
def test_mcl_graph_and_clusters_unchanged_without_repeats(dataset, tmp_path):
    """MCL keeps one edge for repeated (identically weighted) pairs."""
    dog_lists = _dog_member_lists(dataset["protein_dogs"])
    ref, _ = _reference_pair_listing(dog_lists, 3)
    new, _ = _write_pairs(dog_lists, 3)
    weight = lambda line: str(10 + (zlib.crc32(line.encode()) % 90))  # same per pair
    results = {}
    for name, listing in (("old", ref), ("new", new)):
        abc = tmp_path / f"{name}.abc"
        abc.write_text(''.join(l.rstrip('\n') + '\t' + weight(l) + '\n'
                               for l in listing.splitlines(keepends=True)))
        cls, graph = tmp_path / f"{name}.cls", tmp_path / f"{name}.graph"
        subprocess.run(["mcl", str(abc), "--abc", "-I", "1.2", "-o", str(cls),
                        "-write-graph", str(graph)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        results[name] = (cls.read_text(), graph.read_text())
    assert results["old"] == results["new"]
