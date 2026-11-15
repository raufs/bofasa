#!/usr/bin/env python3
"""
Test fixtures for bofasa's phylogenetic splitting, stored as two archives:

  dog_trees.tar.gz    - 129 synthetic coarse domain ortholog group (DOG) trees
                        (dog_trees/OG_*.tre) and their expected splits
                        (expected_dog_splits.json).
  protein_njt.tar.gz  - a synthetic protein dataset (protein_njt/dataset.json:
                        protein -> {DOG: count} and coarse protein clusters), the
                        FastME trees built for each cluster (protein_njt/trees/),
                        and the expected splits (protein_njt/expected_splits.json).

The committed expected results were produced with the bofasa v1.3.0
implementation (commit b6d4e71), before the splitting code was optimised, so the
tests check that the current code reproduces them exactly. Storing the FastME
trees lets the protein tests run without FastME installed.

The archives are extracted automatically into tests/data/.extracted/ when pytest
starts (see tests/conftest.py).

Regenerate (only when splitting behaviour changes on purpose):
    python tests/data/make_test_fixtures.py [--package PKG]
PKG is the importable bofasa package used to compute expected results
(default: bofasa). Protein fixtures require ``fastme`` on PATH.
"""
import argparse
import importlib
import json
import logging
import os
import random
import shutil
import sys
import tarfile
import tempfile
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('ignore', SyntaxWarning)
    from ete3 import Tree

HERE = os.path.dirname(os.path.abspath(__file__))
DOG_ARCHIVE = os.path.join(HERE, 'dog_trees.tar.gz')
PROTEIN_ARCHIVE = os.path.join(HERE, 'protein_njt.tar.gz')
EXTRACT_DIR = os.path.join(HERE, '.extracted')

TREE_DIR = os.path.join(EXTRACT_DIR, 'dog_trees')
EXPECTED_FILE = os.path.join(EXTRACT_DIR, 'expected_dog_splits.json')
PROTEIN_DIR = os.path.join(EXTRACT_DIR, 'protein_njt')

# DOG trees: n_samples:n_paralogs:count
TREE_SPEC = "4:1:20,4:3:20,8:2:20,10:5:20,20:4:20,30:8:15,50:5:10,60:10:5"
TREE_SEED = 21
# Protein families: n_samples:n_paralogs:count
PROTEIN_SPEC = "4:2:20,8:3:20,15:4:20,30:5:15,60:4:8,100:5:3"
PROTEIN_SEED = 5

# (skip_merge_back, rooting_seeds, fixation_index_cutoff)
CONFIGS = [(False, 1, 0.1), (True, 1, 0.1), (False, 1, 0.5), (False, 5, 0.1), (False, 4, 0.5)]


def config_tag(skip, rs, fic):
    return f"skip{int(skip)}_rs{rs}_fic{fic}"


# --------------------------------------------------------------------------- #
# Archives
# --------------------------------------------------------------------------- #
def extract(force=False):
    """Extract both fixture archives into EXTRACT_DIR (if missing or out of date)."""
    os.makedirs(EXTRACT_DIR, exist_ok=True)
    for archive in (DOG_ARCHIVE, PROTEIN_ARCHIVE):
        stamp = os.path.join(EXTRACT_DIR, os.path.basename(archive) + '.extracted')
        if (not force and os.path.isfile(stamp)
                and os.path.getmtime(stamp) >= os.path.getmtime(archive)):
            continue
        with tarfile.open(archive, 'r:gz') as tf:
            for member in tf.getmembers():
                if member.name.startswith(('/', '..')) or '/../' in member.name:
                    raise ValueError(f"Unsafe path in {archive}: {member.name}")
            tf.extractall(EXTRACT_DIR)
        with open(stamp, 'w') as fh:
            fh.write('ok\n')


def _make_archive(archive, base_dir, names):
    with tarfile.open(archive, 'w:gz') as tf:
        for name in names:
            tf.add(os.path.join(base_dir, name), arcname=name)


# --------------------------------------------------------------------------- #
# DOG trees
# --------------------------------------------------------------------------- #
def _rand_topology(names, rng, zero_frac):
    nodes = []
    for nm in names:
        t = Tree()
        t.name = nm
        nodes.append(t)
    while len(nodes) > 1:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        p = Tree()
        for c in (a, b):
            c.dist = 0.0005 if rng.random() < zero_frac else rng.expovariate(10)
            p.add_child(c)
        nodes.append(p)
    return nodes[0]


def _make_family(rng, n_samples, n_paralogs, loss, multi_dom, zero_frac):
    samples = [f"S{i}" for i in range(n_samples)]
    copies = []
    lt_counter = 0
    for _ in range(n_paralogs):
        leaves = []
        for s in samples:
            if rng.random() < loss:
                continue
            lt_counter += 1
            lt = f"{s}_{lt_counter:05d}"
            leaves.append(f"{s}|{lt}|1|100")
            if rng.random() < multi_dom:  # second domain from the same protein
                leaves.append(f"{s}|{lt}|150|250")
        if len(leaves) >= 2:
            copies.append(_rand_topology(leaves, rng, zero_frac))
    nodes = copies
    while len(nodes) > 1:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        p = Tree()
        for c in (a, b):
            c.dist = rng.expovariate(2)
            p.add_child(c)
        nodes.append(p)
    t = nodes[0]
    if len(t) > 2:
        t.unroot()
    return t


def write_trees(outdir, spec=TREE_SPEC, seed=TREE_SEED):
    """Synthetic gene-family trees with bofasa domain leaf names (sample|lt|start|end)."""
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(seed)
    k = 0
    for item in spec.split(','):
        ns, npar, cnt = map(int, item.split(':'))
        for _ in range(cnt):
            t = _make_family(rng, ns, npar, loss=rng.choice([0.0, 0.1, 0.3]),
                             multi_dom=rng.choice([0.0, 0.05]),
                             zero_frac=rng.choice([0.0, 0.2, 0.6]))
            if len(t) < 3:
                continue
            k += 1
            t.write(outfile=os.path.join(outdir, f"OG_{k:04d}.tre"), format=1)
    return k


def read_split_file(path):
    """DOG split file -> ordered list of member lists (index i == cluster <og>_<i>)."""
    if not os.path.isfile(path):
        return None
    clusters = {}
    for line in open(path):
        cid, dom = line.rstrip('\n').split('\t')
        clusters.setdefault(cid, set()).add(dom)
    ordered = sorted(clusters, key=lambda c: int(c.rsplit('_', 1)[1]))
    return [sorted(clusters[c]) for c in ordered]


def protein_style_tree(tre_file):
    """Midpoint-rooted copy of a DOG tree with protein-style leaf names (sample|lt)."""
    t = Tree(tre_file)
    seen = set()
    for i, leaf in enumerate(t.get_leaves()):
        nm = '|'.join(leaf.name.split('|')[:2])
        if nm in seen:
            nm = f"{nm}_dup{i}"
        seen.add(nm)
        leaf.name = nm
    t.set_outgroup(t.get_midpoint_outgroup())
    return t


def compute_dog_expected(core, tree_dir):
    log = logging.getLogger("bofasa_fixture")
    trees = sorted(f for f in os.listdir(tree_dir) if f.endswith('.tre'))
    expected = {"dog_splits": {}, "protein_splits": {}}
    with tempfile.TemporaryDirectory() as tmp:
        for skip, rs, fic in CONFIGS:
            tag = config_tag(skip, rs, fic)
            expected["dog_splits"][tag] = {}
            for f in trees:
                og = f[:-len('.tre')]
                out = os.path.join(tmp, f"{tag}_{og}.txt")
                core.split_dogs([og, os.path.join(tree_dir, f), out, skip, rs, fic, 1, log])
                expected["dog_splits"][tag][og] = read_split_file(out)
    for f in trees:
        og = f[:-len('.tre')]
        t = protein_style_tree(os.path.join(tree_dir, f))
        samples = set(l.name.split('|')[0] for l in t.get_leaves())
        sp = core.recursive_splitting(t.copy(), samples, are_proteins=True)
        sp = core.further_split_disjoint_dog_partitions(t, sp)
        expected["protein_splits"][og] = [sorted(x) for x in sp]
    return expected


# --------------------------------------------------------------------------- #
# Protein (split_njt) dataset
# --------------------------------------------------------------------------- #
def make_protein_dataset(spec=PROTEIN_SPEC, seed=PROTEIN_SEED):
    """
    Synthetic proteins described by DOG counts, grouped into coarse clusters (one
    per gene family): paralogs with accessory domains, fragments (single domain),
    occasional domain fusions and double-copy domains.
    """
    rng = random.Random(seed)
    prot_dogs = {}
    clusters = []
    dog_id = 0
    fam = 0
    for item in spec.split(','):
        ns, npar, cnt = map(int, item.split(':'))
        for _ in range(cnt):
            fam += 1
            pool = [f"DOG_{dog_id + j:06d}" for j in range(rng.randint(3, 12))]
            dog_id += len(pool)
            core = pool[:rng.randint(1, 3)]
            members = []
            for par in range(npar):
                acc = rng.sample(pool, k=min(len(pool), rng.randint(0, 3)))
                for s in range(ns):
                    if rng.random() < 0.15:
                        continue
                    p = f"S{s}|F{fam}_P{par}_{s}"
                    dogs = {}
                    kind = rng.random()
                    if kind < 0.1:
                        dogs[rng.choice(pool)] = 1
                    else:
                        for d in core:
                            if rng.random() < 0.9:
                                dogs[d] = 2 if rng.random() < 0.1 else 1
                        for d in acc:
                            if rng.random() < 0.7:
                                dogs[d] = 1
                        if kind > 0.95:
                            dogs[f"DOG_{rng.randint(0, max(dog_id - 1, 0)):06d}"] = 1
                    if not dogs:
                        dogs[core[0]] = 1
                    prot_dogs[p] = dogs
                    members.append(p)
            if len(members) >= 4:
                clusters.append(sorted(members))
    return prot_dogs, clusters


def njt_inputs(work_dir, ci, members, log_object, prot_dogs=None):
    """Write the cluster listing and return the split_njt input list."""
    og = f"CoarseOG_{ci}"
    listing = os.path.join(work_dir, og + '.txt')
    with open(listing, 'w') as fh:
        fh.write('\n'.join(members) + '\n')
    inp = [og, listing, os.path.join(work_dir, og + '.phylip'),
           os.path.join(work_dir, og + '.tre'), os.path.join(work_dir, og + '.split'),
           log_object, 1]
    if prot_dogs is not None:
        inp.append(prot_dogs)
    return inp


def read_njt_split(path):
    if not os.path.isfile(path):
        return None
    return [line.split() for line in open(path) if line.strip()]


def compute_protein_fixture(core, out_dir):
    """Run split_njt (with real FastME) on every cluster; keep trees and splits."""
    if shutil.which('fastme') is None:
        raise RuntimeError("fastme must be on PATH to regenerate protein fixtures")
    prot_dogs, clusters = make_protein_dataset()
    core.protein_dogs.clear()
    core.protein_dogs.update({p: dict(d) for p, d in prot_dogs.items()})
    log = logging.getLogger("bofasa_fixture")
    tree_out = os.path.join(out_dir, 'trees')
    os.makedirs(tree_out, exist_ok=True)
    expected = {}
    with tempfile.TemporaryDirectory() as tmp:
        for ci, members in enumerate(clusters):
            inp = njt_inputs(tmp, ci, members, log)
            core.split_njt(inp)
            shutil.copy(inp[3], os.path.join(tree_out, f"CoarseOG_{ci}.tre"))
            expected[f"CoarseOG_{ci}"] = read_njt_split(inp[4])
    with open(os.path.join(out_dir, 'dataset.json'), 'w') as fh:
        json.dump({"protein_dogs": prot_dogs, "clusters": clusters}, fh,
                  separators=(',', ':'), sort_keys=True)
    with open(os.path.join(out_dir, 'expected_splits.json'), 'w') as fh:
        json.dump(expected, fh, separators=(',', ':'), sort_keys=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--package', default='bofasa',
                    help='Package whose core.py produces the expected results.')
    ap.add_argument('--skip-protein', action='store_true',
                    help='Only regenerate the DOG tree fixtures.')
    args = ap.parse_args()
    core = importlib.import_module(f"{args.package}.core")
    with tempfile.TemporaryDirectory() as build:
        n = write_trees(os.path.join(build, 'dog_trees'))
        with open(os.path.join(build, 'expected_dog_splits.json'), 'w') as fh:
            json.dump(compute_dog_expected(core, os.path.join(build, 'dog_trees')), fh,
                      separators=(',', ':'), sort_keys=True)
        _make_archive(DOG_ARCHIVE, build, ['dog_trees', 'expected_dog_splits.json'])
        print(f"wrote {DOG_ARCHIVE} ({n} trees)")
        if not args.skip_protein:
            compute_protein_fixture(core, os.path.join(build, 'protein_njt'))
            _make_archive(PROTEIN_ARCHIVE, build, ['protein_njt'])
            print(f"wrote {PROTEIN_ARCHIVE}")
    extract(force=True)


if __name__ == '__main__':
    sys.exit(main())
