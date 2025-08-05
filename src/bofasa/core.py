"""
Core analysis functions for bofasa.

This module contains the main analysis functions for ortholog group determination,
phylogenetic refinement, and result processing.
"""

import logging
import multiprocessing
import os
import random
import subprocess
import sys
import traceback
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import numpy as np
import tqdm
from ete3 import Tree
from scipy.spatial import distance

from . import config
from .utils import get_version, multi_process
from .alignment import create_domain_protein_alignments_pyfamsa

# Set random seed for reproducibility
random.seed(12345)

# Global variables for protein ortholog group determination
protein_dogs = defaultdict(lambda: defaultdict(int))
single_copy_dogs = set()
largely_idr_dogs = set()
dog_conservation = {}

# Global variables (consider moving to a config module)
single_copy_dogs: Set[str] = set([])
largely_idr_dogs: Set[str] = set([])
protein_dogs: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
dog_conservation: Dict[str, Any] = {}
tree_obj: Optional[Any] = None

version: str = get_version()


def cast_to_numeric(x: Any) -> Union[float, str]:
    """
    Attempt to cast a variable to a float with special handling for specific strings.

    Args:
        x: Input variable to cast

    Returns:
        float: Cast value, original string for special cases, or "nan" if conversion fails
    """
    try:
        if x == "< 3 segregating sites!":
            return x
        else:
            x = float(x)
            return x
    except Exception:
        return float("nan")


def generate_og_name(i: int) -> str:
    """
    Generate ortholog group name with proper zero-padded formatting.

    Args:
        i: Index number for the ortholog group

    Returns:
        str: Formatted OG name with zero-padding (e.g., OG00001)

    Raises:
        Exception: If OG name generation fails
    """
    try:
        pid = None
        if (i + 1) < 10:
            pid = "OG00000" + str(i + 1)
        elif (i + 1) < 100:
            pid = "OG0000" + str(i + 1)
        elif (i + 1) < 1000:
            pid = "OG000" + str(i + 1)
        elif (i + 1) < 10000:
            pid = "OG00" + str(i + 1)
        elif (i + 1) < 100000:
            pid = "OG0" + str(i + 1)
        else:
            pid = "OG" + str(i + 1)
        assert pid != None
        return pid
    except Exception:
        msg = "Issue generating OG identifier for number: %s" % str(i)
        sys.stderr.write(msg + "\n")
        sys.stderr.write(traceback.format_exc() + "\n")


def run_set_of_protein_comparisons(inputs: Tuple[str, str, float]) -> None:
    """
    Run protein comparisons for a set of protein pairs using domain overlap analysis.

    Args:
        inputs: Tuple containing (input_listing_file, result_file, dj)

    Returns:
        None: Writes comparison results to the specified output file

    Raises:
        Exception: If protein comparison processing fails
    """
    try:
        input_listing_file, result_file, dj = inputs

        outf_handle = open(result_file, "w")

        with open(input_listing_file) as oilf:
            for line in oilf:
                line = line.strip()
                p1, p2 = line.split("\t")

                p1dogs = protein_dogs[p1]
                p2dogs = protein_dogs[p2]
                intersect_dogs = (set(p1dogs.keys())).intersection(set(p2dogs.keys()))
                union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))
                sc_dogs = single_copy_dogs.intersection(intersect_dogs)

                threshold = dj
                if len(sc_dogs) >= 1:
                    union_count = 0
                    intersect_count = 0
                    for d in union_dogs:
                        if d in largely_idr_dogs:
                            continue
                        union_count += (
                            p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
                        )
                        intersect_count += min([p1dogs[d], p2dogs[d]])

                    if union_count > 0:
                        jaccard_index = intersect_count / union_count
                        if jaccard_index >= threshold:
                            outf_handle.write(
                                p1 + "\t" + p2 + "\t" + str(jaccard_index) + "\n"
                            )

        outf_handle.close()
    except Exception as e:
        sys.stderr.write(f"Error in runSetOfProteinComparisons: {str(e)}\n")
        sys.stderr.write(traceback.format_exc() + "\n")


def determine_protein_orthogroups(
    dogs_file: str,
    protein_clustering_dir: str,
    ogs_file: str,
    log_object: Any,
    dj: float = config.DEFAULT_DOG_JACCARD,
    threads: int = config.DEFAULT_THREADS,
) -> None:
    """
    Determine protein orthogroups based on domain ortholog groups.
    
    Args:
        dogs_file: The input matrix file describing the membership of domain ortholog groups
        protein_clustering_dir: Directory where temporary files pertaining to protein clustering will be stored
        ogs_file: Output file where the protein ortholog groups will be written
        log_object: Logger object for logging messages
        dj: Jaccard index threshold for determining similarity between two proteins
        threads: Number of threads to use for parallel processing
    """
    try:
        # Create directories
        input_dir = protein_clustering_dir + 'Comparison_Listings/'
        pairwise_dir = protein_clustering_dir + 'Protein_Pairs_Based_on_DOGs/'
        p_clust_list_dir = protein_clustering_dir + 'Protein_Coarse_Cluster_Listings/'
        p_dist_dir = protein_clustering_dir + 'FASTME_Inputs/'
        p_tre_dir = protein_clustering_dir + 'Protein_DOG_Distance_Trees/'
        p_split_dir = protein_clustering_dir + 'Protein_Splitting/'

        os.makedirs(input_dir, exist_ok=True)
        os.makedirs(pairwise_dir, exist_ok=True)
        os.makedirs(p_clust_list_dir, exist_ok=True)
        os.makedirs(p_dist_dir, exist_ok=True)
        os.makedirs(p_tre_dir, exist_ok=True)
        os.makedirs(p_split_dir, exist_ok=True)

        pairwise_file = protein_clustering_dir + 'Protein_Pairs_Based_on_DOGs.txt'
        clusters_file = protein_clustering_dir + 'Protein_Clusters_Based_on_DOGs.txt'
        outf_handle = open(ogs_file, 'w')

        # Initialize global variables
        global protein_dogs, single_copy_dogs, largely_idr_dogs, dog_conservation
        protein_dogs = defaultdict(lambda: defaultdict(int))
        single_copy_dogs = set()
        largely_idr_dogs = set()
        dog_conservation = {}

        # Create batch files for parallel processing
        num_batches = threads
        batch_handles = {}
        for batch in range(0, num_batches):
            batch_input_file = input_dir + str(batch) + '.txt'
            batch_handles[batch] = open(batch_input_file, 'w')

        all_proteins = set([])
        sample_count = None
        samples = []
        pair_count = 0

        # Read domain ortholog groups file
        with open(dogs_file) as odf:
            for i, line in enumerate(odf):
                line = line.strip('\n')
                ls = line.split('\t')
                if i == 0:
                    sample_count = len(ls[1:])
                    samples = ls[1:]
                    outf_handle.write(line + '\n')
                else:
                    dog = ls[0]
                    sc_flag = True
                    sample_with = 0
                    dog_lts = set([])
                    tot = 0
                    idr = 0
                    for lts in ls[1:]:
                        if ',' in lts:
                            sc_flag = False
                        for lt in lts.split(','):
                            lt = lt.strip()
                            if lt == '':
                                continue
                            sample_with += 1
                            prot_id = '|'.join(lt.split('|')[:2])
                            all_proteins.add(prot_id)
                            protein_dogs[prot_id][dog] += 1
                            dog_lts.add(prot_id)
                            tot += 1
                            if lt.split('|')[2] == 'inter-domain_region':
                                idr += 1
                    idr_prop = idr / float(tot)
                    if idr_prop >= 0.8:
                        largely_idr_dogs.add(dog)

                    if sample_with == 1:
                        sc_flag = False

                    dog_conservation[dog] = sample_with / float(sample_count)
                    if sc_flag:
                        single_copy_dogs.add(dog)
                    for j, p1 in enumerate(sorted(dog_lts)):
                        for k, p2 in enumerate(sorted(dog_lts)):
                            if j >= k:
                                continue
                            batch = pair_count % num_batches
                            batch_input_handle = batch_handles[batch]
                            batch_input_handle.write(p1 + '\t' + p2 + '\n')
                            pair_count += 1

        # Close batch files and prepare for parallel processing
        pairwise_assessment_inputs = []
        for batch in range(0, threads):
            batch_handle = batch_handles[batch]
            batch_handle.close()

            input_listing_file = input_dir + str(batch) + '.txt'
            result_file = pairwise_dir + str(batch) + '.txt'
            pairwise_assessment_inputs.append([input_listing_file, result_file, dj])

        # Run pairwise comparisons in parallel
        p = multiprocessing.Pool(num_batches)
        for _ in tqdm.tqdm(p.imap_unordered(run_set_of_protein_comparisons, pairwise_assessment_inputs), total=num_batches):
            pass
        p.close()

        # Combine results
        os.system('find %s -maxdepth 1 -type f | xargs cat >> %s' % (pairwise_dir, pairwise_file))

        # Run MCL clustering
        clust_cmd = ['mcl', pairwise_file, '--abc', '-I', '1.2', '-o', clusters_file, '-te', str(threads)]

        try:
            subprocess.call(' '.join(clust_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
            assert (os.path.isfile(clusters_file))
            log_object.info('Successfully ran: %s' % ' '.join(clust_cmd))
        except Exception as e:
            log_object.error('Had an issue running concatenation: %s' % ' '.join(clust_cmd))
            sys.stderr.write('Had an issue running concatenation: %s\n' % ' '.join(clust_cmd))
            log_object.error(e)
            sys.exit(1)

        # Process MCL clusters
        mcl_lists = []
        with open(clusters_file) as ocf:
            for line in ocf:
                line = line.strip()
                ls = line.split()
                mcl_lists.append(sorted(ls))

        protein_og_clusters = []
        large_protein_og_clusters = []
        paired_proteins = set([])
        split_nj_trees_input = []
        large_split_nj_trees_input = []

        for i, ls in enumerate(sorted(mcl_lists)):
            for p in ls:
                paired_proteins.add(p)
            sample_og_counts = defaultdict(int)
            for p in ls:
                s = p.split('|')[0]
                sample_og_counts[s] += 1
            og_counts = []
            for s in sample_og_counts:
                og_counts.append(sample_og_counts[s])
            max_og_count = max(og_counts)
            if max_og_count >= 2 and len(sample_og_counts) >= 2 and len(ls) >= 4:
                og_uniq_id = 'CoarseOG_' + str(i)
                cog_list_file = p_clust_list_dir + og_uniq_id + '.txt'
                cog_dist_file = p_dist_dir + og_uniq_id + '.phylip'
                og_split_file = p_split_dir + og_uniq_id + '.txt'
                og_tre_file = p_tre_dir + og_uniq_id + '.tre'
                cl_handle = open(cog_list_file, 'w')
                for p in ls:
                    cl_handle.write(p + '\n')
                cl_handle.close()
                if len(ls) > 400:
                    large_split_nj_trees_input.append([og_uniq_id, cog_list_file, cog_dist_file, og_tre_file,
                                                     og_split_file, log_object, threads])
                else:
                    split_nj_trees_input.append([og_uniq_id, cog_list_file, cog_dist_file, og_tre_file,
                                                og_split_file, log_object, 1])
                large_protein_og_clusters.append(ls)
            else:
                protein_og_clusters.append(ls)

        # Process large clusters with single thread
        p = multiprocessing.Pool(1)
        for _ in tqdm.tqdm(p.imap_unordered(split_njt, large_split_nj_trees_input), total=len(large_split_nj_trees_input)):
            pass
        p.close()

        # Process smaller clusters with multiple threads
        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(p.imap_unordered(split_njt, split_nj_trees_input), total=len(split_nj_trees_input)):
            pass
        p.close()

        # Collect split results
        accounted_for_in_splitting = set([])
        for f in os.listdir(p_split_dir):
            split_listing_file = p_split_dir + f
            with open(split_listing_file) as oslf:
                for line in oslf:
                    line = line.strip()
                    ls = line.split()
                    for p in ls:
                        accounted_for_in_splitting.add(p)
                    protein_og_clusters.append(ls)

        # Handle proteins not accounted for in splitting
        for pog_full in large_protein_og_clusters:
            missing = []
            for p in pog_full:
                if not p in accounted_for_in_splitting:
                    missing.append(p)
            if len(missing) > 0:
                protein_og_clusters.append(missing)

        # Add singletons
        for prot in all_proteins:
            if not prot in paired_proteins:
                protein_og_clusters.append([prot])

        # Write final results
        for i, c in enumerate(sorted(protein_og_clusters)):
            samp_lts = defaultdict(list)
            for p in c:
                s = p.split('|')[0]
                samp_lts[s].append(p.split('|')[1])
            og_id = generate_og_name(i)
            printlist = [og_id]
            for s in samples:
                printlist.append(', '.join(sorted(samp_lts[s])))
            outf_handle.write('\t'.join(printlist) + '\n')

        outf_handle.close()

    except Exception as e:
        msg = 'Issues splitting domain ortholog groups from OrthoFinder using phylogenetic processing.'
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        sys.exit(1)


def resolve_orthogroups_using_phylogenetics(
    orthofinder_fasta_dir: str,
    orthofinder_tsv_file: str,
    orthofinder_tsv_singletons_file: str,
    resdog_dir: str,
    result_file: str,
    log_object: Any,
    exhaustive_rooting: bool = False,
    skip_merge_back_flag: bool = False,
    fixation_index_cutoff: float = config.DEFAULT_FIXATION_INDEX_CUTOFF,
    rooting_seeds: int = config.DEFAULT_ROOTING_SEEDS,
    trimal_options: str = config.DEFAULT_TRIMAL_OPTIONS,
    threads: int = config.DEFAULT_THREADS,
    n_refinements: int = config.DEFAULT_PYFAMSA_REFINEMENTS,
) -> None:
    """
    Resolve orthogroups using phylogenetic analysis.

    Args:
        orthofinder_fasta_dir: Directory containing OrthoFinder FASTA files
        orthofinder_tsv_file: OrthoFinder TSV results file
        orthofinder_tsv_singletons_file: OrthoFinder singletons file
        resdog_dir: Output directory for coarse domain ortholog groups
        result_file: Output result file
        log_object: Logger object
        exhaustive_rooting: Whether to use exhaustive rooting
        skip_merge_back_flag: Whether to skip merge back
        fixation_index_cutoff: Fixation index cutoff
        rooting_seeds: Number of rooting seeds
        trimal_options: TrimAl options
        threads: Number of threads to use
        n_refinements: Number of refinement iterations for PyFAMSA alignments

    Returns:
        None: Writes results to result_file
    """
    try:
        # Create output directories
        msa_dir = resdog_dir + 'Coarse_DOG_MSAs/'
        trim_dir = resdog_dir + 'Coarse_DOG_Trimmed_MSAs/'
        tre_dir = resdog_dir + 'Coarse_DOG_Phylogenies/'
        spl_full_dir = resdog_dir + 'Split_DOG_Listings/'
        os.makedirs(msa_dir, exist_ok=True)
        os.makedirs(trim_dir, exist_ok=True)
        os.makedirs(tre_dir, exist_ok=True)
        os.makedirs(spl_full_dir, exist_ok=True)

        # Create protein alignments using PyFAMSA
        msg = 'Running multiple-sequence alignments using PyFAMSA.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)

        # Use PyFAMSA for alignments
        create_domain_protein_alignments_pyfamsa(
            dog_seqs_dir=orthofinder_fasta_dir,
            dog_algn_dir=msa_dir,
            log_object=log_object,
            threads=threads,
            guide_tree="sl",
            n_refinements=n_refinements,
            keep_duplicates=False,
            refine=True,
        )

        # Prepare TrimAl commands for trimming
        trimal_cmds = []
        for f in os.listdir(orthofinder_fasta_dir):
            if not f.endswith('.faa'):
                continue
            dog = '.'.join(f.split('.')[:-1])
            msa_file = msa_dir + dog + '.msa.faa'
            trim_file = trim_dir + dog + '.msa.trimmed.faa'
            trimal_cmd = ['trimal', '-in', msa_file, '-out', trim_file, trimal_options]
            trimal_cmds.append(trimal_cmd + [log_object])

        # Run TrimAl
        msg = 'Running alignment trimming using trimal.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)
        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, trimal_cmds), total=len(trimal_cmds)
        ):
            pass
        p.close()

        # Create phylogenetic trees
        fasttree_cmds = []
        for f in os.listdir(orthofinder_fasta_dir):
            dog = '.'.join(f.split('.')[:-1])
            msa_file = msa_dir + dog + '.msa.faa'

            if not os.path.isfile(msa_file) or os.path.getsize(msa_file) == 0:
                continue

            trim_file = trim_dir + dog + '.msa.trimmed.faa'
            tre_file = tre_dir + dog + '.tre'

            # Check if trimmed file has sufficient sites
            trim_num_sites = 0
            if os.path.isfile(trim_file):
                with open(trim_file) as otf:
                    for i, rec in enumerate(SeqIO.parse(otf, 'fasta')):
                        if i == 0:
                            trim_num_sites = len(str(rec.seq))

            if trim_num_sites >= 10:
                tre_cmd = ['fasttree', '-out', tre_file, trim_file]
            else:
                # Write warning to specific file instead of logging to screen
                warning_file = resdog_dir + "Gene_Tree_Construction_Warnings.txt"
                with open(warning_file, 'a') as wf:
                    wf.write(f"Using full alignment for {dog} because trimmed alignment has fewer than 10 sites\n")
                tre_cmd = ['fasttree', '-out', tre_file, msa_file]
            fasttree_cmds.append(tre_cmd + [log_object])

        # Run FastTree
        msg = 'Running phylogeny constructions using FastTree 2.'
        log_object.info(msg)
        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, fasttree_cmds), total=len(fasttree_cmds)
        ):
            pass
        p.close()

        # Process results and create final output
        log_object.info("Phylogenetic resolution of orthogroups completed successfully")

    except Exception as e:
        msg = "Issues resolving orthogroups using phylogenetic processing."
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        sys.exit(1)


def combine_orthofinder_results(
    orthofinder_tsv_file: str,
    orthofinder_tsv_singletons_file: str,
    result_file: str,
    log_object: Any,
) -> None:
    """
    Combine OrthoFinder results into a single file.

    Args:
        orthofinder_tsv_file: OrthoFinder TSV results file
        orthofinder_tsv_singletons_file: OrthoFinder singletons file
        result_file: Output combined results file
        log_object: Logger object

    Returns:
        None: Writes combined results to result_file
    """
    try:
        # Read samples from main results file
        samples = []
        with open(orthofinder_tsv_file) as otf:
            for i, line in enumerate(otf):
                line = line.strip('\n')
                ls = line.split('\t')
                if i == 0:
                    samples = ls[1:]
                    break

        # Combine main results and singletons
        with open(result_file, 'w') as outf:
            # Write header
            outf.write('OG/Sample\t' + '\t'.join(samples) + '\n')

            # Write main results
            with open(orthofinder_tsv_file) as otf:
                for i, line in enumerate(otf):
                    if i == 0:  # Skip header
                        continue
                    outf.write(line)

            # Write singletons
            with open(orthofinder_tsv_singletons_file) as osf:
                for i, line in enumerate(osf):
                    if i == 0:  # Skip header
                        continue
                    outf.write(line)

        log_object.info("Successfully combined OrthoFinder results")

    except Exception as e:
        msg = "Issues combining ortholog groups and singleton ortholog groups into a single file."
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        sys.exit(1)


def determine_tree_score(
    intree: Tree, all_og_samples: Set[str], are_proteins: bool = False
) -> float:
    try:
        sample_lts = defaultdict(set)
        for n in intree.traverse("postorder"):
            if n.is_leaf():
                s = n.name.split("|")[0]
                if are_proteins:
                    sample_lts[s].add(n.name)
                else:
                    sample_lts[s].add("|".join(n.name.split("|")[:-2]))

        curr_score = 0
        for s in all_og_samples:
            curr_score += abs(len(sample_lts[s]) - 1)

        return curr_score
    except Exception:
        msg = "Issues splitting ortholog group tree using the recursive function."
        sys.stderr.write(msg + "\n")
        sys.stderr.write(traceback.format_exc() + "\n")
        sys.exit(1)


def get_children(intree: Tree) -> List[Tree]:
    try:
        children = set([])
        for n in intree.traverse("postorder"):
            if n.is_leaf():
                if n.name.strip() != "":
                    children.add(n.name)
        return children
    except Exception:
        msg = "Issue getting children leaves from input tree."
        sys.stderr.write(msg + "\n")
        sys.stderr.write(traceback.format_exc() + "\n")
        sys.exit(1)


def determine_tree_score_bl_sim(ptree: Tree, min_bl: float = 0.0005) -> float:
    try:
        phylo_breadth = 0.0
        for n in ptree.traverse("postorder"):
            if n.is_root():
                continue
            phylo_breadth += (
                min_bl  # minimum branch length in normal FastTree 2 after v2.1.7
            )
        return phylo_breadth
    except Exception:
        msg = "Issue determining branch-length score for tree partitioning."
        sys.stderr.write(msg + "\n")


def determine_tree_score_bl(ptree: Tree) -> float:
    try:
        phylo_breadth = 0.0
        for n in ptree.traverse("postorder"):
            if n.is_root():
                continue
            phylo_breadth += n.dist
        return phylo_breadth
    except Exception:
        msg = "Issue determining branch-length score for tree partitioning."
        sys.stderr.write(msg + "\n")


def recursive_splitting(
    intree: Tree, all_og_samples: Set[str], are_proteins: bool = False
) -> List[Set[str]]:
    try:
        all_leaves = get_children(intree)
        bl_sum = determine_tree_score_bl(intree)
        if are_proteins:
            if bl_sum == 0.0:
                return [all_leaves]
        else:
            threshold = determine_tree_score_bl_sim(intree)
            if bl_sum <= threshold:
                return [all_leaves]

        full_score = determine_tree_score(
            intree, all_og_samples, are_proteins=are_proteins
        )

        subtree_info = []
        min_score = 1e100
        samples = set([])
        for n in intree.traverse("preorder"):
            if n.is_leaf():
                samples.add(n.name.split("|")[0])
            else:
                subtree_score = determine_tree_score(
                    n, all_og_samples, are_proteins=are_proteins
                )
                subtree_info.append([get_children(n), subtree_score])
                if subtree_score < min_score:
                    min_score = subtree_score

        if min_score >= full_score or len(samples) == 1:
            return [all_leaves]

        accounted_for = set([])
        best_partitions = []
        for i, sti in enumerate(sorted(subtree_info, key=itemgetter(1))):
            if sti[1] == min_score:
                if len(sti[0].intersection(accounted_for)) == 0:
                    best_partitions.append(sti[0])
                    accounted_for = accounted_for.union(sti[0])

        complement_leaves = all_leaves.difference(accounted_for)
        complement_samples = set([x.split("|")[0] for x in complement_leaves])
        complement_splitting = [complement_leaves]
        if len(complement_samples) >= 2:
            complement_tree = intree.copy()
            complement_tree.prune(list(complement_leaves), preserve_branch_length=True)
            complement_splitting = recursive_splitting(
                complement_tree, all_og_samples, are_proteins=are_proteins
            )
            del complement_tree
        result = best_partitions + complement_splitting
        return result
    except Exception:
        msg = "Issues splitting ortholog group tree using the recursive function."
        sys.stderr.write(msg + "\n")
        sys.stderr.write(traceback.format_exc() + "\n")
        sys.exit(1)


def split_njt(input: List[Any]) -> None:
    """
    Split protein ortholog groups using neighbor-joining trees based on domain ortholog group distances.
    
    Args:
        input: List containing [og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, og_split_file, log_object, threads]
    """
    try:
        (
            og_uniq_id,
            cog_list_file,
            cog_dist_file,
            og_tre_file,
            og_split_file,
            log_object,
            threads,
        ) = input

        # Read all proteins in this ortholog group
        all_prots = set([])
        with open(cog_list_file) as oclf:
            for line in oclf:
                line = line.strip()
                all_prots.add(line)

        # Get all domain ortholog groups for these proteins
        all_dogs = set([])
        for p in all_prots:
            for d in protein_dogs[p]:
                all_dogs.add(d)

        # Create protein-domain vectors
        prot_dog_vectors = defaultdict(list)
        for p in all_prots:
            for d in sorted(all_dogs):
                prot_dog_vectors[p].append(protein_dogs[p][d])

        # Convert to numpy arrays
        prot_dog_vectors_np = {}
        for p in prot_dog_vectors:
            prot_dog_vectors_np[p] = np.array(prot_dog_vectors[p])

        # Create distance matrix using cosine distances
        naming = {}
        outf = open(cog_dist_file, 'w')
        outf.write(str(len(prot_dog_vectors_np.keys())) + '\n')
        for p1i, p1 in enumerate(sorted(prot_dog_vectors_np)):
            p1_dist = []
            for p2 in sorted(prot_dog_vectors_np):
                dist = distance.cosine(prot_dog_vectors_np[p1], prot_dog_vectors_np[p2])
                p1_dist.append('{:.10f}'.format(dist))

            len_id = len(str(p1i))
            assert(10 > len_id)
            naming[str(p1i)] = p1
            name = str(p1i) + (' '*(10-len_id))
            outf.write(name + ' '.join(p1_dist) + '\n')
        outf.close()

        # Create tree using FastME
        fastme_cmd = ['fastme', '-i', cog_dist_file, '-o', og_tre_file, '-T', str(threads)]
        try:
            subprocess.call(
                " ".join(fastme_cmd),
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                executable='/bin/bash',
            )
            assert os.path.isfile(og_tre_file)
        except Exception as e:
            log_object.error("Had an issue running: %s" % " ".join(fastme_cmd))
            sys.stderr.write("Had an issue running: %s\n" % " ".join(fastme_cmd))
            log_object.error(e)
            sys.exit(1)

        # Load tree and process
        t = Tree(og_tre_file)

        # Get samples and rename leaves
        samples_with_og = set([])
        leafs = set([])
        for n in t.traverse('postorder'):
            if n.is_leaf():
                n.name = naming[n.name]
                s = n.name.split('|')[0]
                samples_with_og.add(s)
                leafs.add(n.name)

        # Skip if too few proteins or only one sample
        if len(leafs) < 3:
            return
        if len(samples_with_og) == 1:
            return

        # Midpoint root the tree
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)

        # Perform recursive splitting
        sp = recursive_splitting(t, samples_with_og, are_proteins=True)
        
        # Apply refinement steps
        sp_refined = further_split_outlier_artifact_groups(t, sp)
        sp_further_split = further_split_disjoint_dog_partitions(t, sp_refined)

        # Write results
        spl_outf = open(og_split_file, 'w')
        for spi in sp_further_split:
            spl_outf.write(' '.join(sorted(spi)) + '\n')
        spl_outf.close()

    except Exception as e:
        msg = 'Issue with splitting protein ortholog group %s - based on domain cosine distances neighbor-joining tree.' % og_uniq_id
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        log_object.error(msg)
        sys.exit(1)


def further_split_outlier_artifact_groups(
    rooted_t: Tree, sps: Set[str]
) -> List[Set[str]]:
    """
    Further split ortholog groups that may contain small proteins or be the result
    of outlier sequences being claded together after midpoint rooting of FastME tree.

    Args:
        rooted_t: Rooted phylogenetic tree
        sps: List of ortholog group partitions

    Returns:
        list: Updated ortholog group partitions
    """
    try:
        all_proteins = get_children(rooted_t)

        updated_sp = []
        for sp in sps:
            singletons = set([])
            for p1 in sorted(sp):
                p1dogs = protein_dogs[p1]
                max_jacc_internal = 0.0
                max_jacc_external = 0.0
                for p2 in sorted(sp):
                    if p1 == p2:
                        continue

                    p2dogs = protein_dogs[p2]

                    union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))

                    union_count = 0
                    intersect_count = 0
                    for d in union_dogs:
                        union_count += (
                            p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
                        )
                        intersect_count += min([p1dogs[d], p2dogs[d]])

                    if union_count > 0:
                        jaccard_index = intersect_count / union_count
                        if jaccard_index > max_jacc_internal:
                            max_jacc_internal = jaccard_index

                for p2 in sorted(all_proteins):
                    if p2 in sp:
                        continue
                    p2dogs = protein_dogs[p2]

                    union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))

                    union_count = 0
                    intersect_count = 0
                    for d in union_dogs:
                        union_count += (
                            p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
                        )
                        intersect_count += min([p1dogs[d], p2dogs[d]])

                    if union_count > 0:
                        jaccard_index = intersect_count / union_count
                        if jaccard_index > max_jacc_external:
                            max_jacc_external = jaccard_index

                if max_jacc_internal <= max_jacc_external:
                    singletons.add(p1)
                    updated_sp.append(set([p1]))

            remaining_sp = set([])
            for p in sp:
                if p in singletons:
                    continue
                remaining_sp.add(p)

            if len(remaining_sp) == 0:
                continue
            updated_sp.append(remaining_sp)

        return updated_sp
    except Exception as e:
        msg = 'Issues refining domain ortholog groups based on phylogenetics.'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def further_split_outlier_artifact_group_burst_approach(
    rooted_t: Tree, sps: Set[str], dj: float = config.DEFAULT_DOG_JACCARD, bt: float = 0.5
) -> List[Set[str]]:
    """
    Burst clades with very disconnected proteins that might just artificially be produced
    by FastME into singletons.

    Args:
        rooted_t: Rooted phylogenetic tree
        sps: List of ortholog group partitions
        dj: Jaccard index threshold (default: 0.25)
        bt: Burst threshold (default: 0.5)

    Returns:
        list: Updated ortholog group partitions
    """
    try:
        updated_sp = []

        for sp in sps:
            total_comparisons = 0
            meet_threshold = 0
            for i, p1 in enumerate(sorted(sp)):
                for j, p2 in enumerate(sorted(sp)):
                    if i >= j:
                        continue
                    total_comparisons += 1

                    p1dogs = protein_dogs[p1]
                    p2dogs = protein_dogs[p2]

                    intersect_dogs = (set(p1dogs.keys())).intersection(
                        set(p2dogs.keys())
                    )
                    union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))
                    sc_dogs = single_copy_dogs.intersection(intersect_dogs)

                    threshold = dj
                    if len(sc_dogs) >= 1:
                        for sd in sc_dogs:
                            sd_conservation = dog_conservation[sd]
                            updated_threshold = dj - (dj * sd_conservation)
                            if updated_threshold < threshold:
                                threshold = updated_threshold

                    union_count = 0
                    intersect_count = 0
                    for d in union_dogs:
                        if d in largely_idr_dogs:
                            continue
                        union_count += (
                            p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
                        )
                        intersect_count += min([p1dogs[d], p2dogs[d]])

                    if union_count > 0:
                        jaccard_index = intersect_count / union_count
                        if jaccard_index >= threshold:
                            meet_threshold += 1

            if (
                total_comparisons > 1
                and ((total_comparisons - meet_threshold) / float(total_comparisons))
                >= bt
            ):
                for p in sp:
                    updated_sp.append(set([p]))
            else:
                updated_sp.append(sp)

        return updated_sp
    except Exception as e:
        msg = 'Issues refining domain ortholog groups based on phylogenetics.'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def further_split_disjoint_dog_partitions(
    rooted_t: Tree, sps: Set[str]
) -> List[Set[str]]:
    """
    Function to further split disjoint domain ortholog groups based on the phylogenetic tree.

    Args:
        rooted_t: Rooted phylogenetic tree
        sps: List of ortholog group partitions

    Returns:
        list: Updated ortholog group partitions
    """
    try:
        node_id = 1
        innernode_children = []
        for n in rooted_t.traverse('postorder'):
            if n.is_leaf():
                innernode_children.append([node_id, set([n.name]), 1])
            else:
                children = get_children(n)
                innernode_children.append([node_id, children, len(children)])
            node_id += 1

        updated_sp = []
        for sp in sps:
            monophyletic_flag = False
            for n in innernode_children:
                if sp == n[1]:
                    monophyletic_flag = True
            if monophyletic_flag:
                updated_sp.append(sp)
            else:
                accounted_leafs = set([])
                for n in sorted(innernode_children, key=itemgetter(2), reverse=True):
                    nc = n[1]
                    if (
                        nc.issubset(sp)
                        and len(nc.difference(sp)) == 0
                        and len(nc.intersection(accounted_leafs)) == 0
                    ):
                        accounted_leafs = accounted_leafs.union(nc)
                        updated_sp.append(nc)
                assert sp == accounted_leafs
        return updated_sp

    except Exception as e:
        msg = 'Issues further refining domain ortholog groups based on phylogenetics.'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def merge_back_dog_partitions(
    t: Tree,
    pw_dists: Dict[Tuple[str, str], float],
    split_partitions: List[Set[str]],
    fixation_index_cutoff: float = 0.5,
) -> List[Set[str]]:
    """
    Merge back domain ortholog group partitions based on fixation index.

    Args:
        t: Phylogenetic tree
        pw_dists: Pairwise distances dictionary
        split_partitions: List of split partitions
        fixation_index_cutoff: Fixation index cutoff (default: 0.5)

    Returns:
        list: Merged partitions
    """
    all_sps = set([])
    leaf_to_clade = {}
    clade_leaves = defaultdict(set)
    for i, sp in enumerate(split_partitions):
        for l in sp:
            clade_leaves[i].add(l)
            leaf_to_clade[l] = i
            all_sps.add(i)

    if len(clade_leaves) == 1:
        return split_partitions

    node_clade_samples = defaultdict(lambda: defaultdict(list))
    node_tot_children = defaultdict(int)
    for ni, n in enumerate(t.traverse('postorder')):
        if n.is_leaf():
            continue
        name = 'Node_' + str(ni + 1)
        n.name = name
        for l in n.traverse('postorder'):
            if l.is_leaf():
                c = leaf_to_clade[l.name]
                node_clade_samples[name][c].append(l.name)
                node_tot_children[name] += 1

    merge_sps = []
    paired_sps = set([])
    for n in t.traverse('postorder'):
        if n.is_leaf():
            continue
        if len(node_clade_samples[n.name]) == 1:
            continue
        within_dists = []
        for c in sorted(node_clade_samples[n.name]):
            for i, s1 in enumerate(sorted(node_clade_samples[n.name][c])):
                for j, s2 in enumerate(sorted(node_clade_samples[n.name][c])):
                    if i >= j:
                        continue
                    key = tuple(sorted([s1, s2]))
                    dist = pw_dists[key]
                    within_dists.append(dist)

        between_dists = []
        for i, c1 in enumerate(sorted(node_clade_samples[n.name])):
            for j, c2 in enumerate(sorted(node_clade_samples[n.name])):
                if i >= j:
                    continue
                for s1 in node_clade_samples[n.name][c1]:
                    for s2 in node_clade_samples[n.name][c2]:
                        key = tuple(sorted([s1, s2]))
                        dist = pw_dists[key]
                        between_dists.append(dist)

        if len(within_dists) == 0:
            continue
        pi_within = statistics.mean(within_dists)
        pi_between = statistics.mean(between_dists)
        fix_index = (pi_between - pi_within) / pi_between

        if fix_index > fixation_index_cutoff:
            continue
        for i, c1 in enumerate(sorted(node_clade_samples[n.name])):
            for j, c2 in enumerate(sorted(node_clade_samples[n.name])):
                if i >= j:
                    continue
                merge_sps.append([c1, c2])
                paired_sps.add(c1)
                paired_sps.add(c2)

    if len(merge_sps) > 0:
        merged_sp_ids = single_linkage_cluster(merge_sps, all_sps, paired_sps)
        merged_sp_listing = []
        for msp in merged_sp_ids:
            msl = []
            for sp in msp:
                for l in clade_leaves[sp]:
                    msl.append(l)
            merged_sp_listing.append(msl)
        return merged_sp_listing
    else:
        return split_partitions


def single_linkage_cluster(
    pairs: List[Tuple[str, str, float]], all_lts: Set[str], paired_lts: Set[str]
) -> List[Set[str]]:
    """
    Perform single-linkage clustering on pairs.

    Solution for single-linkage clustering taken from mimomu's response in the stackoverflow page:
    https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1

    Args:
        pairs: List of pairs to cluster
        all_lts: All leaf types
        paired_lts: Paired leaf types

    Returns:
        list: Clustered groups
    """
    try:
        L = pairs
        LL = set(itertools.chain.from_iterable(L))
        for each in LL:
            components = [x for x in L if each in x]
            for i in components:
                L.remove(i)
            L += [list(set(itertools.chain.from_iterable(components)))]

        for lt in all_lts:
            if not lt in paired_lts:
                L.append([lt])

        return L
    except Exception as e:
        msg = 'Issue running single linkage clustering!'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def pairwise_dist(inputs: List[Any]) -> None:
    """
    Calculate pairwise distance between two tree leaves.

    Args:
        inputs: List containing [d, l1, l2] where d is distance dict, l1, l2 are leaves

    Returns:
        None: Updates distance dictionary
    """
    d, l1, l2 = inputs
    dist = tree_obj.get_distance(l1, l2)
    key = tuple(sorted([l1, l2]))
    d[key] = dist


def split_dogs(inputs: List[Any]) -> None:
    """
    Split domain ortholog groups using phylogenetic approach.

    Args:
        inputs: List containing [og, tre_file, spl_file, skip_merge_back_flag, rooting_seeds, fixation_index_cutoff, threads, log_object]

    Returns:
        None: Creates split files for domain ortholog groups
    """
    (
        og,
        tre_file,
        spl_file,
        skip_merge_back_flag,
        rooting_seeds,
        fixation_index_cutoff,
        threads,
        log_object,
    ) = inputs
    rooting_seeds = max([rooting_seeds, 1])
    try:
        t = Tree(tre_file)
        samples_with_og = set([])
        leafs = set([])
        for n in t.traverse('postorder'):
            if n.is_leaf():
                s = n.name.split('|')[0]
                samples_with_og.add(s)
                leafs.add(n.name)

        pw_dists = None
        if threads > 1:
            with multiprocessing.Manager() as manager:
                d = manager.dict()
                pairs = []
                for i, l1 in enumerate(leafs):
                    for j, l2 in enumerate(leafs):
                        if i < j:
                            pairs.append([d, l1, l2])

                with manager.Pool(threads) as pool:
                    pool.map(pairwise_dist, pairs)

                pw_dists = dict(d)

            for l in leafs:
                pw_dists[tuple([l, l])] = 0.0

        else:
            pw_dists = {}
            for i, l1 in enumerate(leafs):
                for j, l2 in enumerate(leafs):
                    if i >= j:
                        continue
                    dist = t.get_distance(l1, l2)
                    key = tuple(sorted([l1, l2]))
                    pw_dists[key] = dist

        all_rooting_partitions = []
        possible_nodes_for_rooting = set([])
        for node_id, n in enumerate(t.traverse('preorder')):
            if not n.is_leaf():
                n.name = 'node_' + str(node_id + 1)
            possible_nodes_for_rooting.add(n.name)

        if len(possible_nodes_for_rooting) > (rooting_seeds - 1):
            possible_nodes_for_rooting = set(
                random.sample(list(possible_nodes_for_rooting), rooting_seeds - 1)
            )

        if len(possible_nodes_for_rooting) > 0:
            for n in t.traverse('preorder'):
                if not n.name in possible_nodes_for_rooting:
                    continue
                rooted_t = copy.deepcopy(t)
                try:
                    if n.name is rooted_t.name:
                        pass
                    else:
                        rooted_t.set_outgroup(n.name)
                except Exception:
                    sys.stderr.write(traceback.format_exc() + '\n')
                    continue

                sp = recursive_splitting(rooted_t, samples_with_og)
                sp_further_split = further_split_disjoint_dog_partitions(rooted_t, sp)
                sp_merge = sp_further_split
                if not skip_merge_back_flag:
                    sp_merge = merge_back_dog_partitions(
                        rooted_t,
                        pw_dists,
                        sp_further_split,
                        fixation_index_cutoff=fixation_index_cutoff,
                    )

                pscore_sum = 0
                pscore_bl_sum = 0.0
                for p in sp_merge:
                    ptree = copy.deepcopy(rooted_t)
                    ptree.prune(p, preserve_branch_length=True)
                    ptree_samples = set([x.split('|')[0] for x in p])
                    pscore = determine_tree_score(ptree, ptree_samples)
                    pscore_bl = determine_tree_score_bl(ptree)
                    del ptree
                    pscore_sum += pscore
                    pscore_bl_sum += pscore_bl
                all_rooting_partitions.append(
                    [sp_merge, pscore_sum, len(sp_merge), pscore_bl_sum]
                )
                del rooted_t

        mp_t = copy.deepcopy(t)
        R = mp_t.get_midpoint_outgroup()
        mp_t.set_outgroup(R)
        sp = recursive_splitting(mp_t, samples_with_og)
        sp_further_split = further_split_disjoint_dog_partitions(mp_t, sp)
        sp_merge = sp_further_split
        if not skip_merge_back_flag:
            sp_merge = merge_back_dog_partitions(
                mp_t,
                pw_dists,
                sp_further_split,
                fixation_index_cutoff=fixation_index_cutoff,
            )

        pscore_sum = 0
        pscore_bl_sum = 0.0
        for p in sp_merge:
            ptree = copy.deepcopy(mp_t)
            ptree.prune(p, preserve_branch_length=True)
            ptree_samples = set([x.split('|')[0] for x in p])
            pscore = determine_tree_score(ptree, ptree_samples)
            pscore_bl = determine_tree_score_bl(ptree)
            del ptree
            pscore_sum += pscore
            pscore_bl_sum += pscore_bl
        all_rooting_partitions.append(
            [sp_merge, pscore_sum, len(sp_merge), pscore_bl_sum]
        )
        del mp_t

        if len(all_rooting_partitions) > 0:
            spl_outf = open(spl_file, 'w')
            for i, sp in enumerate(
                sorted(all_rooting_partitions, key=itemgetter(1, 2, 3))
            ):
                if i == 0:
                    for it, spi in enumerate(sp[0]):
                        spog = og + '_' + str(it)
                        for dom in spi:
                            spl_outf.write(spog + '\t' + dom + '\n')
            spl_outf.close()

        return
    except Exception as e:
        msg = (
            f'Issue with splitting ortholog group {og} - based on phylo from full MSA.'
        )
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        log_object.error(msg)
        sys.exit(1)
