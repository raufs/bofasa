"""
Core analysis functions for bofasa.

This module contains the main analysis functions for ortholog group determination,
phylogenetic refinement, and result processing.
"""

import os
import random
import subprocess
import sys
import traceback
import warnings
from collections import defaultdict
from operator import itemgetter
from typing import Any, Dict, List, Optional, Set, Tuple, Union
import concurrent.futures
import numpy as np
from Bio import SeqIO
with warnings.catch_warnings():
    # ete3 3.1.3 emits SyntaxWarnings for invalid escape sequences on Python 3.12+
    warnings.simplefilter('ignore', SyntaxWarning)
    from ete3 import Tree
from scipy.spatial import distance
from . import config
from .utils import get_version, setup_ready_directory, run_cmd, _iter_progress, create_pool
from .alignment import create_domain_alignments
import itertools

# Set random seed for reproducibility
random.seed(12345)

single_copy_dogs: Set[str] = set([])
largely_idr_dogs: Set[str] = set([])
protein_dogs: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
dog_conservation: Dict[str, Any] = {}

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
        if (i) < 10:
            pid = "OG00000" + str(i)
        elif (i) < 100:
            pid = "OG0000" + str(i)
        elif (i) < 1000:
            pid = "OG000" + str(i)
        elif (i) < 10000:
            pid = "OG00" + str(i)
        elif (i) < 100000:
            pid = "OG0" + str(i)
        else:
            pid = "OG" + str(i)
        assert pid != None
        return pid
    except Exception:
        msg = "Issue generating OG identifier for number: %s" % str(i)
        sys.stderr.write(msg + "\n")
        sys.stderr.write(traceback.format_exc() + "\n")


def run_set_of_protein_comparisons(inputs: Tuple[str, str, float, Dict[str, Dict[str, int]], Set[str], Set[str], Dict[str, float]]) -> str:
    """
    Run pairwise protein comparisons based on domain ortholog groups.
    
    Args:
        inputs: Tuple containing (input_listing_file, result_file, dj, protein_dogs, single_copy_dogs, 
        largely_idr_dogs, dog_conservation)
    
    Returns:
        str: Path to the result file
    """
    try:
        input_listing_file, result_file, dj = inputs

        with open(input_listing_file) as oilf, open(result_file, "w") as outf_handle:
            for line in oilf:
                line = line.strip()
                p1, p2 = line.split('\t')

                p1dogs = {}
                for dog, count in protein_dogs[p1].items():
                    if count != 0:
                        p1dogs[dog] = count

                p2dogs = {}
                for dog, count in protein_dogs[p2].items():
                    if count != 0:
                        p2dogs[dog] = count

                p1dogs_keys = set(p1dogs.keys())
                p2dogs_keys = set(p2dogs.keys())
                
                intersect_dogs = p1dogs_keys.intersection(p2dogs_keys)
                union_dogs = p1dogs_keys.union(p2dogs_keys)
                sc_dogs = single_copy_dogs.intersection(intersect_dogs)

                threshold = dj
                if len(sc_dogs) >= 1:
                    for sd in sc_dogs:
                        sd_conservation = dog_conservation[sd]
                        updated_threshold = dj - (dj*sd_conservation)
                        if updated_threshold < threshold:
                            threshold = updated_threshold

                union_count = 0
                intersect_count = 0
                for d in union_dogs:
                    if d in largely_idr_dogs: continue
                    p1dc = 0
                    p2dc = 0
                    if d in p1dogs:
                        p1dc = p1dogs[d]
                    if d in p2dogs:
                        p2dc = p2dogs[d]

                    union_count += p1dc + p2dc - min([p1dc, p2dc])
                    intersect_count += min([p1dc, p2dc])

                if union_count > 0:
                    jaccard_index = intersect_count/float(union_count)
                    if jaccard_index >= threshold:
                        outf_handle.write(p1 + '\t' + p2 + '\t' + str(jaccard_index*100.0) + '\n')

        return result_file
    except Exception as e:
        msg = (
            'Issue performing pairwise assessment between proteins based on DOGs '
            'to determine protein-resolution ortholog groups.'
            f'Error: {e}'
        )
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        raise


def write_unique_protein_pairs(
    dog_member_lists: List[List[str]], batch_handles: Dict[int, Any]
) -> Tuple[int, int]:
    """
    Write every protein pair that shares at least one DOG exactly once.

    Pairs are enumerated DOG by DOG (in input order) as (p1, p2) with p1 < p2
    within each DOG's sorted member list, and the n-th enumerated pair (counting
    repeats) belongs to batch n % len(batch_handles). Previously every
    occurrence was written; now a pair is written only at the occurrence that
    came first in the concatenation of the batch files (lowest batch, then
    earliest within it). The concatenated batch output is therefore exactly the
    previous listing with later repeats removed - same lines, same order - so
    downstream MCL sees the same graph with the same node order.

    Occurrences of a pair are found from the DOGs its two proteins share (each
    protein's short list of DOG indices), so no set of all pairs is kept.

    Args:
        dog_member_lists: For each DOG (in input order), its sorted protein ids.
        batch_handles: Batch index -> open file handle.

    Returns:
        (pairs enumerated including repeats, unique pairs written)
    """
    num_batches = len(batch_handles)
    prot_dog_idx: Dict[str, List[int]] = defaultdict(list)
    for i, members in enumerate(dog_member_lists):
        for p in members:
            prot_dog_idx[p].append(i)
    prot_dog_set = {p: set(v) for p, v in prot_dog_idx.items()}
    sizes = [len(m) for m in dog_member_lists]
    base = []
    total = 0
    for n in sizes:
        base.append(total)
        total += n * (n - 1) // 2
    positions: Dict[int, Dict[str, int]] = {}

    def occurrence_count(j: int, p1: str, p2: str) -> int:
        # Index of pair (p1, p2) in the enumeration, when it comes from DOG j.
        pos = positions.get(j)
        if pos is None:
            pos = {p: k for k, p in enumerate(dog_member_lists[j])}
            positions[j] = pos
        n, x, y = sizes[j], pos[p1], pos[p2]
        return base[j] + x * n - x * (x + 1) // 2 + (y - x - 1)

    pair_count = 0
    written = 0
    for i, members in enumerate(dog_member_lists):
        n = len(members)
        for a in range(n):
            p1 = members[a]
            p1_dogs = prot_dog_idx[p1]
            for b in range(a + 1, n):
                p2 = members[b]
                write = True
                if len(p1_dogs) > 1:
                    s2 = prot_dog_set[p2]
                    shared = [j for j in p1_dogs if j in s2]
                    if len(shared) > 1:
                        best = min(
                            ((c % num_batches, c) for c in
                             (occurrence_count(j, p1, p2) for j in shared))
                        )
                        write = best[1] == pair_count
                if write:
                    batch_handles[pair_count % num_batches].write(p1 + '\t' + p2 + '\n')
                    written += 1
                pair_count += 1
    return pair_count, written


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
        input_dir = protein_clustering_dir + 'Comparison_Listings/'
        pairwise_dir = protein_clustering_dir + 'Protein_Pairs_Based_on_DOGs/'
        p_clust_list_dir = protein_clustering_dir + 'Protein_Coarse_Cluster_Listings/'
        p_dist_dir = protein_clustering_dir + 'FASTME_Inputs/'
        p_tre_dir = protein_clustering_dir + 'Protein_DOG_Distance_Trees/'
        p_split_dir = protein_clustering_dir + 'Protein_Splitting/'

        from .utils import setup_ready_directory
        setup_ready_directory([input_dir, pairwise_dir, p_clust_list_dir, p_dist_dir, p_tre_dir, p_split_dir])
        
        pairwise_file = protein_clustering_dir + 'Protein_Pairs_Based_on_DOGs.txt'
        clusters_file = protein_clustering_dir + 'Protein_Clusters_Based_on_DOGs.txt'
        outf_handle = open(ogs_file, 'w')

        # Initialize local variables instead of using globals
        global protein_dogs, single_copy_dogs, largely_idr_dogs, dog_conservation
                
        num_batches = threads
        batch_handles = {}
        for batch in range(0, num_batches):
            batch_input_file = input_dir + str(batch) + '.txt'
            batch_handles[batch] = open(batch_input_file, 'w')

        all_proteins = set([])
        sample_count = None
        samples = []
        dog_member_lists: List[List[str]] = []

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
                            if prot_id not in protein_dogs:
                                protein_dogs[prot_id] = {}
                            if dog not in protein_dogs[prot_id]:
                                protein_dogs[prot_id][dog] = 0
                            protein_dogs[prot_id][dog] += 1

                            dog_lts.add(prot_id)
                            tot += 1
                            if lt.split('|')[2] == 'inter-domain_region':
                                idr += 1
                                
                    idr_prop = idr/float(tot)
                    if idr_prop >= 0.8:
                        largely_idr_dogs.add(dog)

                    if sample_with == 1:
                        sc_flag = False
                    
                    dog_conservation[dog] = sample_with/float(sample_count)
                    if sc_flag:
                        single_copy_dogs.add(dog)
                    dog_member_lists.append(sorted(dog_lts))

        # Protein pairs sharing at least one DOG, each written once (proteins sharing
        # several DOGs used to be listed - and compared - once per shared DOG; the
        # repeats carry identical weights and MCL keeps a single edge, so dropping
        # them leaves the MCL graph unchanged).
        pair_count, unique_pair_count = write_unique_protein_pairs(dog_member_lists, batch_handles)
        log_object.info(
            f"{unique_pair_count} unique protein pairs sharing DOGs to compare "
            f"({pair_count - unique_pair_count} repeated pairs skipped)"
        )

        pairwise_assessment_inputs = []
        for batch in range(0, threads):
            batch_handle = batch_handles[batch]
            batch_handle.close()

            input_listing_file = input_dir + str(batch) + '.txt'
            result_file = pairwise_dir + str(batch) + '.txt'
            pairwise_assessment_inputs.append([input_listing_file, result_file, dj])
        
        # Configure threading for protein comparisons
        try:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=num_batches
            ) as executor:
                results = {}
                future_to_index = {
                    executor.submit(run_set_of_protein_comparisons, input_data): idx
                    for idx, input_data in enumerate(pairwise_assessment_inputs)
                }
                
                for future in _iter_progress(
                    concurrent.futures.as_completed(future_to_index),
                    total=num_batches,
                    description="Protein comparisons",
                ):
                    try:
                        result = future.result()
                    except Exception as e:
                        result = None
                        msg = f'Failed during protein comparison job: {e}'
                        log_object.error(msg)
                        sys.stderr.write(msg + '\n')
                    if result is None:
                        msg = 'Protein comparison job failed for batch %d; exiting.' % future_to_index[future]
                        log_object.error(msg)
                        sys.stderr.write(msg + '\n')
                        sys.exit(1)
                    results[future_to_index[future]] = result
                
                # Write results in batch (input) order for deterministic MCL input
                with open(pairwise_file, 'w') as out_f:
                    for idx in sorted(results):
                        with open(results[idx], 'r') as f:
                            out_f.write(f.read())
        except Exception as e:
            msg = 'Failed during threading protein comparisons'
            log_object.error(msg)
            sys.stderr.write(msg + '\n')
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)

        clust_cmd = ['mcl', pairwise_file, '--abc', '-I', '1.2', '-o', clusters_file, '-te', str(threads)]
            
        try:
            subprocess.call(' '.join(clust_cmd), shell=True, stdout=subprocess.DEVNULL, 
                            stderr=subprocess.DEVNULL, executable='/bin/bash')
            assert (os.path.isfile(clusters_file))
            log_object.info('Successfully ran: %s' % ' '.join(clust_cmd))
        except Exception as e:
            msg = 'Had an issue running concatenation: %s' % ' '.join(clust_cmd)
            log_object.error(msg)
            sys.stderr.write(msg + '\n')
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)
        
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
                    # Regular groups run in worker processes, which do not share this
                    # process' protein_dogs, so each job carries its proteins' DOG counts.
                    split_nj_trees_input.append([og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, 
                                                og_split_file, log_object, 1,
                                                {p: protein_dogs[p] for p in ls}])
                large_protein_og_clusters.append(ls)
            else:
                protein_og_clusters.append(ls)

        # Run splitting of large protein coarse ortholog groups one at a time
        try:
            for input_data in _iter_progress(
                large_split_nj_trees_input,
                total=len(large_split_nj_trees_input),
                description="Large OG splits",
            ):
                try:
                    split_njt(input_data)
                except Exception as e:
                    msg = f'Failed during large protein ortholog group processing: {e}'
                    log_object.error(msg)
                    sys.stderr.write(msg + '\n')
                    sys.stderr.write(traceback.format_exc() + '\n')
                    # Continue with other jobs instead of exiting
        except Exception as e:
            msg = 'Failed during large protein ortholog group processing'
            log_object.error(msg)
            sys.stderr.write(msg + '\n')
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)

        # Run splitting of regular-sized coarse protein ortholog groups in parallel
        # (processes rather than threads: the splitting is pure Python and GIL-bound),
        # largest groups first.
        split_nj_trees_input.sort(key=lambda x: len(x[7]), reverse=True)
        if split_nj_trees_input:
            p = create_pool(threads)
            try:
                for og_uniq_id, err in _iter_progress(
                    p.imap_unordered(_split_njt_job, split_nj_trees_input),
                    total=len(split_nj_trees_input),
                    description="Regular OG splits",
                ):
                    if err is not None:
                        msg = f'Failed during regular protein ortholog group job {og_uniq_id}: {err}'
                        log_object.error(msg)
                        sys.stderr.write(msg + '\n')
                        # Continue with other jobs instead of exiting
            except Exception as e:
                p.terminate()
                msg = 'Failed during parallel processing of regular protein ortholog groups'
                log_object.error(msg)
                sys.stderr.write(msg + '\n')
                sys.stderr.write(traceback.format_exc() + '\n')
                sys.exit(1)
            else:
                p.close()
            finally:
                p.join()

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

        for pog_full in large_protein_og_clusters:
            missing = []
            for p in pog_full:
                if not p in accounted_for_in_splitting:
                    missing.append(p)
            if len(missing) > 0:
                protein_og_clusters.append(missing)

        for prot in all_proteins:
            if not prot in paired_proteins:
                protein_og_clusters.append([prot])

        for i, c in enumerate(sorted(protein_og_clusters)):
            samp_lts = defaultdict(list)
            for p in c: 
                s = p.split('|')[0]
                samp_lts[s].append(p.split('|')[1])
            og_id = generate_og_name(i+1)
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
    more_deterministic: bool = False,
    aligner: str = config.DEFAULT_ALIGNER,
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
        more_deterministic: Whether to use more deterministic settings
        aligner: Multiple sequence aligner to use ("famsa" or "muscle")

    Returns:
        None: Writes results to result_file
    """
    try:
        # Create output directories
        msa_dir = resdog_dir + 'Coarse_DOG_MSAs/'
        trim_dir = resdog_dir + 'Coarse_DOG_Trimmed_MSAs/'
        tre_dir = resdog_dir + 'Coarse_DOG_Phylogenies/'
        spl_full_dir = resdog_dir + 'Split_DOG_Listings/'
        setup_ready_directory([msa_dir, trim_dir, tre_dir, spl_full_dir], overwrite_mode="overwrite")

        # Create protein alignments
        aligner_name = 'MUSCLE super5' if aligner == 'muscle' else 'FAMSA'
        msg = f'Running multiple-sequence alignments using {aligner_name}.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)

        create_domain_alignments(
            dog_seqs_dir=orthofinder_fasta_dir,
            dog_algn_dir=msa_dir,
            log_object=log_object,
            threads=threads,
            more_deterministic=more_deterministic,
            aligner=aligner,
        )

        # Prepare TrimAl commands for trimming
        trimal_cmds = []
        for f in os.listdir(orthofinder_fasta_dir):
            if not f.endswith('.fa'):
                continue
            dog = '.'.join(f.split('.')[:-1])
            msa_file = msa_dir + dog + '.msa.faa'
            trim_file = trim_dir + dog + '.msa.trimmed.faa'
            
            # Check if MSA file exists and is not empty
            if not os.path.isfile(msa_file) or os.path.getsize(msa_file) == 0:
                continue
                
            # Split trimal_options into individual arguments
            trimal_args = trimal_options.split()
            trimal_cmd = ['trimal', '-in', msa_file, '-out', trim_file] + trimal_args
            trimal_cmds.append(trimal_cmd)

        # Run TrimAl
        msg = 'Running alignment trimming using trimal.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)
        
        if not trimal_cmds:
            log_object.warning("No TrimAl commands to run - no valid MSA files found")
        else:
            log_object.info(f"Running {len(trimal_cmds)} TrimAl commands with {threads} processes")
        
            # Run TrimAl commands using threading with reduced logging
            successful_trimal = 0
            failed_trimal = 0
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                    future_to_cmd = {
                        executor.submit(run_cmd, cmd, None): cmd for cmd in trimal_cmds
                    }
                    for future in _iter_progress(
                        concurrent.futures.as_completed(future_to_cmd),
                        total=len(trimal_cmds),
                        description="Running TrimAl",
                    ):
                        cmd = future_to_cmd[future]
                        try:
                            future.result()
                            successful_trimal += 1
                        except Exception:
                            failed_trimal += 1
                            cmd_str = ' '.join(cmd)
                            log_object.error(f"TrimAl command failed: {cmd_str}")
                            continue
            finally:
                # Log summary of TrimAl execution
                log_object.info(
                    f"TrimAl summary: {successful_trimal} successful, "
                    f"{failed_trimal} failed out of {len(trimal_cmds)} total commands"
                )

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
            fasttree_cmds.append(tre_cmd)

        # Run FastTree
        msg = 'Running phylogeny constructions using FastTree 2.'
        log_object.info(msg)
        
        if not fasttree_cmds:
            log_object.warning("No FastTree commands to run - no valid alignment files found")
        else:
            log_object.info(f"Running {len(fasttree_cmds)} FastTree commands")
        
            # Run FastTree commands using threading with reduced logging
            successful_fasttree = 0
            failed_fasttree = 0
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                    future_to_cmd = {
                        executor.submit(run_cmd, cmd, None): cmd for cmd in fasttree_cmds
                    }
                    for future in _iter_progress(
                        concurrent.futures.as_completed(future_to_cmd),
                        total=len(fasttree_cmds),
                        description="Running FastTree",
                    ):
                        cmd = future_to_cmd[future]
                        try:
                            future.result()
                            successful_fasttree += 1
                        except Exception:
                            failed_fasttree += 1
                            cmd_str = ' '.join(cmd)
                            log_object.error(f"FastTree command failed: {cmd_str}")
                            continue
            finally:
                # Log summary of FastTree execution
                log_object.info(
                    f"FastTree summary: {successful_fasttree} successful, "
                    f"{failed_fasttree} failed out of {len(fasttree_cmds)} total commands"
                )

        # Read samples from OrthoFinder TSV file
        samples = []
        with open(orthofinder_tsv_file) as otf:
            for i, line in enumerate(otf):
                line = line.strip('\n')
                ls = line.split('\t')
                if i == 0:
                    samples = [x[:-len('.ccds')] if x.endswith('.ccds') else x for x in ls[1:]]
                    break

        # Run phylogenetic splitting of orthogroups
        msg = 'Using phylogenetics to split coarse domain resolution ortholog groups.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)

        try:
            split_dogs_in_parallel(
                tre_dir, spl_full_dir, log_object,
                skip_merge_back_flag=skip_merge_back_flag,
                rooting_seeds=rooting_seeds,
                fixation_index_cutoff=fixation_index_cutoff,
                threads=threads,
            )
        except Exception as e:
            msg = f"Phylogenetic splitting of a domain ortholog group failed: {e}"
            log_object.error(msg)
            log_object.error(traceback.format_exc())
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        # Write final results file
        with open(result_file, 'w') as outf_handle:
            # Write header
            outf_handle.write('OG/Sample\t' + '\t'.join(samples) + '\n')
            
            # Process main orthogroups
            with open(orthofinder_tsv_file) as otf:
                for i, line in enumerate(otf):
                    if i == 0:  # Skip header
                        continue
                    line = line.strip('\n')
                    ls = line.split('\t')
                    dog = ls[0]
                    dog_spl_file = spl_full_dir + dog + '.txt'
                    
                    if os.path.isfile(dog_spl_file):
                        # Read split clusters
                        cluster_prots = defaultdict(set)
                        with open(dog_spl_file) as odsf:
                            for split_line in odsf:
                                split_line = split_line.strip()
                                cluster_id, prot = split_line.split('\t')
                                cluster_prots[cluster_id].add(prot)
                        
                        # Write each cluster as separate orthogroup
                        for cid in cluster_prots:
                            printlist = [cid]
                            for lts in ls[1:]:
                                lts_retained = []
                                for lt in lts.split(','):
                                    lt = lt.strip()
                                    if lt in cluster_prots[cid]:
                                        lts_retained.append(lt)
                                printlist.append(', '.join(lts_retained))
                            outf_handle.write('\t'.join(printlist) + '\n')
                    else:
                        # Write original orthogroup if no splitting occurred
                        outf_handle.write(line + '\n')
            
            # Add singletons
            with open(orthofinder_tsv_singletons_file) as ootsf:
                for i, line in enumerate(ootsf):
                    if i == 0:  # Skip header
                        continue
                    line = line.strip('\n')
                    outf_handle.write(line + '\n')

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
                    # Strip trailing '.ccds' from OrthoFinder sample names (as in resolve_orthogroups_using_phylogenetics)
                    samples = [x[:-len('.ccds')] if x.endswith('.ccds') else x for x in ls[1:]]
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


def _leaf_keys(name: str, are_proteins: bool) -> Tuple[str, str]:
    """Return the (sample, protein) keys used for duplication scoring of a leaf."""
    sample = name.split("|")[0]
    if are_proteins:
        return sample, name
    return sample, "|".join(name.split("|")[:-2])


def _subtree_duplication_scores(
    intree: Tree, all_og_samples: Set[str], are_proteins: bool = False
) -> Tuple[Dict[Any, int], Dict[Any, int]]:
    """
    Compute determine_tree_score() for every node of intree in one postorder pass.

    determine_tree_score(n) = sum over s in all_og_samples of |#proteins_s(n) - 1|
                            = |all_og_samples| + sum over s in (samples(n) & all_og_samples)
                              of (#proteins_s(n) - 2)
    Per-node sample->protein sets are merged small-to-large, so the whole pass is
    O(n log n) instead of re-traversing every subtree (O(n^2)).

    Returns:
        (scores, named_leaf_counts) keyed by node object.
    """
    n_all = len(all_og_samples)
    scores: Dict[Any, int] = {}
    named_counts: Dict[Any, int] = {}
    # node -> [sample -> set(protein), total_proteins, acc]
    state: Dict[Any, List[Any]] = {}
    for n in intree.traverse("postorder"):
        if n.is_leaf():
            s, p = _leaf_keys(n.name, are_proteins)
            state[n] = [{s: {p}}, 1, (-1 if s in all_og_samples else 0)]
            named_counts[n] = 1 if n.name.strip() != "" else 0
            scores[n] = n_all + state[n][2]
            continue
        acc = None
        named = 0
        for ch in n.children:
            named += named_counts[ch]
            st = state.pop(ch)
            if acc is None:
                acc = st
                continue
            if st[1] > acc[1]:
                acc, st = st, acc
            big = acc[0]
            for s, prots in st[0].items():
                in_all = s in all_og_samples
                bs = big.get(s)
                if bs is None:
                    big[s] = prots
                    acc[1] += len(prots)
                    if in_all:
                        acc[2] += len(prots) - 2
                else:
                    for p in prots:
                        if p not in bs:
                            bs.add(p)
                            acc[1] += 1
                            if in_all:
                                acc[2] += 1
        state[n] = acc
        named_counts[n] = named
        scores[n] = n_all + acc[2]
    return scores, named_counts


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

        scores, named_counts = _subtree_duplication_scores(
            intree, all_og_samples, are_proteins=are_proteins
        )
        full_score = scores[intree]

        # Preorder pass (same order as ete3's traverse('preorder')) to find the
        # minimum internal-node score and the sample set.
        min_score = 1e100
        samples = set([])
        for n in intree.traverse("preorder"):
            if n.is_leaf():
                samples.add(n.name.split("|")[0])
            elif scores[n] < min_score:
                min_score = scores[n]

        if min_score >= full_score or len(samples) == 1:
            return [all_leaves]

        # Greedily take min-score subtrees in preorder whose leaves do not overlap
        # an already-taken subtree. Non-empty subtree leaf sets only overlap when
        # one node is an ancestor of the other, and ancestors come first in
        # preorder, so a node is taken iff no ancestor was taken (an empty leaf
        # set never overlaps and is always taken, as before).
        best_nodes = []
        stack = [(intree, False)]
        while stack:
            n, under_taken = stack.pop()
            if n.is_leaf():
                continue
            if scores[n] == min_score and (not under_taken or named_counts[n] == 0):
                best_nodes.append(n)
                if named_counts[n] > 0:
                    under_taken = True
            for ch in reversed(n.children):
                stack.append((ch, under_taken))

        accounted_for = set([])
        best_partitions = []
        for n in best_nodes:
            leaves = get_children(n)
            best_partitions.append(leaves)
            accounted_for |= leaves

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


def _protein_dog_count_matrix(
    proteins: List[str], prot_dogs: Dict[str, Dict[str, int]], nonzero_only: bool
) -> np.ndarray:
    """
    Protein x DOG count matrix (rows in the order of ``proteins``; columns are the
    sorted DOGs with a count > 0 (``nonzero_only=False``) or != 0
    (``nonzero_only=True``) in at least one protein).
    """
    dogs = set()
    for p in proteins:
        for d, c in prot_dogs.get(p, {}).items():
            if (c != 0) if nonzero_only else (c > 0):
                dogs.add(d)
    dog_index = {d: j for j, d in enumerate(sorted(dogs))}
    m = np.zeros((len(proteins), len(dog_index)), dtype=np.int64)
    for i, p in enumerate(proteins):
        for d, c in prot_dogs.get(p, {}).items():
            j = dog_index.get(d)
            if j is not None and (c != 0 or not nonzero_only):
                m[i, j] = c
    return m


def write_protein_cosine_distance_matrix(
    all_prots: Set[str],
    prot_dogs: Dict[str, Dict[str, int]],
    cog_dist_file: str,
    block_size: int = 1024,
) -> Dict[str, str]:
    """
    Write the PHYLIP matrix of cosine distances between proteins' DOG count vectors
    (FastME input) and return the PHYLIP-name -> protein mapping.

    Matches scipy.spatial.distance.cosine applied pair by pair: for integer count
    vectors the dot products are exact, and 1 - uv / sqrt(uu * vv) is evaluated
    with the same float operations, then clipped to [0, 2] as scipy does. Computed
    with matrix products in row blocks instead of one Python-level call per pair.
    """
    prots = sorted(all_prots)
    # float64 BLAS products are exact here: counts and their sums of products are
    # small integers (< 2**53), so uv, uu, vv and uu * vv equal scipy's int64 values.
    m = _protein_dog_count_matrix(prots, prot_dogs, nonzero_only=False).astype(np.float64)
    norms = np.einsum('ij,ij->i', m, m)
    naming: Dict[str, str] = {}
    with open(cog_dist_file, 'w') as outf:
        outf.write(str(len(prots)) + '\n')
        for start in range(0, len(prots), block_size):
            stop = min(start + block_size, len(prots))
            uv = m[start:stop] @ m.T
            with np.errstate(divide='ignore', invalid='ignore'):
                denom = np.sqrt(norms[start:stop, None] * norms[None, :])
                dist = 1.0 - uv / denom
            np.clip(dist, 0.0, 2.0, out=dist)
            for bi, row in enumerate(dist.tolist()):
                p1i = start + bi
                len_id = len(str(p1i))
                assert(10 > len_id)
                naming[str(p1i)] = prots[p1i]
                name = str(p1i) + (' ' * (10 - len_id))
                outf.write(name + ' '.join(['{:.10f}'.format(x) for x in row]) + '\n')
    return naming


def split_njt(input: List[Any]) -> None:
    """
    Split protein ortholog groups using neighbor-joining trees based on domain ortholog group distances.

    Args:
        input: List containing [og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, og_split_file,
        log_object, threads] and optionally an 8th element: the protein -> {DOG: count}
        mapping for (at least) the proteins of this group. Without it the module-level
        ``protein_dogs`` is used (only available in the process that loaded it).
    """
    try:
        og, cog_list_file, cog_dist_file, tre_file, spl_file, log_object, threads = input[:7]
        prot_dogs = input[7] if len(input) > 7 else protein_dogs

        all_prots = set([])
        with open(cog_list_file) as oclf:
            for line in oclf:
                line = line.strip()
                all_prots.add(line)

        naming = write_protein_cosine_distance_matrix(all_prots, prot_dogs, cog_dist_file)

        fastme_cmd = ['fastme', '-i', cog_dist_file, '-o', tre_file, '-T', str(threads)]
        run_cmd(fastme_cmd, None, check_files=[tre_file])

        t = Tree(tre_file)

        samples_with_og = set([])
        leafs = set([])
        for n in t.traverse('postorder'):
            if n.is_leaf():
                n.name = naming[n.name]
                s = n.name.split('|')[0]
                samples_with_og.add(s)
                leafs.add(n.name)

        if len(leafs) < 3:
            return
        if len(samples_with_og) == 1:
            return

        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
        sp = recursive_splitting(t, samples_with_og, are_proteins=True)
        sp_further_split = further_split_disjoint_dog_partitions(t, sp)
        sp_refined = further_split_outlier_artifact_groups(t, sp_further_split, prot_dogs=prot_dogs)
        spl_outf = open(spl_file, 'w')
        for spi in sp_refined:
            spl_outf.write(' '.join(sorted(spi)) + '\n')
        spl_outf.close()
        return
    except Exception as e:
        msg = (
            f'Issue with splitting protein ortholog group {og} - based on domain '
            'cosine distances neighbor-joining tree.'
        )
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        log_object.error(msg)
        sys.exit(1)


def _split_njt_job(input: List[Any]) -> Tuple[str, Any]:
    """
    Process-pool wrapper around split_njt: returns (og, None) on success or
    (og, error message) on failure, so one failing group neither kills the
    worker (split_njt raises SystemExit) nor stops the remaining groups.
    """
    try:
        split_njt(input)
        return input[0], None
    except (Exception, SystemExit):
        return input[0], traceback.format_exc()


def _node_leaf_set(n: Tree) -> Set[str]:
    """Leaf set of a node as used by further_split_disjoint_dog_partitions."""
    if n.is_leaf():
        return set([n.name])
    return get_children(n)


def further_split_disjoint_dog_partitions(rooted_t: Tree, sps: List[Set[str]]) -> List[Set[str]]:
    """
    Function to further split disjoint domain ortholog groups based on the phylogenetic tree.

    A partition that is exactly the leaf set of some node is kept as is; otherwise
    it is covered greedily by the largest subtrees fully contained in it.

    Node leaf sets are summarised by (size, hash-sum) so monophyly is checked with
    a dict lookup (verified exactly on the candidate) instead of comparing the
    partition against the materialised leaf set of every node.
    """
    try:
        mask = (1 << 64) - 1
        post = list(rooted_t.traverse('postorder'))
        size: Dict[Any, int] = {}
        hsum: Dict[Any, int] = {}
        index: Dict[Tuple[int, int], List[Any]] = defaultdict(list)
        for n in post:
            if n.is_leaf():
                size[n] = 1
                hsum[n] = hash(n.name) & mask
            else:
                s = 0
                h = 0
                for ch in n.children:
                    if ch.is_leaf():
                        if ch.name.strip() != "":
                            s += 1
                            h += hash(ch.name) & mask
                    else:
                        s += size[ch]
                        h += hsum[ch]
                size[n] = s
                hsum[n] = h & mask
            index[(size[n], hsum[n])].append(n)

        updated_sp = []
        for sp in sps:
            key = (len(sp), sum(hash(x) & mask for x in sp) & mask)
            monophyletic_flag = False
            for n in index.get(key, []):
                if _node_leaf_set(n) == sp:
                    monophyletic_flag = True
                    break
            if monophyletic_flag:
                updated_sp.append(sp)
                continue

            # Nodes whose whole leaf set lies inside sp.
            inside: Dict[Any, bool] = {}
            for n in post:
                if n.is_leaf():
                    inside[n] = n.name in sp
                else:
                    ok = True
                    for ch in n.children:
                        if ch.is_leaf():
                            if ch.name.strip() != "" and ch.name not in sp:
                                ok = False
                                break
                        elif not inside[ch]:
                            ok = False
                            break
                    inside[n] = ok

            # A contained node below a strictly larger contained parent can never
            # be chosen (the parent, or a node with the parent's leaf set, is taken
            # first and overlaps it), so it is skipped. Postorder order is kept so
            # the stable size-descending sort matches the original tie-breaking.
            candidates = []
            for n in post:
                if not inside[n]:
                    continue
                par = n.up
                if size[n] == 0 or par is None or not inside[par] or size[par] == size[n]:
                    candidates.append(n)

            accounted_leafs = set([])
            for n in sorted(candidates, key=lambda x: size[x], reverse=True):
                nc = _node_leaf_set(n)
                if len(nc.intersection(accounted_leafs)) == 0:
                    accounted_leafs = accounted_leafs.union(nc)
                    updated_sp.append(nc)
            assert(sp == accounted_leafs)
        return updated_sp

    except Exception as e:
        msg = 'Issues further refining domain ortholog groups based on phylogenetics.'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def further_split_outlier_artifact_groups(
    rooted_t: Tree,
    sps: List[Set[str]],
    prot_dogs: Optional[Dict[str, Dict[str, int]]] = None,
    block_size: int = 512,
) -> List[Set[str]]:
    """
    Further split ortholog groups that may contain small proteins or be the result
    of outlier sequences being claded together after midpoint rooting of FastME tree.

    A protein becomes a singleton when its best weighted-Jaccard similarity of DOG
    counts to another member of its group is <= its best similarity to a protein
    of the tree outside the group. Weighted Jaccard is sum(min) / sum(max) over DOGs
    with non-zero counts; sum(min) over integer counts is computed exactly as
    sum_k [a >= k] . [b >= k] with matrix products (row blocks), and the final
    division matches the original per-pair Python arithmetic.

    Args:
        rooted_t: Rooted phylogenetic tree
        sps: List of ortholog group partitions
        prot_dogs: protein -> {DOG: count}; defaults to the module-level protein_dogs

    Returns:
        list: Updated ortholog group partitions
    """
    try:
        if prot_dogs is None:
            prot_dogs = protein_dogs
        all_proteins = get_children(rooted_t)
        universe = set(all_proteins)
        for sp in sps:
            universe.update(sp)
        universe = sorted(universe)
        index = {p: i for i, p in enumerate(universe)}

        m = _protein_dog_count_matrix(universe, prot_dogs, nonzero_only=True)
        max_count = int(m.max()) if m.size else 0
        min_count = int(m.min()) if m.size else 0
        if min_count < 0:
            raise ValueError("Negative DOG counts are not supported.")
        levels = [(m >= k).astype(np.float64) for k in range(1, max_count + 1)]
        totals = m.sum(axis=1)
        in_tree = np.zeros(len(universe), dtype=bool)
        for p in all_proteins:
            in_tree[index[p]] = True

        updated_sp = []
        for sp in sps:
            members = sorted(sp)
            rows = np.array([index[p] for p in members], dtype=np.int64)
            in_sp = np.zeros(len(universe), dtype=bool)
            in_sp[rows] = True
            external = in_tree & ~in_sp

            singletons = set([])
            for start in range(0, len(members), block_size):
                block = rows[start:start + block_size]
                inter = np.zeros((len(block), len(universe)), dtype=np.float64)
                for lv in levels:
                    inter += lv[block] @ lv.T
                union = totals[block][:, None] + totals[None, :] - inter
                with np.errstate(divide='ignore', invalid='ignore'):
                    jacc = np.where(union > 0, inter / union, 0.0)
                internal_j = np.where(in_sp[None, :], jacc, 0.0)
                internal_j[np.arange(len(block)), block] = 0.0
                max_internal = internal_j.max(axis=1) if internal_j.shape[1] else np.zeros(len(block))
                max_external = np.where(external[None, :], jacc, 0.0).max(axis=1)
                for bi, p1 in enumerate(members[start:start + block_size]):
                    if max_internal[bi] <= max_external[bi]:
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


def single_linkage_cluster(pairs, all_lts, paired_lts):
    """
    Solution for single-linkage clustering taken from mimomu's response in the stackoverflow page:
    https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
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

        return (L)
    except Exception:
        msg = 'Issue running single linkage clustering!'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        sys.exit(1)


def _exact_depths(t: Tree) -> Tuple[Dict[Any, int], int]:
    """
    Root-to-node depths as exact integers in units of 2**-scale.

    Every float branch length is a dyadic rational, so scaling by the largest
    power-of-two denominator makes all depths (and sums of them) exact integers.
    This keeps the aggregated distance sums in merge_back_dog_partitions exact,
    independent of summation order.
    """
    ratios = {}
    scale = 0
    for n in t.traverse('preorder'):
        if n.is_root():
            continue
        num, den = float(n.dist).as_integer_ratio()
        ratios[n] = (num, den)
        b = den.bit_length() - 1
        if b > scale:
            scale = b
    depth: Dict[Any, int] = {}
    for n in t.traverse('preorder'):
        if n.is_root():
            depth[n] = 0
        else:
            num, den = ratios[n]
            depth[n] = depth[n.up] + num * ((1 << scale) // den)
    return depth, scale


def merge_back_dog_partitions(
    t: Tree,
    split_partitions: List[Set[str]],
    fixation_index_cutoff: float = 0.5,
) -> List[Set[str]]:
    """
    Merge back domain ortholog group partitions based on fixation index.

    For each internal node spanning leaves from >= 2 partitions, the mean patristic
    distance between leaves of the same partition (pi_within) and of different
    partitions (pi_between) under that node give
    fix_index = (pi_between - pi_within) / pi_between; partitions under nodes with
    fix_index <= cutoff are merged by single-linkage.

    The distance sums are aggregated bottom-up using
    dist(a, b) = depth(a) + depth(b) - 2 * depth(LCA(a, b)), merging per-partition
    (count, depth-sum) tables small-to-large, so no pairwise distance matrix is
    built (previously O(n^2) memory and roughly O(n^3) time on large trees).

    Args:
        t: Rooted phylogenetic tree
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

    depth, scale = _exact_depths(t)

    # node -> [clade -> [count, depth_sum], n_leaves, depth_sum,
    #          all_pairs_dist_sum, within_pairs_dist_sum, within_pairs_count]
    state: Dict[Any, List[Any]] = {}
    merge_sps = []
    paired_sps = set([])
    for n in t.traverse('postorder'):
        if n.is_leaf():
            d = depth[n]
            state[n] = [{leaf_to_clade[n.name]: [1, d]}, 1, d, 0, 0, 0]
            continue

        two_d = 2 * depth[n]
        acc = None
        for ch in n.children:
            st = state.pop(ch)
            if acc is None:
                acc = st
                continue
            if len(st[0]) > len(acc[0]):
                acc, st = st, acc
            k1, s1 = acc[1], acc[2]
            k2, s2 = st[1], st[2]
            acc[3] += st[3] + k1 * s2 + k2 * s1 - two_d * k1 * k2
            acc[1] = k1 + k2
            acc[2] = s1 + s2
            acc[4] += st[4]
            acc[5] += st[5]
            big = acc[0]
            for c, (ck2, cs2) in st[0].items():
                e = big.get(c)
                if e is None:
                    big[c] = [ck2, cs2]
                else:
                    ck1, cs1 = e
                    acc[4] += ck1 * cs2 + ck2 * cs1 - two_d * ck1 * ck2
                    acc[5] += ck1 * ck2
                    e[0] = ck1 + ck2
                    e[1] = cs1 + cs2
        state[n] = acc

        clades = acc[0]
        if len(clades) == 1:
            continue
        within_count = acc[5]
        if within_count == 0:
            continue
        n_leaves = acc[1]
        between_count = n_leaves * (n_leaves - 1) // 2 - within_count
        within_sum = acc[4]
        between_sum = acc[3] - within_sum
        # Exact rational means, correctly rounded to float.
        pi_within = within_sum / (within_count << scale)
        pi_between = between_sum / (between_count << scale)
        if pi_between == 0:
            # Identical sequences (zero branch lengths): partitions are not differentiated.
            fix_index = 0.0
        else:
            fix_index = (pi_between - pi_within) / pi_between

        if fix_index > fixation_index_cutoff:
            continue
        sorted_clades = sorted(clades)
        for i, c1 in enumerate(sorted_clades):
            for c2 in sorted_clades[i + 1:]:
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


def score_partitions(rooted_t: Tree, partitions: List[Set[str]]) -> Tuple[int, float]:
    """
    Score a partitioning of rooted_t's leaves without copying/pruning the tree.

    Equivalent to, for each partition p, pruning a copy of rooted_t to p (with
    preserve_branch_length=True) and summing determine_tree_score() and
    determine_tree_score_bl() over the pruned trees:
      - duplication score of p = #distinct proteins in p - #distinct samples in p
      - branch-length score of p = total length of the subtree spanning p (for a
        single leaf, ete3 keeps the original root, so it is the root-to-leaf depth)

    Returns:
        (summed duplication score, summed branch-length score)
    """
    leaf_part: Dict[str, int] = {}
    for i, p in enumerate(partitions):
        for l in p:
            leaf_part[l] = i
    sizes = [len(p) for p in partitions]
    bl = [0.0] * len(partitions)

    float_depth: Dict[Any, float] = {}
    for n in rooted_t.traverse('preorder'):
        float_depth[n] = 0.0 if n.is_root() else float_depth[n.up] + n.dist

    # node -> {partition: number of its leaves below node}, only for partitions
    # not yet complete; an edge lies on p's spanning subtree iff 0 < count < |p|.
    state: Dict[Any, Dict[int, int]] = {}
    for n in rooted_t.traverse('postorder'):
        if n.is_leaf():
            cur: Dict[int, int] = {}
            pi = leaf_part.get(n.name)
            if pi is not None:
                if sizes[pi] == 1:
                    bl[pi] = float_depth[n]
                else:
                    cur[pi] = 1
        else:
            cur = None
            for ch in n.children:
                st = state.pop(ch)
                if cur is None:
                    cur = st
                    continue
                if len(st) > len(cur):
                    cur, st = st, cur
                for pi, c in st.items():
                    cur[pi] = cur.get(pi, 0) + c
            for pi in [pi for pi, c in cur.items() if c == sizes[pi]]:
                del cur[pi]
        if not n.is_root():
            for pi in cur:
                bl[pi] += n.dist
        state[n] = cur

    pscore_sum = 0
    pscore_bl_sum = 0.0
    for i, p in enumerate(partitions):
        prots = set()
        samples = set()
        for x in p:
            s, prot = _leaf_keys(x, False)
            samples.add(s)
            prots.add(prot)
        pscore_sum += len(prots) - len(samples)
        pscore_bl_sum += bl[i]
    return pscore_sum, pscore_bl_sum


def _load_named_dog_tree(tre_file: str) -> Tree:
    """Load a coarse DOG tree and name internal nodes node_<preorder index>."""
    t = Tree(tre_file)
    for node_id, n in enumerate(t.traverse('preorder')):
        if not n.is_leaf():
            n.name = 'node_' + str(node_id + 1)
    return t


def plan_dog_rootings(og: str, t: Tree, rooting_seeds: int) -> List[Any]:
    """
    Rootings to evaluate for a DOG tree, in the order they are scored.

    Returns the sampled outgroup node names (preorder order) followed by None,
    which stands for midpoint rooting. With rooting_seeds <= 1 this is [None].
    """
    rooting_seeds = max([rooting_seeds, 1])
    if rooting_seeds <= 1:
        return [None]
    possible_nodes_for_rooting = set([n.name for n in t.traverse('preorder')])
    if len(possible_nodes_for_rooting) > (rooting_seeds - 1):
        # Per-OG RNG over a sorted population so results do not depend on
        # PYTHONHASHSEED or on the order tasks are processed in worker processes.
        rng = random.Random(f"12345:{og}")
        possible_nodes_for_rooting = set(
            rng.sample(sorted(possible_nodes_for_rooting), rooting_seeds - 1)
        )
    plan: List[Any] = [
        n.name for n in t.traverse('preorder') if n.name in possible_nodes_for_rooting
    ]
    plan.append(None)
    return plan


def _plan_dog_rootings_job(inputs: List[Any]) -> Tuple[str, List[Any]]:
    """Pool worker: (og, tre_file, rooting_seeds) -> (og, rooting plan)."""
    og, tre_file, rooting_seeds = inputs
    return og, plan_dog_rootings(og, _load_named_dog_tree(tre_file), rooting_seeds)


def split_dog_rooting(inputs: List[Any]) -> Tuple[str, int, Any]:
    """
    Split one coarse DOG tree under one rooting.

    Args:
        inputs: [og, tre_file, rooting_index, root (outgroup node name, or None for
            midpoint rooting), skip_merge_back_flag, fixation_index_cutoff, log_object]

    Returns:
        (og, rooting_index, result) where result is
        [partitions, duplication_score, n_partitions, branch_length_score], or None
        if the tree could not be rooted at the requested node.
    """
    (og, tre_file, rooting_index, root, skip_merge_back_flag,
     fixation_index_cutoff, log_object) = inputs
    try:
        rooted_t = _load_named_dog_tree(tre_file)
        samples_with_og = set([])
        for n in rooted_t.get_leaves():
            samples_with_og.add(n.name.split('|')[0])

        if root is None:
            R = rooted_t.get_midpoint_outgroup()
            rooted_t.set_outgroup(R)
        elif root != rooted_t.name:
            try:
                rooted_t.set_outgroup(root)
            except Exception:
                log_object.error(f"Error setting outgroup for node {root}: {traceback.format_exc()}")
                return og, rooting_index, None

        sp = recursive_splitting(rooted_t, samples_with_og)
        sp_merge = further_split_disjoint_dog_partitions(rooted_t, sp)
        if not skip_merge_back_flag:
            sp_merge = merge_back_dog_partitions(
                rooted_t,
                sp_merge,
                fixation_index_cutoff=fixation_index_cutoff,
            )
        pscore_sum, pscore_bl_sum = score_partitions(rooted_t, sp_merge)
        return og, rooting_index, [sp_merge, pscore_sum, len(sp_merge), pscore_bl_sum]
    except (Exception, SystemExit) as e:
        # Helpers call sys.exit(1) on error; inside a multiprocessing pool worker a
        # SystemExit kills the worker and the parent's imap_unordered hangs forever.
        # Convert to a regular exception so it propagates to the parent.
        msg = (
            f'Issue with splitting ortholog group {og} - based on phylo from full MSA.'
        )
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        raise RuntimeError(msg + '\n' + traceback.format_exc()) from None


def write_best_dog_split(og: str, rooting_results: List[Any], spl_file: str) -> None:
    """
    Write the best-scoring partitioning across rootings for one DOG.

    rooting_results must be in rooting-plan order (sampled roots, then midpoint);
    ties on (duplication score, n partitions, branch length) keep that order.
    """
    all_rooting_partitions = [r for r in rooting_results if r is not None]
    if len(all_rooting_partitions) == 0:
        return
    best = sorted(all_rooting_partitions, key=itemgetter(1, 2, 3))[0]
    with open(spl_file, 'w') as spl_outf:
        for it, spi in enumerate(best[0]):
            spog = og + '_' + str(it)
            for dom in spi:
                spl_outf.write(spog + '\t' + dom + '\n')


def split_dogs(inputs: List[Any]) -> None:
    """
    Split domain ortholog groups using phylogenetic approach (serially).

    resolve_orthogroups_using_phylogenetics() runs the individual rootings of all
    trees in parallel via split_dog_rooting(); this helper does the same for a
    single tree in the current process.

    Args:
        inputs: List containing [og, tre_file, spl_file, skip_merge_back_flag, rooting_seeds,
        fixation_index_cutoff, threads, log_object]

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
    try:
        plan = plan_dog_rootings(og, _load_named_dog_tree(tre_file), rooting_seeds)
    except (Exception, SystemExit):
        msg = (
            f'Issue with splitting ortholog group {og} - based on phylo from full MSA.'
        )
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        raise RuntimeError(msg + '\n' + traceback.format_exc()) from None
    results = []
    for i, root in enumerate(plan):
        _, _, res = split_dog_rooting(
            [og, tre_file, i, root, skip_merge_back_flag, fixation_index_cutoff, log_object]
        )
        results.append(res)
    write_best_dog_split(og, results, spl_file)


def _count_newick_leaves(tre_file: str) -> int:
    """Cheap leaf count of a Newick file (commas + 1) used for job scheduling."""
    with open(tre_file) as f:
        return f.read().count(',') + 1


def split_dogs_in_parallel(
    tre_dir: str,
    spl_full_dir: str,
    log_object: Any,
    skip_merge_back_flag: bool = False,
    rooting_seeds: int = config.DEFAULT_ROOTING_SEEDS,
    fixation_index_cutoff: float = config.DEFAULT_FIXATION_INDEX_CUTOFF,
    threads: int = config.DEFAULT_THREADS,
) -> None:
    """
    Phylogenetically split every coarse DOG tree in tre_dir, in parallel.

    Each (tree, rooting) pair is an independent job: one job per tree by default
    (midpoint rooting), or rooting_seeds jobs per tree with -rs > 1. Jobs are
    dispatched largest tree first so big trees do not end up running alone at the
    end. Results are gathered per tree, and once all rootings of a tree are in,
    the best partitioning is chosen in rooting-plan order (so ties resolve exactly
    as in a serial run) and written to <spl_full_dir>/<og>.txt.
    """
    trees = []
    for f in sorted(os.listdir(tre_dir)):
        if not f.endswith('.tre'):
            continue
        dog = '.tre'.join(f.split('.tre')[:-1])
        tre_file = os.path.join(tre_dir, f)
        trees.append((dog, tre_file, _count_newick_leaves(tre_file)))
    if not trees:
        return
    # Largest trees first (longest-processing-time-first scheduling).
    trees.sort(key=lambda x: x[2], reverse=True)

    p = create_pool(threads)
    try:
        # Rooting plans: trivially [midpoint] unless extra rooting seeds are requested.
        if max([rooting_seeds, 1]) > 1:
            plans = dict(
                _iter_progress(
                    p.imap_unordered(
                        _plan_dog_rootings_job,
                        [[dog, tre_file, rooting_seeds] for dog, tre_file, _ in trees],
                    ),
                    total=len(trees),
                    description="Planning DOG rootings",
                )
            )
        else:
            plans = {dog: [None] for dog, _, _ in trees}

        jobs = []
        for dog, tre_file, _ in trees:
            for i, root in enumerate(plans[dog]):
                jobs.append([dog, tre_file, i, root, skip_merge_back_flag,
                             fixation_index_cutoff, log_object])
        log_object.info(
            f"Splitting {len(trees)} coarse DOG trees ({len(jobs)} rooting jobs) "
            f"with {threads} processes"
        )

        pending: Dict[str, Dict[int, Any]] = defaultdict(dict)
        for dog, idx, res in _iter_progress(
            p.imap_unordered(split_dog_rooting, jobs),
            total=len(jobs),
            description="Domain OG splits",
        ):
            pending[dog][idx] = res
            if len(pending[dog]) == len(plans[dog]):
                results = [pending[dog][i] for i in range(len(plans[dog]))]
                write_best_dog_split(dog, results, os.path.join(spl_full_dir, dog + '.txt'))
                del pending[dog]
    except BaseException:
        p.terminate()
        raise
    else:
        p.close()
    finally:
        p.join()
