"""
Core analysis functions for bofasa.

This module contains the main analysis functions for ortholog group determination,
phylogenetic refinement, and result processing.
"""

import copy
import multiprocessing
import os
import random
import statistics   
import subprocess
import sys
import traceback
from collections import defaultdict
from operator import itemgetter
from typing import Any, Dict, List, Set, Tuple, Union
import concurrent.futures
import numpy as np
from Bio import SeqIO
from ete3 import Tree
from scipy.spatial import distance
from . import config
from .utils import get_version, setup_ready_directory, run_cmd, _iter_progress
from .alignment import create_domain_alignments

# Set random seed for reproducibility
random.seed(12345)

single_copy_dogs: Set[str] = set([])
largely_idr_dogs: Set[str] = set([])
protein_dogs: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
dog_conservation: Dict[str, Any] = {}
tree_obj: Any = None

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
        outf_handle = open(result_file, "w")

        with open(input_listing_file) as oilf:
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

        outf_handle.close()
        return result_file
    except Exception as e:
        msg = (
            'Issue performing pairwise assessment between proteins based on DOGs '
            'to determine protein-resolution ortholog groups.'
            f'Error: {e}'
        )
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        return None


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
        pair_count = 0

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
                    for j, p1 in enumerate(sorted(dog_lts)):
                        for k, p2 in enumerate(sorted(dog_lts)):
                            if j >= k: continue
                            batch = pair_count % num_batches
                            batch_input_handle = batch_handles[batch]
                            batch_input_handle.write(p1 + '\t' + p2 + '\n')
                            pair_count += 1

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
                results = []
                future_to_input = {
                    executor.submit(run_set_of_protein_comparisons, input_data): input_data
                    for input_data in pairwise_assessment_inputs
                }
                
                for future in _iter_progress(
                    concurrent.futures.as_completed(future_to_input),
                    total=num_batches,
                    description="Protein comparisons",
                ):
                    try:
                        result = future.result()
                        if result is not None:
                            results.append(result)
                    except Exception as e:
                        input_data = future_to_input[future]
                        msg = f'Failed during protein comparison job: {e}'
                        log_object.error(msg)
                        sys.stderr.write(msg + '\n')
                        sys.stderr.write(traceback.format_exc() + '\n')
                        # Continue with other jobs instead of exiting
                
                # Write results in sorted order
                for result in results:
                    if result is not None:
                        with open(result, 'r') as f:
                            content = f.read()
                            with open(pairwise_file, 'a') as out_f:
                                out_f.write(content)
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
                    split_nj_trees_input.append([og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, 
                                                og_split_file, log_object, 1])
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

        # Run splitting of regular-sized caorse protein ortholog groups in parallel
        try:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=threads
            ) as executor:
                future_to_input = {
                    executor.submit(split_njt, input_data): input_data
                    for input_data in split_nj_trees_input
                }
                
                for future in _iter_progress(
                    concurrent.futures.as_completed(future_to_input),
                    total=len(split_nj_trees_input),
                    description="Regular OG splits",
                ):
                    try:
                        future.result()
                    except Exception as e:
                        input_data = future_to_input[future]
                        msg = f'Failed during regular protein ortholog group job: {e}'
                        log_object.error(msg)
                        sys.stderr.write(msg + '\n')
                        sys.stderr.write(traceback.format_exc() + '\n')
                        # Continue with other jobs instead of exiting
        except Exception as e:
            msg = 'Failed during threading regular protein ortholog groups'
            log_object.error(msg)
            sys.stderr.write(msg + '\n')
            sys.stderr.write(traceback.format_exc() + '\n')
            sys.exit(1)

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
    more_deterministic: bool = False,
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

        # Create protein alignments using MUSCLE super5
        msg = 'Running multiple-sequence alignments using MUSCLE super5.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)

        # Use MUSCLE super5 for alignments
        create_domain_alignments(
            dog_seqs_dir=orthofinder_fasta_dir,
            dog_algn_dir=msa_dir,
            log_object=log_object,
            threads=threads,
            more_deterministic=more_deterministic,
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
            return
            
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
            return
            
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
                    samples = ['.ccds'.join(x.split('.ccds')[:-1]) for x in ls[1:]]
                    break

        # Set up inputs for splitting DOGs using phylogenetics
        split_inputs = []
        for f in os.listdir(tre_dir):
            dog = '.tre'.join(f.split('.tre')[:-1])
            tre_file = tre_dir + f
            spl_full_file = spl_full_dir + dog + '.txt'
            t = Tree(tre_file)
            if len(t.get_leaves()) < 500:
                split_inputs.append([dog, tre_file, spl_full_file, skip_merge_back_flag, rooting_seeds, fixation_index_cutoff, 1, log_object])
            else:
                # Handle large trees with global tree object for efficiency
                global tree_obj
                tree_obj = Tree(tre_file)
                split_dogs([dog, tre_file, spl_full_file, skip_merge_back_flag, rooting_seeds, fixation_index_cutoff, threads, log_object])
                tree_obj = None

        # Run phylogenetic splitting of orthogroups
        msg = 'Using phylogenetics to split coarse domain resolution ortholog groups.'
        sys.stderr.write(msg + '\n')
        log_object.info(msg)
        
        # Configure threading for domain ortholog group splitting
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
                future_to_input = {executor.submit(split_dogs, input_data): input_data 
                                 for input_data in split_inputs}
                
                for future in _iter_progress(
                    concurrent.futures.as_completed(future_to_input),
                    total=len(split_inputs),
                    description="Domain OG splits",
                ):
                    try:
                        future.result()
                    except Exception as e:
                        input_data = future_to_input[future]
                        msg = f'Failed during domain ortholog group splitting job: {e}'
                        log_object.error(msg)
                        sys.stderr.write(msg + '\n')
                        sys.stderr.write(traceback.format_exc() + '\n')
                        # Continue with other jobs instead of exiting
        except Exception as e:
            msg = 'Failed during threading domain ortholog group splitting'
            log_object.error(msg)
            sys.stderr.write(msg + '\n')
            sys.stderr.write(traceback.format_exc() + '\n')
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
        input: List containing [og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, og_split_file, 
        log_object, threads]
    """
    try:
        og, cog_list_file, cog_dist_file, tre_file, spl_file, log_object, threads = input

        all_prots = set([])
        with open(cog_list_file) as oclf:
            for line in oclf:
                line = line.strip()
                all_prots.add(line)

        all_dogs = set([])
        for p in all_prots:
            for d in protein_dogs[p]:
                if protein_dogs[p][d] > 0:
                    all_dogs.add(d)

        prot_dog_vectors = defaultdict(list)
        for p in all_prots:
            for d in sorted(all_dogs):
                if p in protein_dogs and d in protein_dogs[p]:
                    prot_dog_vectors[p].append(protein_dogs[p][d])
                else:
                    prot_dog_vectors[p].append(0)

        prot_dog_vectors_np = {}
        for p in prot_dog_vectors:
            prot_dog_vectors_np[p] = np.array(prot_dog_vectors[p])

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
        sp_refined = further_split_outlier_artifact_groups(t, sp)
        sp_further_split = further_split_disjoint_dog_partitions(t, sp_refined)
        spl_outf = open(spl_file, 'w')
        for spi in sp_further_split:
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


def further_split_outlier_artifact_groups(rooted_t: Tree, sps: List[Set[str]]) -> List[Set[str]]:
    """
    Further split ortholog groups that may contain small proteins or be the result
    of outlier sequences being claded together after midpoint rooting of FastME tree.
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
                        union_count += p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
                        intersect_count += min([p1dogs[d], p2dogs[d]])

                    if union_count > 0:
                        jaccard_index = intersect_count/union_count
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
                        union_count += p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
                        intersect_count += min([p1dogs[d], p2dogs[d]])

                    if union_count > 0:
                        jaccard_index = intersect_count/union_count
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


def further_split_disjoint_dog_partitions(rooted_t: Tree, sps: List[Set[str]]) -> List[Set[str]]:
    """
    Function to further split disjoint domain ortholog groups based on the phylogenetic tree.
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
                    if nc.issubset(sp) and len(nc.difference(sp)) == 0 and len(nc.intersection(accounted_leafs)) == 0:
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
                p1dogs = {}
                for dog, count in protein_dogs[p1].items():
                    if count != 0:
                        p1dogs[dog] = count
                p1dogs_keys = set(p1dogs.keys())

                max_jacc_internal = 0.0
                max_jacc_external = 0.0
                for p2 in sorted(sp):
                    if p1 == p2:
                        continue

                    p2dogs = {}
                    for dog, count in protein_dogs[p2].items():
                        if count != 0:
                            p2dogs[dog] = count

                    p2dogs_keys = set(p2dogs.keys())
                    union_dogs = p1dogs_keys.union(p2dogs_keys)

                    union_count = 0
                    intersect_count = 0
                    for d in union_dogs:
                        p1dc = 0
                        p2dc = 0
                        if d in p1dogs:
                            p1dc = p1dogs[d]
                        if d in p2dogs:
                            p2dc = p2dogs[d]
                        union_count += p1dc + p2dc - min([p1dc, p2dc])
                        intersect_count += min([p1dc, p2dc])

                    if union_count > 0:
                        jaccard_index = intersect_count / union_count
                        if jaccard_index > max_jacc_internal:
                            max_jacc_internal = jaccard_index

                for p2 in sorted(all_proteins):
                    if p2 in sp:
                        continue
                    p2dogs = {}
                    for dog, count in protein_dogs[p2].items():
                        if count != 0:
                            p2dogs[dog] = count

                    p2dogs_keys = set(p2dogs.keys())
                    union_dogs = p1dogs_keys.union(p2dogs_keys)

                    union_count = 0
                    intersect_count = 0
                    for d in union_dogs:
                        p1dc = 0
                        p2dc = 0
                        if d in p1dogs:
                            p1dc = p1dogs[d]
                        if d in p2dogs:
                            p2dc = p2dogs[d]
                        union_count += p1dc + p2dc - min([p1dc, p2dc])
                        intersect_count += min([p1dc, p2dc])

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
        from .utils import single_linkage_cluster
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


def pairwise_dist(inputs: List[Any]) -> None:
    """
    Calculate pairwise distance between two tree leaves.

    Args:
        inputs: List containing [d, l1, l2] where d is distance dict, 
        l1, l2 are leaves

    Returns:
        None: Updates distance dictionary
    """
    global tree_obj
    d, l1, l2 = inputs
    dist = tree_obj.get_distance(l1, l2)
    key = tuple(sorted([l1, l2]))
    d[key] = dist


def split_dogs(inputs: List[Any]) -> None:
    """
    Split domain ortholog groups using phylogenetic approach.

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
    rooting_seeds = max([rooting_seeds, 1])
    try:
        # Use global tree object if available, otherwise load from file
        if tree_obj is not None:
            t = tree_obj
        else:
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
                    log_object.error(f"Error setting outgroup for node {n.name}: {traceback.format_exc()}")
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
            for i, sp in enumerate(sorted(all_rooting_partitions, key=itemgetter(1, 2, 3))):
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
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        sys.exit(1)
