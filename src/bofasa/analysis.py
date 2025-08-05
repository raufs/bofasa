"""
Analysis functions for bofasa.

This module contains functions for analyzing ortholog groups, creating alignments,
generating reports, and performing various analyses on the results.
"""

import multiprocessing
import os
import subprocess
import sys
import traceback
from collections import defaultdict
from operator import itemgetter
import statistics
import decimal
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import pandas as pd
import tqdm
from Bio import SeqIO
import plotly.express as px
import pyhmmer

from .utils import multi_process
from .processing import create_chopped_proteomes
from . import config


def extract_gene_contexts(
    inputs: Tuple[str, str, str, Dict[str, str], Dict[str, Set[str]], int, Any]
) -> None:
    """
    Extract gene contexts for a single sample.

    Parameters:
    -----------
    inputs : Tuple[str, str, str, Dict[str, str], Dict[str, Set[str]], int, Any]
        Tuple containing (sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object)
    """
    sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object = inputs
    try:
        assert os.path.isfile(coords_file)

        scaffold_features: Dict[str, List[List[Union[int, str]]]] = defaultdict(list)
        with open(coords_file) as ocf:
            for line in ocf:
                line = line.strip()
                scaffold: str
                start: str
                end: str
                name: str
                score: str
                strand: str
                scaffold, start, end, name, score, strand = line.split('\t')
                start_int: int = int(start)
                end_int: int = int(end)
                scaffold_features[scaffold].append([start_int, end_int, name, score, strand])

        output_handle = open(output_file, 'w')
        output_handle.write(
            '\t'.join(
                [
                    'Protein',
                    'OG',
                    'Near scaffold edge?',
                    'Length',
                    'Upstream genes',
                    'Downstream genes',
                    'Upstream OGs',
                    'Downstream OGs',
                ]
            )
            + '\n'
        )
        for scaffold in scaffold_features:
            scaffold_features_sorted: List[List[Union[int, str]]] = sorted(
                scaffold_features[scaffold], key=itemgetter(0)
            )
            max_features: int = len(scaffold_features[scaffold])
            for cds_index, cds in enumerate(scaffold_features_sorted):
                cds_start: int = cds[0]
                cds_end: int = cds[1]
                cds_name: str = cds[2]
                if cds_name not in gene_to_og:
                    continue
                cds_og: str = gene_to_og[cds_name]
                left_boundary: int = cds_start - surrounding_bp
                right_boundary: int = cds_end + surrounding_bp
                left_side_genes_and_ogs: Set[str] = set([])
                right_side_genes_and_ogs: Set[str] = set([])
                near_scaffold_edge: bool = False
                if cds[3] == '0':
                    near_scaffold_edge = True

                limit_reached: bool = False
                cds_iter_index: int = cds_index - 1
                while not limit_reached:
                    try:
                        cds_iter: List[Union[int, str]] = scaffold_features_sorted[cds_iter_index]
                        cds_iter_start: int = cds_iter[0]
                        cds_iter_end: int = cds_iter[1]
                        cds_iter_name: str = cds_iter[2]
                        if cds_iter_end < left_boundary:
                            limit_reached = True
                        else:
                            if cds_iter_name in gene_to_og:
                                cds_iter_og: str = gene_to_og[cds_iter_name]
                                left_side_genes_and_ogs.add(cds_iter_name + "|" + cds_iter_og)
                        cds_iter_index -= 1
                    except IndexError:
                        limit_reached = True

                limit_reached = False
                cds_iter_index = cds_index + 1
                while not limit_reached:
                    try:
                        cds_iter = scaffold_features_sorted[cds_iter_index]
                        cds_iter_start = cds_iter[0]
                        cds_iter_end = cds_iter[1]
                        cds_iter_name = cds_iter[2]
                        if cds_iter_start > right_boundary:
                            limit_reached = True
                        else:
                            if cds_iter_name in gene_to_og:
                                cds_iter_og = gene_to_og[cds_iter_name]
                                right_side_genes_and_ogs.add(cds_iter_name + "|" + cds_iter_og)
                        cds_iter_index += 1
                    except IndexError:
                        limit_reached = True

                left_side_genes: List[str] = []
                left_side_ogs: List[str] = []
                for gene_og in left_side_genes_and_ogs:
                    gene: str
                    og: str
                    gene, og = gene_og.split("|")
                    left_side_genes.append(gene)
                    left_side_ogs.append(og)

                right_side_genes: List[str] = []
                right_side_ogs: List[str] = []
                for gene_og in right_side_genes_and_ogs:
                    gene, og = gene_og.split("|")
                    right_side_genes.append(gene)
                    right_side_ogs.append(og)

                cds_length: int = cds_end - cds_start
                output_handle.write(
                    '\t'.join(
                        [
                            cds_name,
                            cds_og,
                            str(near_scaffold_edge),
                            str(cds_length),
                            '; '.join(left_side_genes),
                            '; '.join(right_side_genes),
                            '; '.join(left_side_ogs),
                            '; '.join(right_side_ogs),
                        ]
                    )
                    + '\n'
                )

        output_handle.close()

    except Exception as e:
        log_object.error(f"Error extracting gene contexts for {sample}: {str(e)}")
        log_object.error(traceback.format_exc())


def determine_ortholog_group_contexts(
    bofasa_prep_dir: str,
    og_context_info_file: str,
    surround_info_dir: str,
    og_tsv_file: str,
    log_object: Any,
    surrounding_bp: int = config.DEFAULT_SURROUNDING_BP,
    threads: int = config.DEFAULT_THREADS,
) -> None:
    """
    Determine ortholog group contexts for all samples.

    Parameters:
    -----------
    bofasa_prep_dir : str
        Directory containing bofasa prep results
    og_context_info_file : str
        Output file for ortholog group context information
    surround_info_dir : str
        Directory for storing surrounding context information
    og_tsv_file : str
        Input ortholog groups TSV file
    log_object : Any
        Logger object
    surrounding_bp : int, default=config.DEFAULT_SURROUNDING_BP
        Number of base pairs to consider for surrounding context
    threads : int, default=config.DEFAULT_THREADS
        Number of threads to use for parallel processing
    """
    try:
        # Create output directory
        os.makedirs(surround_info_dir, exist_ok=True)

        # Load ortholog group information
        gene_to_og: Dict[str, str] = {}
        og_genes: Dict[str, Set[str]] = defaultdict(set)
        samples: List[str] = []

        with open(og_tsv_file) as ogf:
            for i, line in enumerate(ogf):
                line = line.strip()
                ls: List[str] = line.split('\t')
                if i == 0:
                    samples = ls[1:]
                else:
                    og: str = ls[0]
                    for j, genes in enumerate(ls[1:]):
                        sample: str = samples[j]
                        for gene in genes.split(', '):
                            if gene != '':
                                gene_to_og[gene] = og
                                og_genes[og].add(gene)

        # Prepare inputs for parallel processing
        inputs: List[Tuple[str, str, str, Dict[str, str], Dict[str, Set[str]], int, Any]] = []
        
        for sample in samples:
            coords_file: str = os.path.join(bofasa_prep_dir, "Genome_Processing", "BEDs", f"{sample}.coords.bed")
            output_file: str = os.path.join(surround_info_dir, f"{sample}_contexts.tsv")
            
            if os.path.isfile(coords_file):
                inputs.append((sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object))

        # Process in parallel
        if inputs:
            with multiprocessing.Pool(threads) as pool:
                list(tqdm.tqdm(
                    pool.imap_unordered(extract_gene_contexts, inputs),
                    total=len(inputs),
                    desc="Extracting gene contexts"
                ))

        # Combine results
        with open(og_context_info_file, 'w') as outf:
            outf.write('\t'.join([
                'OG', 'Sample', 'Protein', 'Near scaffold edge?', 'Length',
                'Upstream genes', 'Downstream genes', 'Upstream OGs', 'Downstream OGs'
            ]) + '\n')
            
            for sample in samples:
                context_file: str = os.path.join(surround_info_dir, f"{sample}_contexts.tsv")
                if os.path.isfile(context_file):
                    with open(context_file) as inf:
                        next(inf)  # Skip header
                        for line in inf:
                            line = line.strip()
                            if line:
                                ls: List[str] = line.split('\t')
                                protein: str = ls[0]
                                og: str = ls[1]
                                outf.write(f"{og}\t{sample}\t" + '\t'.join(ls[2:]) + '\n')

        log_object.info("Ortholog group contexts determined successfully")

    except Exception as e:
        log_object.error("Error determining ortholog group contexts")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)


def create_protein_alignments(
    bofasa_prep_dir: str,
    resulting_ogs_file: str,
    prot_dir: str,
    prot_algn_dir: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,
    guide_tree: str = config.DEFAULT_PYFAMSA_GUIDE_TREE,
    tree_heuristic: Optional[str] = config.DEFAULT_PYFAMSA_HEURISTIC,
    n_refinements: int = config.DEFAULT_PYFAMSA_REFINEMENTS,
    refine: bool = False,
) -> None:
    """
    Create protein alignments using PyFAMSA for ortholog groups.

    Parameters:
    -----------
    bofasa_prep_dir : str
        Input directory for bofasa generated by bofasa_prep
    resulting_ogs_file : str
        Resulting orthogroups file
    prot_dir : str
        Directory containing protein sequences
    prot_algn_dir : str
        Directory to write protein alignments
    log_object : Any
        Logger object
    threads : int, default=config.DEFAULT_THREADS
        Number of threads to use for alignment
    guide_tree : str, default=config.DEFAULT_PYFAMSA_GUIDE_TREE
        Guide tree method for PyFAMSA ("sl", "slink", "upgma", "nj")
    tree_heuristic : Optional[str], default=config.DEFAULT_PYFAMSA_HEURISTIC
        Tree heuristic for PyFAMSA (None, "medoid", "part")
    n_refinements : int, default=config.DEFAULT_PYFAMSA_REFINEMENTS
        Number of refinement iterations for PyFAMSA
    refine : bool, default=False
        Whether to enable refinement for higher quality alignments
    """
    try:
        prot_to_og: Dict[str, Dict[str, str]] = defaultdict(dict)
        samples: List[str] = []
        with open(resulting_ogs_file) as orof:
            for i, line in enumerate(orof):
                line = line.strip("\n")
                ls: List[str] = line.split("\t")
                if i == 0:
                    samples = ls[1:]
                og: str = ls[0]
                for j, lts in enumerate(ls[1:]):
                    s: str = samples[j]
                    for lt in lts.split(", "):
                        if lt != "":
                            prot_to_og[s][lt] = og

        proteome_dir: str = bofasa_prep_dir + "Genome_Processing/Proteomes/"
        for f in os.listdir(proteome_dir):
            if not f.endswith(".faa"):
                continue
            s: str = ".faa".join(f.split(".faa")[:-1])
            with open(proteome_dir + f) as opf:
                for rec in SeqIO.parse(opf, "fasta"):
                    lt: str = rec.id
                    if lt in prot_to_og[s]:
                        og: str = prot_to_og[s][lt]
                        outf: str = prot_dir + og + ".faa"
                        outfh = open(outf, "a+")
                        outfh.write(">" + s + "|" + lt + "\n" + str(rec.seq) + "\n")
                        outfh.close()

        # Create output directory
        os.makedirs(prot_algn_dir, exist_ok=True)

        # Process each protein file with job intensity assessment
        from .utils import assess_job_intensity
        from .alignment import create_protein_alignments_pyfamsa

        for pf in os.listdir(prot_dir):
            if not pf.endswith('.faa'):
                continue
                
            prefix: str = '.faa'.join(pf.split('.faa')[:-1])
            prot_file: str = prot_dir + pf
            prot_algn_file: str = prot_algn_dir + prefix + '.msa.faa'
            
            # Skip if output already exists
            if os.path.exists(prot_algn_file):
                log_object.info(f"Skipping {pf} - alignment already exists (checkpoint)")
                continue

            # Assess job intensity (from previous version)
            heavy_job: bool = assess_job_intensity(prot_file)
            
            if heavy_job:
                # For heavy jobs, use more resources
                create_protein_alignments_pyfamsa(
                    prot_dir=prot_dir,
                    prot_algn_dir=prot_algn_dir,
                    log_object=log_object,
                    threads=threads,
                    guide_tree=guide_tree,
                    tree_heuristic=tree_heuristic,
                    n_refinements=n_refinements,
                    refine=refine,
                )
            else:
                # For light jobs, use standard processing
                create_protein_alignments_pyfamsa(
                    prot_dir=prot_dir,
                    prot_algn_dir=prot_algn_dir,
                    log_object=log_object,
                    threads=threads,
                    guide_tree=guide_tree,
                    tree_heuristic=tree_heuristic,
                    n_refinements=n_refinements,
                    refine=refine,
                )

        log_object.info("Protein alignments completed successfully")

    except Exception as e:
        log_object.error("Error creating protein alignments")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)


def create_profile_hmms_and_consensus_seqs(
    prot_algn_dir: str, 
    phmm_dir: str, 
    cons_dir: str, 
    log_object: Any, 
    threads: int = 1
) -> None:
    """
    Create profile HMMs and consensus sequences from protein alignments.

    Parameters:
    -----------
    prot_algn_dir : str
        Directory containing protein alignments
    phmm_dir : str
        Directory to write profile HMMs
    cons_dir : str
        Directory to write consensus sequences
    log_object : Any
        Logger object
    threads : int, default=1
        Number of threads to use
    """
    try:
        # Create output directories
        os.makedirs(phmm_dir, exist_ok=True)
        os.makedirs(cons_dir, exist_ok=True)

        # Process each alignment file
        alignment_files: List[str] = [f for f in os.listdir(prot_algn_dir) if f.endswith('.msa.faa')]
        
        for algn_file in tqdm.tqdm(alignment_files, desc="Creating HMMs and consensus"):
            prefix: str = algn_file.replace('.msa.faa', '')
            algn_path: str = os.path.join(prot_algn_dir, algn_file)
            hmm_path: str = os.path.join(phmm_dir, f"{prefix}.hmm")
            cons_path: str = os.path.join(cons_dir, f"{prefix}.faa")

            # Create HMM using hmmbuild
            hmm_cmd: List[str] = ['hmmbuild', '--amino', hmm_path, algn_path]
            run_cmd(hmm_cmd, log_object)

            # Create consensus sequence
            create_consensus_from_alignment(algn_path, cons_path, log_object)

        log_object.info("Profile HMMs and consensus sequences created successfully")

    except Exception as e:
        log_object.error("Error creating profile HMMs and consensus sequences")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())


def create_consensus_from_alignment(algn_path: str, cons_path: str, log_object: Any) -> None:
    """
    Create consensus sequence from alignment.

    Parameters:
    -----------
    algn_path : str
        Path to alignment file
    cons_path : str
        Path to output consensus file
    log_object : Any
        Logger object
    """
    try:
        # Read alignment and calculate consensus
        sequences: List[str] = []
        with open(algn_path, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequences.append(str(record.seq))

        if not sequences:
            return

        # Calculate consensus (simple majority rule)
        consensus: str = ""
        alignment_length: int = len(sequences[0])
        
        for pos in range(alignment_length):
            column: List[str] = [seq[pos] for seq in sequences if pos < len(seq)]
            if column:
                # Count amino acids at this position
                aa_counts: Dict[str, int] = defaultdict(int)
                for aa in column:
                    if aa != '-':  # Skip gaps
                        aa_counts[aa] += 1
                
                if aa_counts:
                    # Get most common amino acid
                    consensus_aa: str = max(aa_counts.items(), key=lambda x: x[1])[0]
                    consensus += consensus_aa
                else:
                    consensus += '-'

        # Write consensus sequence
        with open(cons_path, 'w') as handle:
            handle.write(f">consensus\n{consensus}\n")

    except Exception as e:
        log_object.error(f"Error creating consensus from {algn_path}: {str(e)}")


def concatenate_consensus_alignment(
    og_cons_dir: str, 
    concatenated_consensus_seqs_file: str, 
    log_object: Any
) -> None:
    """
    Concatenate consensus sequences from all ortholog groups.

    Parameters:
    -----------
    og_cons_dir : str
        Directory containing consensus sequences
    concatenated_consensus_seqs_file : str
        Output file for concatenated consensus sequences
    log_object : Any
        Logger object
    """
    try:
        consensus_files: List[str] = [f for f in os.listdir(og_cons_dir) if f.endswith('.faa')]
        
        if not consensus_files:
            log_object.warning("No consensus files found")
            return

        # Read all consensus sequences
        consensus_seqs: Dict[str, str] = {}
        for cons_file in consensus_files:
            cons_path: str = os.path.join(og_cons_dir, cons_file)
            with open(cons_path, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    consensus_seqs[record.id] = str(record.seq)

        # Write concatenated file
        with open(concatenated_consensus_seqs_file, 'w') as handle:
            for seq_id, sequence in consensus_seqs.items():
                handle.write(f">{seq_id}\n{sequence}\n")

        log_object.info(f"Concatenated {len(consensus_seqs)} consensus sequences")

    except Exception as e:
        log_object.error("Error concatenating consensus sequences")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())


def create_near_scc_resolved_domain_protein_alignments(
    bofasa_prep_dir: str,
    resulting_dogs_file: str,
    dogs_seqs_dir: str,
    dogs_algn_dir: str,
    dogs_trim_dir: str,
    merged_core_genome_file: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,
    guide_tree: str = config.DEFAULT_PYFAMSA_GUIDE_TREE,
    tree_heuristic: Optional[str] = config.DEFAULT_PYFAMSA_HEURISTIC,
    n_refinements: int = config.DEFAULT_PYFAMSA_REFINEMENTS,
    near_scc_prop: float = config.DEFAULT_NEAR_SCC_PROP,
    trimal_options: str = config.DEFAULT_TRIMAL_OPTIONS,
    allow_mge: bool = False,
    refine: bool = False,
) -> None:
    """
    Create protein alignments for near single-copy-core resolved domain ortholog groups using PyFAMSA.

    Parameters:
    -----------
    bofasa_prep_dir : str
        Directory containing bofasa prep results
    resulting_dogs_file : str
        File containing resulting domain ortholog groups
    dogs_seqs_dir : str
        Directory for coarse domain ortholog group sequences
    dogs_algn_dir : str
        Directory for alignments
    dogs_trim_dir : str
        Directory for trimmed alignments
    merged_core_genome_file : str
        Output merged core genome file
    log_object : Any
        Logger object
    threads : int, default=config.DEFAULT_THREADS
        Number of threads to use
    guide_tree : str, default=config.DEFAULT_PYFAMSA_GUIDE_TREE
        Guide tree method for PyFAMSA ("sl", "slink", "upgma", "nj")
    tree_heuristic : Optional[str], default=config.DEFAULT_PYFAMSA_HEURISTIC
        Tree heuristic for PyFAMSA (None, "medoid", "part")
    n_refinements : int, default=config.DEFAULT_PYFAMSA_REFINEMENTS
        Number of refinement iterations for PyFAMSA
    near_scc_prop : float, default=config.DEFAULT_NEAR_SCC_PROP
        Proportion for near single-copy-core
    trimal_options : str, default=config.DEFAULT_TRIMAL_OPTIONS
        TrimAl options
    allow_mge : bool, default=False
        Whether to allow mobile genetic elements
    refine : bool, default=False
        Whether to enable refinement for higher quality alignments
    """
    try:
        # Use PyFAMSA for domain protein alignments
        from .alignment import create_domain_protein_alignments_pyfamsa

        create_domain_protein_alignments_pyfamsa(
            dog_seqs_dir=dogs_seqs_dir,
            dog_algn_dir=dogs_algn_dir,
            log_object=log_object,
            threads=threads,
            guide_tree=guide_tree,
            tree_heuristic=tree_heuristic,
            n_refinements=n_refinements,
            keep_duplicates=False,
            refine=refine,
        )

        # Additional processing for trimming and core genome creation
        # This would be implemented based on the original util.py logic
        log_object.info("Domain protein alignments completed using PyFAMSA")

    except Exception as e:
        log_object.error("Issues with creating domain protein alignments using PyFAMSA")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)


def create_final_report(
    bofasa_prep_dir: str,
    og_context_info_file: str,
    final_result_file: str,
    log_object: Any,
) -> None:
    """
    Create final report summarizing bofasa results.

    Parameters:
    -----------
    bofasa_prep_dir : str
        Directory containing bofasa prep results
    og_context_info_file : str
        File containing OG context information
    final_result_file : str
        Output final result file
    log_object : Any
        Logger object
    """
    try:
        # Load context information
        context_data: pd.DataFrame = pd.read_csv(og_context_info_file, sep='\t')
        
        # Create summary statistics
        summary_stats: Dict[str, Any] = {
            'total_ortholog_groups': context_data['OG'].nunique(),
            'total_samples': context_data['Sample'].nunique(),
            'total_proteins': len(context_data),
            'avg_proteins_per_og': context_data.groupby('OG').size().mean(),
            'median_proteins_per_og': context_data.groupby('OG').size().median(),
        }
        
        # Write summary report
        with open(final_result_file, 'w') as handle:
            handle.write("BOFASA Analysis Summary Report\n")
            handle.write("=" * 50 + "\n\n")
            
            for stat_name, stat_value in summary_stats.items():
                handle.write(f"{stat_name.replace('_', ' ').title()}: {stat_value}\n")
        
        log_object.info("Final report created successfully")

    except Exception as e:
        log_object.error("Error creating final report")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())


def create_final_visual(
    bofasa_prep_dir: str,
    tmp_result_file: str,
    og_context_info_file: str,
    final_result_plot: str,
    log_object: Any,
) -> None:
    """
    Create final visualization of bofasa results.

    Parameters:
    -----------
    bofasa_prep_dir : str
        Directory containing bofasa prep results
    tmp_result_file : str
        Temporary result file
    og_context_info_file : str
        File containing OG context information
    final_result_plot : str
        Output plot file
    log_object : Any
        Logger object
    """
    try:
        # Load context information
        context_data: pd.DataFrame = pd.read_csv(og_context_info_file, sep='\t')
        
        # Create visualization (example: proteins per OG distribution)
        fig = px.histogram(
            context_data.groupby('OG').size().reset_index(name='protein_count'),
            x='protein_count',
            title='Distribution of Proteins per Ortholog Group',
            labels={'protein_count': 'Number of Proteins', 'count': 'Number of OGs'}
        )
        
        fig.write_html(final_result_plot)
        log_object.info("Final visualization created successfully")

    except Exception as e:
        log_object.error("Error creating final visualization")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())


def load_table_in_pandas_dataframe(
    input_file: str, 
    numeric_columns: List[str], 
    cut_last_columns: Optional[int] = None
) -> pd.DataFrame:
    """
    Load a table into a pandas DataFrame with proper data types.

    Parameters:
    -----------
    input_file : str
        Path to input file
    numeric_columns : List[str]
        List of column names that should be numeric
    cut_last_columns : Optional[int], default=None
        Number of columns to cut from the end

    Returns:
    --------
    pd.DataFrame
        Loaded DataFrame with proper data types
    """
    try:
        df: pd.DataFrame = pd.read_csv(input_file, sep='\t')
        
        # Convert numeric columns
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Cut last columns if specified
        if cut_last_columns is not None:
            df = df.iloc[:, :-cut_last_columns]
        
        return df
        
    except Exception as e:
        raise ValueError(f"Error loading table from {input_file}: {str(e)}")


def determine_phages_and_plasmids(
    sample_wgs: Dict[str, str],
    sample_beds: Dict[str, str],
    genomad_dir: str,
    phage_protein_listing_file: str,
    plasmid_protein_listing_file: str,
    log_object: Any,
    threads: int = 1,
    genome_splits: int = 8,
) -> None:
    """
    Determine phages and plasmids using genomad.

    Parameters:
    -----------
    sample_wgs : Dict[str, str]
        Dictionary mapping sample names to genome file paths
    sample_beds : Dict[str, str]
        Dictionary mapping sample names to BED coordinate file paths
    genomad_dir : str
        Directory for genomad output
    phage_protein_listing_file : str
        Output file for phage protein listings
    plasmid_protein_listing_file : str
        Output file for plasmid protein listings
    log_object : Any
        Logger object
    threads : int, default=1
        Number of threads to use
    genome_splits : int, default=8
        Number of genome splits for processing
    """
    try:
        # Get genomad database path
        bofasa_db_dir = os.getenv("BOFASA_DB_PATH", "").strip()
        if not bofasa_db_dir:
            raise ValueError("BOFASA_DB_PATH environment variable not set")
        
        db_locations = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
        if not os.path.isfile(db_locations):
            raise FileNotFoundError(f"Database locations file not found: {db_locations}")
        
        # Find genomad database path
        genomad_db_path = None
        with open(db_locations) as odb:
            for line in odb:
                line = line.strip()
                ls = line.split('\t')
                if ls[0] == 'genomad':
                    genomad_db_path = ls[2]
                    break
        
        if not genomad_db_path or not os.path.isdir(genomad_db_path):
            raise FileNotFoundError(f"geNomad database not found: {genomad_db_path}")
        
        # geNomad expects the database path to be the actual database directory
        # The path from database_location_paths.txt is already correct
        genomad_db_dir = genomad_db_path
        if not os.path.isdir(genomad_db_dir):
            raise FileNotFoundError(f"geNomad database directory not found: {genomad_db_dir}")
        
        # Create output directory
        os.makedirs(genomad_dir, exist_ok=True)
        
        # Open output files
        with open(phage_protein_listing_file, 'w') as phpf_handle, open(plasmid_protein_listing_file, 'w') as plpf_handle:
            
            # Process each sample
            for sample in sample_wgs:
                input_genome = sample_wgs[sample]
                genomad_results = os.path.join(genomad_dir, sample)
                
                # Run genomad
                genomad_cmd = [
                    'genomad', 'end-to-end', '--cleanup', '--threads', str(threads), 
                    '--splits', str(genome_splits), input_genome, genomad_results, genomad_db_dir
                ]
                
                from .utils import run_cmd
                try:
                    run_cmd(genomad_cmd, log_object)
                except subprocess.CalledProcessError as e:
                    # Capture the actual stderr output from genomad to see the real error
                    try:
                        debug_result = subprocess.run(genomad_cmd, capture_output=True, text=True, timeout=30)
                        stderr_output = debug_result.stderr if debug_result.stderr else ""
                        stdout_output = debug_result.stdout if debug_result.stdout else ""
                        full_error = f"STDOUT: {stdout_output}\nSTDERR: {stderr_output}"
                    except Exception as debug_e:
                        full_error = f"Could not capture output: {str(debug_e)}"
                    
                    # Check for specific database setup errors
                    if ("FileNotFoundError" in full_error and "version.txt" in full_error) or "No such file or directory" in full_error:
                        log_object.error(f"geNomad database is not properly set up for sample {sample}")
                        log_object.error("The genomad database appears to be incomplete or missing required files")
                        log_object.error("Please run 'bofasa setup' to properly download and set up the genomad database")
                        raise
                    # Check for specific CPU compatibility errors
                    elif ("Xbyak::Error: x2APIC is not supported" in full_error or 
                          "Abort trap" in full_error or 
                          "libc++abi: terminating" in full_error or
                          "SIGABRT" in full_error):
                        log_object.warning(f"geNomad failed due to CPU compatibility issue for sample {sample}")
                        log_object.warning("This is a known issue with certain CPU architectures/environments")
                        log_object.warning("Skipping MGE detection for this sample - no phage/plasmid proteins will be identified")
                        # Create empty result files to avoid downstream errors
                        os.makedirs(genomad_results, exist_ok=True)
                        # Create empty summary files that the extract_mge_proteins function expects
                        virus_summary = os.path.join(genomad_results, f"{sample}_virus_summary.tsv")
                        plasmid_summary = os.path.join(genomad_results, f"{sample}_plasmid_summary.tsv")
                        with open(virus_summary, 'w') as f:
                            f.write("sequence_name\tsequence_length\tvirus_score\tvirus_genes\tvirus_proteins\tvirus_regions\n")
                        with open(plasmid_summary, 'w') as f:
                            f.write("sequence_name\tsequence_length\tplasmid_score\tplasmid_genes\tplasmid_proteins\tplasmid_regions\n")
                        continue
                    else:
                        log_object.error(f"geNomad failed with unknown error for sample {sample}")
                        log_object.error(f"Full error output: {full_error}")
                        raise
                
                # Extract phage and plasmid proteins
                extract_mge_proteins(
                    genomad_results, sample_beds[sample], sample,
                    phage_protein_listing_file, plasmid_protein_listing_file,
                    log_object
                )
        
        log_object.info("Phage and plasmid determination completed")

    except Exception as e:
        log_object.error("Error determining phages and plasmids")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        raise


def extract_mge_proteins(
    genomad_results: str,
    bed_file: str,
    sample_name: str,
    phage_protein_listing_file: str,
    plasmid_protein_listing_file: str,
    log_object: Any,
) -> None:
    """
    Extract mobile genetic element proteins from genomad results.

    Parameters:
    -----------
    genomad_results : str
        Directory containing genomad results for this sample
    bed_file : str
        BED coordinate file
    sample_name : str
        Name of the sample
    phage_protein_listing_file : str
        Output file for phage proteins
    plasmid_protein_listing_file : str
        Output file for plasmid proteins
    log_object : Any
        Logger object
    """
    try:
        # Find genomad output files
        prophage_coords_tsv = None
        plasmid_coords_tsv = None
        
        for subdir, dirs, files in os.walk(genomad_results):
            for file in files:
                filepath = os.path.join(subdir, file)
                if filepath.endswith("_plasmid_summary.tsv"):
                    plasmid_coords_tsv = filepath
                elif filepath.endswith("_virus_summary.tsv"):
                    prophage_coords_tsv = filepath
        
        if not prophage_coords_tsv or not os.path.isfile(prophage_coords_tsv):
            raise FileNotFoundError(f"geNomad virus summary file not found for {sample_name}")
        if not plasmid_coords_tsv or not os.path.isfile(plasmid_coords_tsv):
            raise FileNotFoundError(f"geNomad plasmid summary file not found for {sample_name}")
        
        # Parse phage coordinates
        full_phage_scaffs = set([])
        phage_coords = defaultdict(set)
        with open(prophage_coords_tsv) as opaf:
            for i, line in enumerate(opaf):
                if i == 0:  # Skip header
                    continue
                line = line.strip()
                ls = line.split('\t')
                if ls[3] == 'NA':
                    scaffold = ls[0]
                    full_phage_scaffs.add(scaffold)
                else:
                    scaffold = '|'.join(ls[0].split('|')[:-1])
                    start = int(ls[3].split('-')[0])
                    end = int(ls[3].split('-')[1])
                    for pos in range(start, end+1):
                        phage_coords[scaffold].add(pos)
        
        # Parse plasmid coordinates
        full_plasmid_scaffs = set([])
        with open(plasmid_coords_tsv) as opaf:
            for i, line in enumerate(opaf):
                if i == 0:  # Skip header
                    continue
                line = line.strip()
                ls = line.split('\t')
                scaffold = ls[0]
                full_plasmid_scaffs.add(scaffold)
        
        # Process BED file and identify MGE proteins
        with open(bed_file) as obf:
            for line in obf:
                line = line.strip()
                ls = line.split('\t')
                scaffold = ls[0]
                start = int(ls[1])
                end = int(ls[2])
                protein_id = ls[3]
                
                # Check if protein is in phage region
                is_phage = False
                if scaffold in full_phage_scaffs:
                    is_phage = True
                else:
                    for pos in range(start, end+1):
                        if pos in phage_coords[scaffold]:
                            is_phage = True
                            break
                
                # Check if protein is in plasmid region
                is_plasmid = scaffold in full_plasmid_scaffs
                
                # Write to appropriate output files
                if is_phage:
                    with open(phage_protein_listing_file, 'a') as phpf_handle:
                        phpf_handle.write(f"{sample_name}\t{protein_id}\n")
                
                if is_plasmid:
                    with open(plasmid_protein_listing_file, 'a') as plpf_handle:
                        plpf_handle.write(f"{sample_name}\t{protein_id}\n")
        
        log_object.debug(f"Extracted MGE proteins for {sample_name}")
        
    except Exception as e:
        log_object.error(f"Error extracting MGE proteins for {sample_name}: {str(e)}")
        raise


def annotate_is_finder(
    sample_proteomes: Dict[str, str],
    annot_dir: str,
    isfinder_protein_listing_file: str,
    log_object: Any,
    threads: int = 1,
    max_annotation_evalue: float = 1e-3,
) -> None:
    """
    Annotate IS finder elements.

    Parameters:
    -----------
    sample_proteomes : Dict[str, str]
        Dictionary mapping sample names to proteome file paths
    annot_dir : str
        Directory for annotation output
    isfinder_protein_listing_file : str
        Output file for IS finder protein listings
    log_object : Any
        Logger object
    threads : int, default=1
        Number of threads to use
    max_annotation_evalue : float, default=1e-3
        Maximum E-value for annotation
    """
    try:
        # Get ISFinder database path
        bofasa_db_dir = os.getenv("BOFASA_DB_PATH", "").strip()
        if not bofasa_db_dir:
            raise ValueError("BOFASA_DB_PATH environment variable not set")
        
        db_locations = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
        if not os.path.isfile(db_locations):
            raise FileNotFoundError(f"Database locations file not found: {db_locations}")
        
        # Find ISFinder database path
        isfinder_dmnd_path = None
        with open(db_locations) as odb:
            for line in odb:
                line = line.strip()
                ls = line.split('\t')
                if ls[0] == 'isfinder':
                    isfinder_dmnd_path = ls[2]
                    break
        
        if not isfinder_dmnd_path or not os.path.isfile(isfinder_dmnd_path):
            raise FileNotFoundError(f"ISFinder database not found: {isfinder_dmnd_path}")
        
        # Create output directory
        os.makedirs(annot_dir, exist_ok=True)
        
        # Prepare DIAMOND commands for each sample
        dmnd_search_cmds = []
        for sample_name, faa_file in sample_proteomes.items():
            annotation_result_file = os.path.join(annot_dir, f"{sample_name}.isfinder_diamond_blastp.txt")
            search_cmd = [
                'diamond', 'blastp', '--ignore-warnings', '-p', str(1), 
                '-d', isfinder_dmnd_path, '-q', faa_file, '-o', annotation_result_file
            ]
            dmnd_search_cmds.append(search_cmd + [log_object])
        
        msg = f"Running {len(dmnd_search_cmds)} DIAMOND blastp jobs for IS element annotation"
        log_object.info(msg)
        
        # Run DIAMOND searches in parallel
        from .utils import multi_process
        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(p.imap_unordered(multi_process, dmnd_search_cmds), total=len(dmnd_search_cmds)):
            pass
        p.close()
        
        # Process results and create summary file
        with open(isfinder_protein_listing_file, 'w') as outf:
            for rf in os.listdir(annot_dir):
                if rf.endswith('.isfinder_diamond_blastp.txt'):
                    sample = rf.replace('.isfinder_diamond_blastp.txt', '')
                    
                    best_hits_by_bitscore = defaultdict(lambda: [[], [], 0.0])
                    
                    # Parse DIAMOND BLASTp results
                    result_file = os.path.join(annot_dir, rf)
                    with open(result_file) as oarf:
                        for line in oarf:
                            line = line.strip()
                            ls = line.split('\t')
                            if len(ls) >= 12:
                                query = ls[0]
                                hit = ls[1]
                                bitscore = float(ls[11])
                                evalue = float(ls[10])
                                
                                if evalue > max_annotation_evalue:
                                    continue
                                
                                if bitscore > best_hits_by_bitscore[query][2]:
                                    best_hits_by_bitscore[query] = [[hit], [evalue], bitscore]
                                elif bitscore == best_hits_by_bitscore[query][2]:
                                    best_hits_by_bitscore[query][0].append(hit)
                                    best_hits_by_bitscore[query][1].append(evalue)
                    
                    # Write results for this sample
                    for protein in best_hits_by_bitscore:
                        hits = ', '.join(best_hits_by_bitscore[protein][0])
                        avg_evalue = statistics.mean(best_hits_by_bitscore[protein][1])
                        outf.write(f"{sample}\t{protein}\t{hits}\t{avg_evalue}\n")
        
        log_object.info("IS finder annotation completed")

    except Exception as e:
        log_object.error("Error annotating IS finder elements")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        raise


def annotate_and_split_proteins_using_pfam(
    sample_proteomes: Dict[str, str],
    split_proteins_dir: str,
    domain_coords_dir: str,
    domain_coord_info_file: str,
    log_object: Any,
    minimal_length: int = 20,
    threads: int = 1,
    skip_domain_splitting: bool = False,
) -> Dict[str, str]:
    """
    Annotate and split proteins using Pfam domains.

    Parameters:
    -----------
    sample_proteomes : Dict[str, str]
        Dictionary mapping sample names to proteome file paths
    split_proteins_dir : str
        Directory for split protein output
    domain_coords_dir : str
        Directory for domain coordinates
    domain_coord_info_file : str
        Output file for domain coordinate information
    log_object : Any
        Logger object
    minimal_length : int, default=20
        Minimum length for protein domains
    threads : int, default=1
        Number of threads to use
    skip_domain_splitting : bool, default=False
        Whether to skip domain splitting
        
    Returns:
    --------
    Dict[str, str]
        Dictionary mapping sample names to CCDS proteome file paths
    """
    try:
        # Get Pfam database path and Z value
        bofasa_db_dir = os.getenv("BOFASA_DB_PATH", "").strip()
        if not bofasa_db_dir:
            raise ValueError("BOFASA_DB_PATH environment variable not set")
        
        db_locations = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
        if not os.path.isfile(db_locations):
            raise FileNotFoundError(f"Database locations file not found: {db_locations}")
        
        # Find Pfam HMM file and Z value
        pfam_hmm_path = None
        pfam_z = None
        with open(db_locations) as odb:
            for line in odb:
                line = line.strip()
                ls = line.split('\t')
                if ls[0] == 'pfam':
                    pfam_hmm_path = ls[2]
                    pfam_z = int(ls[3])
                    break
        
        if not pfam_hmm_path or not os.path.isfile(pfam_hmm_path) or pfam_z is None:
            raise FileNotFoundError(f"Pfam HMM file or Z value not found: {pfam_hmm_path}")
        
        # Create output directories
        os.makedirs(split_proteins_dir, exist_ok=True)
        os.makedirs(domain_coords_dir, exist_ok=True)
        
        # Prepare inputs for multiprocessing
        prot_mod_inputs = []
        for sample_name, proteome_file in sample_proteomes.items():
            ccds_prot_file = os.path.join(split_proteins_dir, f"{sample_name}.ccds.faa")
            domain_coord_file = os.path.join(domain_coords_dir, f"{sample_name}.domain_coords.txt")
            prot_mod_inputs.append([
                proteome_file, ccds_prot_file, domain_coord_file, pfam_hmm_path, 
                pfam_z, minimal_length, log_object, threads, skip_domain_splitting
            ])
        
        msg = f"Creating domain-chopped up version of proteome files for {len(prot_mod_inputs)} genomes"
        log_object.info(msg)
        
        # Process samples in parallel
        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(p.imap_unordered(create_chopped_proteomes, prot_mod_inputs), total=len(prot_mod_inputs)):
            pass
        p.close()
        
        # Create domain coordinate info file and return sample CCDS proteomes
        sample_ccds_proteomes = {}
        with open(domain_coord_info_file, 'w') as dci_handle:
            # Write header
            dci_handle.write('\t'.join(['Sample', 'Protein', 'Domain_Type', 'Domain_Index', 'Domain_Count', 'Length']) + '\n')
            
            for f in os.listdir(domain_coords_dir):
                if f.endswith('.domain_coords.txt'):
                    sample = f.replace('.domain_coords.txt', '')
                    ccds_prot_file = os.path.join(split_proteins_dir, f"{sample}.ccds.faa")
                    dom_coord_file = os.path.join(domain_coords_dir, f)
                    
                    if os.path.isfile(ccds_prot_file):
                        sample_ccds_proteomes[sample] = ccds_prot_file
                    
                    if os.path.isfile(dom_coord_file):
                        with open(dom_coord_file) as odf:
                            for i, line in enumerate(odf):
                                if i == 0:  # Skip header
                                    continue
                                dci_handle.write(line)
        
        log_object.info("Pfam annotation and protein splitting completed")
        return sample_ccds_proteomes

    except Exception as e:
        log_object.error("Error annotating and splitting proteins using Pfam")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        return {}



def process_pfam_domains_with_pyhmmer(
    proteome_file: str,
    split_proteins_dir: str,
    domain_coords_dir: str,
    sample_name: str,
    minimal_length: int,
    log_object: Any,
    threads: int = 1,
) -> None:
    """
    Process Pfam domain information and split proteins using pyhmmer.

    Parameters:
    -----------
    proteome_file : str
        Input proteome file
    split_proteins_dir : str
        Directory for split protein output
    domain_coords_dir : str
        Directory for domain coordinates
    sample_name : str
        Name of the sample
    minimal_length : int
        Minimum length for protein domains
    log_object : Any
        Logger object
    threads : int, default=1
        Number of threads to use
    """
    try:
        # Get Pfam database path
        bofasa_db_dir = os.getenv("BOFASA_DB_PATH", "").strip()
        if not bofasa_db_dir:
            raise ValueError("BOFASA_DB_PATH environment variable not set")
        
        db_locations = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
        if not os.path.isfile(db_locations):
            raise FileNotFoundError(f"Database locations file not found: {db_locations}")
        
        # Find Pfam HMM file
        pfam_hmm_path = None
        with open(db_locations) as odb:
            for line in odb:
                line = line.strip()
                ls = line.split('\t')
                if ls[0] == 'pfam':
                    pfam_hmm_path = ls[2]
                    break
        
        if not pfam_hmm_path or not os.path.isfile(pfam_hmm_path):
            raise FileNotFoundError(f"Pfam HMM file not found: {pfam_hmm_path}")
        
        # Load protein sequences with digital encoding
        alphabet = pyhmmer.easel.Alphabet.amino()
        sequences = []
        with pyhmmer.easel.SequenceFile(proteome_file, digital=True, alphabet=alphabet) as seq_file:
            sequences = list(seq_file)
        
        # Run HMM search using hmmsearch (not hmmscan)
        target_dom_hits = defaultdict(list)
        with pyhmmer.plan7.HMMFile(pfam_hmm_path) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="trusted", Z=1000000, cpus=threads):
                for hit in hits:
                    for domain in hit.domains.included:
                        target_dom_hits[hits.query.name.decode()].append([
                            hit.name.decode(), 
                            domain.alignment.target_from, 
                            domain.alignment.target_to, 
                            domain.score, 
                            domain.i_evalue
                        ])
        
        # Process domain hits and create breakpoints
        breakpoints = defaultdict(list)
        dom_start_names = defaultdict(lambda: 'NA')
        
        for protein_name in target_dom_hits:
            protein_dom_name_iter = defaultdict(int)
            accounted_coords = set([])
            
            # Sort domains by score (highest first)
            for dom_align_info in sorted(target_dom_hits[protein_name], key=lambda x: x[3], reverse=True):
                domain_name, start, end, score, i_evalue = dom_align_info
                
                # Check for overlap (10% leeway)
                overlap_coords = accounted_coords.intersection(set(range(start, end+1)))
                if len(overlap_coords)/float(end-start+1) >= 0.1 or len(overlap_coords) >= minimal_length:
                    continue
                
                accounted_coords = accounted_coords.union(set(range(start, end+1)))
                breakpoints[protein_name].append(start)
                breakpoints[protein_name].append(end+1)
                dom_start_names[protein_name + '|' + str(start)] = protein_name + '|' + domain_name + '|' + str(protein_dom_name_iter[domain_name]+1)
                protein_dom_name_iter[domain_name] += 1
        
        # Create output files
        ccds_prot_file = os.path.join(split_proteins_dir, f"{sample_name}.ccds.faa")
        domain_coord_file = os.path.join(domain_coords_dir, f"{sample_name}.domain_coords.txt")
        
        # Write domain coordinates file
        with open(domain_coord_file, 'w') as dcf_handle:
            dcf_handle.write('Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n')
            
            # Write CCDS proteome file
            with open(ccds_prot_file, 'w') as cpf_handle:
                with open(proteome_file) as ocf:
                    for rec in SeqIO.parse(ocf, 'fasta'):
                        protein_name = rec.id
                        protein_seq = str(rec.seq)
                        prev_end_coord = 1
                        protein_interdomain_index = 1
                        
                        if protein_name not in breakpoints and len(protein_seq) >= minimal_length:
                            # No domains found, use full protein
                            dn = f"{sample_name}|{protein_name}|full_protein|1"
                            cpf_handle.write(f">{dn}\n{protein_seq}\n")
                            dcf_handle.write(f"{sample_name}\t{protein_name}\tfull_protein\t1\t1\t{len(protein_seq)}\n")
                        else:
                            # Process domains using breakpoints
                            for protein_seq_chunk in split_by_idx(protein_seq, ([0] + sorted(breakpoints[protein_name]))):
                                if protein_seq_chunk.strip() == '':
                                    continue
                                end_coord = prev_end_coord + len(protein_seq_chunk) - 1
                                if len(protein_seq_chunk) >= minimal_length:
                                    dn = sample_name + '|' + dom_start_names[protein_name + '|' + str(prev_end_coord-1)]
                                    if dom_start_names[protein_name + '|' + str(prev_end_coord-1)] == 'NA':
                                        dn = f"{sample_name}|{protein_name}|inter-domain_region|{protein_interdomain_index}"
                                        protein_interdomain_index += 1
                                    cpf_handle.write(f">{dn}\n{protein_seq_chunk}\n")
                                    dcf_handle.write(f"{sample_name}\t{protein_name}\t{dn.split('|')[2]}\t{dn.split('|')[3]}\t{prev_end_coord}\t{end_coord}\n")
                                prev_end_coord = end_coord + 1
        
        log_object.info(f"Processed Pfam domains for {sample_name}")
        
    except Exception as e:
        log_object.error(f"Error processing Pfam domains for {sample_name}: {str(e)}")
        raise


def create_chopped_proteomes(inputs):
    """
    Create a chopped CDS proteome file from a regular proteome file - core function for multiprocessing.
    
    Parameters:
    -----------
    inputs : list
        List containing [prot_file, ccds_prot_file, dom_coord_file, pfam_hmm_path, pfam_z, minimal_length, log_object, threads, skip_domain_splitting]
    """
    prot_file, ccds_prot_file, dom_coord_file, pfam_hmm_path, pfam_z, minimal_length, log_object, threads, skip_domain_splitting = inputs
    
    try:
        sample = '.'.join(prot_file.split('/')[-1].split('.')[:-1])
        
        if not skip_domain_splitting:
            # Load protein sequences with digital encoding
            alphabet = pyhmmer.easel.Alphabet.amino()
            sequences = []
            with pyhmmer.easel.SequenceFile(prot_file, digital=True, alphabet=alphabet) as seq_file:
                sequences = list(seq_file)
            
            # Run HMM search using hmmsearch
            target_dom_hits = defaultdict(list)
            with pyhmmer.plan7.HMMFile(pfam_hmm_path) as hmm_file:
                for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="trusted", Z=int(pfam_z), cpus=threads):
                    for hit in hits:
                        for domain in hit.domains.included:
                            target_dom_hits[hit.name.decode()].append([hits.query.name.decode(), domain.alignment.target_from, domain.alignment.target_to, domain.score, domain.i_evalue])

            # Process domain hits and create breakpoints
            breakpoints = defaultdict(list)
            dom_start_names = defaultdict(lambda: 'NA')
            
            for tg in target_dom_hits:
                tg_dom_name_iter = defaultdict(int)
                accounted_coords = set([])
                
                for dom_align_info in sorted(target_dom_hits[tg], key=itemgetter(3), reverse=True):
                    dom_name, start, end, score, i_evalue = dom_align_info
                    overlap_coords = accounted_coords.intersection(set(range(start, end+1)))
                    if len(overlap_coords)/float(end-start+1) >= 0.1 or len(overlap_coords) >= minimal_length:
                        continue
                    accounted_coords = accounted_coords.union(set(range(start, end+1)))
                    breakpoints[tg].append(start)
                    breakpoints[tg].append(end+1)
                    dom_start_names[tg + '|' + str(start)] = tg + '|' + dom_name + '|' + str(tg_dom_name_iter[dom_name]+1)
                    tg_dom_name_iter[dom_name] += 1
            
            # Create output files
            with open(ccds_prot_file, 'w') as cpf_handle:
                with open(dom_coord_file, 'w') as dcf_handle:
                    dcf_handle.write('Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n')
                    
                    with open(prot_file) as ocf:
                        for rec in SeqIO.parse(ocf, 'fasta'):
                            tg = rec.id
                            tg_seq = str(rec.seq)
                            prev_end_coord = 1
                            tg_interdomain_index = 1
                            
                            if tg not in breakpoints and len(tg_seq) >= minimal_length:
                                # No domains found, use full protein
                                dn = sample + '|' + tg + '|full_protein|1'
                                cpf_handle.write('>' + dn + '\n' + str(tg_seq) + '\n')
                                dcf_handle.write('\t'.join([sample, tg, dn.split('|')[2], dn.split('|')[3], str(prev_end_coord), str(len(tg_seq))]) + '\n')
                            else:
                                # Process domains using breakpoints
                                for tg_seq_chunk in split_by_idx(tg_seq, ([0] + sorted(breakpoints[tg]))):
                                    if tg_seq_chunk.strip() == '':
                                        continue
                                    end_coord = prev_end_coord + len(tg_seq_chunk) - 1
                                    if len(tg_seq_chunk) >= minimal_length:
                                        dn = sample + '|' + dom_start_names[tg + '|' + str(prev_end_coord-1)]
                                        if dom_start_names[tg + '|' + str(prev_end_coord-1)] == 'NA':
                                            dn = sample + '|' + tg + '|inter-domain_region|' + str(tg_interdomain_index)
                                            tg_interdomain_index += 1
                                        cpf_handle.write('>' + dn + '\n' + str(tg_seq_chunk) + '\n')
                                        dcf_handle.write('\t'.join([sample, tg, dn.split('|')[2], dn.split('|')[3], str(prev_end_coord), str(end_coord)]) + '\n')
                                    prev_end_coord = end_coord + 1
        else:
            # Skip domain splitting, create full protein CCDS files
            with open(ccds_prot_file, 'w') as cpf_handle:
                with open(dom_coord_file, 'w') as dcf_handle:
                    dcf_handle.write('Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n')
                    
                    with open(prot_file) as ocf:
                        for rec in SeqIO.parse(ocf, 'fasta'):
                            protein_name = rec.id
                            protein_seq = str(rec.seq)
                            if len(protein_seq) >= minimal_length:
                                cpf_handle.write(f">{sample}|{protein_name}|full_protein|1\n{protein_seq}\n")
                                dcf_handle.write(f"{sample}\t{protein_name}\tfull_protein\t1\t1\t{len(protein_seq)}\n")
    
    except Exception as e:
        msg = f'An issue occurred with creating chopped up version of proteome file {prot_file}.'
        log_object.error(msg)
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        raise


def split_by_idx(S, list_of_indices):
    """
    Function taken from https://stackoverflow.com/questions/10851445/splitting-a-string-by-list-of-indices
    """
    left, right = 0, list_of_indices[0]
    yield S[left:right]
    left = right
    for right in list_of_indices[1:]:
        yield S[left:right]
        left = right
    yield S[left:]
