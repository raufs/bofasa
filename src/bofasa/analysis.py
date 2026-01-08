"""
Analysis functions for bofasa.

This module contains functions for analyzing ortholog groups, creating alignments,
generating reports, and performing various analyses on the results.
"""

from logging import warn
from math import log
import multiprocessing
import os
import platform
import subprocess
import sys
import traceback
from collections import defaultdict
from operator import itemgetter
import statistics
from typing import Any, Dict, List, Set
import pandas as pd
from Bio import SeqIO
import plotly.express as px
import pyhmmer
from .utils import _iter_progress, multi_process, setup_ready_directory, load_table_in_pandas_dataframe, run_cmd
from .processing import create_chopped_proteomes, extract_gene_contexts
from .core import generate_og_name
from . import config
from scipy import stats 

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

    This function loads coordinate information of CDSs for each input genome into a dictionary and also
    calculates context entropy scores and MGE annotations.

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
        # Load MGE annotation files
        isfinder_file = os.path.join(bofasa_prep_dir, 'Sample_IS_Element_Proteins.txt')
        plasmid_file = os.path.join(bofasa_prep_dir, 'Sample_Plasmid_Proteins.txt')
        phage_file = os.path.join(bofasa_prep_dir, 'Sample_Phage_Proteins.txt')

        ise_set: Set[str] = set([])
        plasmid_set: Set[str] = set([])
        phage_set: Set[str] = set([])

        genomad_flag: bool = False
        
        if os.path.isfile(isfinder_file):
            with open(isfinder_file) as oif:
                for line in oif:
                    line = line.strip()
                    ls = line.split('\t')
                    ise_set.add(ls[1])
        
        if os.path.isfile(plasmid_file):
            with open(plasmid_file) as opf:
                for line in opf:
                    line = line.strip()
                    ls = line.split('\t')
                    plasmid_set.add(ls[1])
            genomad_flag = True
        
        if os.path.isfile(phage_file):
            with open(phage_file) as opf:
                for line in opf:
                    line = line.strip()
                    ls = line.split('\t')
                    phage_set.add(ls[1])
            genomad_flag = True

        # Load ortholog group information
        og_genes: Dict[str, Set[str]] = defaultdict(set)
        og_samples: Dict[str, Set[str]] = defaultdict(set)
        og_is_single_copy: Dict[str, bool] = {}
        gene_to_og: Dict[str, str] = {}
        samples: List[str] = []

        if os.path.isfile(og_tsv_file):
            with open(og_tsv_file) as oot:
                for i, line in enumerate(oot):
                    line = line.strip('\n')
                    ls = line.split('\t')
                    if i == 0:
                        samples = ls[1:]
                    else:
                        og = ls[0]
                        is_single_copy = True
                        for j, gs in enumerate(ls[1:]):
                            sample = samples[j]
                            gene_count = 0
                            for g in gs.split(','):
                                g = g.strip()
                                if g != '':
                                    og_genes[og].add(g)
                                    og_samples[og].add(sample)
                                    gene_to_og[g] = og
                                    gene_count += 1
                            # If any sample has > 1 gene, it's not single-copy
                            if gene_count > 1:
                                is_single_copy = False
                        og_is_single_copy[og] = is_single_copy

        try:
            assert(len(og_genes) > 0)
        except Exception as e:
            msg = f"Difficulties parsing input orthogroup results in the file: {og_tsv_file}"
            log_object.error(msg)
            log_object.error(traceback.format_exc())
            sys.exit(1)

        # Load genome parameters
        listing_file = os.path.join(bofasa_prep_dir, 'Info_on_Input_Genome_Files.txt')
        assert(os.path.isfile(listing_file))
        
        genome_params: List[List[Any]] = []
        with open(listing_file) as olf:
            for i, line in enumerate(olf):
                if i == 0: 
                    continue
                line = line.strip()
                sample, ccds_proteome_file, proteome_file, coords_file, genome_file = line.split('\t')
                coords_file = os.path.join(bofasa_prep_dir, coords_file)
                output_file = os.path.join(surround_info_dir, f"{sample}.tsv")
                if os.path.isfile(coords_file):
                    genome_params.append([sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object])
                else:
                    msg = f"Coordinates file {coords_file} not found for {sample}."
                    log_object.error(msg)
                    raise RuntimeError(msg)

        # Process in parallel
        with multiprocessing.Pool(threads) as pool:
            list(
                _iter_progress(
                    pool.imap_unordered(extract_gene_contexts, genome_params),
                    total=len(genome_params),
                    description="Extracting gene contexts",
                )
            )

        # Process results and calculate entropy
        og_surrounding_nogs: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        og_completed_surrounding_nogs: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
        og_contexts: Dict[str, List[str]] = defaultdict(list)
        og_contexts_nses: Dict[str, int] = defaultdict(int)
        og_gene_lengths: Dict[str, List[float]] = defaultdict(list)
        context_nogs: Dict[str, List[int]] = defaultdict(list)
        context_nogs_complete: Dict[str, List[int]] = defaultdict(list)
        og_proteins: Dict[str, List[str]] = defaultdict(list)

        for genome_info in genome_params:
            sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object = genome_info
            with open(output_file) as oof:
                for i, line in enumerate(oof):
                    if i == 0: 
                        continue
                    line = line.strip('\n')
                    protein, og, nse, length, upstream_genes, downstream_genes, upstream_ogs, downstream_ogs = line.split('\t')
                    og_gene_lengths[og].append(float(length))
                    og_proteins[og].append(protein)
                    complete_status = 'C:'
                    if nse == 'True':
                        og_contexts_nses[og] += 1
                        complete_status = 'I:'
                    context = complete_status
                    if len(upstream_ogs) > 0:
                        context += upstream_ogs + ', '
                    context += og
                    if len(downstream_ogs) > 0:
                        context += ', ' + downstream_ogs
                    
                    og_contexts[og].append(context)
                    context_ogs = set([])
                    complete_context_ogs = set([])
                    for ogc in upstream_ogs.split(', '):
                        if nse == 'False':
                            complete_context_ogs.add(ogc)
                        if not ogc in context_ogs:
                            og_surrounding_nogs[og][ogc] += 1
                            if nse == 'False':
                                og_completed_surrounding_nogs[og][ogc] += 1
                        context_ogs.add(ogc)
                    for ogc in downstream_ogs.split(', '):
                        if nse == 'False':
                            complete_context_ogs.add(ogc)
                        if not ogc in context_ogs:
                            og_surrounding_nogs[og][ogc] += 1
                            if nse == 'False':
                                og_completed_surrounding_nogs[og][ogc] += 1
                        context_ogs.add(ogc)
                    context_nogs[og].append(len(context_ogs))
                    context_nogs_complete[og].append(len(complete_context_ogs))

        # Write detailed context information
        with open(og_context_info_file, 'w') as og_context_info_handle:
            og_context_info_handle.write('\t'.join([
                'OG', 'Median OG length (bp)', 'Percentage contexts near scaffold edge', 'Number of genomes with OG', 
                'Number of protein in OG', 'Single-copy', 'Context conservation score', 'Context conservation score - complete contexts', 
                'Number of distinct neighbor OGs', 
                'Number of distinct OGs from complete contexts', 'Avg. number of distinct neighbor OGs', 
                'Avg. number of distinct neighbor OGs from complete contexts', 
                'Percentage instances on plasmid (based on geNomad annotation)', 
                'Percentage instances on phage (based on geNomad annotation)', 
                'Percentage homologous to IS-elements (based on ISfinder database)', 'Instances', 'Contexts'
            ]) + '\n')
            
            for og in sorted(og_contexts):
                median_length = og_gene_lengths[og][0]
                if len(og_gene_lengths[og]) > 1:
                    median_length = statistics.median(og_gene_lengths[og])
                num_samples = len(og_samples[og])
                num_contexts = len(og_contexts[og])
                nse_perc = round(100.0*(og_contexts_nses[og]/float(num_contexts)), 2)
                
                nog_freqs = []
                total_nog = 0
                nog_freqs_complete = []
                total_nog_complete = 0
                
                for nog in og_surrounding_nogs[og]:
                    nog_freqs.append(og_surrounding_nogs[og][nog])
                    total_nog += 1
                    if nog in og_completed_surrounding_nogs[og]:
                        nog_freqs_complete.append(og_completed_surrounding_nogs[og][nog])
                        total_nog_complete += 1
                
                sum_nog_freqs = sum(nog_freqs)
                sum_nog_freqs_complete = sum(nog_freqs_complete)

                avg_nog = round(statistics.mean(context_nogs[og]), 2)
                avg_nog_complete = round(statistics.mean(context_nogs_complete[og]), 2)
                
                context_var_score = 'NA'
                context_var_score_complete = 'NA'
                
                if total_nog > 0:
                    context_var_score = round(avg_nog/total_nog, 2)
                
                if total_nog_complete > 0:
                    context_var_score_complete = round(avg_nog_complete/total_nog_complete, 2)
                
                # Calculate MGE percentages
                plasmid_count = 0
                phage_count = 0
                is_count = 0
                for p in og_proteins[og]:
                    if p in phage_set:
                        phage_count += 1
                    if p in plasmid_set:
                        plasmid_count += 1
                    if p in ise_set:
                        is_count += 1
                    
                plasmid_per = 100.0*(plasmid_count / float(num_contexts))
                phage_per = 100.0*(phage_count / float(num_contexts))
                ise_per = 100.0*(is_count / float(num_contexts))

                if not genomad_flag:
                    plasmid_per = 'NA'
                    phage_per = 'NA'
                
                # Determine single-copy status
                single_copy_status = 'Yes' if og_is_single_copy.get(og, False) else 'No'
                        
                og_context_info_handle.write('\t'.join([str(x) for x in [
                    og, round(median_length, 2), nse_perc, num_samples, num_contexts, single_copy_status, context_var_score, 
                    context_var_score_complete, total_nog, total_nog_complete,
                    avg_nog, avg_nog_complete, plasmid_per, phage_per, ise_per, '; '.join(og_proteins[og]), '; '.join(og_contexts[og])
                ]]) + '\n')

        # Create simplified plotting file
        log_object.info("Ortholog group contexts determined successfully")

    except Exception as e:
        log_object.error("Problem determining contexts of ortholog groups. Exiting now...")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)
    

def create_profile_hmms_and_consensus_seqs(
    prot_algn_dir: str, 
    phmm_dir: str, 
    cons_dir: str, 
    concatenate_consensus_faa: str,
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
        setup_ready_directory([phmm_dir, cons_dir], overwrite_mode="overwrite")

        # Process each alignment file
        alignment_files: List[str] = [f for f in os.listdir(prot_algn_dir) if f.endswith('.msa.faa')]
        
        hmmbuild_cmds = []
        hmmemit_cmds = []
        for algn_file in _iter_progress(alignment_files, description="Creating HMMs and consensus"):
            prefix: str = algn_file.replace('.msa.faa', '')
            algn_path: str = os.path.join(prot_algn_dir, algn_file)
            hmm_path: str = os.path.join(phmm_dir, f"{prefix}.hmm")
            cons_path: str = os.path.join(cons_dir, f"{prefix}.faa")

            # Create HMM using hmmbuild
            hmmbuild_cmd: List[str] = ['hmmbuild', '--amino', hmm_path, algn_path, log_object]
            hmmemit_cmd: List[str] = ['hmmemit', '-c', '-o', cons_path, hmm_path, log_object]
            hmmbuild_cmds.append(hmmbuild_cmd)
            hmmemit_cmds.append(hmmemit_cmd)

        # Run hmmbuild commands in parallel
        if hmmbuild_cmds:
            msg = f"Running {len(hmmbuild_cmds)} hmmbuild jobs"
            log_object.info(msg)
            
            p = multiprocessing.Pool(threads)
            try:
                for _ in _iter_progress(
                    p.imap_unordered(multi_process, hmmbuild_cmds),
                    total=len(hmmbuild_cmds),
                    description="Running hmmbuild",
                ):
                    pass
            except Exception as e:
                log_object.error("Error in hmmbuild multiprocessing")
                log_object.error(str(e))
                log_object.error(traceback.format_exc())
            finally:
                p.close()

        # Run hmmemit commands in parallel
        if hmmemit_cmds:
            msg = f"Running {len(hmmemit_cmds)} hmmemit jobs"
            log_object.info(msg)
            
            p = multiprocessing.Pool(threads)
            try:
                for _ in _iter_progress(
                    p.imap_unordered(multi_process, hmmemit_cmds),
                    total=len(hmmemit_cmds),
                    description="Running hmmemit",
                ):
                    pass
            except Exception as e:
                log_object.error("Error in hmmemit multiprocessing")
                log_object.error(str(e))
                log_object.error(traceback.format_exc())
            finally:
                p.close()

        log_object.info("Profile HMMs and consensus sequences created successfully")

        concatenate_consensus_sequences(cons_dir, concatenate_consensus_faa, log_object)

    except Exception as e:
        log_object.error("Error creating profile HMMs and consensus sequences")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())

def concatenate_consensus_sequences(
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





def create_final_report(
    bofasa_prep_dir: str,
    og_context_info_file: str,
    final_result_file: str,
    log_object: Any,
) -> None:
    """
    Create the formatted Excel report summarizing BOFASA results.

    Parameters:
    - bofasa_prep_dir: Directory containing bofasa prep results
    - og_context_info_file: Path to OG context information TSV
    - final_result_file: Output .xlsx file path
    - log_object: Logger instance
    """
    try:
        # Detect presence of geNomad annotations to decide formatting of certain columns
        plasmid_file: str = os.path.join(bofasa_prep_dir, "Sample_Plasmid_Proteins.txt")
        phage_file: str = os.path.join(bofasa_prep_dir, "Sample_Phage_Proteins.txt")
        genomad_flag: bool = os.path.isfile(plasmid_file) or os.path.isfile(phage_file)

        # Columns to treat as numeric in the context table
        numeric_columns: List[str] = [
            "Median OG length (bp)",
            "Percentage contexts near scaffold edge",
            "Number of genomes with OG",
            "Number of protein in OG",
            "Context conservation score",
            "Context conservation score - complete contexts",
            "Number of distinct neighbor OGs",
            "Number of distinct OGs from complete contexts",
            "Avg. number of distinct neighbor OGs",
            "Avg. number of distinct neighbor OGs from complete contexts",
            "Percentage instances on plasmid (based on geNomad annotation)",
            "Percentage instances on phage (based on geNomad annotation)",
            "Percentage homologous to IS-elements (based on ISfinder database)",
        ]

        # Load the context table, dropping the last two columns (Instances, Contexts)
        results_df: pd.DataFrame = load_table_in_pandas_dataframe(
            og_context_info_file,
            numeric_columns=numeric_columns,
            cut_last_columns=2,
        )

        # Create Excel writer and add data dictionary sheet
        writer = pd.ExcelWriter(final_result_file, engine="xlsxwriter")
        workbook = writer.book
        dd_sheet = workbook.add_worksheet("Data Dictionary")
        dd_sheet.write(
            0,
            0,
            "Data Dictionary describing columns of \"bofasa Results\" spreadsheet can be found below and on bofasa's Wiki page at:",
        )
        dd_sheet.write(1, 0, "https://github.com/raufs/bofasa")

        # Formats
        na_format = workbook.add_format({"font_color": "#a6a6a6", "bg_color": "#FFFFFF", "italic": True})
        header_format = workbook.add_format(
            {"bold": True, "text_wrap": True, "valign": "top", "fg_color": "#D7E4BC", "border": 1}
        )

        # Write main sheet
        sheet_name = "bofasa Results"
        results_df.to_excel(writer, sheet_name=sheet_name, index=False, na_rep="NA")
        worksheet = writer.sheets[sheet_name]
        num_rows: int = results_df.shape[0] + 1  # +1 for header row

        # Apply conditional formatting for NA cells and header row
        worksheet.conditional_format(
            f"A2:BA{num_rows}", {"type": "cell", "criteria": "==", "value": '"NA"', "format": na_format}
        )
        worksheet.conditional_format("A1:BA1", {"type": "cell", "criteria": "!=", "value": "NA", "format": header_format})

        # Compute max values for color scales from numeric columns that exist in the DataFrame
        def max_or_zero(col: str) -> float:
            if col in results_df.columns:
                series = results_df[col]
                try:
                    return float(series.max(skipna=True)) if series.size > 0 else 0.0
                except Exception:
                    return 0.0
            return 0.0

        max_num_genomes = max_or_zero("Number of genomes with OG")
        max_num_proteins = max_or_zero("Number of protein in OG")
        max_context_var = max_or_zero("Context conservation score")
        max_context_var_comp = max_or_zero("Context conservation score - complete contexts")
        # Entropy score disabled - see TODO in determine_ortholog_group_contexts()
        # max_context_ent = max_or_zero("Context entropy score")
        # max_context_ent_comp = max_or_zero("Context entropy score - complete contexts")

        # Column color scales (ranges mirror v1.1.1)
        worksheet.conditional_format(
            f"B2:B{num_rows}",
            {"type": "2_color_scale", "min_color": "#a9cafc", "max_color": "#736991", "min_value": 0, "max_value": 2500, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"C2:C{num_rows}",
            {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"D2:D{num_rows}",
            {"type": "2_color_scale", "min_color": "#f2c6f7", "max_color": "#b07fb5", "min_value": 0.0, "max_value": max_num_genomes, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"E2:E{num_rows}",
            {"type": "2_color_scale", "min_color": "#e7cdf7", "max_color": "#a186b3", "min_value": 0.0, "max_value": max_num_proteins, "min_type": "num", "max_type": "num"},
        )
        # F is Single-copy (Yes/No) - no color scale needed
        
        worksheet.conditional_format(
            f"G2:G{num_rows}",
            {"type": "2_color_scale", "min_color": "#e6f5ab", "max_color": "#a4b36b", "min_value": 0.0, "max_value": max_context_var, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"H2:H{num_rows}",
            {"type": "2_color_scale", "min_color": "#e6f5ab", "max_color": "#a4b36b", "min_value": 0.0, "max_value": max_context_var_comp, "min_type": "num", "max_type": "num"},
        )

        if genomad_flag:
            worksheet.conditional_format(
                f"M2:M{num_rows}",
                {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 100.0, "min_type": "num", "max_type": "num"},
            )
            worksheet.conditional_format(
                f"N2:N{num_rows}",
                {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 100.0, "min_type": "num", "max_type": "num"},
            )

        worksheet.conditional_format(
            f"O2:O{num_rows}",
            {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 100.0, "min_type": "num", "max_type": "num"},
        )

        # Autofilter and finalize
        worksheet.autofilter(f"A1:BA{num_rows}")
        workbook.close()
        if log_object:
            log_object.info("Final Excel report created successfully")
    except Exception as e:
        if log_object:
            log_object.error("Error creating final Excel report")
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
        isfinder_file = os.path.join(bofasa_prep_dir, 'Sample_IS_Element_Proteins.txt')
        plasmid_file = os.path.join(bofasa_prep_dir, 'Sample_Plasmid_Proteins.txt')
        phage_file = os.path.join(bofasa_prep_dir, 'Sample_Phage_Proteins.txt')

        mge_set = set([])

        if os.path.isfile(isfinder_file):
            with open(isfinder_file) as oif:
                for line in oif:
                    line = line.strip()
                    ls = line.split('\t')
                    mge_set.add(ls[1])
        
        if os.path.isfile(plasmid_file):
            with open(plasmid_file) as opf:
                for line in opf:
                    line = line.strip()
                    ls = line.split('\t')
                    mge_set.add(ls[1])
        
        if os.path.isfile(phage_file):
            with open(phage_file) as opf:
                for line in opf:
                    line = line.strip()
                    ls = line.split('\t')
                    mge_set.add(ls[1])
            
        tmp_file_header = [
            'OG', 
            'Number of protein in OG', 
            'Context conservation score', 
            'Majority of protein instances homologous to IS-element or on plasmid or phage'
        ]
        outf_handle = open(tmp_result_file, 'w')
        outf_handle.write('\t'.join(tmp_file_header) + '\n')
        with open(og_context_info_file) as oocif:
            for i, line in enumerate(oocif):
                if i == 0: 
                    continue
                line = line.strip()
                ls = line.split('\t')
                mge_related = 'No'
                tot = 0
                mge = 0
                for p in ls[-2].split('; '):
                    tot += 1
                    if p in mge_set:
                        mge += 1
                if mge/tot > 0.5:
                    mge_related = 'Yes'
                outf_handle.write('\t'.join([ls[0], ls[4], ls[5], mge_related]) + '\n')
        outf_handle.close()
                
        numeric_columns = set(['Number of protein in OG', 'Context conservation score'])
        simple_df = load_table_in_pandas_dataframe(tmp_result_file, numeric_columns)
        fig = px.scatter(
            simple_df, 
            x="Number of protein in OG", 
            y="Context conservation score",
            color="Majority of protein instances homologous to IS-element or on plasmid or phage", 
            color_discrete_map={"No": "grey", "Yes": "red"},
            opacity=0.4,
            marginal_x="histogram", 
            marginal_y="histogram"
        )
        fig.write_html(final_result_plot)
    except Exception as e:
        msg = 'Issues with create final HTML report.'
        log_object.error(msg)
        log_object.error(str(e))
        log_object.error(traceback.format_exc())


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
        setup_ready_directory([genomad_dir], overwrite_mode="overwrite")
        
        # Open output files
        with open(phage_protein_listing_file, 'w') as phpf_handle, open(plasmid_protein_listing_file, 'w') as plpf_handle:
            
            # Process each sample
            for sample in sample_wgs:
                input_genome = sample_wgs[sample]
                genomad_results = os.path.join(genomad_dir, sample)
                
                # Run genomad
                # Build base command
                genomad_cmd = ['genomad', 'end-to-end', '--cleanup']
                
                # On macOS, disable neural network classification to avoid TensorFlow crashes
                # (especially on Apple Silicon). On Linux, neural network classification works fine.
                if platform.system() == 'Darwin':
                    genomad_cmd.append('--disable-nn-classification')
                
                # Add remaining parameters
                genomad_cmd.extend([
                    '--threads', str(threads), '--splits', str(genome_splits), 
                    input_genome, genomad_results, genomad_db_dir
                ])
                
                try:
                    run_cmd(genomad_cmd, log_object)
                except subprocess.CalledProcessError as e:
                    log_object.error(f"geNomad failed with unknown error for sample {sample}")
                    log_object.error(f"Command: {' '.join(genomad_cmd)}")
                    log_object.error(f"Full error output: {e.output}")
                    raise RuntimeError(f"geNomad failed with unknown error for sample {sample}")
                
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


def extract_mge_genome_sequences(
    genomad_dir: str,
    sample_wgs: Dict[str, str],
    output_dir: str,
    log_object: Any,
) -> Dict[str, Dict[str, str]]:
    """
    Extract plasmid and phage genome sequences from genomad results as separate genome files.
    
    Parameters:
    -----------
    genomad_dir : str
        Directory containing genomad results
    sample_wgs : Dict[str, str]
        Dictionary mapping sample names to genome file paths
    output_dir : str
        Output directory for extracted MGE genomes
    log_object : Any
        Logger object
        
    Returns:
    --------
    Dict[str, Dict[str, str]]
        Dictionary with keys 'plasmids' and 'phages', each containing a dict mapping
        MGE identifiers to their FASTA file paths
    """
    from Bio import SeqIO
    import os
    
    try:
        mge_genomes = {
            'plasmids': {},
            'phages': {}
        }
        
        # Create output directories
        plasmid_genomes_dir = os.path.join(output_dir, "Extracted_Plasmid_Genomes/")
        phage_genomes_dir = os.path.join(output_dir, "Extracted_Phage_Genomes/")
        setup_ready_directory([plasmid_genomes_dir, phage_genomes_dir], overwrite_mode="overwrite")
        
        msg = "Extracting plasmid and phage genome sequences from geNomad results..."
        log_object.info(msg)
        
        total_plasmids = 0
        total_phages = 0
        
        # Process each sample
        for sample in sample_wgs:
            genomad_results = os.path.join(genomad_dir, sample)
            
            if not os.path.isdir(genomad_results):
                log_object.warning(f"geNomad results directory not found for {sample}: {genomad_results}")
                continue
            
            # Find genomad output files
            prophage_summary_tsv = None
            plasmid_summary_tsv = None
            plasmid_fasta = None
            virus_fasta = None
            
            for subdir, dirs, files in os.walk(genomad_results):
                for file in files:
                    filepath = os.path.join(subdir, file)
                    if filepath.endswith("_plasmid_summary.tsv"):
                        plasmid_summary_tsv = filepath
                    elif filepath.endswith("_virus_summary.tsv"):
                        prophage_summary_tsv = filepath
                    elif filepath.endswith("_plasmid.fna"):
                        plasmid_fasta = filepath
                    elif filepath.endswith("_virus.fna"):
                        virus_fasta = filepath
            
            # Extract plasmids
            if plasmid_summary_tsv and plasmid_fasta and os.path.isfile(plasmid_fasta):
                try:
                    plasmid_count = 0
                    with open(plasmid_fasta, 'r') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            # Create unique identifier for this plasmid
                            plasmid_id = f"{sample}_plasmid_{plasmid_count + 1}"
                            output_file = os.path.join(plasmid_genomes_dir, f"{plasmid_id}.fna")
                            
                            # Write plasmid sequence
                            with open(output_file, 'w') as out_handle:
                                out_handle.write(f">{plasmid_id} {record.description}\n")
                                out_handle.write(f"{str(record.seq)}\n")
                            
                            mge_genomes['plasmids'][plasmid_id] = output_file
                            plasmid_count += 1
                            total_plasmids += 1
                    
                    if plasmid_count > 0:
                        log_object.info(f"Extracted {plasmid_count} plasmid(s) from {sample}")
                except Exception as e:
                    log_object.warning(f"Could not extract plasmids from {sample}: {str(e)}")
            
            # Extract phages/viruses
            if prophage_summary_tsv and virus_fasta and os.path.isfile(virus_fasta):
                try:
                    phage_count = 0
                    with open(virus_fasta, 'r') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            # Create unique identifier for this phage
                            phage_id = f"{sample}_phage_{phage_count + 1}"
                            output_file = os.path.join(phage_genomes_dir, f"{phage_id}.fna")
                            
                            # Write phage sequence
                            with open(output_file, 'w') as out_handle:
                                out_handle.write(f">{phage_id} {record.description}\n")
                                out_handle.write(f"{str(record.seq)}\n")
                            
                            mge_genomes['phages'][phage_id] = output_file
                            phage_count += 1
                            total_phages += 1
                    
                    if phage_count > 0:
                        log_object.info(f"Extracted {phage_count} phage(s) from {sample}")
                except Exception as e:
                    log_object.warning(f"Could not extract phages from {sample}: {str(e)}")
        
        msg = f"Extraction complete: {total_plasmids} plasmids and {total_phages} phages extracted as separate genomes"
        log_object.info(msg)
        
        return mge_genomes
        
    except Exception as e:
        log_object.error(f"Error extracting MGE genome sequences: {str(e)}")
        log_object.error(traceback.format_exc())
        raise


def extract_mge_annotations_from_parent_genomes(
    genomad_dir: str,
    sample_wgs: Dict[str, str],
    sample_proteomes: Dict[str, str],
    sample_beds: Dict[str, str],
    output_dir: str,
    log_object: Any,
) -> Dict[str, Dict[str, str]]:
    """
    Extract gene annotations for MGEs from their parent genomes and remove them from parent files.
    
    This function extracts phage and plasmid sequences along with their gene annotations
    from the parent genome's existing gene calls, avoiding the need to re-run gene calling
    on small MGE sequences which often fail due to size constraints.
    
    IMPORTANT: This function also filters out MGE genes from the parent genome files to
    avoid double-counting when MGEs are treated as separate entities.
    
    Parameters:
    -----------
    genomad_dir : str
        Directory containing genomad results
    sample_wgs : Dict[str, str]
        Dictionary mapping sample names to genome file paths
    sample_proteomes : Dict[str, str]
        Dictionary mapping sample names to proteome file paths
    sample_beds : Dict[str, str]
        Dictionary mapping sample names to BED coordinate file paths
    output_dir : str
        Output directory for extracted MGE genomes and annotations
    log_object : Any
        Logger object
        
    Returns:
    --------
    Dict[str, Dict[str, str]]
        Dictionary with keys 'mge_wgs', 'mge_proteomes', 'mge_beds' containing paths to
        MGE genome, proteome, and BED files
    """
    from Bio import SeqIO
    import os
    
    try:
        mge_data = {
            'mge_wgs': {},
            'mge_proteomes': {},
            'mge_beds': {}
        }
        
        # Create output directories
        plasmid_genomes_dir = os.path.join(output_dir, "Extracted_Plasmid_Genomes/")
        phage_genomes_dir = os.path.join(output_dir, "Extracted_Phage_Genomes/")
        setup_ready_directory([plasmid_genomes_dir, phage_genomes_dir], overwrite_mode="overwrite")
        
        # Also create directories in Genome_Processing
        gp_dir = os.path.join(output_dir, "Genome_Processing/")
        wgs_dir = os.path.join(gp_dir, "Genomes/")
        faa_dir = os.path.join(gp_dir, "Proteomes/")
        bed_dir = os.path.join(gp_dir, "BEDs/")
        
        msg = "Extracting plasmid and phage sequences with annotations from parent genomes..."
        log_object.info(msg)
        
        total_plasmids = 0
        total_phages = 0
        
        # Track MGE genes per sample for filtering parent genomes
        sample_mge_genes = {}  # sample -> set of gene_ids that belong to MGEs
        
        # Process each sample
        for sample in sample_wgs:
            genomad_results = os.path.join(genomad_dir, sample)
            
            if not os.path.isdir(genomad_results):
                log_object.warning(f"geNomad results directory not found for {sample}: {genomad_results}")
                continue
            
            # Initialize MGE gene tracking for this sample
            sample_mge_genes[sample] = set()
            
            # Get parent genome files
            # Note: sample_proteomes and sample_beds already contain full paths
            parent_proteome = sample_proteomes[sample]
            parent_bed = sample_beds[sample]
            
            if not os.path.isfile(parent_proteome) or not os.path.isfile(parent_bed):
                log_object.warning(f"Parent genome annotations not found for {sample}")
                log_object.warning(f"  Proteome path: {parent_proteome}")
                log_object.warning(f"  BED path: {parent_bed}")
                continue
            
            # Load parent genome proteins
            parent_proteins = {}
            with open(parent_proteome, 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    parent_proteins[record.id] = record
            
            # Load parent genome gene coordinates
            parent_genes = {}  # gene_id -> (scaffold, start, end, score, strand)
            with open(parent_bed) as obf:
                for line in obf:
                    line = line.strip()
                    if not line:
                        continue
                    ls = line.split('\t')
                    scaffold = ls[0]
                    start = int(ls[1])
                    end = int(ls[2])
                    gene_id = ls[3]
                    score = ls[4] if len(ls) > 4 else '1'
                    strand = ls[5] if len(ls) > 5 else '+'
                    parent_genes[gene_id] = (scaffold, start, end, score, strand)
            
            # Find genomad output files
            plasmid_fasta = None
            virus_fasta = None
            
            for subdir, dirs, files in os.walk(genomad_results):
                for file in files:
                    filepath = os.path.join(subdir, file)
                    if filepath.endswith("_plasmid.fna"):
                        plasmid_fasta = filepath
                    elif filepath.endswith("_virus.fna"):
                        virus_fasta = filepath
            
            # Process plasmids
            if plasmid_fasta and os.path.isfile(plasmid_fasta):
                try:
                    plasmid_count = 0
                    with open(plasmid_fasta, 'r') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            plasmid_id = f"{sample}_plasmid_{plasmid_count + 1}"
                            seq_name = record.id
                            
                            # Write plasmid genome sequence
                            output_fna = os.path.join(wgs_dir, f"{plasmid_id}.fna")
                            with open(output_fna, 'w') as out_handle:
                                out_handle.write(f">{plasmid_id} {record.description}\n")
                                out_handle.write(f"{str(record.seq)}\n")
                            
                            # Extract genes for this plasmid
                            plasmid_genes = []
                            for gene_id, (gene_scaffold, gene_start, gene_end, gene_score, gene_strand) in parent_genes.items():
                                if gene_scaffold == seq_name:
                                    plasmid_genes.append(gene_id)
                                    # Track that this gene belongs to an MGE
                                    sample_mge_genes[sample].add(gene_id)
                            
                            # Debug: warn if no genes found
                            if len(plasmid_genes) == 0:
                                log_object.warning(f"No genes found for plasmid {plasmid_id} (scaffold: {seq_name})")
                                log_object.warning(f"  Available scaffolds in parent BED: {list(set([g[0] for g in parent_genes.values()]))[:10]}")
                                log_object.warning(f"  This may indicate a scaffold name mismatch between geNomad output and parent genome")
                                # Skip creating empty files
                                continue
                            
                            # Write proteome and BED files
                            output_faa = os.path.join(faa_dir, f"{plasmid_id}.faa")
                            output_bed = os.path.join(bed_dir, f"{plasmid_id}.bed")
                            
                            with open(output_faa, 'w') as faa_handle:
                                for gene_id in plasmid_genes:
                                    if gene_id in parent_proteins:
                                        prot_record = parent_proteins[gene_id]
                                        faa_handle.write(f">{prot_record.id} {prot_record.description}\n")
                                        faa_handle.write(f"{str(prot_record.seq)}\n")
                            
                            with open(output_bed, 'w') as bed_handle:
                                for gene_id in plasmid_genes:
                                    if gene_id in parent_genes:
                                        gene_scaffold, gene_start, gene_end, gene_score, gene_strand = parent_genes[gene_id]
                                        bed_handle.write(f"{gene_scaffold}\t{gene_start}\t{gene_end}\t{gene_id}\t{gene_score}\t{gene_strand}\n")
                            
                            # Track MGE data with absolute paths
                            mge_data['mge_wgs'][plasmid_id] = os.path.join(output_dir, "Genome_Processing/Genomes", f"{plasmid_id}.fna")
                            mge_data['mge_proteomes'][plasmid_id] = os.path.join(output_dir, "Genome_Processing/Proteomes", f"{plasmid_id}.faa")
                            mge_data['mge_beds'][plasmid_id] = os.path.join(output_dir, "Genome_Processing/BEDs", f"{plasmid_id}.bed")
                            
                            plasmid_count += 1
                            total_plasmids += 1
                    
                    if plasmid_count > 0:
                        log_object.info(f"Extracted {plasmid_count} plasmid(s) from {sample}")
                except Exception as e:
                    log_object.warning(f"Could not extract plasmids from {sample}: {str(e)}")
                    log_object.warning(traceback.format_exc())
            
            # Process phages/viruses
            if virus_fasta and os.path.isfile(virus_fasta):
                try:
                    phage_count = 0
                    with open(virus_fasta, 'r') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            phage_id = f"{sample}_phage_{phage_count + 1}"
                            seq_name = record.id
                            
                            # Write phage genome sequence
                            output_fna = os.path.join(wgs_dir, f"{phage_id}.fna")
                            with open(output_fna, 'w') as out_handle:
                                out_handle.write(f">{phage_id} {record.description}\n")
                                out_handle.write(f"{str(record.seq)}\n")
                            
                            # Parse coordinates from geNomad FASTA header
                            # geNomad provirus format: scaffold|provirus_START_END
                            # geNomad full virus format: scaffold_name
                            if '|provirus_' in seq_name:
                                # Provirus - extract scaffold and coordinates
                                scaffold = seq_name.split('|')[0]
                                coord_part = seq_name.split('|')[1].replace('provirus_', '')
                                coord_parts = coord_part.split('_')
                                v_start = int(coord_parts[0])
                                v_end = int(coord_parts[1])
                            else:
                                # Full viral scaffold (not a provirus)
                                scaffold = seq_name
                                v_start = None
                                v_end = None
                            
                            # Extract genes for this phage
                            phage_genes = []
                            for gene_id, (gene_scaffold, gene_start, gene_end, gene_score, gene_strand) in parent_genes.items():
                                if gene_scaffold == scaffold:
                                    # If we have specific coordinates, check overlap
                                    if v_start is not None and v_end is not None:
                                        # Gene overlaps if it starts or ends within the virus region
                                        if (v_start <= gene_start <= v_end) or (v_start <= gene_end <= v_end):
                                            phage_genes.append(gene_id)
                                            # Track that this gene belongs to an MGE
                                            sample_mge_genes[sample].add(gene_id)
                                    else:
                                        # Full scaffold - include all genes from that scaffold
                                        phage_genes.append(gene_id)
                                        # Track that this gene belongs to an MGE
                                        sample_mge_genes[sample].add(gene_id)
                            
                            # Debug: warn if no genes found
                            if len(phage_genes) == 0:
                                log_object.warning(f"No genes found for phage {phage_id}")
                                log_object.warning(f"  geNomad sequence ID: {seq_name}")
                                log_object.warning(f"  Parsed scaffold: {scaffold}")
                                if v_start is not None:
                                    log_object.warning(f"  Parsed coordinates: {v_start}-{v_end}")
                                    # Count genes in this region
                                    genes_in_region = 0
                                    for gene_id, (gene_scaffold, gene_start, gene_end, gene_score, gene_strand) in parent_genes.items():
                                        if gene_scaffold == scaffold:
                                            if (v_start <= gene_start <= v_end) or (v_start <= gene_end <= v_end):
                                                genes_in_region += 1
                                    log_object.warning(f"  Genes found in region {v_start}-{v_end} on {scaffold}: {genes_in_region}")
                                log_object.warning(f"  Available scaffolds in parent BED: {list(set([g[0] for g in parent_genes.values()]))[:10]}")
                                # Skip creating empty files
                                continue
                            
                            # Write proteome and BED files
                            output_faa = os.path.join(faa_dir, f"{phage_id}.faa")
                            output_bed = os.path.join(bed_dir, f"{phage_id}.bed")
                            
                            with open(output_faa, 'w') as faa_handle:
                                for gene_id in phage_genes:
                                    if gene_id in parent_proteins:
                                        prot_record = parent_proteins[gene_id]
                                        faa_handle.write(f">{prot_record.id} {prot_record.description}\n")
                                        faa_handle.write(f"{str(prot_record.seq)}\n")
                            
                            with open(output_bed, 'w') as bed_handle:
                                for gene_id in phage_genes:
                                    if gene_id in parent_genes:
                                        gene_scaffold, gene_start, gene_end, gene_score, gene_strand = parent_genes[gene_id]
                                        bed_handle.write(f"{gene_scaffold}\t{gene_start}\t{gene_end}\t{gene_id}\t{gene_score}\t{gene_strand}\n")
                            
                            # Track MGE data with absolute paths
                            mge_data['mge_wgs'][phage_id] = os.path.join(output_dir, "Genome_Processing/Genomes", f"{phage_id}.fna")
                            mge_data['mge_proteomes'][phage_id] = os.path.join(output_dir, "Genome_Processing/Proteomes", f"{phage_id}.faa")
                            mge_data['mge_beds'][phage_id] = os.path.join(output_dir, "Genome_Processing/BEDs", f"{phage_id}.bed")
                            
                            phage_count += 1
                            total_phages += 1
                    
                    if phage_count > 0:
                        log_object.info(f"Extracted {phage_count} phage(s) from {sample}")
                except Exception as e:
                    log_object.warning(f"Could not extract phages from {sample}: {str(e)}")
                    log_object.warning(traceback.format_exc())
        
        msg = f"Extraction complete: {total_plasmids} plasmids and {total_phages} phages extracted with annotations from parent genomes"
        log_object.info(msg)
        
        # Now filter MGE genes from parent genome files
        log_object.info("Filtering MGE genes from parent genome files to avoid double-counting...")
        total_filtered_genes = 0
        
        for sample, mge_gene_set in sample_mge_genes.items():
            if len(mge_gene_set) == 0:
                continue
            
            try:
                # Note: sample_proteomes and sample_beds already contain full paths
                parent_proteome = sample_proteomes[sample]
                parent_bed = sample_beds[sample]
                
                if not os.path.isfile(parent_proteome) or not os.path.isfile(parent_bed):
                    log_object.warning(f"Cannot filter parent genome for {sample}: files not found")
                    continue
                
                # Read parent proteome and filter out MGE genes
                filtered_proteins = []
                with open(parent_proteome, 'r') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        if record.id not in mge_gene_set:
                            filtered_proteins.append(record)
                
                # Write filtered proteome
                with open(parent_proteome, 'w') as out_handle:
                    SeqIO.write(filtered_proteins, out_handle, 'fasta')
                
                # Read parent BED and filter out MGE genes
                filtered_bed_lines = []
                with open(parent_bed, 'r') as bed_handle:
                    for line in bed_handle:
                        line = line.strip()
                        if not line:
                            continue
                        ls = line.split('\t')
                        gene_id = ls[3]
                        if gene_id not in mge_gene_set:
                            filtered_bed_lines.append(line)
                
                # Write filtered BED file
                with open(parent_bed, 'w') as bed_handle:
                    for line in filtered_bed_lines:
                        bed_handle.write(line + '\n')
                
                genes_filtered = len(mge_gene_set)
                total_filtered_genes += genes_filtered
                log_object.info(f"Filtered {genes_filtered} MGE gene(s) from {sample}")
                
            except Exception as e:
                log_object.warning(f"Could not filter MGE genes from {sample}: {str(e)}")
                log_object.warning(traceback.format_exc())
        
        if total_filtered_genes > 0:
            log_object.info(f"Total MGE genes filtered from parent genomes: {total_filtered_genes}")
        
        return mge_data
        
    except Exception as e:
        log_object.error(f"Error extracting MGE annotations from parent genomes: {str(e)}")
        log_object.error(traceback.format_exc())
        raise


def update_mge_protein_listings_from_proteomes(
    output_dir: str,
    phage_protein_listing_file: str,
    plasmid_protein_listing_file: str,
    log_object: Any,
) -> None:
    """
    Update MGE protein listing files by reading extracted MGE proteome files.
    
    After MGE extraction, this reads all proteome files with _phage_ or _plasmid_
    in their names and rewrites the listing files accordingly.
    
    Parameters:
    -----------
    output_dir : str
        Output directory containing Genome_Processing/Proteomes/
    phage_protein_listing_file : str
        Path to Sample_Phage_Proteins.txt
    plasmid_protein_listing_file : str
        Path to Sample_Plasmid_Proteins.txt
    log_object : Any
        Logger object
    """
    from Bio import SeqIO
    
    try:
        proteomes_dir = os.path.join(output_dir, "Genome_Processing/Proteomes/")
        
        if not os.path.isdir(proteomes_dir):
            log_object.warning(f"Proteomes directory not found: {proteomes_dir}")
            return
        
        # Collect phage proteins
        phage_entries = []
        plasmid_entries = []
        
        for filename in os.listdir(proteomes_dir):
            if not filename.endswith('.faa'):
                continue
            
            sample_name = filename.replace('.faa', '')
            proteome_path = os.path.join(proteomes_dir, filename)
            
            # Check if this is a phage or plasmid MGE
            if '_phage_' in sample_name:
                # Read all protein IDs from this phage proteome
                with open(proteome_path, 'r') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        phage_entries.append(f"{sample_name}\t{record.id}\n")
            
            elif '_plasmid_' in sample_name:
                # Read all protein IDs from this plasmid proteome
                with open(proteome_path, 'r') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        plasmid_entries.append(f"{sample_name}\t{record.id}\n")
        
        # Rewrite phage listing file
        if phage_entries:
            with open(phage_protein_listing_file, 'w') as f:
                f.writelines(phage_entries)
            log_object.info(f"Updated {phage_protein_listing_file} with {len(phage_entries)} entries from extracted MGE proteomes")
        
        # Rewrite plasmid listing file
        if plasmid_entries:
            with open(plasmid_protein_listing_file, 'w') as f:
                f.writelines(plasmid_entries)
            log_object.info(f"Updated {plasmid_protein_listing_file} with {len(plasmid_entries)} entries from extracted MGE proteomes")
    
    except Exception as e:
        log_object.error(f"Error updating MGE protein listings: {str(e)}")
        log_object.error(traceback.format_exc())
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
        setup_ready_directory([annot_dir], overwrite_mode="overwrite")
        
        # Prepare DIAMOND commands for each sample
        # Skip empty proteome files (e.g., from MGEs with no extractable genes)
        dmnd_search_cmds = []
        skipped_samples = []
        for sample_name, faa_file in sample_proteomes.items():
            # Check if file exists and is not empty
            if not os.path.isfile(faa_file) or os.path.getsize(faa_file) == 0:
                skipped_samples.append(sample_name)
                continue
            
            annotation_result_file = os.path.join(annot_dir, f"{sample_name}.isfinder_diamond_blastp.txt")
            search_cmd = [
                'diamond', 'blastp', '--ignore-warnings', '-p', str(1), 
                '-d', isfinder_dmnd_path, '-q', faa_file, '-o', annotation_result_file
            ]
            dmnd_search_cmds.append(search_cmd + [log_object])
        
        if skipped_samples:
            log_object.warning(f"Skipping {len(skipped_samples)} samples with empty proteome files: {', '.join(skipped_samples[:5])}{'...' if len(skipped_samples) > 5 else ''}")
        
        msg = f"Running {len(dmnd_search_cmds)} DIAMOND blastp jobs for IS element annotation"
        log_object.info(msg)
        
        # Run DIAMOND searches in parallel
        from .utils import multi_process
        p = multiprocessing.Pool(threads)
        for _ in _iter_progress(
            p.imap_unordered(multi_process, dmnd_search_cmds),
            total=len(dmnd_search_cmds),
            description="Running DIAMOND blastp",
        ):
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
    sample_is_mge: Set[str] = None,
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
    sample_is_mge : Set[str], optional
        Set of sample names that are MGEs (plasmids/phages)
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
        
        # Create output directories (all files in common directories)
        setup_ready_directory([split_proteins_dir, domain_coords_dir], overwrite_mode="overwrite")
        
        # Prepare inputs for multiprocessing
        # All samples (bacterial genomes and MGEs) go to the same directories
        # The sample_is_mge set is used downstream to filter which samples go to OrthoFinder
        if sample_is_mge is None:
            sample_is_mge = set()
        
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
        for _ in _iter_progress(
            p.imap_unordered(create_chopped_proteomes, prot_mod_inputs),
            total=len(prot_mod_inputs),
            description="Creating domain-chopped proteomes",
        ):
            pass
        p.close()
        
        # Create domain coordinate info file and return sample CCDS proteomes
        sample_ccds_proteomes = {}
        with open(domain_coord_info_file, 'w') as dci_handle:
            # Write header
            dci_handle.write('\t'.join(['Sample', 'Protein', 'Domain_Type', 'Domain_Index', 'Domain_Count', 'Length']) + '\n')
            
            # Process all files in the domain_coords_dir
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
            for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="gathering", Z=1000000, cpus=threads):
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
                for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="gathering", Z=int(pfam_z), cpus=threads):
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


def integrate_mge_proteins_into_orthogroups(
    input_dir: str,
    sample_is_mge: Set[str],
    orthofinder_tsv_file: str,
    orthofinder_tsv_singletons_file: str,
    workspace_dir: str,
    orthofinder_mod_tsv_file: str,
    orthofinder_mod_tsv_singletons_file: str,
    mge_og_assignment_file: str,
    orthogroup_sequences_dir: str,
    threads: int = 1,
    ultra_sens: bool = False,
    evalue_cutoff: float = 1e-3,
    log_object: Any = None,
) -> None:
    """
    Integrate MGE (phage/plasmid) proteins into ortholog groups via DIAMOND alignment.
    
    Strategy:
    1. Concatenate all MGE proteins into a single FASTA file
    2. Concatenate all bacterial genome proteins and create a DIAMOND database
    3. Run DIAMOND blastp to align MGE proteins against bacterial genome proteins
    4. Merge MGE proteins into existing OGs based on best bitscore match (e-value <= cutoff)
    5. Create new OGs for MGE proteins with no matches meeting the cutoff
    
    Parameters:
    -----------
    input_dir : str
        Directory containing all ccds.faa files (both bacterial genomes and MGEs)
    sample_is_mge : Set[str]
        Set of sample names that are MGEs (plasmids/phages)
    orthofinder_tsv_file : str
        Path to Orthogroups.tsv from OrthoFinder
    orthofinder_tsv_singletons_file : str
        Path to Orthogroups_UnassignedGenes.tsv from OrthoFinder
    orthofinder_mod_tsv_file : str
        Path to Orthogroups_Modified.tsv from OrthoFinder
    orthofinder_mod_tsv_singletons_file : str
        Path to Orthogroups_UnassignedGenes_Modified.tsv from OrthoFinder
    workspace_dir : str
        Workspace directory for results
    mge_og_assignment_file : str
        Output file for MGE orthogroup assignments
    orthogroup_sequences_dir : str
        Path to OrthoFinder's Orthogroup_Sequences directory (will append MGE proteins to existing files)
    threads : int
        Number of threads for DIAMOND
    ultra_sens : bool
        Whether to use ultra-sensitive mode for DIAMOND
    evalue_cutoff : float
        E-value cutoff for assigning MGE proteins to OGs
    log_object : Any
        Logger object
    """
    try:
        log_object.info("Step 1b.1: Concatenating MGE proteins...")

        # Concatenate all MGE proteins into a single file
        mge_concat_fasta = os.path.join(workspace_dir, "MGE_protein_chunks_concatenated.faa")
        mge_protein_to_sample = {}  # Track which sample each MGE protein comes from
        
        # Concatenate bacterial genome proteins
        bacterial_genome_concat_fasta = os.path.join(workspace_dir, "Bacterial_genome_protein_chunks_concatenated.faa")
        bacterial_genome_protein_to_sample = {}  # Track which sample each bacterial genome protein comes from
        
        mge_count = 0
        bac_count = 0
        all_mge_protein_chunks = set([])
        mge_protein_chunks = defaultdict(set)
        all_mges = set([])
        mge_pchunk_to_seq = defaultdict(dict)
        with open(mge_concat_fasta, 'w') as mge_out, open(bacterial_genome_concat_fasta, 'w') as bac_out:
            for fasta_file in os.listdir(input_dir):
                if not fasta_file.endswith('.ccds.faa'): continue
              
                sample_name = fasta_file.replace('.ccds.faa', '')
                fasta_path = os.path.join(input_dir, fasta_file)
                
                # Determine if this sample is an MGE or bacterial genome
                is_mge = sample_name in sample_is_mge
                if is_mge:
                    all_mges.add(sample_name)
                with open(fasta_path) as in_handle:
                    for rec in SeqIO.parse(in_handle, 'fasta'):
                        if is_mge:
                            mge_protein_to_sample[rec.id] = sample_name
                            all_mge_protein_chunks.add(rec.id)
                            mge_protein_chunks[sample_name].add(rec.id)
                            mge_pchunk_to_seq[sample_name][rec.id] = str(rec.seq)
                            SeqIO.write(rec, mge_out, 'fasta')
                        else:
                            bacterial_genome_protein_to_sample[rec.id] = sample_name
                            SeqIO.write(rec, bac_out, 'fasta')
                
                if is_mge:
                    mge_count += 1
                else:
                    bac_count += 1
        mges_sorted = sorted(list(all_mges))


        log_object.info(f"Concatenated {len(mge_protein_to_sample)} MGE proteins from {mge_count} MGEs")
        log_object.info(f"Concatenated {len(bacterial_genome_protein_to_sample)} bacterial genome proteins from {bac_count} bacteria genomes")
        
        # Create DIAMOND bacterial genome database
        log_object.info("Step 1b.2: Creating DIAMOND databases from bacterial and MGE domain-resolution protein chunks...")
        diamond_bac_db = os.path.join(workspace_dir, "Bacterial_genome_domain-resolution_protein_chunks.dmnd")
        diamond_makedb_cmd = [
            'diamond', 'makedb', '--in', bacterial_genome_concat_fasta, '--db', diamond_bac_db, '--threads', str(threads)
        ]
        run_cmd(diamond_makedb_cmd, log_object, check_files=[diamond_bac_db])
        
        # Create DIAMOND MGE domain-resolution protein chunks database
        diamond_mge_db = os.path.join(workspace_dir, "MGE_domain-resolution_protein_chunks.dmnd")
        diamond_makedb_cmd = [
            'diamond', 'makedb', '--in', mge_concat_fasta, '--db', diamond_mge_db, '--threads', str(threads)
        ]
        run_cmd(diamond_makedb_cmd, log_object, check_files=[diamond_mge_db])

        sensitivity = 'ultra-sensitive' if ultra_sens else 'very-sensitive'
        log_object.info("Step 1b.3: Running reflexive alignment of MGE domain-resolution protein chunks against themselves to cluster them into homologous. groups...")

        mge_reflexive_diamond_output = os.path.join(workspace_dir, "MGE_vs_MGE_diamond.tsv")
        
        # Run reflexive alignment of MGE proteins against themselves
        diamond_blastp_cmd = [
            'diamond', 'blastp',
            '--db', diamond_mge_db,
            '--query', mge_concat_fasta,
            '--out', mge_reflexive_diamond_output,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovhsp',
                        'scovhsp',
            '--threads', str(threads),
            '--' + sensitivity,
            '--evalue', str(evalue_cutoff),
            '--max-target-seqs', '10'  # Only keep best 10 hits per query
        ]
        run_cmd(diamond_blastp_cmd, log_object, check_files=[mge_reflexive_diamond_output])

        pair_listing_file = os.path.join(workspace_dir, "Homologous_MGE_pairs.txt")
        with open(pair_listing_file, 'w') as plf:
            with open(mge_reflexive_diamond_output, 'r') as of:
                for line in of:
                    line = line.strip('\n')
                    ls = line.split('\t')
                    if len(ls) >= 12:
                        query = ls[0]
                        hit = ls[1]
                        qcovhsp = float(ls[12])
                        scovhsp = float(ls[13])
                        if query == hit: continue
                        if qcovhsp < 0.5 or scovhsp < 0.5: continue
                        plf.write(f"{query} {hit}\n")

        # Run slclust with explicit file redirection (shell operators don't work with subprocess.run without shell=True)
        slclust_output = os.path.join(workspace_dir, "Homologous_MGE_slclusters.txt")
        log_object.info(f"Running slclust with input: {pair_listing_file} and output: {slclust_output}")
        
        with open(pair_listing_file, 'r') as stdin_file, open(slclust_output, 'w') as stdout_file:
            result = subprocess.run(
                ['slclust'],
                stdin=stdin_file,
                stdout=stdout_file,
                stderr=subprocess.PIPE,
                check=True,
                text=True
            )
        
        log_object.info("slclust completed successfully")
        
        if not os.path.isfile(slclust_output):
            raise FileNotFoundError(f"slclust output file not created: {slclust_output}")

        cluster_proteins = defaultdict(set)
        with open(slclust_output, 'r') as of:
            for i, line in enumerate(of):
                line = line.strip('\n')
                ls = line.split()
                if len(ls) >= 2:
                    for pchunk in ls:
                        cluster_proteins[i].add(pchunk)

        # Run DIAMOND blastp
        log_object.info("Step 1b.4: Running DIAMOND blastp to align MGE proteins against bacterial genome proteins...")
        
        diamond_output = os.path.join(workspace_dir, "MGE_vs_BacterialGenomes_diamond.tsv")
        sensitivity = 'ultra-sensitive' if ultra_sens else 'very-sensitive'
        
        diamond_blastp_cmd = [
            'diamond', 'blastp',
            '--db', diamond_bac_db,
            '--query', mge_concat_fasta,
            '--out', diamond_output,
            '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
            '--threads', str(threads),
            '--' + sensitivity,
            '--evalue', str(evalue_cutoff),
            '--max-target-seqs', '10'  # Only keep best 10 hits per query
        ]
        
        run_cmd(diamond_blastp_cmd, log_object, check_files=[diamond_output])
        
        # Sort DIAMOND output by bitscore (column 12) in descending order using Unix sort
        log_object.info("Sorting DIAMOND output by bitscore (this may take a while for large files)...")
        diamond_output_sorted = os.path.join(workspace_dir, "MGE_vs_BacterialGenomes_diamond_sorted.tsv")
        
        # Use Unix sort for efficient sorting of large files
        # -t $'\t' : tab delimiter
        # -k12,12 : sort by column 12 (bitscore)
        # -n : numeric sort
        # -r : reverse order (highest first)
        # --parallel : use multiple threads for sorting
        sort_cmd = f"sort -t $'\\t' -k12,12 -n -r --parallel={threads} {diamond_output} > {diamond_output_sorted}"
        
        try:
            subprocess.run(sort_cmd, shell=True, check=True, executable='/bin/bash')
            log_object.info(f"Successfully sorted DIAMOND output by bitscore")
            # Use sorted file for downstream processing
            diamond_output = diamond_output_sorted
        except subprocess.CalledProcessError as e:
            log_object.error(f"Failed to sort DIAMOND output: {str(e)}")
            raise
        
        last_og_id = -1
        pchunk_to_og = dict()  
        bac_sample_count = 0
        with open(orthofinder_tsv_file, 'r') as ogf:
            for i, line in enumerate(ogf):
                line = line.strip('\n')
                ls = line.split('\t')
                if i == 0: 
                    bac_sample_count = len(ls[1:])
                    continue
                og_id = ls[0]
                protein_chunks = ls[1:]
                for pchunks in protein_chunks:
                    for pchunk in pchunks.split(','):
                        pchunk = pchunk.strip()
                        if pchunk == '': continue
                        pchunk_to_og[pchunk] = og_id
                        last_og_id = max([last_og_id, int(og_id[2:])])
                        bac_sample_count = len(ls[1:])

        with open(orthofinder_tsv_singletons_file, 'r') as ogf:
            for i, line in enumerate(ogf):
                if i == 0: continue 
                line = line.strip('\n')
                ls = line.split('\t')                
                og_id = ls[0]
                protein_chunks = ls[1:]
                for pchunks in protein_chunks:
                    for pchunk in pchunks.split(','):
                        pchunk = pchunk.strip()
                        if pchunk == '': continue
                        pchunk_to_og[pchunk] = og_id
                        last_og_id = max([last_og_id, int(og_id[2:])])
        last_og_id += 1

        mge_pchunk_to_og = dict()
        accounted_mge_pchunks = set([])
        ogs_hit = set([])
        with open(diamond_output_sorted, 'r') as of:
            for i, line in enumerate(of):
                line = line.strip('\n')
                ls = line.split('\t')
                if len(ls) != 12: continue
                query = ls[0]
                hit = ls[1]
                bitscore = float(ls[11])
                if query in accounted_mge_pchunks: continue
                mge_pchunk_to_og[query] = pchunk_to_og[hit]
                accounted_mge_pchunks.add(query)
                ogs_hit.add(pchunk_to_og[hit])

        difficult_to_fit = 0
        not_sogs = []
        for clust in cluster_proteins:
            total = len(cluster_proteins[clust])
            accounted = 0
            best_ogs_by_members = set([])
            for pchunk in cluster_proteins[clust]:
                if pchunk in accounted_mge_pchunks:
                    accounted += 1
                    best_ogs_by_members.add(mge_pchunk_to_og[pchunk])
            if accounted == total: continue
            if len(best_ogs_by_members) == 1:
                best_og = list(best_ogs_by_members)[0]
                for pchunk in cluster_proteins[clust]:
                    mge_pchunk_to_og[pchunk] = best_og
                    accounted_mge_pchunks.add(pchunk)
            elif len(best_ogs_by_members) == 0:
                new_og_id = generate_og_name(last_og_id)
                new_og_line = [new_og_id] + (['']*bac_sample_count)
                mge_og_pchunks = defaultdict(set)

                og_fasta_path = os.path.join(orthogroup_sequences_dir, new_og_id + '.fa')
                ofp_handle = open(og_fasta_path, 'w')
                for pchunk in cluster_proteins[clust]:
                    mge = pchunk.split('|')[0]
                    mge_og_pchunks[mge].add(pchunk)
                    accounted_mge_pchunks.add(pchunk)
                    ofp_handle.write('>' + pchunk + '\n' + mge_pchunk_to_seq[mge][pchunk] + '\n')
                ofp_handle.close()

                mge_parts = []
                for mge in mges_sorted:
                    mge_parts.append(', '.join(sorted(mge_og_pchunks[mge])))
                not_sogs.append('\t'.join(new_og_line + mge_parts))
                last_og_id += 1
            else:
                # this can be improved - but for now we prioritize 
                # minimizing false positive orthology prediction
                # and so MGE proteins that do not map to a protein
                # from a non-MGE context are not assigned to an OG
                # but are left as singletons
                difficult_to_fit += len(cluster_proteins[clust])

        warning_msg = f"There were {difficult_to_fit} MGE-derived proteins\n"
        warning_msg += f"which might be homologous/orthologous to each other, "
        warning_msg += f"but were differentially assigned to orthogroups determined\n"
        warning_msg += f"by OrthoFinder based on proteins from non-MGE contexts.\n"
        warning_msg += f"For proteins which didn't map to OrthoFinder orthogroups\n"
        warning_msg += f"directly, we currently leave them as singletons for such\n"
        warning_msg += f"cases.\n"
        log_object.warning(warning_msg)

        mge_og_pchunks = defaultdict(lambda: defaultdict(set))
        for pchunk in mge_pchunk_to_og:
            mge = pchunk.split('|')[0]
            og_id = mge_pchunk_to_og[pchunk]
            mge_og_pchunks[mge][og_id].add(pchunk)

        with open(orthofinder_mod_tsv_file, 'w') as ogf:
            with open(orthofinder_tsv_file, 'r') as ogf_orig:
                for i, line in enumerate(ogf_orig):
                    line = line.strip('\n')
                    if i == 0:
                        ogf.write(line + '\t' + '\t'.join([x + '.ccds' for x in mges_sorted]) + '\n')
                    else:
                        ls = line.split('\t')
                        og_id = ls[0]
                        og_fasta_path = os.path.join(orthogroup_sequences_dir, og_id + '.fa')
                        ofp_handle = open(og_fasta_path, 'a+')
                        mge_pcs = []
                        for mge in mges_sorted:
                            if mge in mge_og_pchunks:
                                if og_id in mge_og_pchunks[mge]:
                                    mge_pcs.append(', '.join(list(mge_og_pchunks[mge][og_id])))
                                    for pchunk in mge_og_pchunks[mge][og_id]:
                                        ofp_handle.write('>' + pchunk + '\n' + mge_pchunk_to_seq[mge][pchunk] + '\n')
                                else:
                                    mge_pcs.append('')
                            else:
                                mge_pcs.append('')
                        ofp_handle.close()
                        ogf.write(line + '\t' + '\t'.join(mge_pcs) + '\n')

        with open(orthofinder_mod_tsv_singletons_file, 'w') as ogf:
            with open(orthofinder_tsv_singletons_file, 'r') as ogf_orig:
                for i, line in enumerate(ogf_orig):
                    line = line.strip('\n')
                    if i == 0:
                        ogf.write(line + '\t' + '\t'.join([x + '.ccds' for x in mges_sorted]) + '\n')
                    else:
                        ls = line.split('\t')
                        og_id = ls[0]
                        lts = [x for x in ls[1:] if x.strip() != '' and ',' not in x]
                        if len(lts) != 1:
                            raise ValueError(f"Expected 1 singleton protein for {og_id}, got {len(lts)}")                            
                        if og_id in ogs_hit:
                            mge_matches = []
                            og_fasta_path = os.path.join(orthogroup_sequences_dir, og_id + '.fa')
                            ofp_handle = open(og_fasta_path, 'a+')
                            for mge in mges_sorted:
                                if mge in mge_og_pchunks:
                                    if og_id in mge_og_pchunks[mge]:
                                        mge_matches.append(', '.join(list(mge_og_pchunks[mge][og_id])))
                                        for pchunk in mge_og_pchunks[mge][og_id]:
                                            ofp_handle.write('>' + pchunk + '\n' + mge_pchunk_to_seq[mge][pchunk] + '\n')
                                    else:
                                        mge_matches.append('')
                                else:
                                    mge_matches.append('')
                            ofp_handle.close()

                            if len([x for x in mge_matches if x != '']) == 0:
                                ogf.write(line + '\t' + '\t'.join(['']*len(mges_sorted)) + '\n')
                            else:
                                updated_line = line + '\t' + '\t'.join(mge_matches)
                                not_sogs.append(updated_line)
                        else:
                            ogf.write(line + '\t' + '\t'.join(['']*len(mges_sorted)) + '\n')

            for i, mge in enumerate(mges_sorted):
                for pchunk in mge_protein_chunks[mge]:
                    if pchunk in accounted_mge_pchunks: continue
                    new_sog_id = generate_og_name(last_og_id)
                    line_part_1 = [new_sog_id] + (['']*bac_sample_count) + (['']*i)
                    remainder_i = len(mges_sorted) - i - 1
                    line_part_2 = ['']*remainder_i
                    ogf.write('\t'.join(line_part_1 + [pchunk] + line_part_2) + '\n')
                    if (len(line_part_1[1:]) + len(line_part_2) + 1) != (len(mges_sorted) + bac_sample_count):
                        raise ValueError(f"Expected {len(mges_sorted) + bac_sample_count} columns, got {len(line_part_1[1:]) + len(line_part_2) + 1}")

                    og_fasta_path = os.path.join(orthogroup_sequences_dir, new_sog_id + '.fa')
                    ofp_handle = open(og_fasta_path, 'a+')
                    ofp_handle.write('>' + pchunk + '\n' + str(mge_pchunk_to_seq[mge][pchunk]) + '\n')
                    ofp_handle.close()
                    last_og_id += 1
        
        with open(orthofinder_mod_tsv_file, 'a+') as ogf:
            for line in not_sogs:
                ogf.write(line + '\n')

    except Exception as e:
        log_object.error("Error integrating MGE protein chunks into domain-resolution ortholog groups")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        raise
