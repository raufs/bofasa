"""
Analysis functions for bofasa.

This module contains functions for analyzing ortholog groups, creating alignments,
generating reports, and performing various analyses on the results.
"""

from math import log
import multiprocessing
import os
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
                        for j, gs in enumerate(ls[1:]):
                            sample = samples[j]
                            for g in gs.split(','):
                                g = g.strip()
                                if g != '':
                                    og_genes[og].add(g)
                                    og_samples[og].add(sample)
                                    gene_to_og[g] = og

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
                output_file = os.path.join(surround_info_dir, f"{sample}.tsv")
                genome_params.append([sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object])

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
                    context_nogs_complete[og].append(len(context_ogs))

        # Write detailed context information
        with open(og_context_info_file, 'w') as og_context_info_handle:
            og_context_info_handle.write('\t'.join([
                'OG', 'Median OG length (bp)', 'Percentage contexts near scaffold edge', 'Number of genomes with OG', 
                'Number of protein in OG', 'Context conservation score', 'Context conservation score - complete contexts', 
                'Context entropy score', 'Context entropy score - complete contexts', 'Number of distinct neighbor OGs', 
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
                
                context_entropy = 'NA'
                context_var_score = 'NA'
                context_entropy_complete = 'NA'
                context_var_score_complete = 'NA'
                
                if total_nog > 0:
                    context_var_score = round(avg_nog/total_nog, 2)
                    if total_nog > 1:
                        context_entropy = round(stats.entropy([x/sum_nog_freqs for x in nog_freqs]), total_nog)
                
                if total_nog_complete > 0:
                    context_var_score_complete = round(avg_nog_complete/total_nog_complete, 2)
                    if total_nog_complete > 1:
                        context_entropy_complete = round(stats.entropy([x/sum_nog_freqs_complete for x in nog_freqs_complete]), total_nog_complete)
                
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
                        
                og_context_info_handle.write('\t'.join([str(x) for x in [
                    og, round(median_length, 2), nse_perc, num_samples, num_contexts, context_var_score, 
                    context_var_score_complete, context_entropy, context_entropy_complete, total_nog, total_nog_complete,
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
            "Context entropy score",
            "Context entropy score - complete contexts",
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
        max_context_ent = max_or_zero("Context entropy score")
        max_context_ent_comp = max_or_zero("Context entropy score - complete contexts")

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
        worksheet.conditional_format(
            f"F2:F{num_rows}",
            {"type": "2_color_scale", "min_color": "#e6f5ab", "max_color": "#a4b36b", "min_value": 0.0, "max_value": max_context_var, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"G2:G{num_rows}",
            {"type": "2_color_scale", "min_color": "#e6f5ab", "max_color": "#a4b36b", "min_value": 0.0, "max_value": max_context_var_comp, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"H2:H{num_rows}",
            {"type": "2_color_scale", "min_color": "#b3e3d6", "max_color": "#6aa192", "min_value": 0.0, "max_value": max_context_ent, "min_type": "num", "max_type": "num"},
        )
        worksheet.conditional_format(
            f"I2:I{num_rows}",
            {"type": "2_color_scale", "min_color": "#b3e3d6", "max_color": "#6aa192", "min_value": 0.0, "max_value": max_context_ent_comp, "min_type": "num", "max_type": "num"},
        )

        if genomad_flag:
            worksheet.conditional_format(
                f"N2:N{num_rows}",
                {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"},
            )
            worksheet.conditional_format(
                f"O2:O{num_rows}",
                {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"},
            )

        worksheet.conditional_format(
            f"P2:P{num_rows}",
            {"type": "2_color_scale", "min_color": "#ffffff", "max_color": "#ed9393", "min_value": 0.0, "max_value": 1.0, "min_type": "num", "max_type": "num"},
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
            'Context entropy score', 
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
                outf_handle.write('\t'.join([ls[0], ls[4], ls[7], mge_related]) + '\n')
        outf_handle.close()
                
        numeric_columns = set(['Number of protein in OG', 'Context entropy score'])
        simple_df = load_table_in_pandas_dataframe(tmp_result_file, numeric_columns)
        fig = px.scatter(
            simple_df, 
            x="Number of protein in OG", 
            y="Context entropy score",
            color="Majority of protein instances homologous to IS-element or on plasmid or phage", 
            color_discrete_map={"No": "grey", "Yes": "red"},
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
                genomad_cmd = [
                    'genomad', 'end-to-end', '--cleanup', '--threads', str(threads), 
                    '--splits', str(genome_splits), input_genome, genomad_results, genomad_db_dir
                ]
                
                try:
                    run_cmd(genomad_cmd, log_object)
                except subprocess.CalledProcessError as e:
                    log_object.error(f"geNomad failed with unknown error for sample {sample}")
                    log_object.error(f"Command: {' '.join(genomad_cmd)}")
                    log_object.error(f"Full error output: {e.output}")
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
        setup_ready_directory([annot_dir], overwrite_mode="overwrite")
        
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
        setup_ready_directory([split_proteins_dir, domain_coords_dir], overwrite_mode="overwrite")
        
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
