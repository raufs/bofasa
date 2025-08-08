"""
Alignment functions for bofasa.

This module contains functions for creating multiple sequence alignments
using MUSCLE super5 alignment tool.
"""

import os
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from typing import Any, Dict, List, Set
from Bio import SeqIO
from .utils import assess_job_intensity, setup_ready_directory, _iter_progress, multi_process
from . import config


def run_alignment(inputs) -> None:
    """
    Create protein alignments using MUSCLE super5.

    Parameters:
    -----------
    inputs : list[str]
        List of inputs for the MUSCLE super5 alignment job
        prot_file: str,
        prot_algn_file: str,
        threads: int,
    """
    prot_file, prot_algn_file, threads = inputs
    try:
        # Check if input file exists and has sequences
        try:
            with open(prot_file, 'r') as handle:
                sequences = list(SeqIO.parse(handle, 'fasta'))
                if len(sequences) < 2:
                    return True
        except Exception as e:
            sys.stderr.write(f"Error reading input file {prot_file}: {str(e)}\n")
            return False

        cmd = f"muscle -super5 {prot_file} -threads {threads} -output {prot_algn_file}"
        
        # Run MUSCLE super5 command with shell redirection
        try:
            subprocess.call(cmd, shell=True, executable="/bin/bash", 
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"MUSCLE super5 command failed for {prot_file}\n")
            return False

        # Verify alignment file was created and is not empty
        if not os.path.exists(prot_algn_file):
            sys.stderr.write(f"Alignment file was not created: {prot_algn_file}\n")
            return False

        if os.path.getsize(prot_algn_file) == 0:
            sys.stderr.write(f"Alignment file is empty: {prot_algn_file}\n")
            return False

        # Verify alignment content
        try:
            with open(prot_algn_file, 'r') as handle:
                alignment_records = list(SeqIO.parse(handle, 'fasta'))
                if len(alignment_records) == 0:
                    msg = f"Alignment file contains no sequences: {prot_algn_file}\n"
                    raise RuntimeError(msg)
                if len(alignment_records) != len(sequences):
                    msg = (
                        f"Alignment file sequence count mismatch for {prot_algn_file}: "
                        f"expected {len(sequences)}, got {len(alignment_records)}\n"
                    )
                    raise RuntimeError(msg)
        except Exception as e:
            msg = f"Error reading alignment file {prot_algn_file}: {str(e)}\n"
            raise RuntimeError(msg)
        return True

    except Exception as e:
        sys.stderr.write(f"Error aligning {prot_file}: {str(e)}\n")
        return False


def create_domain_alignments(
    dog_seqs_dir: str,
    dog_algn_dir: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,  
    more_deterministic: bool = False,
) -> None:
    """
    Create domain alignments using parallel processing.

    This function aligns domain-specific protein sequences for coarse domain
    ortholog groups. Each alignment job uses 1 thread, and the work
    is parallelized across N processes where N is the number of threads.

    Parameters:
    -----------
    dog_seqs_dir : str
        Directory containing DOG protein sequence files
    dog_algn_dir : str
        Directory where aligned DOG protein files will be saved
    log_object : logging.Logger
        Logger object for recording progress and errors
    threads : int, default=4
        Number of processes to use for parallel processing (each uses 1 thread)
    more_deterministic : bool, default=False
        Whether to use more deterministic settings for reproducible results
    """
    try:
        # Get list of DOG sequence files
        dog_files: List[str] = [f for f in os.listdir(dog_seqs_dir) if f.endswith('.fa')]

        if not dog_files:
            log_object.warning(f"No DOG sequence files found in {dog_seqs_dir}")
            return

        # Prepare tasks for parallel processing
        heavy_jobs = []
        light_jobs = []
        all_jobs = []
        skipped_files = 0
        for dog_file in dog_files:
            seq_file: str = os.path.join(dog_seqs_dir, dog_file)
            prefix: str = '.fa'.join(dog_file.split('.fa')[:-1])
            algn_file: str = os.path.join(dog_algn_dir, f"{prefix}.msa.faa")

            try:
                if more_deterministic:
                    all_jobs.append([seq_file, algn_file, 1])
                else:
                    heavy_job = assess_job_intensity(seq_file)
                    if heavy_job:
                        heavy_jobs.append([seq_file, algn_file, threads])
                    else:
                        light_jobs.append([seq_file, algn_file, 1])
            except Exception as e:
                skipped_files += 1
                log_object.error(
                    f"Error assessing job intensity for {seq_file}: {str(e)}"
                )
                continue

        _perform_parallelization(
            all_jobs,
            heavy_jobs,
            light_jobs,
            log_object,
            skipped_files,
            more_deterministic,
            threads,
            "domain",
        )

    except Exception as e:
        log_object.error(
            "Issues with creating domain protein alignments using MUSCLE super5"
        )
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
    more_deterministic: bool = False,
) -> None:
    """
    Create protein alignments using MUSCLE super5 for full 
    protein ortholog groups.

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
    more_deterministic : bool, default=False
        Whether to use more deterministic settings for reproducible results
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
        proteins_extracted = 0
        try:
            for f in os.listdir(proteome_dir):
                if not f.endswith(".faa"):
                    continue
                s: str = ".faa".join(f.split(".faa")[:-1])
                try:
                    with open(os.path.join(proteome_dir, f)) as opf:
                        for rec in SeqIO.parse(opf, "fasta"):
                            lt: str = rec.id
                            if lt in prot_to_og[s]:
                                og: str = prot_to_og[s][lt]
                                outf: str = os.path.join(prot_dir, og + ".faa")
                                with open(outf, "a+") as outfh:
                                    outfh.write(">" + s + "|" + lt + "\n" + str(rec.seq) + "\n")
                                proteins_extracted += 1
                except Exception as e:
                    log_object.error(
                        f"Error processing proteome file {f}: {str(e)}"
                    )
                    continue
        except Exception as e:
            log_object.error(
                f"Error accessing proteome directory {proteome_dir}: {str(e)}"
            )
            raise
        
        log_object.info(f"Extracted {proteins_extracted} proteins for alignment")

        # Create output directory
        setup_ready_directory([prot_algn_dir], overwrite_mode="overwrite")

        # Process each protein file with job intensity assessment
        from .utils import assess_job_intensity

        heavy_jobs = []
        light_jobs = []
        all_jobs = []
        skipped_files = 0        
        try:
            for pf in os.listdir(prot_dir):
                if not pf.endswith('.faa'):
                    continue
                    
                prefix: str = '.faa'.join(pf.split('.faa')[:-1])
                prot_file: str = os.path.join(prot_dir, pf)
                prot_algn_file: str = os.path.join(prot_algn_dir, prefix + '.msa.faa')

                if more_deterministic:
                    all_jobs.append([prot_file, prot_algn_file, 1])
                else:
                    try:
                        # Assess job intensity
                        heavy_job: bool = assess_job_intensity(prot_file)

                        if heavy_job:
                            heavy_jobs.append([prot_file, prot_algn_file, threads])
                        else:
                            light_jobs.append([prot_file, prot_algn_file, 1])
                    except Exception as e:
                        log_object.error(
                            f"Error assessing job intensity for {pf}: {str(e)}"
                        )
                        skipped_files += 1
                        continue
                
        except Exception as e:
            log_object.error(
                f"Error processing protein directory {prot_dir}: {str(e)}"
            )
            raise
        
        _perform_parallelization(
            all_jobs,
            heavy_jobs,
            light_jobs,
            log_object,
            skipped_files,
            more_deterministic,
            threads,
            "protein",
        )
    except Exception as e:
        log_object.error("Error creating protein alignments")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)
        

def _perform_parallelization(
    all_jobs: List[List[str]],
    heavy_jobs: List[List[str]],
    light_jobs: List[List[str]],
    log_object: Any,
    skipped_files: int,
    more_deterministic: bool,
    threads: int,
    type: str = "protein",
) -> None:
    """
    Perform parallelization of protein alignments.
    """

    if more_deterministic:
        log_object.info(f"Prepared {len(all_jobs)} jobs for alignment")
        if skipped_files > 0:
            log_object.info(f"Skipped {skipped_files} files (already aligned or errors)")
            
        # process all jobs
        try:
            with ThreadPoolExecutor(max_workers=1) as executor:
                results: List[bool] = list(
                    _iter_progress(
                        executor.map(run_alignment, all_jobs),
                        total=len(all_jobs),
                        description=f"Creating {type} alignments",
                    )
                )
        except Exception as e:
            log_object.error(
                f"Error during threaded protein alignment: {str(e)}"
            )
            log_object.error(traceback.format_exc())
            raise

        # Count successful alignments
        successful_alignments = sum(results)
        failed_alignments = len(all_jobs) - successful_alignments

        log_object.info(
            f"{type} alignments completed: {successful_alignments} successful, "
            f"{failed_alignments} failed"
        )

    else:
        log_object.info(
            f"Prepared {len(heavy_jobs)} heavy jobs and {len(light_jobs)} light jobs for alignment"
        )
        if skipped_files > 0:
            log_object.info(
                f"Skipped {skipped_files} files (already aligned or errors)"
            )
        
        if not heavy_jobs:
            log_object.info(f"No heavy {type} alignment jobs")
            
        # process heavy jobs
        if heavy_jobs:
            try:
                with ThreadPoolExecutor(max_workers=1) as executor:
                    results: List[bool] = list(
                        _iter_progress(
                            executor.map(run_alignment, heavy_jobs),
                            total=len(heavy_jobs),
                            description=f"Creating heavy {type} alignments",
                        )
                    )
            except Exception as e:
                log_object.error(
                    f"Error during threaded protein alignment: {str(e)}"
                )
                log_object.error(traceback.format_exc())
                raise

            # Count successful alignments
            successful_alignments = sum(results)
            failed_alignments = len(heavy_jobs) - successful_alignments

            log_object.info(
                f"Heavy {type} alignments completed: {successful_alignments} successful, {failed_alignments} failed"
            )

        # Process light jobs
        try:
            with ThreadPoolExecutor(max_workers=threads) as executor:
                results: List[bool] = list(
                    _iter_progress(
                        executor.map(run_alignment, light_jobs),
                        total=len(light_jobs),
                        description=f"Creating light {type} alignments",
                    )
                )
        except Exception as e:
            log_object.error(
                f"Error during threaded light protein alignment: {str(e)}"
            )
            log_object.error(traceback.format_exc())
            raise

        # Count successful alignments
        successful_alignments = sum(results)
        failed_alignments = len(light_jobs) - successful_alignments

        log_object.info(
            f"Light {type} alignments completed: {successful_alignments} successful, "
            f"{failed_alignments} failed"
        )


def create_near_scc_resolved_domain_protein_alignments(
    bofasa_prep_dir: str,
    resulting_dogs_file: str,
    dogs_seqs_dir: str,
    dogs_algn_dir: str,
    dogs_trim_dir: str,
    merged_core_genome_file: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,
    near_scc_prop: float = config.DEFAULT_NEAR_SCC_PROP,
    trimal_options: str = config.DEFAULT_TRIMAL_OPTIONS,
    allow_mge: bool = False,
    more_deterministic: bool = False,
) -> None:
    """
    Create protein alignments for near single-copy-core resolved domain ortholog groups using MUSCLE super5.

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
    near_scc_prop : float, default=config.DEFAULT_NEAR_SCC_PROP
        Proportion for near single-copy-core
    trimal_options : str, default=config.DEFAULT_TRIMAL_OPTIONS
        TrimAl options
    allow_mge : bool, default=False
        Whether to allow mobile genetic elements
    more_deterministic : bool, default=False
        Whether to use more deterministic settings for reproducible results
    """
    try:
        # Create output directories
        setup_ready_directory([dogs_seqs_dir, dogs_algn_dir, dogs_trim_dir], overwrite_mode="overwrite")

        # Load MGE protein sets
        isfinder_file: str = os.path.join(bofasa_prep_dir, 'Sample_IS_Element_Proteins.txt')
        plasmid_file: str = os.path.join(bofasa_prep_dir, 'Sample_Plasmid_Proteins.txt')
        phage_file: str = os.path.join(bofasa_prep_dir, 'Sample_Phage_Proteins.txt')

        mge_set: Set[str] = set()

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

        # Parse domain ortholog groups and filter by MGE content and single-copy status
        prot_to_og: Dict[str, str] = {}
        samples: List[str] = []

        with open(resulting_dogs_file) as orof:
            for i, line in enumerate(orof):
                line = line.strip()
                ls = line.split('\t')
                if i == 0:
                    samples = ls[1:]
                    continue

                og: str = ls[0]
                og_has_mge_prot: bool = False

                # Check if OG contains MGE proteins
                for lts in ls[1:]:
                    for lt in lts.split(','):
                        lt = lt.strip()
                        if lt == '':
                            continue
                        if lt.split('|')[1] in mge_set:
                            og_has_mge_prot = True

                if not allow_mge and og_has_mge_prot:
                    continue

                # Check single-copy status across samples
                samples_with_sc: Set[str] = set()
                for j, lts in enumerate(ls[1:]):
                    s: str = samples[j]
                    if ',' in lts:
                        continue
                    for lt in lts.split(', '):
                        lt = lt.strip()
                        if lt != '':
                            samples_with_sc.add(lt)

                # Only include OGs that meet the near single-copy-core threshold
                if float(len(samples_with_sc) / len(samples)) >= near_scc_prop:
                    for lt in samples_with_sc:
                        prot_to_og[lt] = og

        # Extract protein sequences for selected OGs
        proteome_dir: str = os.path.join(bofasa_prep_dir, 'Domain_and_Interdomain_FASTAs')
        for f in os.listdir(proteome_dir):
            if not f.endswith('.faa'):
                continue
            s: str = '.faa'.join(f.split('.faa')[:-1])
            with open(os.path.join(proteome_dir, f)) as opf:
                for rec in SeqIO.parse(opf, 'fasta'):
                    lt: str = rec.id
                    if lt in prot_to_og:
                        og: str = prot_to_og[lt]
                        outf: str = os.path.join(dogs_seqs_dir, og + '.faa')
                        with open(outf, 'a+') as outfh:
                            outfh.write(f'>{lt}\n{str(rec.seq)}\n')

        # Create alignment and trimming commands
        msa_trim_cmds: List[List[str]] = []
        heavy_jobs = []
        light_jobs = []
        all_jobs = []
        for pf in os.listdir(dogs_seqs_dir):
            prefix: str = '.faa'.join(pf.split('.faa')[:-1])
            prot_file: str = os.path.join(dogs_seqs_dir, pf)
            prot_algn_file: str = os.path.join(dogs_algn_dir, prefix + '.msa.faa')
            prot_algn_trim_file: str = os.path.join(dogs_trim_dir, prefix + '.trimmed.msa.faa')

            if more_deterministic:
                all_jobs.append([prot_file, prot_algn_file, 1])
            else:   
                # Assess job intensity to determine threading strategy
                heavy_job: bool = assess_job_intensity(prot_file)

                if heavy_job:
                    # For heavy jobs, run MUSCLE separately
                    heavy_jobs.append([prot_file, prot_algn_file, threads])
                    trimal_cmd: List[str] = ['trimal', '-in', prot_algn_file, '-out', prot_algn_trim_file] + trimal_options.split() + [log_object]
                    msa_trim_cmds.append(trimal_cmd)
                else:
                    # For light jobs, run MUSCLE in parallel
                    light_jobs.append([prot_file, prot_algn_file, 1])
                    trimal_cmd: List[str] = ['trimal', '-in', prot_algn_file, '-out', prot_algn_trim_file] + trimal_options.split() + [log_object]
                    msa_trim_cmds.append(trimal_cmd)
        
        _perform_parallelization(
            all_jobs,
            heavy_jobs,
            light_jobs,
            log_object,
            0,
            more_deterministic,
            threads,
            "domain",
        )

        # Run alignment and trimming commands in parallel
        if msa_trim_cmds:
            msg: str = f"Running {len(msa_trim_cmds)} alignment and trimming jobs"
            log_object.info(msg)

            import multiprocessing
            p = multiprocessing.Pool(threads)
            try:
                for _ in _iter_progress(
                    p.imap_unordered(multi_process, msa_trim_cmds),
                    total=len(msa_trim_cmds),
                    description="Running alignments and trimming",
                ):
                    pass
            except Exception as e:
                log_object.error("Error in alignment and trimming multiprocessing")
                log_object.error(str(e))
                log_object.error(traceback.format_exc())
            finally:
                p.close()

        # Create merged core genome alignment
        sample_seqs: Dict[str, str] = defaultdict(lambda: "")

        for f in os.listdir(dogs_trim_dir):
            samples_accounted: Set[str] = set()
            seqlen: int = 0

            with open(os.path.join(dogs_trim_dir, f)) as ortf:
                for rec in SeqIO.parse(ortf, 'fasta'):
                    seq: str = str(rec.seq)
                    seqlen = len(seq)
                    if seqlen == 0:
                        continue
                    s: str = rec.id.split('|')[0]
                    samples_accounted.add(s)
                    sample_seqs[s] += seq

            # Add gaps for samples not present in this alignment
            for s in samples:
                if s not in samples_accounted:
                    sample_seqs[s] += ('-' * seqlen)

        # Write merged core genome alignment
        with open(merged_core_genome_file, 'w') as aln_handle:
            for s in sample_seqs:
                aln_handle.write(f'>{s}\n{sample_seqs[s]}\n')

        log_object.info("Near single-copy-core resolved domain protein alignments completed successfully")

    except Exception as e:
        log_object.error("Issues with creating near single-copy-core resolved domain protein alignments")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        raise
