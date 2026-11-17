"""
Pfam domain landscape (domscape) functionality for bofasa.

This module annotates proteomes with Pfam HMMs (using gathering thresholds) and
produces a Pfam-by-genome matrix of the number of proteins carrying each domain.
It is designed to scale to thousands of genomes: Pfam HMMs are loaded once and
proteomes are searched in batches, with per-genome results written to disk so
that interrupted runs can be resumed.
"""

import os
import sys
import traceback
from collections import defaultdict
from typing import Any, Dict, List, Tuple

import numpy as np
import pyhmmer

from .utils import _iter_progress, pyhmmer_text as _to_str


def get_pfam_db_info() -> Tuple[str, int]:
    """
    Get the Pfam HMM file path and Z value from the bofasa database manifest.

    Returns:
    --------
    Tuple[str, int]
        Path to the Pfam HMM file and the Z value to use for hmmsearch
    """
    bofasa_db_dir = os.getenv("BOFASA_DB_PATH", "").strip()
    if not bofasa_db_dir:
        raise ValueError("BOFASA_DB_PATH environment variable not set")

    db_locations = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
    if not os.path.isfile(db_locations):
        raise FileNotFoundError(f"Database locations file not found: {db_locations}")

    pfam_hmm_path = None
    pfam_z = None
    with open(db_locations) as odb:
        for line in odb:
            ls = line.strip().split('\t')
            if ls[0] == 'pfam':
                pfam_hmm_path = ls[2]
                pfam_z = int(ls[3])
                break

    if not pfam_hmm_path or not os.path.isfile(pfam_hmm_path) or pfam_z is None:
        raise FileNotFoundError(f"Pfam HMM file or Z value not found: {pfam_hmm_path}")

    return pfam_hmm_path, pfam_z


def annotate_pfam_domains(
    sample_proteomes: Dict[str, str],
    hits_dir: str,
    log_object: Any,
    threads: int = 1,
    batch_size: int = 100,
) -> Dict[str, str]:
    """
    Annotate proteomes with Pfam HMMs using gathering thresholds.

    For each sample a tab-delimited file is written to hits_dir listing every
    protein-domain pair where at least one domain hit met the Pfam gathering
    threshold. Samples whose hits file already exists are skipped, allowing
    interrupted runs to resume.

    Parameters:
    -----------
    sample_proteomes : Dict[str, str]
        Dictionary mapping sample names to proteome FASTA file paths
    hits_dir : str
        Directory for per-sample Pfam hit files
    log_object : Any
        Logger object
    threads : int, default=1
        Number of threads to use for hmmsearch
    batch_size : int, default=100
        Number of proteomes searched together in a single hmmsearch call

    Returns:
    --------
    Dict[str, str]
        Dictionary mapping sample names to their Pfam hit files (only samples
        that were annotated successfully)
    """
    pfam_hmm_path, pfam_z = get_pfam_db_info()

    sample_hits_files = {}
    pending_samples = []
    for sample in sorted(sample_proteomes):
        hits_file = os.path.join(hits_dir, f"{sample}.pfam_hits.tsv")
        sample_hits_files[sample] = hits_file
        if not os.path.isfile(hits_file):
            pending_samples.append(sample)

    num_done = len(sample_proteomes) - len(pending_samples)
    if num_done > 0:
        log_object.info(f"Found existing Pfam annotations for {num_done} samples - skipping them")

    if not pending_samples:
        return sample_hits_files

    log_object.info(f"Loading Pfam HMMs from {pfam_hmm_path}")
    with pyhmmer.plan7.HMMFile(pfam_hmm_path) as hmm_file:
        pfam_hmms = list(hmm_file)
    log_object.info(f"Loaded {len(pfam_hmms)} Pfam HMMs")

    alphabet = pyhmmer.easel.Alphabet.amino()
    batches = [pending_samples[i:i + batch_size] for i in range(0, len(pending_samples), batch_size)]
    log_object.info(
        f"Annotating {len(pending_samples)} proteomes with Pfam in {len(batches)} batch(es) "
        f"of up to {batch_size} proteomes using {threads} threads"
    )

    failed_samples = []
    for batch in _iter_progress(batches, total=len(batches), description="Annotating Pfam domains"):
        # Load proteins, renaming each to "<batch sample index>_<protein index>" so hits
        # can be traced back to the right sample and protein regardless of input IDs
        sequences = []
        batch_protein_ids: List[List[str]] = []
        batch_samples: List[str] = []
        for sample in batch:
            sample_idx = len(batch_samples)
            protein_ids = []
            try:
                with pyhmmer.easel.SequenceFile(sample_proteomes[sample], digital=True, alphabet=alphabet) as seq_file:
                    sample_seqs = list(seq_file)
            except Exception as e:
                log_object.error(f"Unable to read proteome for sample {sample}: {str(e)}")
                failed_samples.append(sample)
                continue
            for seq in sample_seqs:
                protein_ids.append(_to_str(seq.name))
                seq.name = f"{sample_idx}_{len(protein_ids) - 1}".encode()
                sequences.append(seq)
            batch_samples.append(sample)
            batch_protein_ids.append(protein_ids)

        # (sample index, protein index) -> {pfam accession: [domains meeting GA, best domain score]}
        protein_pfam_hits: Dict[Tuple[int, int], Dict[str, List[Any]]] = defaultdict(dict)
        pfam_names: Dict[str, str] = {}
        if sequences:
            for hits in pyhmmer.hmmsearch(pfam_hmms, sequences, bit_cutoffs="gathering", Z=int(pfam_z), cpus=threads):
                pfam_name = _to_str(hits.query.name)
                pfam_acc = _to_str(hits.query.accession) or pfam_name
                pfam_names[pfam_acc] = pfam_name
                for hit in hits:
                    included_domains = hit.domains.included
                    domain_scores = [domain.score for domain in included_domains]
                    if not domain_scores:
                        continue
                    sample_idx, protein_idx = (int(x) for x in _to_str(hit.name).split('_'))
                    protein_pfam_hits[(sample_idx, protein_idx)][pfam_acc] = [len(domain_scores), max(domain_scores)]

        # Write per-sample hit files (via a temporary file so partial files are never treated as done)
        sample_rows: Dict[int, List[str]] = defaultdict(list)
        for (sample_idx, protein_idx) in sorted(protein_pfam_hits):
            protein_id = batch_protein_ids[sample_idx][protein_idx]
            for pfam_acc, (num_domains, best_score) in sorted(protein_pfam_hits[(sample_idx, protein_idx)].items()):
                sample_rows[sample_idx].append(
                    f"{protein_id}\t{pfam_acc}\t{pfam_names[pfam_acc]}\t{num_domains}\t{best_score:.1f}\n"
                )

        for sample_idx, sample in enumerate(batch_samples):
            hits_file = sample_hits_files[sample]
            tmp_file = hits_file + '.tmp'
            with open(tmp_file, 'w') as out_handle:
                out_handle.write('Protein\tPfam_Accession\tPfam_Name\tDomains_Meeting_GA\tBest_Domain_Score\n')
                out_handle.writelines(sample_rows[sample_idx])
            os.replace(tmp_file, hits_file)

    for sample in failed_samples:
        del sample_hits_files[sample]

    if failed_samples:
        log_object.warning(f"Pfam annotation failed for {len(failed_samples)} samples: {', '.join(failed_samples)}")

    return sample_hits_files


def create_pfam_count_matrix(
    sample_hits_files: Dict[str, str],
    matrix_file: str,
    pfam_info_file: str,
    protein_listing_file: str,
    log_object: Any,
) -> None:
    """
    Create a Pfam-by-sample matrix of protein counts.

    Each cell is the number of proteins in a sample with at least one domain hit
    meeting the gathering threshold for that Pfam domain - a protein carrying
    multiple copies of a domain is only counted once, and overlapping hits to
    different domains are each counted.

    Parameters:
    -----------
    sample_hits_files : Dict[str, str]
        Dictionary mapping sample names to their Pfam hit files
    matrix_file : str
        Output path for the Pfam (rows) by sample (columns) count matrix
    pfam_info_file : str
        Output path for Pfam domain information (names and prevalence)
    protein_listing_file : str
        Output path for the listing of proteins carrying each Pfam domain
    log_object : Any
        Logger object
    """
    try:
        samples = sorted(sample_hits_files)

        # First pass - gather all Pfam domains observed
        pfam_names: Dict[str, str] = {}
        for sample in samples:
            with open(sample_hits_files[sample]) as ohf:
                next(ohf)
                for line in ohf:
                    ls = line.rstrip('\n').split('\t')
                    pfam_names[ls[1]] = ls[2]

        pfams = sorted(pfam_names)
        pfam_index = {pfam: i for i, pfam in enumerate(pfams)}
        log_object.info(f"Building count matrix for {len(pfams)} Pfam domains across {len(samples)} samples")

        # Second pass - fill counts and list proteins; each line of a hits file is a
        # unique protein-domain pair
        counts = np.zeros((len(pfams), len(samples)), dtype=np.uint32)
        with open(protein_listing_file, 'w') as pl_handle:
            pl_handle.write('Pfam\tGenome\tProtein\n')
            for sample_idx, sample in enumerate(samples):
                with open(sample_hits_files[sample]) as ohf:
                    next(ohf)
                    for line in ohf:
                        protein, pfam = line.split('\t', 2)[:2]
                        counts[pfam_index[pfam], sample_idx] += 1
                        pl_handle.write(f"{pfam}\t{sample}\t{protein}\n")

        with open(matrix_file, 'w') as mf_handle:
            mf_handle.write('\t'.join(['Pfam'] + samples) + '\n')
            for i, pfam in enumerate(pfams):
                mf_handle.write(pfam + '\t' + '\t'.join(map(str, counts[i].tolist())) + '\n')

        samples_with_pfam = (counts > 0).sum(axis=1)
        total_proteins = counts.sum(axis=1)
        with open(pfam_info_file, 'w') as pi_handle:
            pi_handle.write('Pfam\tPfam_Name\tSamples_With_Domain\tTotal_Proteins_With_Domain\n')
            for i, pfam in enumerate(pfams):
                pi_handle.write(f"{pfam}\t{pfam_names[pfam]}\t{samples_with_pfam[i]}\t{total_proteins[i]}\n")

        log_object.info(f"Pfam count matrix written to {matrix_file}")

    except Exception as e:
        log_object.error("Error creating Pfam count matrix")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        raise
