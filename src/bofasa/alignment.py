"""
Alignment functions for bofasa.

This module contains functions for creating multiple sequence alignments
using various alignment tools and methods.
"""

import os
import subprocess
import sys
import traceback
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import tqdm
from Bio import SeqIO
from pyfamsa import Aligner, Sequence

from .utils import run_cmd, multi_process
from . import config


def create_protein_alignments_pyfamsa(
    prot_dir: str,
    prot_algn_dir: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,
    guide_tree: str = config.DEFAULT_PYFAMSA_GUIDE_TREE,
    tree_heuristic: Optional[str] = config.DEFAULT_PYFAMSA_HEURISTIC,
    n_refinements: int = config.DEFAULT_PYFAMSA_REFINEMENTS,
    keep_duplicates: bool = False,
    refine: bool = False,
) -> None:
    """
    Create protein alignments using PyFAMSA.

    Parameters:
    -----------
    prot_dir : str
        Directory containing protein files in FASTA format
    prot_algn_dir : str
        Directory where aligned protein files will be saved
    log_object : logging.Logger
        Logger object for recording progress and errors
    threads : int, default=1
        Number of threads to use for parallel processing
    guide_tree : str, default="sl"
        Guide tree method: "sl" (single linkage), "slink", "upgma", "nj" (neighbor joining)
    tree_heuristic : Optional[str], default=None
        Tree heuristic: None, "medoid", or "part"
    n_refinements : int, default=100
        Number of refinement iterations
    keep_duplicates : bool, default=False
        Whether to keep duplicate sequences
    """
    try:
        # Create output directory if it doesn't exist
        os.makedirs(prot_algn_dir, exist_ok=True)

        # Get list of protein files
        protein_files: List[str] = [f for f in os.listdir(prot_dir) if f.endswith('.faa')]

        if not protein_files:
            log_object.warning(f"No protein files found in {prot_dir}")
            return

        msg: str = (
            f"Creating protein alignments using PyFAMSA for {len(protein_files)} files"
        )
        log_object.info(msg)

        # Create aligner with specified parameters
        aligner: Aligner = Aligner(
            threads=threads,
            guide_tree=guide_tree,
            tree_heuristic=tree_heuristic,
            n_refinements=n_refinements,
            keep_duplicates=keep_duplicates,
            refine=refine,
        )

        # Track statistics
        skipped_single_sequences = 0
        successful_alignments = 0
        
        # Process each protein file
        for protein_file in tqdm.tqdm(protein_files, desc="Aligning proteins"):
            try:
                # Input and output file paths
                prot_file: str = os.path.join(prot_dir, protein_file)
                prefix: str = protein_file.replace('.faa', '')
                prot_algn_file: str = os.path.join(prot_algn_dir, f"{prefix}.msa.faa")

                # Skip if output already exists
                if os.path.exists(prot_algn_file):
                    log_object.info(
                        f"Skipping {protein_file} - alignment already exists (checkpoint)"
                    )
                    continue

                # Load sequences from FASTA file
                sequences: List[Sequence] = []
                with open(prot_file, 'r') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        seq_id: bytes = record.id.encode('utf-8')
                        seq_data: bytes = str(record.seq).encode('utf-8')
                        sequences.append(Sequence(seq_id, seq_data))

                if len(sequences) < 2:
                    skipped_single_sequences += 1
                    continue

                # Perform alignment
                alignment: List[Sequence] = aligner.align(sequences)

                # Write alignment
                with open(prot_algn_file, 'w') as handle:
                    for gapped_seq in alignment:
                        seq_id: str = gapped_seq.id.decode('utf-8')
                        seq_data: str = gapped_seq.sequence.decode('utf-8')
                        handle.write(f">{seq_id}\n{seq_data}\n")

                successful_alignments += 1
                log_object.debug(
                    f"Successfully aligned {protein_file} ({len(sequences)} sequences)"
                )

            except Exception as e:
                log_object.error(f"Error aligning {protein_file}: {str(e)}")
                log_object.error(traceback.format_exc())
                continue

        # Report statistics
        total_processed = len(protein_files)
        log_object.info(f"Protein alignments completed using PyFAMSA: {successful_alignments}/{total_processed} files successfully aligned, {skipped_single_sequences} skipped (single sequences)")

    except Exception as e:
        log_object.error("Issues with creating protein alignments using PyFAMSA")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)


def create_protein_alignments_parallel_pyfamsa(
    prot_dir: str,
    prot_algn_dir: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,
    guide_tree: str = config.DEFAULT_PYFAMSA_GUIDE_TREE,
    tree_heuristic: Optional[str] = config.DEFAULT_PYFAMSA_HEURISTIC,
    n_refinements: int = config.DEFAULT_PYFAMSA_REFINEMENTS,
    keep_duplicates: bool = False,
    refine: bool = False,
) -> None:
    """
    Create protein alignments using PyFAMSA with parallel processing.

    Parameters:
    -----------
    prot_dir : str
        Directory containing protein files in FASTA format
    prot_algn_dir : str
        Directory where aligned protein files will be saved
    log_object : logging.Logger
        Logger object for recording progress and errors
    threads : int, default=1
        Number of threads to use for parallel processing
    guide_tree : str, default="sl"
        Guide tree method: "sl" (single linkage), "slink", "upgma", "nj" (neighbor joining)
    tree_heuristic : Optional[str], default=None
        Tree heuristic: None, "medoid", or "part"
    n_refinements : int, default=100
        Number of refinement iterations
    keep_duplicates : bool, default=False
        Whether to keep duplicate sequences
    refine : bool, default=False
        Whether to enable refinement for higher quality alignments
    """
    try:
        # Create output directory if it doesn't exist
        os.makedirs(prot_algn_dir, exist_ok=True)

        # Get list of protein files
        protein_files: List[str] = [f for f in os.listdir(prot_dir) if f.endswith('.faa')]

        if not protein_files:
            log_object.warning(f"No protein files found in {prot_dir}")
            return

        msg: str = (
            f"Creating protein alignments using PyFAMSA (parallel) for {len(protein_files)} files"
        )
        log_object.info(msg)

        # Prepare tasks for parallel processing
        tasks: List[Dict[str, Any]] = []
        for protein_file in protein_files:
            prot_file: str = os.path.join(prot_dir, protein_file)
            prefix: str = protein_file.replace('.faa', '')
            prot_algn_file: str = os.path.join(prot_algn_dir, f"{prefix}.msa.faa")

            # Skip if output already exists
            if os.path.exists(prot_algn_file):
                log_object.info(
                    f"Skipping {protein_file} - alignment already exists (checkpoint)"
                )
                continue

            task: Dict[str, Any] = {
                'input_file': prot_file,
                'output_file': prot_algn_file,
                'threads': threads,
                'guide_tree': guide_tree,
                'tree_heuristic': tree_heuristic,
                'n_refinements': n_refinements,
                'keep_duplicates': keep_duplicates,
                'refine': refine,
                'log_object': log_object
            }
            tasks.append(task)

        if not tasks:
            log_object.info("No new alignments to create")
            return

        # Process tasks in parallel
        with multi_process(threads) as pool:
            results: List[bool] = list(
                tqdm.tqdm(
                    pool.imap_unordered(_align_single_file, tasks),
                    total=len(tasks),
                    desc="Aligning proteins (parallel)"
                )
            )

        successful: int = sum(results)
        skipped_single_sequences = len(tasks) - successful
        log_object.info(f"Successfully aligned {successful}/{len(tasks)} protein files, {skipped_single_sequences} skipped (single sequences)")

    except Exception as e:
        log_object.error("Issues with creating protein alignments using PyFAMSA (parallel)")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)


def _align_single_file(task: Dict[str, Any]) -> bool:
    """
    Align a single protein file using PyFAMSA.

    Parameters:
    -----------
    task : Dict[str, Any]
        Dictionary containing alignment task parameters

    Returns:
    --------
    bool
        True if alignment was successful, False otherwise
    """
    try:
        input_file: str = task['input_file']
        output_file: str = task['output_file']
        threads: int = task['threads']
        guide_tree: str = task['guide_tree']
        tree_heuristic: Optional[str] = task['tree_heuristic']
        n_refinements: int = task['n_refinements']
        keep_duplicates: bool = task['keep_duplicates']
        refine: bool = task['refine']
        log_object: Any = task['log_object']

        # Load sequences from FASTA file
        sequences: List[Sequence] = []
        with open(input_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                seq_id: bytes = record.id.encode('utf-8')
                seq_data: bytes = str(record.seq).encode('utf-8')
                sequences.append(Sequence(seq_id, seq_data))

        if len(sequences) < 2:
            return False

        # Create aligner
        aligner: Aligner = Aligner(
            threads=threads,
            guide_tree=guide_tree,
            tree_heuristic=tree_heuristic,
            n_refinements=n_refinements,
            keep_duplicates=keep_duplicates,
            refine=refine,
        )

        # Perform alignment
        alignment: List[Sequence] = aligner.align(sequences)

        # Write alignment
        with open(output_file, 'w') as handle:
            for gapped_seq in alignment:
                seq_id: str = gapped_seq.id.decode('utf-8')
                seq_data: str = gapped_seq.sequence.decode('utf-8')
                handle.write(f">{seq_id}\n{seq_data}\n")

        return True

    except Exception as e:
        log_object.error(f"Error aligning {os.path.basename(input_file)}: {str(e)}")
        return False


def create_domain_protein_alignments_pyfamsa(
    dog_seqs_dir: str,
    dog_algn_dir: str,
    log_object: Any,
    threads: int = config.DEFAULT_THREADS,
    guide_tree: str = config.DEFAULT_PYFAMSA_GUIDE_TREE,
    tree_heuristic: Optional[str] = config.DEFAULT_PYFAMSA_HEURISTIC,
    n_refinements: int = config.DEFAULT_PYFAMSA_REFINEMENTS,
    keep_duplicates: bool = False,
    refine: bool = False,
) -> None:
    """
    Create domain protein alignments using PyFAMSA.

    This function aligns domain-specific protein sequences for coarse domain
    ortholog groups.

    Parameters:
    -----------
    dog_seqs_dir : str
        Directory containing DOG protein sequence files
    dog_algn_dir : str
        Directory where aligned DOG protein files will be saved
    log_object : logging.Logger
        Logger object for recording progress and errors
    threads : int, default=1
        Number of threads to use for parallel processing
    guide_tree : str, default="sl"
        Guide tree method: "sl" (single linkage), "slink", "upgma", "nj" (neighbor joining)
    tree_heuristic : Optional[str], default=None
        Tree heuristic: None, "medoid", or "part"
    n_refinements : int, default=100
        Number of refinement iterations
    keep_duplicates : bool, default=False
        Whether to keep duplicate sequences
    refine : bool, default=False
        Whether to enable refinement for higher quality alignments
    """
    try:
        # Create output directory if it doesn't exist
        os.makedirs(dog_algn_dir, exist_ok=True)

        # Get list of DOG sequence files
        dog_files: List[str] = [f for f in os.listdir(dog_seqs_dir) if f.endswith('.fa')]

        if not dog_files:
            log_object.warning(f"No DOG sequence files found in {dog_seqs_dir}")
            return

        msg: str = f"Creating domain protein alignments using PyFAMSA for {len(dog_files)} files"
        log_object.info(msg)

        # Create aligner
        aligner: Aligner = Aligner(
            threads=threads,
            guide_tree=guide_tree,
            tree_heuristic=tree_heuristic,
            n_refinements=n_refinements,
            keep_duplicates=keep_duplicates,
            refine=refine,
        )

        # Track statistics
        skipped_single_sequences = 0
        successful_alignments = 0
        
        # Process each DOG file
        for dog_file in tqdm.tqdm(dog_files, desc="Aligning domain proteins"):
            try:
                # Input and output file paths
                seq_file: str = os.path.join(dog_seqs_dir, dog_file)
                prefix: str = dog_file.replace('.fa', '')
                algn_file: str = os.path.join(dog_algn_dir, f"{prefix}.msa.faa")

                # Skip if output already exists
                if os.path.exists(algn_file):
                    log_object.info(f"Skipping {dog_file} - alignment already exists (checkpoint)")
                    continue

                # Load sequences
                sequences: List[Sequence] = []
                with open(seq_file, 'r') as handle:
                    for record in SeqIO.parse(handle, 'fasta'):
                        seq_id: bytes = record.id.encode('utf-8')
                        seq_data: bytes = str(record.seq).encode('utf-8')
                        sequences.append(Sequence(seq_id, seq_data))

                if len(sequences) < 2:
                    skipped_single_sequences += 1
                    continue

                # Perform alignment
                alignment: List[Sequence] = aligner.align(sequences)

                # Check if alignment was successful
                if not alignment or len(alignment) == 0:
                    log_object.error(f"PyFAMSA failed to create alignment for {dog_file}")
                    continue

                if len(alignment) != len(sequences):
                    log_object.warning(f"Alignment length mismatch for {dog_file}: expected {len(sequences)}, got {len(alignment)}")

                # Write alignment
                with open(algn_file, 'w') as handle:
                    for gapped_seq in alignment:
                        seq_id: str = gapped_seq.id.decode('utf-8')
                        seq_data: str = gapped_seq.sequence.decode('utf-8')
                        handle.write(f">{seq_id}\n{seq_data}\n")

                # Verify file was created and is not empty
                if not os.path.exists(algn_file):
                    log_object.error(f"Alignment file was not created: {algn_file}")
                    continue

                if os.path.getsize(algn_file) == 0:
                    log_object.error(f"Alignment file is empty: {algn_file}")
                    continue

                # Verify alignment content
                try:
                    with open(algn_file, 'r') as handle:
                        alignment_records = list(SeqIO.parse(handle, 'fasta'))
                        if len(alignment_records) == 0:
                            log_object.error(f"Alignment file contains no sequences: {algn_file}")
                            continue
                        if len(alignment_records) != len(sequences):
                            log_object.warning(f"Alignment file sequence count mismatch for {dog_file}: expected {len(sequences)}, got {len(alignment_records)}")
                except Exception as e:
                    log_object.error(f"Error reading alignment file {algn_file}: {str(e)}")
                    continue

                successful_alignments += 1
                log_object.debug(
                    f"Successfully aligned {dog_file} ({len(sequences)} sequences)"
                )

            except Exception as e:
                log_object.error(f"Error aligning {dog_file}: {str(e)}")
                log_object.error(traceback.format_exc())
                continue

        # Report statistics
        total_processed = len(dog_files)
        log_object.info(f"Domain protein alignments completed using PyFAMSA: {successful_alignments}/{total_processed} files successfully aligned, {skipped_single_sequences} skipped (single sequences)")

        if successful_alignments == 0:
            log_object.error("No alignments were successfully created!")
            raise RuntimeError("PyFAMSA failed to create any valid alignments")

    except Exception as e:
        log_object.error("Issues with creating domain protein alignments using PyFAMSA")
        log_object.error(str(e))
        log_object.error(traceback.format_exc())
        sys.exit(1)


def assess_alignment_complexity(protein_file: str) -> bool:
    """
    Assess the complexity of a protein alignment job.

    Parameters:
    -----------
    protein_file : str
        Path to the protein FASTA file

    Returns:
    --------
    bool
        True if the job is considered complex/heavy, False otherwise
    """
    try:
        sequence_count: int = 0
        total_length: int = 0
        
        with open(protein_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                sequence_count += 1
                total_length += len(record.seq)
        
        # Consider job complex if:
        # - More than 100 sequences, or
        # - Average sequence length > 500 amino acids, or
        # - Total alignment size > 50,000 amino acids
        avg_length: float = total_length / sequence_count if sequence_count > 0 else 0
        
        return (sequence_count > 100 or 
                avg_length > 500 or 
                total_length > 50000)
                
    except Exception:
        # Default to simple job if assessment fails
        return False
