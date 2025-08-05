"""
Input/Output functions for bofasa.

This module contains functions for processing genomes, proteomes, and other
data files, including file format validation and data loading.
"""

import gzip
import itertools
import logging
import multiprocessing
import os
import subprocess
import sys
import traceback
import tempfile
import time
from collections import defaultdict
from operator import itemgetter
from typing import Dict, List, Any, Optional

import tqdm
from Bio import SeqIO

from .utils import multi_process, create_locus_tag_options
from . import config

# No global variables needed for per-sample logging

def _get_cds_log_file(outdir: str, sample_name: str) -> str:
    """
    Get the path to the CDS logging file for a specific sample.
    
    Args:
        outdir: Output directory
        sample_name: Name of the sample
        
    Returns:
        Path to the CDS log file for the sample
    """
    return os.path.join(outdir, f"cds_processing_issues_{sample_name}.log")

def _log_cds_issue(outdir: str, sample_name: str, scaffold: str, feature_info: str, issue_type: str, details: str = ""):
    """
    Log CDS processing issues to a sample-specific file.
    
    Args:
        outdir: Output directory
        sample_name: Name of the sample being processed
        scaffold: Scaffold/contig name
        feature_info: Information about the feature (e.g., locus tag, coordinates)
        issue_type: Type of issue (e.g., "skipped", "coordinate_error", "translation_error")
        details: Additional details about the issue
    """
    cds_log_file = _get_cds_log_file(outdir, sample_name)
    
    try:
        # Write directly to the sample-specific log file
        with open(cds_log_file, 'a') as log_handle:
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
            log_entry = f"{timestamp}\t{scaffold}\t{feature_info}\t{issue_type}\t{details}\n"
            log_handle.write(log_entry)
            
    except Exception:
        # If logging fails, just continue - don't let logging break the main process
        pass


def run_individual_gene_calling(
    input_file: str,
    outdir: str,
    sample_name: str,
    gene_calling_method: str = "pyrodigal",
    meta_mode: bool = False,
) -> None:
    """
    Run gene calling using prodigal or pyrodigal.
    
    Args:
        input_file: Path to input FASTA file
        outdir: Output directory
        sample_name: Sample name for output files
        gene_calling_method: Gene calling method ("pyrodigal" or "prodigal")
        meta_mode: Whether to run in metagenomics mode
    """
    if gene_calling_method == "pyrodigal":
        run_pyrodigal_gene_calling(input_file, outdir, sample_name, meta_mode)
    elif gene_calling_method == "prodigal":
        run_prodigal_command(input_file, outdir, sample_name, meta_mode)
    else:
        raise ValueError(f"Unsupported gene calling method: {gene_calling_method}")


def run_pyrodigal_gene_calling(
    input_file: str, outdir: str, sample_name: str, meta_mode: bool = False
) -> None:
    """Run pyrodigal gene calling."""
    try:
        import pyrodigal
    except ImportError:
        raise ImportError(
            "pyrodigal is not installed. Please install it with: pip install pyrodigal"
        )

    # Create pyrodigal gene finder with parameters matching v1.1.1 command line usage
    gene_finder = pyrodigal.GeneFinder(
        meta=meta_mode,
        closed=False,  # v1.1.1 didn't specify -c, so default is False (allows genes to run off edges)
        mask=False,  # v1.1.1 didn't specify -m, so default is False
        min_gene=90,  # Default minimum gene length
        min_edge_gene=60,  # Default minimum edge gene length
        max_overlap=60,  # Default maximum overlap
        backend="detect"  # Use fastest available backend
    )

    # Read sequences
    sequences = []
    with open(input_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequences.append(record)

    # Predict genes with sequence tracking
    genes_with_contigs = []
    if meta_mode:
        # In metagenomic mode, we can directly find genes
        for i, record in enumerate(sequences):
            genes = gene_finder.find_genes(bytes(record.seq))
            for gene in genes:
                genes_with_contigs.append((gene, record.id, i))
    else:
        # In single mode, train on all sequences first, then predict
        if len(sequences) > 0:
            # Train on all sequences to build the model (pass all sequences as separate arguments)
            gene_finder.train(*(bytes(record.seq) for record in sequences))
            
            # Now predict genes on all sequences using the trained model
            for i, record in enumerate(sequences):
                genes = gene_finder.find_genes(bytes(record.seq))
                for gene in genes:
                    genes_with_contigs.append((gene, record.id, i))

    # Write results
    write_gene_predictions(genes_with_contigs, sequences, outdir, sample_name)


def run_prodigal_command(
    input_file: str, outdir: str, sample_name: str, meta_mode: bool = False
) -> None:
    """Run prodigal command-line tool."""
    # Prepare command
    cmd = ['prodigal']
    if meta_mode:
        cmd.append('-p')
        cmd.append('meta')
    else:
        cmd.append('-p')
        cmd.append('single')

    cmd.extend(
        [
            '-i',
            input_file,
            '-o',
            os.path.join(outdir, f"{sample_name}.gff"),
            '-a',
            os.path.join(outdir, f"{sample_name}.faa"),
            '-d',
            os.path.join(outdir, f"{sample_name}.fna"),
            '-f',
            'gff',
        ]
    )

    # Run prodigal
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Prodigal failed: {result.stderr}")


def write_gene_predictions(
    genes_with_contigs: List, sequences: List, outdir: str, sample_name: str
) -> None:
    """Write gene predictions to files."""
    # Create output files
    gff_file = os.path.join(outdir, f"{sample_name}.gff")
    faa_file = os.path.join(outdir, f"{sample_name}.faa")
    fna_file = os.path.join(outdir, f"{sample_name}.fna")

    with open(gff_file, 'w') as gff_handle, open(faa_file, 'w') as faa_handle, open(
        fna_file, 'w'
    ) as fna_handle:

        # Write GFF header
        gff_handle.write("##gff-version 3\n")

        for i, (gene, contig_id, contig_idx) in enumerate(genes_with_contigs):
            # Write GFF entry
            strand = '+' if gene.strand == 1 else '-'
            gff_handle.write(
                f"{contig_id}\tpyrodigal\tCDS\t{gene.begin}\t{gene.end}\t{gene.score}\t{strand}\t0\tID=gene_{i};locus_tag=gene_{i}\n"
            )

            # Write protein sequence
            faa_handle.write(f">gene_{i}\n{gene.translate()}\n")

            # Write DNA sequence
            fna_handle.write(f">gene_{i}\n{gene.sequence()}\n")


def create_formatted_annotation_files(
    outdir: str, sample_name: str, locus_tag: str
) -> None:
    """Create formatted BED and proteome files with custom locus tags from Prodigal output."""
    # Read GFF file
    gff_file = os.path.join(outdir, f"{sample_name}.gff")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")

    # Create output files
    bed_file = os.path.join(outdir, f"{sample_name}.coords.bed")
    proteome_file = os.path.join(outdir, f"{sample_name}.faa")
    name_map_file = os.path.join(outdir, f"{sample_name}.name_map.txt")

    # Process GFF and create formatted files
    locus_tag_counter = 1

    with open(gff_file, 'r') as gff_handle, open(bed_file, 'w') as bed_handle, open(
        name_map_file, 'w'
    ) as map_handle:

        for line in gff_handle:
            if line.startswith('#'):
                continue

            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'CDS':
                continue

            scaffold = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]

            # Create new locus tag
            new_locus_tag = f"{locus_tag}_{locus_tag_counter:06d}"
            locus_tag_counter += 1

            # Write BED entry
            bed_handle.write(
                f"{scaffold}\t{start}\t{end}\t{new_locus_tag}\t1\t{strand}\n"
            )

            # Write name mapping
            old_name = f"gene_{locus_tag_counter-2}"  # Adjust for 0-based indexing
            map_handle.write(f"{old_name}\t{new_locus_tag}\n")

    # Update proteome file with new locus tags
    update_proteome_locus_tags(proteome_file, name_map_file)


def update_proteome_locus_tags(proteome_file: str, name_map_file: str) -> None:
    """Update proteome file with new locus tags with coordinate information."""
    # Read name mapping
    name_map = {}
    with open(name_map_file, 'r') as handle:
        for line in handle:
            old_name, new_name = line.strip().split('\t')
            name_map[old_name] = new_name

    # Get the GFF file path (same directory as proteome file)
    gff_file = proteome_file.replace('.faa', '.gff')
    
    # Read coordinate information from GFF file
    coord_map = {}
    if os.path.exists(gff_file):
        with open(gff_file, 'r') as gff_handle:
            for line in gff_handle:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue
                
                # Extract information from GFF
                scaffold = parts[0]
                start = parts[3]
                end = parts[4]
                strand = parts[6]
                
                # Extract gene ID from attributes
                attributes = parts[8]
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
                
                if gene_id:
                    coord_map[gene_id] = (scaffold, start, end, strand)

    # Create temporary file
    temp_file = proteome_file + '.tmp'

    with open(proteome_file, 'r') as input_handle, open(
        temp_file, 'w'
    ) as output_handle:

        for line in input_handle:
            if line.startswith('>'):
                # Extract old name and replace with new locus tag
                old_name = line.strip()[1:].split()[0]  # Remove '>' and get first part
                if old_name in name_map:
                    new_name = name_map[old_name]
                    # Add coordinate information if available
                    if old_name in coord_map:
                        scaffold, start, end, strand = coord_map[old_name]
                        output_handle.write(f">{new_name} {scaffold} {start} {end} {strand}\n")
                    else:
                        output_handle.write(f">{new_name}\n")
                else:
                    output_handle.write(line)
            else:
                output_handle.write(line)

    # Replace original file
    os.replace(temp_file, proteome_file)


def process_genbank_file(
    input_file: str,
    outdir: str,
    sample_name: str,
    locus_tag: Optional[str] = None,
    min_length: int = 20,
) -> None:
    """
    Process GenBank file and create formatted output files.

    This function silently logs CDS features that are skipped or have processing issues
    to a sample-specific file named 'cds_processing_issues_{sample_name}.log' in the output 
    directory. The log file contains tab-separated entries with timestamp, scaffold, 
    feature info, issue type, and details.
    
    Complex location formats (join/order) are parsed using Biopython functionality,
    extracting the first location segment for processing. CDS features without
    translation qualifiers are logged as issues and skipped.

    Args:
        input_file: Path to input GenBank file
        outdir: Output directory
        sample_name: Sample name for output files
        locus_tag: Locus tag prefix (optional)
        min_length: Minimum protein length
    """
    # Create output file paths
    bed_file = os.path.join(outdir, f"{sample_name}.coords.bed")
    proteome_file = os.path.join(outdir, f"{sample_name}.faa")
    genome_file = os.path.join(outdir, f"{sample_name}.fna")
    name_map_file = os.path.join(outdir, f"{sample_name}.name_map.txt")

    # Initialize counters
    locus_tag_counter = 1
    protein_count = 0

    # Open output files
    with open(bed_file, 'w') as bed_handle, open(
        proteome_file, 'w'
    ) as proteome_handle, open(genome_file, 'w') as genome_handle, open(
        name_map_file, 'w'
    ) as map_handle:

        # Initialize CDS log file with header if it doesn't exist
        cds_log_file = _get_cds_log_file(outdir, sample_name)
        if not os.path.exists(cds_log_file):
            try:
                with open(cds_log_file, 'w') as log_handle:
                    log_handle.write("Timestamp\tScaffold\tFeature_Info\tIssue_Type\tDetails\n")
            except Exception:
                pass  # Don't let logging initialization break the main process

        # Read GenBank file
        genbank_handle = None
        if input_file.endswith('.gz'):
            genbank_handle = gzip.open(input_file, 'rt')
        else:
            genbank_handle = open(input_file, 'r')

        try:
            record_count = 0
            for record in SeqIO.parse(genbank_handle, 'genbank'):
                record_count += 1
                scaffold = record.id
                scaffold_length = len(str(record.seq))
                
                # Count CDS features
                cds_count = sum(1 for feature in record.features if feature.type == "CDS")
                if cds_count == 0:
                    print(f"Warning: No CDS features found in record {scaffold} of {input_file}")
                    _log_cds_issue(outdir, sample_name, scaffold, "N/A", "no_cds_features", f"No CDS features found in record")

                # Write genome sequence
                genome_handle.write(f">{scaffold}\n{str(record.seq)}\n")

                # Process features
                for feature in record.features:
                    if feature.type == "CDS":
                        # Extract coordinates using Biopython functionality
                        start, end, direction, success = _parse_complex_location(feature, record, outdir, sample_name, scaffold)
                        if not success:
                            continue

                        # Extract locus tag and protein sequence
                        old_locus_tag = None
                        protein_sequence = None

                        try:
                            locus_tag_list = feature.qualifiers.get('locus_tag')
                            old_locus_tag = locus_tag_list[0] if locus_tag_list else None
                        except (KeyError, IndexError):
                            old_locus_tag = None

                        # Extract translation from qualifiers only
                        translation_list = feature.qualifiers.get('translation')
                        if translation_list:
                            protein_sequence = str(translation_list[0]).replace('*', '')
                        else:
                            # No translation qualifier - log as issue and skip
                            feature_info = f"no_translation_{old_locus_tag or 'unknown'}"
                            _log_cds_issue(outdir, sample_name, scaffold, feature_info, "no_translation", "CDS feature has no translation qualifier")
                            continue

                        if not protein_sequence or len(protein_sequence) < min_length:
                            feature_info = f"short_protein_{old_locus_tag or 'unknown'}"
                            protein_len = len(protein_sequence) if protein_sequence else 0
                            _log_cds_issue(outdir, sample_name, scaffold, feature_info, "short_protein", f"Protein too short: {protein_len} < {min_length}")
                            continue

                        # Create new locus tag if requested
                        if locus_tag is not None:
                            new_locus_tag = f"{locus_tag}_{locus_tag_counter:06d}"
                            locus_tag_counter += 1
                        else:
                            new_locus_tag = (
                                old_locus_tag or f"gene_{locus_tag_counter:06d}"
                            )
                            locus_tag_counter += 1

                        # Determine protein score based on scaffold position
                        protein_score = '1'
                        if (scaffold_length - end) < 1000 or start < 1000:
                            protein_score = '0'

                        # Write BED entry
                        bed_handle.write(
                            f"{scaffold}\t{start}\t{end}\t{new_locus_tag}\t{protein_score}\t{direction}\n"
                        )

                        # Write protein sequence
                        proteome_handle.write(
                            f">{new_locus_tag} {scaffold} {start} {end} {direction}\n{protein_sequence}\n"
                        )

                        # Write name mapping if locus tag was changed
                        if old_locus_tag and old_locus_tag != new_locus_tag:
                            map_handle.write(f"{old_locus_tag}\t{new_locus_tag}\n")

                        protein_count += 1

        finally:
            genbank_handle.close()

    print(f"Processed {protein_count} proteins from {record_count} records in GenBank file.")
    
    # Check if any CDS issues were logged
    cds_log_file = _get_cds_log_file(outdir, sample_name)
    if os.path.exists(cds_log_file) and os.path.getsize(cds_log_file) > 0:
        # Count lines in log file (excluding header)
        try:
            with open(cds_log_file, 'r') as log_handle:
                lines = log_handle.readlines()
                if len(lines) > 1:  # More than just header
                    issue_count = len(lines) - 1  # Exclude header
                    print(f"CDS processing issues logged to: {cds_log_file} ({issue_count} issues)")
        except Exception:
            pass  # Don't let logging summary break the main process
    
    # Validate that we processed at least some proteins
    if protein_count == 0:
        raise ValueError(f"No proteins were extracted from GenBank file: {input_file}")

    # No cleanup needed for per-sample logging


def validate_genbank_file(input_file: str) -> bool:
    """
    Validate that the input file is a valid GenBank file.

    Args:
        input_file: Path to input file

    Returns:
        bool: True if valid GenBank file, False otherwise
    """
    try:
        if not is_genbank(input_file):
            return False

        # Try to parse the file
        handle = None
        if input_file.endswith('.gz'):
            handle = gzip.open(input_file, 'rt')
        else:
            handle = open(input_file, 'r')

        try:
            records = list(SeqIO.parse(handle, 'genbank'))
            return len(records) > 0
        finally:
            handle.close()
    except Exception:
        return False


def extract_genbank_metadata(input_file: str) -> Dict[str, Any]:
    """
    Extract metadata from GenBank file.

    Args:
        input_file: Path to GenBank file

    Returns:
        dict: Dictionary containing metadata
    """
    metadata = {
        'num_records': 0,
        'num_cds': 0,
        'total_length': 0,
        'sources': set(),
        'organisms': set(),
    }

    handle = None
    if input_file.endswith('.gz'):
        handle = gzip.open(input_file, 'rt')
    else:
        handle = open(input_file, 'r')

    try:
        for record in SeqIO.parse(handle, 'genbank'):
            metadata['num_records'] += 1
            metadata['total_length'] += len(record.seq)

            # Extract organism and source information
            if 'organism' in record.annotations:
                metadata['organisms'].add(record.annotations['organism'])
            if 'source' in record.annotations:
                metadata['sources'].add(record.annotations['source'])

            # Count CDS features
            for feature in record.features:
                if feature.type == "CDS":
                    metadata['num_cds'] += 1
    finally:
        handle.close()

    return metadata


def _process_single_genome(args):
    """
    Process a single genome for gene calling (multiprocessing worker function).
    
    Args:
        args: Tuple containing (sample, sample_assembly, prodigal_outdir, sample_locus_tag, 
              gene_calling_method, meta_mode, log_object)
    
    Returns:
        tuple: (sample, success, error_message)
    """
    sample, sample_assembly, prodigal_outdir, sample_locus_tag, gene_calling_method, meta_mode, log_object = args
    
    try:
        # Run gene calling
        run_individual_gene_calling(
            input_file=sample_assembly,
            outdir=prodigal_outdir,
            sample_name=sample,
            gene_calling_method=gene_calling_method,
            meta_mode=meta_mode
        )
        
        # Create formatted files with locus tags
        create_formatted_annotation_files(
            outdir=prodigal_outdir,
            sample_name=sample,
            locus_tag=sample_locus_tag
        )
        
        return (sample, True, None)
        
    except Exception as e:
        error_msg = f"Error processing genome {sample}: {str(e)}"
        return (sample, False, error_msg)


def run_gene_calling(
    sample_genomes: Dict[str, str],
    prodigal_outdir: str,
    log_object: logging.Logger,
    threads: int = config.DEFAULT_THREADS,
    locus_tag_length: int = config.DEFAULT_LOCUS_TAG_LENGTH,
    gene_calling_method: str = "pyrodigal",
    meta_mode: bool = False,
) -> None:
    """
    Process input genomes using prodigal/pyrodigal for gene prediction and annotation.

    This function oversees processing of input genomes to create proteome and GenBank files
    using p(y)rodigal for gene calling and annotation.

    Args:
        sample_genomes: Dictionary mapping sample identifiers to genome FASTA file paths
        prodigal_outdir: Output directory for prodigal/pyrodigal intermediate and final results
        log_object: Logger object for logging operations
        threads: Number of threads to use for parallel processing (default: 1)
        locus_tag_length: Length of locus tags to generate (default: 3)
        gene_calling_method: Gene calling method to use - "pyrodigal", "prodigal", or "prodigal-gv" (default: "pyrodigal")
        meta_mode: Whether to run in metagenomics mode (default: False)

    Returns:
        None: Creates proteome and GenBank files in the output directory

    Raises:
        Exception: If prodigal/pyrodigal processing fails or required files are not     
                   created
    """
    try:
        possible_locus_tags = create_locus_tag_options(locus_tag_length)

        # Create output directory if it doesn't exist
        os.makedirs(prodigal_outdir, exist_ok=True)

        # Prepare arguments for multiprocessing
        process_args = []
        for i, sample in enumerate(sorted(sample_genomes)):
            sample_assembly = sample_genomes[sample]
            sample_locus_tag = possible_locus_tags[i]
            process_args.append((
                sample, sample_assembly, prodigal_outdir, sample_locus_tag,
                gene_calling_method, meta_mode, log_object
            ))

        # Process genomes in parallel
        if threads > 1 and len(sample_genomes) > 1:
            with multiprocessing.Pool(processes=threads) as pool:
                results = list(tqdm.tqdm(
                    pool.imap(_process_single_genome, process_args),
                    total=len(process_args),
                    desc="Processing genomes"
                ))
        else:
            # Process sequentially if single thread or single genome
            results = []
            for args in tqdm.tqdm(process_args, desc="Processing genomes"):
                results.append(_process_single_genome(args))

        # Check results and report
        successfully_processed = 0
        failed_samples = []
        for sample, success, error_msg in results:
            if success:
                successfully_processed += 1
            else:
                log_object.error(error_msg)
                failed_samples.append(sample)

        msg = f"Successfully processed {successfully_processed} out of {len(sample_genomes)} genomes."
        log_object.info(msg)

        if failed_samples:
            raise Exception(f"Failed to process genomes: {', '.join(failed_samples)}")

    except Exception as e:
        log_object.error(
            f"Problem with running {gene_calling_method} gene calling. Exiting now ..."
        )
        log_object.error(traceback.format_exc())
        sys.exit(1)


def _process_single_genbank(args):
    """
    Process a single GenBank file (multiprocessing worker function).
    
    Args:
        args: Tuple containing (sample, sample_genbank, gp_dir, sample_locus_tag, log_object)
    
    Returns:
        tuple: (sample, success, error_message)
    """
    sample, sample_genbank, gp_dir, sample_locus_tag, log_object = args
    
    try:
        # Process GenBank file
        process_genbank_file(
            input_file=sample_genbank,
            outdir=gp_dir,
            sample_name=sample,
            locus_tag=sample_locus_tag,
            min_length=20
        )
        
        # Verify output files
        faa_file = os.path.join(gp_dir, f"{sample}.faa")
        bed_file = os.path.join(gp_dir, f"{sample}.coords.bed")
        fna_file = os.path.join(gp_dir, f"{sample}.fna")
        
        if (os.path.isfile(faa_file) and os.path.isfile(bed_file) and os.path.isfile(fna_file) and
            os.path.getsize(faa_file) > 0 and os.path.getsize(bed_file) > 0 and os.path.getsize(fna_file) > 0):
            return (sample, True, None)
        else:
            error_msg = f"Sample {sample} was not processed successfully."
            return (sample, False, error_msg)
            
    except Exception as e:
        error_msg = f"Error processing sample {sample}: {str(e)}"
        return (sample, False, error_msg)


def process_genomes_as_genbanks(
    sample_genomes: Dict[str, str],
    gp_dir: str,
    log_object: logging.Logger,
    threads: int = config.DEFAULT_THREADS,
    locus_tag_length: int = config.DEFAULT_LOCUS_TAG_LENGTH,
    rename_locus_tags: bool = False,
) -> None:
    """
    Process input genomes provided as GenBank files with existing CDS features.

    This function oversees processing of input genomes as GenBanks with CDS features
    already available, extracting proteomes and creating standardized formats.
    
    CDS processing issues are silently logged to individual files 'cds_processing_issues_{sample_name}.log' 
    in the output directory, with a summary provided at the end of processing.

    Args:
        sample_genomes: Dictionary mapping sample identifiers to GenBank file paths
        gp_dir: Output directory for processed genome files
        log_object: Logger object for logging operations
        threads: Number of threads to use for parallel processing (default: 1)
        locus_tag_length: Length of locus tags to generate (default: 3)
        rename_locus_tags: Whether to rename existing locus tags (default: False)

    Returns:
        dict: Dictionary mapping sample names to paths of processed sample files

    Raises:
        Exception: If GenBank processing fails or required files are not created
    """

    try:
        # Create output directory if it doesn't exist
        os.makedirs(gp_dir, exist_ok=True)

        # Generate locus tags if needed
        locus_tags = None
        if rename_locus_tags:
            locus_tags = create_locus_tag_options(locus_tag_length)

        msg = f"Attempting to process/re-format {len(sample_genomes)} genomes provided as GenBank files using {threads} threads"
        log_object.info(msg)

        # Prepare arguments for multiprocessing
        process_args = []
        for i, sample in enumerate(sorted(sample_genomes)):
            sample_genbank = sample_genomes[sample]
            sample_locus_tag = locus_tags[i] if rename_locus_tags else None
            process_args.append((sample, sample_genbank, gp_dir, sample_locus_tag, log_object))

        # Process genomes in parallel
        if threads > 1 and len(sample_genomes) > 1:
            with multiprocessing.Pool(processes=threads) as pool:
                results = list(tqdm.tqdm(
                    pool.imap(_process_single_genbank, process_args),
                    total=len(process_args),
                    desc="Processing GenBank files"
                ))
        else:
            # Process sequentially if single thread or single genome
            results = []
            for args in tqdm.tqdm(process_args, desc="Processing GenBank files"):
                results.append(_process_single_genbank(args))

        # Check results and report
        successfully_processed = 0
        failed_samples = []
        for sample, success, error_msg in results:
            if success:
                successfully_processed += 1
            else:
                log_object.error(error_msg)
                failed_samples.append(sample)

        msg = f"Successfully processed {successfully_processed} out of {len(sample_genomes)} genomes."
        log_object.info(msg)

        if failed_samples:
            log_object.warning(f"Failed to process samples: {', '.join(failed_samples)}")

    except Exception as e:
        log_object.error(
            "Problem with creating commands for processing GenBank files. Exiting now ..."
        )
        log_object.error(traceback.format_exc())
        sys.exit(1)

    _summarize_cds_issues(gp_dir, log_object)


def extract_gene_contexts(inputs: List[Any]) -> None:
    """
    Extract gene contexts for a single sample from coordinate information.

    Args:
        inputs: Tuple containing (sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, log_object)

    Returns:
        None: Writes gene context information to output file

    Raises:
        Exception: If gene context extraction fails
    """
    (
        sample,
        coords_file,
        output_file,
        gene_to_og,
        og_genes,
        surrounding_bp,
        log_object,
    ) = inputs
    try:
        assert os.path.isfile(coords_file)

        scaffold_features = defaultdict(list)
        with open(coords_file) as ocf:
            for line in ocf:
                line = line.strip()
                scaffold, start, end, name, score, strand = line.split('\t')
                start = int(start)
                end = int(end)
                scaffold_features[scaffold].append([start, end, name, score, strand])

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
            scaffold_features_sorted = sorted(
                scaffold_features[scaffold], key=itemgetter(0)
            )
            max_features = len(scaffold_features[scaffold])
            for cds_index, cds in enumerate(scaffold_features_sorted):
                cds_start = cds[0]
                cds_end = cds[1]
                cds_name = cds[2]
                if not cds_name in gene_to_og:
                    continue
                cds_og = gene_to_og[cds_name]
                left_boundary = cds_start - surrounding_bp
                right_boundary = cds_end + surrounding_bp
                left_side_genes_and_ogs = set([])
                right_side_genes_and_ogs = set([])
                near_scaffold_edge = False
                if cds[3] == '0':
                    near_scaffold_edge = True

                limit_reached = False
                cds_iter_index = cds_index - 1
                while not limit_reached:
                    try:
                        cds_iter = scaffold_features_sorted[cds_iter_index]
                        if cds_iter_index < 0:
                            near_scaffold_edge = True
                            limit_reached = True
                        elif cds_iter[1] < left_boundary:
                            if cds_iter[3] == '0':
                                near_scaffold_edge = True
                            limit_reached = True
                        elif cds_iter[0] < left_boundary:
                            if cds_iter[3] == '0':
                                near_scaffold_edge = True
                            left_side_genes_and_ogs.add(
                                tuple(
                                    [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                )
                            )
                            limit_reached = True
                        else:
                            left_side_genes_and_ogs.add(
                                tuple(
                                    [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                )
                            )
                    except Exception:
                        near_scaffold_edge = True
                        limit_reached = True
                    cds_iter_index -= 1

                limit_reached = False
                cds_iter_index = cds_index + 1
                while not limit_reached:
                    try:
                        cds_iter = scaffold_features_sorted[cds_iter_index]
                        if cds_iter_index > max_features - 1:
                            near_scaffold_edge = True
                            limit_reached = True
                        elif cds_iter[0] > right_boundary:
                            if cds_iter[3] == '0':
                                near_scaffold_edge = True
                            limit_reached = True
                        elif cds_iter[1] > right_boundary:
                            if cds_iter[3] == '0':
                                near_scaffold_edge = True
                            right_side_genes_and_ogs.add(
                                tuple(
                                    [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                )
                            )
                            limit_reached = True
                        else:
                            right_side_genes_and_ogs.add(
                                tuple(
                                    [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                )
                            )
                    except Exception:
                        near_scaffold_edge = True
                        limit_reached = True
                    cds_iter_index += 1

                upstream_genes = []
                downstream_genes = []
                upstream_ogs = []
                downstream_ogs = []
                if cds[2] == '1':
                    for go in sorted(left_side_genes_and_ogs, key=itemgetter(2)):
                        upstream_genes.append(go[0])
                        upstream_ogs.append(go[1])
                    for go in sorted(right_side_genes_and_ogs, key=itemgetter(2)):
                        downstream_genes.append(go[0])
                        downstream_ogs.append(go[1])
                else:
                    for go in sorted(right_side_genes_and_ogs, key=itemgetter(2)):
                        upstream_genes.append(go[0])
                        upstream_ogs.append(go[1])
                    for go in sorted(left_side_genes_and_ogs, key=itemgetter(2)):
                        downstream_genes.append(go[0])
                        downstream_ogs.append(go[1])
                output_handle.write(
                    '\t'.join(
                        [
                            cds_name,
                            cds_og,
                            str(near_scaffold_edge),
                            str(abs(cds_end - cds_start + 1)),
                            ', '.join(upstream_genes),
                            ', '.join(downstream_genes),
                            ', '.join(upstream_ogs),
                            ', '.join(downstream_ogs),
                        ]
                    )
                    + '\n'
                )
        output_handle.close()

    except Exception as e:
        log_object.error(
            "Problem with determining context of genes for sample %s" % sample
        )
        log_object.error(traceback.format_exc())
        sys.exit(1)


def is_fasta(fasta: str) -> bool:
    """
    Check if a file is in FASTA format.

    Args:
        fasta: Path to the file to check

    Returns:
        bool: True if file is in FASTA format, False otherwise
    """
    try:
        with open(fasta) as of:
            for i, line in enumerate(of):
                if i == 0:
                    if not line.startswith('>'):
                        return False
                    break
        return True
    except Exception:
        return False


def is_genbank(gbk: str, check_for_cds: bool = False) -> bool:
    """
    Check if a file is in GenBank format.

    Args:
        gbk: Path to the file to check
        check_for_cds: Whether to verify CDS features are present (default: False)

    Returns:
        bool: True if file is in GenBank format, False otherwise
    """
    try:
        recs = 0
        cds_flag = False
        if (
            gbk.endswith(".gbk")
            or gbk.endswith(".gb")
            or gbk.endswith(".gbff")
            or gbk.endswith(".gbk.gz")
            or gbk.endswith(".gb.gz")
            or gbk.endswith(".gbff.gz")
        ):
            if gbk.endswith(".gz"):
                with gzip.open(gbk, "rt") as ogf:
                    for rec in SeqIO.parse(ogf, "genbank"):
                        if check_for_cds:
                            for feature in rec.features:
                                if feature.type == "CDS":
                                    cds_flag = True
                        if not check_for_cds or cds_flag:
                            recs += 1
                            break
            else:
                with open(gbk) as ogf:
                    for rec in SeqIO.parse(ogf, "genbank"):
                        if check_for_cds:
                            for feature in rec.features:
                                if feature.type == "CDS":
                                    cds_flag = True
                        if not check_for_cds or cds_flag:
                            recs += 1
                            break
            if recs > 0:
                return True
            else:
                return False
        else:
            return False
    except Exception:
        return False


def create_chopped_proteomes(inputs: List[Any]) -> None:
    """
    Create chopped proteome files based on Pfam domain annotations.

    This function creates chopped CDS proteome files from regular proteome files
    by identifying Pfam domains and splitting proteins at domain boundaries.

    Args:
        inputs: Tuple containing (prot_file, ccds_prot_file, dom_coord_file, pfam_db_file, pfam_z, minimal_length, log_object, threads, skip_domain_splitting)

    Returns:
        None: Creates chopped proteome files and domain coordinate information

    Raises:
        Exception: If proteome processing fails
    """
    (
        prot_file,
        ccds_prot_file,
        dom_coord_file,
        pfam_db_file,
        pfam_z,
        minimal_length,
        log_object,
        threads,
        skip_domain_splitting,
    ) = inputs

    try:
        sample = '.'.join(prot_file.split('/')[-1].split('.')[:-1])

        if not skip_domain_splitting:
            # Check if pyhmmer is available
            try:
                import pyhmmer
            except ImportError:
                msg = f"pyhmmer is required for domain splitting but not available for {sample}"
                log_object.error(msg)
                raise ImportError(msg)

            # Align Pfam domains and remove overlap similar to BiG-SCAPE
            alphabet = pyhmmer.easel.Alphabet.amino()
            sequences = []
            with pyhmmer.easel.SequenceFile(
                prot_file, digital=True, alphabet=alphabet
            ) as seq_file:
                sequences = list(seq_file)

            target_dom_hits = defaultdict(list)
            with pyhmmer.plan7.HMMFile(pfam_db_file) as hmm_file:
                for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="trusted", Z=int(pfam_z), cpus=threads):
                    for hit in hits:
                        for domain in hit.domains.included:
                            target_dom_hits[hit.name.decode()].append(
                                [
                                    hits.query.name.decode(),
                                    domain.alignment.target_from,
                                    domain.alignment.target_to,
                                    domain.score,
                                    domain.i_evalue,
                                ]
                            )

            # Chop up FASTA based on mostly non-overlapping domains, 10% leeway is given
            breakpoints = defaultdict(list)
            dom_start_names = defaultdict(lambda: 'NA')
            for tg in target_dom_hits:
                tg_dom_name_iter = defaultdict(int)
                accounted_coords = set([])
                for dom_align_info in sorted(
                    target_dom_hits[tg], key=itemgetter(3), reverse=True
                ):
                    dom_name, start, end, score, i_evalue = dom_align_info
                    overlap_coords = accounted_coords.intersection(
                        set(range(start, end + 1))
                    )
                    if (
                        len(overlap_coords) / float(end - start + 1) >= 0.1
                        or len(overlap_coords) >= minimal_length
                    ):
                        continue
                    accounted_coords = accounted_coords.union(
                        set(range(start, end + 1))
                    )
                    breakpoints[tg].append(start)
                    breakpoints[tg].append(end + 1)
                    dom_start_names[tg + '|' + str(start)] = (
                        tg + '|' + dom_name + '|' + str(tg_dom_name_iter[dom_name] + 1)
                    )
                    tg_dom_name_iter[dom_name] += 1

            cpf_handle = open(ccds_prot_file, 'w')
            dcf_handle = open(dom_coord_file, 'w')
            dcf_handle.write(
                'Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n'
            )
            with open(prot_file) as ocf:
                for rec in SeqIO.parse(ocf, 'fasta'):
                    tg = rec.id
                    tg_seq = str(rec.seq)
                    prev_end_coord = 1
                    tg_interdomain_index = 1
                    if not tg in breakpoints and len(tg_seq) >= minimal_length:
                        cpf_handle.write('>' + sample + '|' + tg + '|full_protein|1' + '\n' + str(tg_seq) + '\n')
                        dcf_handle.write(
                            '\t'.join(
                                [
                                    sample,
                                    tg,
                                    'full_protein',
                                    '1',
                                    str(prev_end_coord),
                                    str(len(tg_seq)),
                                ]
                            )
                            + '\n'
                        )
                    else:
                        for tg_seq_chunk in split_by_idx(
                            tg_seq, ([0] + sorted(breakpoints[tg]))
                        ):
                            if tg_seq_chunk.strip() == '':
                                continue
                            end_coord = prev_end_coord + len(tg_seq_chunk) - 1
                            if len(tg_seq_chunk) >= minimal_length:
                                dn = (
                                    sample
                                    + '|'
                                    + dom_start_names[
                                        tg + '|' + str(prev_end_coord - 1)
                                    ]
                                )
                                if (
                                    dom_start_names[tg + '|' + str(prev_end_coord - 1)]
                                    == 'NA'
                                ):
                                    dn = (
                                        sample
                                        + '|'
                                        + tg
                                        + '|inter-domain_region|'
                                        + str(tg_interdomain_index)
                                    )
                                    tg_interdomain_index += 1
                                if (
                                    dom_start_names[tg + '|' + str(prev_end_coord - 1)]
                                    == 'NA'
                                ):
                                    domain_name = f"inter-domain_region"
                                    tg_interdomain_index += 1
                                else:
                                    domain_info = dom_start_names[tg + '|' + str(prev_end_coord - 1)]
                                    domain_name = domain_info.split('|')[1]
                                cpf_handle.write(
                                    '>' + sample + '|' + tg + '|' + domain_name + '|' + (str(tg_interdomain_index) if dom_start_names[tg + '|' + str(prev_end_coord - 1)] == 'NA' else domain_info.split('|')[2]) + '\n' + str(tg_seq_chunk) + '\n'
                                )
                                dcf_handle.write(
                                    '\t'.join(
                                        [
                                            sample,
                                            tg,
                                            domain_name,
                                            '1',
                                            str(prev_end_coord),
                                            str(end_coord),
                                        ]
                                    )
                                    + '\n'
                                )
                            prev_end_coord = end_coord + 1
            cpf_handle.close()
            dcf_handle.close()
        else:
            # Skip domain splitting - create full protein entries
            cpf_handle = open(ccds_prot_file, 'w')
            dcf_handle = open(dom_coord_file, 'w')
            dcf_handle.write(
                'Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n'
            )
            with open(prot_file) as ocf:
                for rec in SeqIO.parse(ocf, 'fasta'):
                    tg = rec.id
                    tg_seq = str(rec.seq)
                    tg_interdomain_index = 1
                    if len(tg_seq) >= minimal_length:
                        cpf_handle.write('>' + sample + '|' + tg + '|full_protein|1' + '\n' + str(tg_seq) + '\n')
                        dcf_handle.write(
                            '\t'.join(
                                [
                                    sample,
                                    tg,
                                    'full_protein',
                                    '1',
                                    '1',
                                    str(len(tg_seq)),
                                ]
                            )
                            + '\n'
                        )
            cpf_handle.close()
            dcf_handle.close()

    except Exception as e:
        msg = f'An issue occurred with creating chopped up version of proteome file {prot_file}.'
        log_object.error(msg)
        log_object.error(traceback.format_exc())
        raise e


def split_by_idx(S: List[Any], list_of_indices: List[int]) -> List[List[Any]]:
    """
    Split a string by a list of indices.

    Function taken from https://stackoverflow.com/questions/10851445/splitting-a-string-by-list-of-indices

    Args:
        S: String to split
        list_of_indices: List of indices to split at

    Yields:
        str: Substrings split at the specified indices
    """
    left, right = 0, list_of_indices[0]
    yield S[left:right]
    left = right
    for right in list_of_indices[1:]:
        yield S[left:right]
        left = right
    yield S[left:]


def load_table_in_pandas_dataframe(
    input_file: str, numeric_columns: List[str], cut_last_columns: Optional[int] = None
) -> Any:
    """
    Load a table into a pandas DataFrame with numeric column handling.

    Args:
        input_file: Path to the input table file
        numeric_columns: List of column names to treat as numeric
        cut_last_columns: Number of columns to remove from the end (default: None)

    Returns:
        pandas.DataFrame: Loaded data with proper numeric column types, or None if loading fails
    """
    try:
        import pandas as pd

        df = pd.read_csv(input_file, sep="\t")

        if cut_last_columns:
            df = df.iloc[:, :-cut_last_columns]

        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        return df

    except Exception as e:
        # Note: This function doesn't have access to log_object, so we keep stderr for this case
        sys.stderr.write(f"Error loading table: {str(e)}\n")
        return None


def process_prokka_directory(
    prokka_dir: str,
    outdir: str,
    sample_name: str,
    locus_tag: Optional[str] = None,
    min_length: int = 20,
) -> None:
    """
    Process Prokka annotation directory and create formatted output files.

    Args:
        prokka_dir: Path to Prokka output directory
        outdir: Output directory
        sample_name: Sample name for output files
        locus_tag: Locus tag prefix (optional)
        min_length: Minimum protein length
    """
    # Find Prokka output files
    gff_file = None
    faa_file = None
    fna_file = None
    
    # Look for Prokka output files
    for file in os.listdir(prokka_dir):
        if file.endswith('.gff'):
            gff_file = os.path.join(prokka_dir, file)
        elif file.endswith('.faa'):
            faa_file = os.path.join(prokka_dir, file)
        elif file.endswith('.fna'):
            fna_file = os.path.join(prokka_dir, file)
    
    if not gff_file:
        raise FileNotFoundError(f"No GFF file found in Prokka directory: {prokka_dir}")
    
    # Create output file paths
    bed_file = os.path.join(outdir, f"{sample_name}.coords.bed")
    proteome_file = os.path.join(outdir, f"{sample_name}.faa")
    genome_file = os.path.join(outdir, f"{sample_name}.fna")

    # Initialize counters and mappings
    locus_tag_counter = 1
    protein_count = 0
    sample_old_to_new_lts = {}
    prot_id_counts = defaultdict(int)

    # Get scaffold lengths if genome file is available
    scaffold_lengths = {}
    if fna_file:
        with open(fna_file, 'r') as fna_handle:
            for record in SeqIO.parse(fna_handle, 'fasta'):
                scaffold = record.id
                scaffold_lengths[scaffold] = len(str(record.seq))

    # Process GFF file to create BED file
    with open(bed_file, 'w') as bed_handle:
        with open(gff_file, 'r') as gff_handle:
            for line in gff_handle:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue

                scaffold = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]

                # Parse attributes
                attributes = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value

                # Extract locus tag
                old_locus_tag = attributes.get('ID', '').replace('locus_tag=', '')
                if not old_locus_tag:
                    continue

                final_lt = old_locus_tag
                
                # Create new locus tag if requested (matching previous version's format)
                if locus_tag is not None:
                    if (locus_tag_counter + 1) < 10:
                        pid = '00000' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 100:
                        pid = '0000' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 1000:
                        pid = '000' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 10000:
                        pid = '00' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 100000:
                        pid = '0' + str(locus_tag_counter + 1)
                    else:
                        pid = str(locus_tag_counter + 1)
                    
                    new_lt = locus_tag + '_' + pid
                    sample_old_to_new_lts[old_locus_tag] = [new_lt, scaffold, str(start), str(end), strand]
                    final_lt = new_lt
                
                locus_tag_counter += 1
                
                # Determine if near scaffold edge
                prot_score = '1'
                if scaffold in scaffold_lengths:
                    if (scaffold_lengths[scaffold] - end) < 1000 or start < 1000:
                        prot_score = '0'
                elif start < 1000:
                    prot_score = '0'
                
                prot_id_counts[final_lt] += 1
                
                # Write BED entry
                bed_handle.write(
                    f"{scaffold}\t{start}\t{end}\t{final_lt}\t{prot_score}\t{strand}\n"
                )

    # Check for duplicate locus tags
    for pi in prot_id_counts:
        if prot_id_counts[pi] != 1:
            raise ValueError(f"Not all protein locus tags across genomes are unique! Please try running bofasa_prep again with the --rename-locus-tags option.")

    # Process FAA file to create proteome file
    if faa_file:
        with open(proteome_file, 'w') as proteome_handle:
            with open(faa_file, 'r') as faa_handle:
                for record in SeqIO.parse(faa_handle, 'fasta'):
                    old_lt = record.id
                    final_lt = old_lt
                    
                    # Use new locus tag if available (matching previous version's format)
                    if locus_tag is not None and old_lt in sample_old_to_new_lts:
                        final_lt = ' '.join(sample_old_to_new_lts[old_lt])
                    
                    # Check minimum length
                    if len(str(record.seq)) >= min_length:
                        proteome_handle.write(f">{final_lt}\n{str(record.seq)}\n")
                        protein_count += 1

    # Copy genome sequence if available
    if fna_file:
        with open(fna_file, 'r') as fna_handle:
            with open(genome_file, 'w') as genome_handle:
                genome_handle.write(fna_handle.read())

    print(f"Processed {protein_count} proteins from Prokka directory.")


def process_bakta_directory(
    bakta_dir: str,
    outdir: str,
    sample_name: str,
    locus_tag: Optional[str] = None,
    min_length: int = 20,
) -> None:
    """
    Process Bakta annotation directory and create formatted output files.

    Args:
        bakta_dir: Path to Bakta output directory
        outdir: Output directory
        sample_name: Sample name for output files
        locus_tag: Locus tag prefix (optional)
        min_length: Minimum protein length
    """
    # Find Bakta output files
    gff_file = None
    faa_file = None
    fna_file = None
    
    # Look for Bakta output files
    for file in os.listdir(bakta_dir):
        if file.endswith('.gff3'):
            gff_file = os.path.join(bakta_dir, file)
        elif file.endswith('.faa') and not file.endswith('.hypotheticals.faa'):
            faa_file = os.path.join(bakta_dir, file)
        elif file.endswith('.fna'):
            fna_file = os.path.join(bakta_dir, file)
    
    if not gff_file:
        raise FileNotFoundError(f"No GFF3 file found in Bakta directory: {bakta_dir}")
    
    # Create output file paths
    bed_file = os.path.join(outdir, f"{sample_name}.coords.bed")
    proteome_file = os.path.join(outdir, f"{sample_name}.faa")
    genome_file = os.path.join(outdir, f"{sample_name}.fna")

    # Initialize counters and mappings
    locus_tag_counter = 1
    protein_count = 0
    sample_old_to_new_lts = {}
    prot_id_counts = defaultdict(int)

    # Get scaffold lengths if genome file is available
    scaffold_lengths = {}
    if fna_file:
        with open(fna_file, 'r') as fna_handle:
            for record in SeqIO.parse(fna_handle, 'fasta'):
                scaffold = record.id
                scaffold_lengths[scaffold] = len(str(record.seq))

    # Process GFF3 file to create BED file
    with open(bed_file, 'w') as bed_handle:
        with open(gff_file, 'r') as gff_handle:
            for line in gff_handle:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue

                scaffold = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]

                # Parse attributes
                attributes = {}
                for attr in parts[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key] = value

                # Extract locus tag
                old_locus_tag = attributes.get('ID', '').replace('locus_tag=', '')
                if not old_locus_tag:
                    continue

                final_lt = old_locus_tag
                
                # Create new locus tag if requested (matching previous version's format)
                if locus_tag is not None:
                    if (locus_tag_counter + 1) < 10:
                        pid = '00000' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 100:
                        pid = '0000' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 1000:
                        pid = '000' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 10000:
                        pid = '00' + str(locus_tag_counter + 1)
                    elif (locus_tag_counter + 1) < 100000:
                        pid = '0' + str(locus_tag_counter + 1)
                    else:
                        pid = str(locus_tag_counter + 1)
                    
                    new_lt = locus_tag + '_' + pid
                    sample_old_to_new_lts[old_locus_tag] = [new_lt, scaffold, str(start), str(end), strand]
                    final_lt = new_lt
                
                locus_tag_counter += 1
                
                # Determine if near scaffold edge
                prot_score = '1'
                if scaffold in scaffold_lengths:
                    if (scaffold_lengths[scaffold] - end) < 1000 or start < 1000:
                        prot_score = '0'
                elif start < 1000:
                    prot_score = '0'
                
                prot_id_counts[final_lt] += 1
                
                # Write BED entry
                bed_handle.write(
                    f"{scaffold}\t{start}\t{end}\t{final_lt}\t{prot_score}\t{strand}\n"
                )

    # Check for duplicate locus tags
    for pi in prot_id_counts:
        if prot_id_counts[pi] != 1:
            raise ValueError(f"Not all protein locus tags across genomes are unique! Please try running bofasa_prep again with the --rename-locus-tags option.")

    # Process FAA file to create proteome file
    if faa_file:
        with open(proteome_file, 'w') as proteome_handle:
            with open(faa_file, 'r') as faa_handle:
                for record in SeqIO.parse(faa_handle, 'fasta'):
                    old_lt = record.id
                    final_lt = old_lt
                    
                    # Use new locus tag if available (matching previous version's format)
                    if locus_tag is not None and old_lt in sample_old_to_new_lts:
                        final_lt = ' '.join(sample_old_to_new_lts[old_lt])
                    
                    # Check minimum length
                    if len(str(record.seq)) >= min_length:
                        proteome_handle.write(f">{final_lt}\n{str(record.seq)}\n")
                        protein_count += 1

    # Copy genome sequence if available
    if fna_file:
        with open(fna_file, 'r') as fna_handle:
            with open(genome_file, 'w') as genome_handle:
                genome_handle.write(fna_handle.read())

    print(f"Processed {protein_count} proteins from Bakta directory.")


def extract_protein_sequence(faa_file: str, locus_tag: str) -> Optional[str]:
    """
    Extract protein sequence for a specific locus tag from FAA file.
    
    Args:
        faa_file: Path to FAA file
        locus_tag: Locus tag to search for
        
    Returns:
        Protein sequence or None if not found
    """
    try:
        with open(faa_file, 'r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                if locus_tag in record.id:
                    return str(record.seq)
    except Exception:
        pass
    return None


def detect_annotation_type(annotation_dir: str) -> str:
    """
    Detect the type of annotation directory (prokka or bakta).
    
    Args:
        annotation_dir: Path to annotation directory
        
    Returns:
        'prokka' or 'bakta' or raises ValueError if unknown
    """
    if not os.path.isdir(annotation_dir):
        raise ValueError(f"Directory does not exist: {annotation_dir}")
    
    files = os.listdir(annotation_dir)
    
    # Check for Prokka files
    if any(f.endswith('.gff') for f in files):
        return 'prokka'
    
    # Check for Bakta files
    if any(f.endswith('.gff3') for f in files):
        return 'bakta'
    
    raise ValueError(f"Could not determine annotation type for directory: {annotation_dir}")


def _process_single_annotation_dir(args):
    """
    Process a single annotation directory (multiprocessing worker function).
    
    Args:
        args: Tuple containing (annotation_dir, outdir, sample_locus_tag, log_object)
    
    Returns:
        tuple: (sample_name, success, error_message, annotation_type)
    """
    annotation_dir, outdir, sample_locus_tag, log_object = args
    
    try:
        # Detect annotation type
        annotation_type = detect_annotation_type(annotation_dir)
        
        # Generate sample name from directory name
        sample_name = os.path.basename(annotation_dir.rstrip('/'))
        
        if annotation_type == 'prokka':
            process_prokka_directory(
                prokka_dir=annotation_dir,
                outdir=outdir,
                sample_name=sample_name,
                locus_tag=sample_locus_tag,
                min_length=20
            )
        elif annotation_type == 'bakta':
            process_bakta_directory(
                bakta_dir=annotation_dir,
                outdir=outdir,
                sample_name=sample_name,
                locus_tag=sample_locus_tag,
                min_length=20
            )
        
        # Verify output files
        faa_file = os.path.join(outdir, f"{sample_name}.faa")
        bed_file = os.path.join(outdir, f"{sample_name}.coords.bed")
        fna_file = os.path.join(outdir, f"{sample_name}.fna")
        
        if (os.path.isfile(faa_file) and os.path.isfile(bed_file) and os.path.isfile(fna_file) and
            os.path.getsize(faa_file) > 0 and os.path.getsize(bed_file) > 0 and os.path.getsize(fna_file) > 0):
            return (sample_name, True, None, annotation_type)
        else:
            error_msg = f"Annotation directory {sample_name} was not processed successfully."
            return (sample_name, False, error_msg, annotation_type)
            
    except Exception as e:
        error_msg = f"Error processing annotation directory {os.path.basename(annotation_dir.rstrip('/'))}: {str(e)}"
        return (os.path.basename(annotation_dir.rstrip('/')), False, error_msg, "unknown")


def process_annotation_directories(
    annotation_dirs: List[str],
    outdir: str,
    log_object: logging.Logger,
    locus_tag_length: int = config.DEFAULT_LOCUS_TAG_LENGTH,
    rename_locus_tags: bool = False,
    threads: int = config.DEFAULT_THREADS,
) -> Dict[str, Dict[str, str]]:
    """
    Process annotation directories (Prokka/Bakta) and create formatted output files.
    
    Args:
        annotation_dirs: List of paths to annotation directories
        outdir: Output directory
        log_object: Logger object
        locus_tag_length: Length of locus tag prefix
        rename_locus_tags: Whether to rename locus tags
        threads: Number of threads to use
        
    Returns:
        Dictionary containing sample mappings for wgs, proteomes, and beds
    """
    # Create locus tag options if renaming is requested
    possible_locustags = None
    if rename_locus_tags:
        possible_locustags = create_locus_tag_options(locus_tag_length)
    
    # Prepare arguments for multiprocessing
    process_args = []
    for i, annotation_dir in enumerate(annotation_dirs):
        sample_locus_tag = possible_locustags[i] if rename_locus_tags else None
        process_args.append((annotation_dir, outdir, sample_locus_tag, log_object))
    
    # Process annotation directories in parallel
    sample_wgs = {}
    sample_proteomes = {}
    sample_beds = {}
    
    if log_object:
        log_object.info(f"Processing {len(annotation_dirs)} annotation directories...")
    
    # Use multiprocessing if multiple directories
    if len(annotation_dirs) > 1 and threads > 1:
        with multiprocessing.Pool(threads) as pool:
            results = list(tqdm.tqdm(
                pool.imap(_process_single_annotation_dir, process_args),
                total=len(process_args),
                desc="Processing annotation directories"
            ))
    else:
        results = []
        for args in process_args:
            results.append(_process_single_annotation_dir(args))
    
    # Collect results
    for sample_name, success, error_message, annotation_type in results:
        if success:
            # Define expected file paths
            faa_file = os.path.join(outdir, f"{sample_name}.faa")
            bed_file = os.path.join(outdir, f"{sample_name}.coords.bed")
            fna_file = os.path.join(outdir, f"{sample_name}.fna")
            
            # Verify files exist and have content
            if (os.path.isfile(faa_file) and os.path.isfile(bed_file) and 
                os.path.getsize(faa_file) > 0 and os.path.getsize(bed_file) > 0):
                
                sample_proteomes[sample_name] = faa_file
                sample_beds[sample_name] = bed_file
                
                # Add genome file if it exists
                if os.path.isfile(fna_file) and os.path.getsize(fna_file) > 0:
                    sample_wgs[sample_name] = fna_file
                
                if log_object:
                    log_object.info(f"Successfully processed {annotation_type} directory: {sample_name}")
            else:
                if log_object:
                    log_object.error(f"Missing or empty output files for {sample_name}")
        else:
            if log_object:
                log_object.error(f"Failed to process annotation directory {sample_name}: {error_message}")
    
    return {
        'sample_wgs': sample_wgs,
        'sample_proteomes': sample_proteomes,
        'sample_beds': sample_beds
    }


def _summarize_cds_issues(outdir: str, log_object: logging.Logger):
    """
    Provide a summary of CDS processing issues across all processed genomes.
    
    Args:
        outdir: Output directory containing the CDS log files
        log_object: Logger object for logging operations
    """
    try:
        # Find all CDS log files
        cds_log_files = []
        for file in os.listdir(outdir):
            if file.startswith("cds_processing_issues_") and file.endswith(".log"):
                cds_log_files.append(os.path.join(outdir, file))
        
        if not cds_log_files:
            return
        
        # Parse issues from all log files
        issue_counts = defaultdict(int)
        sample_issues = defaultdict(int)
        total_issues = 0
        
        for log_file in cds_log_files:
            sample_name = os.path.basename(log_file).replace("cds_processing_issues_", "").replace(".log", "")
            
            try:
                with open(log_file, 'r') as log_handle:
                    lines = log_handle.readlines()
                    if len(lines) <= 1:  # Only header or empty
                        continue
                    
                    sample_issue_count = len(lines) - 1  # Exclude header
                    sample_issues[sample_name] = sample_issue_count
                    total_issues += sample_issue_count
                    
                    # Parse issue types for this sample
                    for line in lines[1:]:  # Skip header
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            issue_type = parts[3]
                            issue_counts[issue_type] += 1
                            
            except Exception:
                continue  # Skip files that can't be read
        
        # Log summary
        if total_issues > 0:
            msg = f"CDS processing summary: {total_issues} total issues across {len(sample_issues)} samples"
            log_object.info(msg)
            
            # Log issue type breakdown
            for issue_type, count in sorted(issue_counts.items()):
                msg = f"  - {issue_type}: {count} issues"
                log_object.info(msg)
            
            # Log sample breakdown
            for sample, count in sorted(sample_issues.items()):
                if count > 0:
                    msg = f"  - {sample}: {count} issues"
                    log_object.info(msg)
            
            msg = f"Detailed CDS processing issues logged to individual files: cds_processing_issues_*.log"
            log_object.info(msg)
                
    except Exception as e:
        log_object.warning(f"Could not generate CDS processing summary: {str(e)}")


def _parse_complex_location(feature, record, outdir, sample_name, scaffold):
    """
    Parse complex location strings (join/order) using Biopython functionality.
    
    Args:
        feature: Biopython feature object
        record: Biopython record object
        outdir: Output directory for logging
        sample_name: Sample name for logging
        scaffold: Scaffold name for logging
        
    Returns:
        tuple: (start, end, direction, success) or (None, None, None, False) if failed
    """
    try:
        # Try to use Biopython's location parsing for simple locations
        if hasattr(feature.location, 'start') and hasattr(feature.location, 'end'):
            start = feature.location.start.position + 1  # Convert to 1-based
            end = feature.location.end.position
            direction = '+' if feature.location.strand == 1 else '-'
            return start, end, direction, True
    except (AttributeError, TypeError):
        pass
    
    try:
        # For simple locations that don't have start/end attributes, try string parsing
        location_str = str(feature.location)
        
        # Skip complex locations for now - they'll be handled below
        if 'join{' in location_str or 'order{' in location_str:
            pass
        else:
            # Simple location format: [start:end](strand)
            parts = location_str.split(']')[0].replace('[', '').replace(']', '').split(':')
            if len(parts) == 2:
                start = int(parts[0]) + 1  # Convert to 1-based
                end = int(parts[1])
                direction = '+' if '(' not in location_str or '(+)' in location_str else '-'
                return start, end, direction, True
            else:
                feature_info = f"location_parse_error_{location_str[:50]}"
                _log_cds_issue(outdir, sample_name, scaffold, feature_info, "location_parse_error", f"Could not parse location format: {location_str}")
                return None, None, None, False
    except (ValueError, IndexError) as e:
        location_str = str(feature.location)
        feature_info = f"location_parse_exception_{location_str[:50]}"
        _log_cds_issue(outdir, sample_name, scaffold, feature_info, "location_parse_exception", f"Exception parsing location: {location_str} - {str(e)}")
        return None, None, None, False
    
    try:
        # For complex locations, try to extract the first part
        location_str = str(feature.location)
        
        # Handle join{...} format
        if 'join{' in location_str:
            # Extract the first location from join
            join_start = location_str.find('join{') + 5
            join_end = location_str.rfind('}')
            if join_start < join_end:
                join_content = location_str[join_start:join_end]
                # Find the first location in the join
                first_loc_start = join_content.find('[')
                first_loc_end = join_content.find(']', first_loc_start)
                if first_loc_start >= 0 and first_loc_end > first_loc_start:
                    first_location = join_content[first_loc_start:first_loc_end + 1]
                    # Parse the first location
                    parts = first_location.replace('[', '').replace(']', '').split(':')
                    if len(parts) == 2:
                        start = int(parts[0]) + 1  # Convert to 1-based
                        end = int(parts[1])
                        direction = '+' if '(' not in location_str or '(+)' in location_str else '-'
                        return start, end, direction, True
        
        # Handle order{...} format
        elif 'order{' in location_str:
            # Extract the first location from order
            order_start = location_str.find('order{') + 6
            order_end = location_str.rfind('}')
            if order_start < order_end:
                order_content = location_str[order_start:order_end]
                # Find the first location in the order
                first_loc_start = order_content.find('[')
                first_loc_end = order_content.find(']', first_loc_start)
                if first_loc_start >= 0 and first_loc_end > first_loc_start:
                    first_location = order_content[first_loc_start:first_loc_end + 1]
                    # Parse the first location
                    parts = first_location.replace('[', '').replace(']', '').split(':')
                    if len(parts) == 2:
                        start = int(parts[0]) + 1  # Convert to 1-based
                        end = int(parts[1])
                        direction = '+' if '(' not in location_str or '(+)' in location_str else '-'
                        return start, end, direction, True
        
        # If we can't parse it, log as complex location issue
        feature_info = f"complex_location_{location_str[:50]}"
        _log_cds_issue(outdir, sample_name, scaffold, feature_info, "complex_location", f"Complex location format: {location_str}")
        return None, None, None, False
        
    except Exception as e:
        location_str = str(feature.location)
        feature_info = f"location_parse_exception_{location_str[:50]}"
        _log_cds_issue(outdir, sample_name, scaffold, feature_info, "location_parse_exception", f"Exception parsing location: {location_str} - {str(e)}")
        return None, None, None, False



