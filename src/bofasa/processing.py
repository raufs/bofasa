"""
Input/Output functions for bofasa.

This module contains functions for processing genomes, proteomes, and other
data files, including file format validation and data loading.
"""

from curses import raw
import gzip
import logging
import os
from re import M
import shutil
import subprocess
import sys
import traceback
import time
from collections import defaultdict
from operator import itemgetter
from typing import Dict, List, Any, Optional
from Bio import SeqIO
from .utils import create_locus_tag_options, _iter_progress, create_pool, thread_imap
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

    shutil.copy(input_file, os.path.join(outdir, f"{sample_name}.fna"))

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

    cmd.extend(
        [
            '-i',
            input_file,
            '-o',
            os.path.join(outdir, f"{sample_name}.gff"),
            '-a',
            os.path.join(outdir, f"{sample_name}.faa"),
            '-f',
            'gff',
        ]
    )

    # Run prodigal
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Prodigal failed: {result.stderr}")

    # Match the pyrodigal path: copy the original genome (not gene sequences)
    shutil.copy(input_file, os.path.join(outdir, f"{sample_name}.fna"))

    # Rename prodigal gene IDs (GFF ID=1_1, FASTA >contig_1 # ... # ID=1_1;...)
    # into the same gene_i scheme written by the pyrodigal path so downstream
    # formatting maps GFF entries and proteins consistently.
    _normalize_prodigal_output(outdir, sample_name)


def _normalize_prodigal_output(outdir: str, sample_name: str) -> None:
    """Rewrite prodigal GFF/FAA output to use pyrodigal-style gene_i identifiers."""
    gff_file = os.path.join(outdir, f"{sample_name}.gff")
    faa_file = os.path.join(outdir, f"{sample_name}.faa")

    id_map = {}
    gff_lines = []
    gene_index = 0
    with open(gff_file, 'r') as gff_handle:
        for line in gff_handle:
            parts = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(parts) < 9 or parts[2] != 'CDS':
                gff_lines.append(line)
                continue
            new_id = f"gene_{gene_index}"
            gene_index += 1
            other_attrs = []
            for attr in parts[8].split(';'):
                if attr.startswith('ID='):
                    id_map[attr[3:]] = new_id
                elif attr:
                    other_attrs.append(attr)
            parts[8] = ';'.join([f"ID={new_id}", f"locus_tag={new_id}"] + other_attrs)
            gff_lines.append('\t'.join(parts) + '\n')

    with open(gff_file, 'w') as gff_handle:
        gff_handle.writelines(gff_lines)

    records = []
    with open(faa_file, 'r') as faa_handle:
        for record in SeqIO.parse(faa_handle, 'fasta'):
            prodigal_id = None
            for field in record.description.split('#'):
                field = field.strip()
                if field.startswith('ID='):
                    prodigal_id = field.split(';')[0][3:]
                    break
            if prodigal_id not in id_map:
                raise ValueError(f"Could not match prodigal protein {record.id} to a GFF entry for sample {sample_name}")
            records.append((id_map[prodigal_id], str(record.seq)))

    with open(faa_file, 'w') as faa_handle:
        for new_id, seq in records:
            faa_handle.write(f">{new_id}\n{seq}\n")


def write_gene_predictions(
    genes_with_contigs: List, sequences: List, outdir: str, sample_name: str
) -> None:
    """Write gene predictions to files."""
    # Create output files
    gff_file = os.path.join(outdir, f"{sample_name}.gff")
    faa_file = os.path.join(outdir, f"{sample_name}.faa")

    with open(gff_file, 'w') as gff_handle, open(faa_file, 'w') as faa_handle:

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
            
            # Note: We don't write individual gene DNA sequences to FNA
            # The original genome scaffolds should be copied separately


def create_formatted_annotation_files(
    outdir: str, sample_name: str, locus_tag: str, min_length: int = config.DEFAULT_MIN_LENGTH,
    input_genome_file: str = None
) -> None:
    """Create formatted BED and proteome files with custom locus tags from Prodigal output.
    
    Args:
        outdir: Output directory
        sample_name: Sample name
        locus_tag: Locus tag prefix
        min_length: Minimum protein length
        input_genome_file: Path to original input genome file (to copy to Genomes/)
    """
    # Read GFF file
    gff_file = os.path.join(outdir, f"{sample_name}.gff")
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")

    bed_dir = os.path.join(outdir, "BEDs/")
    faa_dir = os.path.join(outdir, "Proteomes/")
    wgs_dir = os.path.join(outdir, "Genomes/")

    # Copy original genome file to Genomes/ directory (not gene sequences!)
    if input_genome_file and os.path.exists(input_genome_file):
        shutil.copy(input_genome_file, os.path.join(wgs_dir, f"{sample_name}.fna"))

    # Create output files
    bed_file = os.path.join(bed_dir, f"{sample_name}.bed")
    proteome_file = os.path.join(faa_dir, f"{sample_name}.faa")
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

            protein_len = (end-start+1)/3
            if protein_len < min_length:
                continue

            # Original gene ID (e.g. gene_0) from the GFF attributes, used to
            # map the raw protein to its new locus tag
            old_name = None
            for attr in parts[8].split(';'):
                if attr.startswith('ID='):
                    old_name = attr[3:]
                    break

            # Create new locus tag
            new_locus_tag = f"{locus_tag}_{locus_tag_counter:06d}"
            locus_tag_counter += 1

            # Write BED entry
            bed_handle.write(
                f"{scaffold}\t{start}\t{end}\t{new_locus_tag}\t1\t{strand}\n"
            )

            # Write name mapping
            if old_name:
                map_handle.write(f"{old_name}\t{new_locus_tag}\n")

    # Update proteome file with new locus tags
    raw_proteome_file = os.path.join(outdir, f"{sample_name}.faa")
    update_proteome_locus_tags(proteome_file, raw_proteome_file, name_map_file)


def update_proteome_locus_tags(proteome_file: str, raw_proteome_file: str, name_map_file: str) -> None:
    """Update proteome file with new locus tags with coordinate information."""
    # Read name mapping
    name_map = {}
    with open(name_map_file, 'r') as handle:
        for line in handle:
            old_name, new_name = line.strip().split('\t')
            name_map[old_name] = new_name

    # Get the GFF file path (same directory as proteome file)
    gff_file = raw_proteome_file.replace('.faa', '.gff')
    
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

    with open(raw_proteome_file, 'r') as input_handle, open(proteome_file, 'w') as output_handle:
        keep = False
        for line in input_handle:
            if line.startswith('>'):
                # Extract old name and replace with new locus tag
                old_name = line.strip()[1:].split()[0]  # Remove '>' and get first part
                # Proteins absent from the name map were filtered out of the
                # BED (e.g. by min_length) and are dropped to keep both in sync
                keep = old_name in name_map
                if keep:
                    new_name = name_map[old_name]
                    # Add coordinate information if available
                    if old_name in coord_map:
                        scaffold, start, end, strand = coord_map[old_name]
                        output_handle.write(f">{new_name} {scaffold} {start} {end} {strand}\n")
                    else:
                        output_handle.write(f">{new_name}\n")
            elif keep:
                output_handle.write(line)


def process_genbank_file(
    input_file: str,
    outdir: str,
    sample_name: str,
    locus_tag: Optional[str] = None,
    min_length: int = config.DEFAULT_MIN_LENGTH,
) -> None:
    """
    Process GenBank file and create formatted output files.

    This function silently logs CDS features that are skipped or have processing issues
    to a sample-specific file named 'cds_processing_issues_{sample_name}.log' in the output 
    directory. The log file contains tab-separated entries with timestamp, scaffold, 
    feature info, issue type, and details.
    
    Complex location formats (join/order) are parsed using Biopython functionality
    and represented by their overall span (see _parse_complex_location). CDS features without
    translation qualifiers are logged as issues and skipped.

    Args:
        input_file: Path to input GenBank file
        outdir: Output directory
        sample_name: Sample name for output files
        locus_tag: Locus tag prefix (optional)
        min_length: Minimum protein length
    """
    # Create output file paths
    bed_file = os.path.join(outdir, f"{sample_name}.bed")
    proteome_file = os.path.join(outdir, f"{sample_name}.faa")
    genome_file = os.path.join(outdir, f"{sample_name}.fna")
    name_map_file = os.path.join(outdir, f"{sample_name}.name_map.txt")

    # Initialize counters
    locus_tag_counter = 1
    protein_count = 0
    prot_id_counts = defaultdict(int)

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
                    sys.stdout.write(f"Warning: No CDS features found in record {scaffold} of {input_file}\n")
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
                            # Fallback includes the sample name so CDSs lacking a
                            # locus_tag do not collide across samples
                            new_locus_tag = (
                                old_locus_tag or f"{sample_name}_gene_{locus_tag_counter:06d}"
                            )
                            locus_tag_counter += 1
                        prot_id_counts[new_locus_tag] += 1

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

    # Check for duplicate locus tags (same handling as the Prokka/Bakta paths)
    for pi in prot_id_counts:
        if prot_id_counts[pi] != 1:
            raise ValueError(f"Not all protein locus tags across genomes are unique! Please try running bofasa_prep again with the --rename-locus-tags option.")

    sys.stdout.write(f"Processed {protein_count} proteins from {record_count} records in GenBank file.\n")
    
    # Check if any CDS issues were logged
    cds_log_file = _get_cds_log_file(outdir, sample_name)
    if os.path.exists(cds_log_file) and os.path.getsize(cds_log_file) > 0:
        # Count lines in log file (excluding header)
        try:
            with open(cds_log_file, 'r') as log_handle:
                lines = log_handle.readlines()
                if len(lines) > 1:  # More than just header
                    issue_count = len(lines) - 1  # Exclude header
                    sys.stdout.write(f"CDS processing issues logged to: {cds_log_file} ({issue_count} issues)\n")
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
              gene_calling_method, meta_mode, log_object, min_length    )
    
    Returns:
        tuple: (sample, success, error_message)
    """
    sample, sample_assembly, prodigal_outdir, sample_locus_tag, gene_calling_method, meta_mode, log_object, min_length = args
    
    try:
        # Run gene calling
        run_individual_gene_calling(
            input_file=sample_assembly,
            outdir=prodigal_outdir,
            sample_name=sample,
            gene_calling_method=gene_calling_method,
            meta_mode=meta_mode
        )
        
        # Create formatted files with locus tags and copy original genome file
        create_formatted_annotation_files(
            outdir=prodigal_outdir,
            sample_name=sample,
            locus_tag=sample_locus_tag,
            min_length=min_length,
            input_genome_file=sample_assembly
        )
        
        return (sample, True, None)
        
    except Exception as e:
        error_msg = f"Error processing genome {sample}: {str(e)}"
        return (sample, False, error_msg)


def _get_available_locus_tags(
    locus_tag_length: int, outdir: str, current_samples: List[str], locus_tag_offset: int = 0
) -> List[str]:
    """
    Get locus tag prefixes not already used by other samples' BED files.

    Prefixes used in BED files already present in outdir or outdir/BEDs/ (from
    samples processed earlier in the same run, e.g. FASTA inputs before GenBank
    inputs) are excluded so mixed input types do not receive duplicate prefixes.
    BED files belonging to current_samples are ignored (they will be rewritten).

    Args:
        locus_tag_length: Length of locus tag prefixes
        outdir: Genome processing output directory
        current_samples: Sample names about to be (re)processed
        locus_tag_offset: Number of prefixes to skip at the start of the options list

    Returns:
        List of available locus tag prefixes, in order
    """
    possible_locus_tags = create_locus_tag_options(locus_tag_length)[locus_tag_offset:]
    current = set(current_samples)
    used_prefixes = set()
    for bed_dir in [outdir, os.path.join(outdir, "BEDs")]:
        if not os.path.isdir(bed_dir):
            continue
        for bed_name in os.listdir(bed_dir):
            if not bed_name.endswith('.bed') or bed_name[:-4] in current:
                continue
            try:
                with open(os.path.join(bed_dir, bed_name)) as bed_handle:
                    for line in bed_handle:
                        fields = line.rstrip('\n').split('\t')
                        if len(fields) >= 4:
                            used_prefixes.add(fields[3].rsplit('_', 1)[0])
            except OSError:
                continue
    available = [lt for lt in possible_locus_tags if lt not in used_prefixes]
    if len(available) < len(current):
        raise ValueError(
            f"Not enough unused locus tag prefixes of length {locus_tag_length} for {len(current)} samples."
        )
    return available


def run_gene_calling(
    sample_genomes: Dict[str, str],
    prodigal_outdir: str,
    log_object: logging.Logger,
    threads: int = config.DEFAULT_THREADS,
    locus_tag_length: int = config.DEFAULT_LOCUS_TAG_LENGTH,
    gene_calling_method: str = "pyrodigal",
    meta_mode: bool = False,
    min_length: int = config.DEFAULT_MIN_LENGTH,
    locus_tag_offset: int = 0,
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
        min_length: Minimum length of protein to include (default: 20)
        locus_tag_offset: Number of locus tag prefixes to skip (default: 0); prefixes
                          used by other samples' BEDs already in prodigal_outdir are
                          always skipped
    Returns:
        None: Creates proteome and GenBank files in the output directory

    Raises:
        Exception: If prodigal/pyrodigal processing fails or required files are not     
                   created
    """
    try:
        possible_locus_tags = _get_available_locus_tags(
            locus_tag_length, prodigal_outdir, list(sample_genomes), locus_tag_offset
        )

        # Prepare arguments for multiprocessing
        process_args = []
        for i, sample in enumerate(sorted(sample_genomes)):
            sample_assembly = sample_genomes[sample]
            sample_locus_tag = possible_locus_tags[i]
            process_args.append((
                sample, sample_assembly, prodigal_outdir, sample_locus_tag,
                gene_calling_method, meta_mode, log_object, min_length
            ))

        # Process genomes in parallel
        if threads > 1 and len(sample_genomes) > 1:
            # Threads: pyrodigal/prodigal run outside the GIL.
            results = list(
                _iter_progress(
                    thread_imap(_process_single_genome, process_args, threads),
                    total=len(process_args),
                    description="Processing genomes",
                )
            )
        else:
            # Process sequentially if single thread or single genome
            results = []
            for args in _iter_progress(process_args, description="Processing genomes"):
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
        args: Tuple containing (sample, sample_genbank, gp_dir, sample_locus_tag, min_length)
    
    Returns:
        tuple: (sample, success, error_message)
    """
    sample, sample_genbank, gp_dir, sample_locus_tag, min_length = args
    
    try:
        # Process GenBank file
        process_genbank_file(
            input_file=sample_genbank,
            outdir=gp_dir,
            sample_name=sample,
            locus_tag=sample_locus_tag,
            min_length=min_length
        )
        
        # Verify output files
        faa_file = os.path.join(gp_dir, f"{sample}.faa")
        bed_file = os.path.join(gp_dir, f"{sample}.bed")
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
    min_length: int = config.DEFAULT_MIN_LENGTH,
    locus_tag_offset: int = 0,
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
        min_length: Minimum length of protein to include
        locus_tag_offset: Number of locus tag prefixes to skip when renaming (default: 0);
                          prefixes used by other samples' BEDs already in gp_dir are
                          always skipped
    Returns:
        dict: Dictionary mapping sample names to paths of processed sample files

    Raises:
        Exception: If GenBank processing fails or required files are not created
    """

    try:
        # Generate locus tags if needed
        locus_tags = None
        if rename_locus_tags:
            locus_tags = _get_available_locus_tags(
                locus_tag_length, gp_dir, list(sample_genomes), locus_tag_offset
            )

        msg = f"Attempting to process/re-format {len(sample_genomes)} genomes provided as GenBank files using {threads} threads"
        log_object.info(msg)

        # Prepare arguments for multiprocessing
        process_args = []
        for i, sample in enumerate(sorted(sample_genomes)):
            sample_genbank = sample_genomes[sample]
            sample_locus_tag = locus_tags[i] if rename_locus_tags else None
            process_args.append((sample, sample_genbank, gp_dir, sample_locus_tag, min_length))

        # Process genomes in parallel
        if threads > 1 and len(sample_genomes) > 1:
            with create_pool(threads) as pool:
                results = list(
                    _iter_progress(
                        pool.imap(_process_single_genbank, process_args),
                        total=len(process_args),
                        description="Processing GenBank files",
                    )
                )
        else:
            # Process sequentially if single thread or single genome
            results = []
            for args in _iter_progress(process_args, description="Processing GenBank files"):
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
                            # Neighbours without an OG assignment are skipped
                            if cds_iter[2] in gene_to_og:
                                left_side_genes_and_ogs.add(
                                    tuple(
                                        [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                    )
                                )
                            limit_reached = True
                        else:
                            # Neighbours without an OG assignment are skipped
                            if cds_iter[2] in gene_to_og:
                                left_side_genes_and_ogs.add(
                                    tuple(
                                        [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                    )
                                )
                    except IndexError:
                        # Walked past the last feature on the scaffold
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
                            if cds_iter[2] in gene_to_og:
                                right_side_genes_and_ogs.add(
                                    tuple(
                                        [cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]
                                    )
                                )
                            limit_reached = True
                        else:
                            if cds_iter[2] in gene_to_og:
                                right_side_genes_and_ogs.add(
                                    tuple([cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]])
                                )
                    except IndexError:
                        # Walked past the last feature on the scaffold
                        near_scaffold_edge = True
                        limit_reached = True
                    cds_iter_index += 1

                upstream_genes = []
                downstream_genes = []
                upstream_ogs = []
                downstream_ogs = []
                # cds = [start, end, name, score, strand]; BED strands are '+'/'-'
                if cds[4] == '+':
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
    bed_file = os.path.join(outdir, f"BEDs/{sample_name}.bed")
    proteome_file = os.path.join(outdir, f"Proteomes/{sample_name}.faa")
    genome_file = os.path.join(outdir, f"Genomes/{sample_name}.fna")

    # Initialize counters and mappings
    locus_tag_counter = 0
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

    # Protein lengths, so the BED applies the same min_length filter as the proteome
    protein_lengths = {}
    if faa_file:
        with open(faa_file, 'r') as faa_handle:
            for record in SeqIO.parse(faa_handle, 'fasta'):
                protein_lengths[record.id] = len(str(record.seq))

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

                if old_locus_tag in protein_lengths and protein_lengths[old_locus_tag] < min_length:
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

    sys.stdout.write(f"Processed {protein_count} proteins from Prokka directory.\n")


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
    bed_file = os.path.join(outdir, f"BEDs/{sample_name}.bed")
    proteome_file = os.path.join(outdir, f"Proteomes/{sample_name}.faa")
    genome_file = os.path.join(outdir, f"Genomes/{sample_name}.fna")

    # Initialize counters and mappings
    locus_tag_counter = 0
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

    # Protein lengths, so the BED applies the same min_length filter as the proteome
    protein_lengths = {}
    if faa_file:
        with open(faa_file, 'r') as faa_handle:
            for record in SeqIO.parse(faa_handle, 'fasta'):
                protein_lengths[record.id] = len(str(record.seq))

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

                if old_locus_tag in protein_lengths and protein_lengths[old_locus_tag] < min_length:
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

    sys.stdout.write(f"Processed {protein_count} proteins from Bakta directory.\n")


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
        args: Tuple containing (annotation_dir, outdir, sample_locus_tag, min_length)
    
    Returns:
        tuple: (sample_name, success, error_message, annotation_type)
    """
    annotation_dir, outdir, sample_locus_tag, min_length = args
    
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
                min_length=min_length
            )
        elif annotation_type == 'bakta':
            process_bakta_directory(
                bakta_dir=annotation_dir,
                outdir=outdir,
                sample_name=sample_name,
                locus_tag=sample_locus_tag,
                min_length=min_length
            )
        
        # Verify output files
        faa_file = os.path.join(outdir, f"Proteomes/{sample_name}.faa")
        bed_file = os.path.join(outdir, f"BEDs/{sample_name}.bed")
        fna_file = os.path.join(outdir, f"Genomes/{sample_name}.fna")
        
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
    min_length: int = config.DEFAULT_MIN_LENGTH,
    locus_tag_offset: int = 0,
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
        min_length: Minimum length of protein to include
        locus_tag_offset: Number of locus tag prefixes to skip when renaming (default: 0);
                          prefixes used by other samples' BEDs already in outdir/BEDs
                          are always skipped
    Returns:
        Dictionary containing sample mappings for wgs, proteomes, and beds
    """
    # Create locus tag options if renaming is requested
    possible_locustags = None
    if rename_locus_tags:
        possible_locustags = _get_available_locus_tags(
            locus_tag_length, outdir,
            [os.path.basename(d.rstrip('/')) for d in annotation_dirs], locus_tag_offset
        )
    
    # Prepare arguments for multiprocessing
    process_args = []
    for i, annotation_dir in enumerate(annotation_dirs):
        sample_locus_tag = possible_locustags[i] if rename_locus_tags else None
        process_args.append((annotation_dir, outdir, sample_locus_tag, min_length))
    
    # Process annotation directories in parallel
    sample_wgs = {}
    sample_proteomes = {}
    sample_beds = {}
    
    if log_object:
        log_object.info(f"Processing {len(annotation_dirs)} annotation directories...")
    
    # Use multiprocessing if multiple directories
    if len(annotation_dirs) > 1 and threads > 1:
        with create_pool(threads) as pool:
            results = list(
                _iter_progress(
                    pool.imap(_process_single_annotation_dir, process_args),
                    total=len(process_args),
                    description="Processing annotation directories",
                )
            )
    else:
        results = []
        for args in process_args:
            results.append(_process_single_annotation_dir(args))
    
    # Collect results
    for sample_name, success, error_message, annotation_type in results:
        if success:
            # Define expected file paths
            faa_file = os.path.join(outdir, f"Proteomes/{sample_name}.faa")
            bed_file = os.path.join(outdir, f"BEDs/{sample_name}.bed")
            fna_file = os.path.join(outdir, f"Genomes/{sample_name}.fna")
            
            # Verify files exist and have content
            if (os.path.isfile(faa_file) and os.path.isfile(bed_file) and 
                os.path.getsize(faa_file) > 0 and os.path.getsize(bed_file) > 0):
                
                sample_proteomes[sample_name] = f"Genome_Processing/Proteomes/{sample_name}.faa"
                sample_beds[sample_name] = f"Genome_Processing/BEDs/{sample_name}.bed"
                
                # Add genome file if it exists
                if os.path.isfile(fna_file) and os.path.getsize(fna_file) > 0:
                    sample_wgs[sample_name] = f"Genome_Processing/Genomes/{sample_name}.fna"
                
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
    Parse CDS locations (simple, fuzzy, and compound join/order) using Biopython.

    Coordinates are taken from integer positions (int() works for exact and
    fuzzy positions such as <1 or >500), converted to 1-based inclusive.
    Compound locations are represented by their overall span (min start to
    max end). Origin-spanning compound locations on circular records (parts
    wrapping from the record end back to position 1) are represented in
    unwrapped circular coordinates: start is on the high side of the origin
    and end = record length + bases after the origin, so start <= end holds
    and end - start + 1 equals the CDS length (end may exceed the record
    length); such cases are logged.

    Args:
        feature: Biopython feature object
        record: Biopython record object
        outdir: Output directory for logging
        sample_name: Sample name for logging
        scaffold: Scaffold name for logging

    Returns:
        tuple: (start, end, direction, success) or (None, None, None, False) if failed
    """
    location = feature.location
    location_str = str(location)
    try:
        direction = '+' if location.strand == 1 else '-'
        parts = list(getattr(location, 'parts', [location]))
        record_length = len(record.seq)
        touches_both_ends = (
            any(int(x.start) == 0 for x in parts)
            and any(int(x.end) == record_length for x in parts)
        )
        if len(parts) > 1 and location.strand in (1, -1) and touches_both_ends:
            # Parts are in biological order; a coordinate reversal between
            # consecutive parts (with parts touching both record ends) means
            # the feature wraps around the origin.
            wrap_index = None
            for i in range(1, len(parts)):
                prev_start, cur_start = int(parts[i - 1].start), int(parts[i].start)
                if (location.strand == 1 and cur_start < prev_start) or (
                    location.strand == -1 and cur_start > prev_start
                ):
                    wrap_index = i
                    break
            if wrap_index is not None:
                first_side = parts[:wrap_index]
                second_side = parts[wrap_index:]
                # Biological order: on + the CDS starts before the origin (high
                # coordinates), on - it starts after it (low coordinates).
                high_side = first_side if location.strand == 1 else second_side
                low_side = second_side if location.strand == 1 else first_side
                # Represent in unwrapped circular coordinates: start on the high
                # side, end extended past the record end by the bases after the
                # origin. BED length then equals the full CDS length and
                # start <= end; end may exceed the record length.
                start = min(int(x.start) for x in high_side) + 1
                end = record_length + max(int(x.end) for x in low_side)
                _log_cds_issue(
                    outdir, sample_name, scaffold, f"origin_spanning_{location_str[:50]}",
                    "origin_spanning_location",
                    f"Origin-spanning location {location_str} represented in unwrapped "
                    f"coordinates {start}-{end} (record length {record_length})",
                )
                return start, end, direction, True
        start = min(int(x.start) for x in parts) + 1  # Convert to 1-based
        end = max(int(x.end) for x in parts)
        if start > end:
            raise ValueError(f"start {start} > end {end}")
        return start, end, direction, True
    except Exception as e:
        feature_info = f"location_parse_exception_{location_str[:50]}"
        _log_cds_issue(outdir, sample_name, scaffold, feature_info, "location_parse_exception", f"Exception parsing location: {location_str} - {str(e)}")
        return None, None, None, False

