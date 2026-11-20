#!/usr/bin/env python3
"""
Main entry point for the unified bofasa interface.

This module provides the main function that handles the unified command-line
interface with 'run' and 'prep' subcommands.
"""

import sys
import argparse
import io
from typing import Any, Optional, TextIO
from . import config
from .utils import print_bofasa_help

try:
    from rich.console import Console
    RICH_AVAILABLE: bool = True
    console: Optional[Console] = Console(force_terminal=True)
except ImportError:
    RICH_AVAILABLE: bool = False
    console: Optional[Any] = None


class BofasaArgumentParser(argparse.ArgumentParser):
    """Custom argument parser that uses our colored help functions."""

    def __init__(
        self, *args: Any, help_type: Optional[str] = None, **kwargs: Any
    ) -> None:
        self.help_type = help_type
        super().__init__(*args, **kwargs)

    def print_help(self, file: Optional[TextIO] = None) -> None:
        """Override print_help to use our colored help functions."""
        help_type = self.help_type
        
        # If this is the main parser (no help_type), use the main help function
        if help_type is None:
            from .utils import print_bofasa_help
            print_bofasa_help()
            return
        
        # For subcommand parsers, use the colored help function
        # Capture the actual argparse help output
        help_buffer = io.StringIO()
        super().print_help(file=help_buffer)
        help_text = help_buffer.getvalue()
        help_buffer.close()
        
        # Color the actual help text
        from .utils import print_colored_help
        
        print_colored_help(help_text, help_type)

    def add_subparsers(self, **kwargs: Any) -> argparse._SubParsersAction:
        """Override add_subparsers to use our custom parser class."""
        subparsers = super().add_subparsers(**kwargs)

        # Override the add_parser method to use our custom parser
        original_add_parser = subparsers.add_parser

        def add_parser_with_custom_help(
            name: str, **kwargs: Any
        ) -> argparse.ArgumentParser:
            parser = original_add_parser(name, **kwargs)
            # Set the help_type based on the subcommand name
            parser.help_type = name
            return parser

        subparsers.add_parser = add_parser_with_custom_help
        return subparsers


class BofasaHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom formatter that preserves newlines in help text."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)

    def _format_action_invocation(self, action: argparse.Action) -> str:
        if not action.option_strings:
            (metavar,) = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            # if the Optional doesn't take a value, format is: -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)
            # if the Optional takes a value, format is: -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append('%s %s' % (option_string, args_string))
            return ', '.join(parts)

    def _split_lines(self, text: Optional[str], width: int) -> list[str]:
        """Split text into lines, preserving newlines."""
        if text is None:
            return []
        lines = []
        for line in text.split('\n'):
            if line.strip():
                # Split long lines at word boundaries
                words = line.split()
                current_line = []
                current_length = 0

                for word in words:
                    word_length = len(word)
                    if current_length + word_length + 1 <= width:
                        current_line.append(word)
                        current_length += word_length + 1
                    else:
                        if current_line:
                            lines.append(' '.join(current_line))
                        current_line = [word]
                        current_length = word_length + 1

                if current_line:
                    lines.append(' '.join(current_line))
            else:
                lines.append('')
        return lines


class RawTextHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Custom formatter that preserves newlines in argument help text."""

    def _format_action_invocation(self, action: argparse.Action) -> str:
        if not action.option_strings:
            (metavar,) = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            # if the Optional doesn't take a value, format is: -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)
            # if the Optional takes a value, format is: -s ARGS, --long ARGS
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append('%s %s' % (option_string, args_string))
            return ', '.join(parts)

    def _split_lines(self, text: Optional[str], width: int) -> list[str]:
        """Split text into lines, preserving newlines."""
        if text is None:
            return []
        lines = []
        for line in text.split('\n'):
            if line.strip():
                # Split long lines at word boundaries
                words = line.split()
                current_line = []
                current_length = 0

                for word in words:
                    word_length = len(word)
                    if current_length + word_length + 1 <= width:
                        current_line.append(word)
                        current_length += word_length + 1
                    else:
                        if current_line:
                            lines.append(' '.join(current_line))
                        current_line = [word]
                        current_length = word_length + 1

                if current_line:
                    lines.append(' '.join(current_line))
            else:
                lines.append('')
        return lines


def create_main_parser() -> BofasaArgumentParser:
    """Create the main argument parser for the bofasa interface."""
    parser = BofasaArgumentParser(
        description="bacterial orthology finding and syntenic analysis (bofasa)",
        formatter_class=BofasaHelpFormatter,
    )

    # Add version option
    parser.add_argument(
        "-v", "--version", action="store_true", help="Get version and exit."
    )

    # Add subparsers
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Add setup_annotation_dbs subcommand
    setup_parser = subparsers.add_parser(
        "setup", 
        help="Setup annotation databases required for bofasa analysis",
        description="Setup annotation databases required for bofasa analysis.",
        formatter_class=BofasaHelpFormatter,
    )
    add_setup_annotation_dbs_arguments(setup_parser)
    setup_parser.epilog = (
        "Examples:\n"
        "  export BOFASA_DB_PATH=/path/to/databases\n"
        "  bofasa setup\n"
        "  bofasa setup --threads 8 --force\n"
        "  bofasa setup --skip-genomad\n"
        "  bofasa setup --skip-genomad --pfam-version 37.4"
    )



    # Add prep subcommand
    prep_parser = subparsers.add_parser(
        "prep",
        help="Prepare input data for bofasa analysis",
        description="Prepare input data for bofasa run.\n\n"
        "It can take in contig/scaffold FASTA files (*.fasta, *.fa, *.fna),\n"
        "GenBank files (*.genbank, *.gbk, *.gbff), with CDS features and\n"
        "translation qualifiers available, or Prokka/Bakta annotation\n"
        "folders. Note, gzipped files are not supported.\n\n"
        "Pfam HMMs will be used to split proteins into domain/inter-domain\n"
        "units which will serve as the basis for orthologroup inference\n"
        "downstream.\n\n"
        "Annotations for transposases and MGEs will also be performed using\n"
        "sequences from the ISFinder database and geNomad.\n\n"
        "If Bakta or Prokka annotations are provided, please make sure\n"
        "such annotations were run using the \"--prefix\" argument as this\n"
        "will be important for parsing out sample names.\n\n"
        "If using --run-genomad, please also cite geNomad:\n"
        "  Camargo AP, Roux S, Schulz F, Babinski M, Xu Y, Hu B, Chain PSG,\n"
        "  Nayfach S, Kyrpides NC. Identification of mobile genetic elements\n"
        "  with geNomad. Nature Biotechnology 42, 1303-1312 (2024).",
        formatter_class=BofasaHelpFormatter,
    )
    add_prep_arguments(prep_parser)
    prep_parser.epilog = (
        "Examples:\n"
        "  bofasa prep -i genome1.fasta genome2.fasta -o prepared_data\n"
        "  bofasa prep -i *.fasta -o prepared_data -c 8 --gene-calling-method prodigal\n"
        "  bofasa prep -i *.fasta -o prepared_data --run-genomad --extract-mge-genomes"
    )

    # Add run subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Execute the main bofasa orthology analysis",
        description="Execute the main bofasa orthology analysis.\n\n"
        "OrthoFinder is used to determine coarse domain ortholog groups,\n"
        "please also cite it:\n"
        "  Emms DM, Kelly S. OrthoFinder: phylogenetic orthology inference\n"
        "  for comparative genomics. Genome Biology 20, 238 (2019).",
        formatter_class=BofasaHelpFormatter,
    )
    add_run_arguments(run_parser)
    run_parser.epilog = (
        "Examples:\n"
        "  bofasa run -i prepared_data -o analysis_results\n"
        "  bofasa run -i prepared_data -o analysis_results --threads 8 --og-consensus --core-genome\n"
        "  bofasa run -i prepared_data -o analysis_results --surrounding-bp 5000 --dog-jaccard 0.3"
    )

    # Add domscape subcommand
    domscape_parser = subparsers.add_parser(
        "domscape",
        help="Count proteins carrying Pfam domains across genomes",
        description="Annotate Pfam domains across a set of genomes and produce a matrix\n"
        "of the number of proteins carrying each Pfam domain (rows) in each\n"
        "genome (columns).\n\n"
        "It can take in contig/scaffold FASTA files (*.fasta, *.fa, *.fna),\n"
        "GenBank files (*.genbank, *.gbk, *.gbff) with CDS features and\n"
        "translation qualifiers available, Prokka/Bakta annotation folders,\n"
        "or protein FASTA files (proteomes). Note, gzipped files are not\n"
        "supported.\n\n"
        "A protein is counted once for a Pfam domain if at least one hit to\n"
        "the domain meets the Pfam gathering threshold - multiple copies of\n"
        "a domain within a protein are not counted separately, and hits to\n"
        "different domains are counted regardless of coordinate overlap.\n\n"
        "There is no upper limit on the number of genomes. Runs can be\n"
        "resumed by re-running the same command with the same output\n"
        "directory.",
        formatter_class=BofasaHelpFormatter,
    )
    add_domscape_arguments(domscape_parser)
    domscape_parser.epilog = (
        "Examples:\n"
        "  bofasa domscape -i genomes/*.fna -o domscape_results -c 16\n"
        "  bofasa domscape -p proteomes/*.faa -o domscape_results -c 16\n"
        "  bofasa domscape -a prokka_dirs/* -i genbanks/*.gbk -o domscape_results"
    )

    return parser


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the run subcommand."""
    parser.add_argument(
        "-i",
        "--bofasa-prep-dir",
        required=True,
        help="Input - a directory produced by bofasa prep.",
    )
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory.")
    parser.add_argument(
        "-sr",
        "--surrounding-bp",
        type=int,
        default=config.DEFAULT_SURROUNDING_BP,
        help=f"The number of basepairs to look up and downstream for syntenic\n"
        f"analysis [Default is {config.DEFAULT_SURROUNDING_BP}].",
    )
    parser.add_argument(
        "-ogc",
        "--og-consensus",
        action="store_true",
        help="Determine consensus sequences of ortholog groups for annotation of\n"
        "new genomes.",
    )
    parser.add_argument(
        "-cg",
        "--core-genome",
        action="store_true",
        help="Construct core genome alignment(s) using (near) single-copy-core\n"
        "resolved domain ortholog groups for phylogenomics.",
    )
    parser.add_argument(
        "-dj",
        "--dog-jaccard",
        type=float,
        default=config.DEFAULT_DOG_JACCARD,
        help="The Jaccard index threshold for domain ortholog group overlap between\n"
        "two pairs of proteins needed to consider them as sharing an edge. Will\n"
        "be automatically lowered for cases where single-copy domain ortholog\n"
        f"groups are observed [Default is {config.DEFAULT_DOG_JACCARD}].",
    )
    parser.add_argument(
        "-fic",
        "--fixation-index-cutoff",
        type=float,
        default=config.DEFAULT_FIXATION_INDEX_CUTOFF,
        help="Fixation index cutoff for domain ortholog group re-merging following\n"
        f"phylogenetic splitting [Default is {config.DEFAULT_FIXATION_INDEX_CUTOFF}].",
    )
    parser.add_argument(
        "-smb",
        "--skip-merge-back",
        action="store_true",
        help="Skip fixation index-based assessment for deciding whether to perform\n"
        "re-joining of split ortholog partitions.",
    )
    parser.add_argument(
        "-spr",
        "--skip-phylo-refine",
        action="store_true",
        help="Skip phylogenetic refinement of domain ortholog groups.",
    )
    parser.add_argument(
        "-rs",
        "--rooting-seeds",
        type=int,
        default=config.DEFAULT_ROOTING_SEEDS,
        help="The maximum number of nodes to try for rooting. A random sampling is\n"
        f"performed [Default is {config.DEFAULT_ROOTING_SEEDS}; uses midpoint rooting].",
    )
    parser.add_argument(
        "-us",
        "--ultra-sens",
        action="store_true",
        help="Use ultra-sensitivity parameters for running OrthoFinder DIAMOND\n"
        "searches.",
    )
    parser.add_argument(
        "-mmc",
        "--mge-min-coverage",
        type=float,
        default=None,
        help="Minimum query and subject coverage (%%, 0-100) of MGE-vs-MGE DIAMOND\n"
        "hits used to cluster MGE proteins before fitting them into ortholog\n"
        "groups (only relevant if MGE genomes were extracted with prep -emg).\n"
        "Higher values reduce chaining of distantly/partially similar MGE\n"
        "proteins but also break links between gene fragments and full-length\n"
        "copies [Default is off - no coverage filter].",
    )
    parser.add_argument(
        "-mi",
        "--mcl-inflation",
        type=float,
        default=config.DEFAULT_MCL_INFLATION,
        help="MCL inflation parameter for determining coarse domain-resolution\n"
        f"ortholog groups via OrthoFinder [Default is {config.DEFAULT_MCL_INFLATION}].",
    )
    parser.add_argument(
        "-ns",
        "--near-scc-prop",
        type=float,
        default=config.DEFAULT_NEAR_SCC_PROP,
        help="Proportion of genomes which single-copy coarse domain ortholog\n"
        "groups need to be found in for constructing core genome alignment\n"
        f"for phylogenomics [Default is {config.DEFAULT_NEAR_SCC_PROP}].",
    )
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        default=config.DEFAULT_THREADS,
        help=f"Total number of cores/threads to use [Default is {config.DEFAULT_THREADS}].",
    )
    parser.add_argument(
        "-mrd",
        "--max-recursion-depth",
        type=int,
        default=config.DEFAULT_MAX_RECURSION_DEPTH,
        help=f"The maximum recursion depth [Default is {config.DEFAULT_MAX_RECURSION_DEPTH}].",
    )
    parser.add_argument(
        "-mm",
        "--max-memory",
        type=int,
        default=config.DEFAULT_MAX_MEMORY,
        help=config.MAX_MEMORY_HELP,
    )
    parser.add_argument(
        "-al",
        "--aligner",
        choices=config.ALIGNMENT_METHODS,
        default=config.DEFAULT_ALIGNER,
        help="Multiple sequence aligner to use: 'famsa' or 'muscle'\n"
        f"(MUSCLE super5) [Default is {config.DEFAULT_ALIGNER}]. Warning: MUSCLE "
        "super5 alignments are not deterministic, so runs using muscle are not "
        "exactly reproducible; FAMSA is deterministic.",
    )
    parser.add_argument(
        "-y",
        "--auto",
        action="store_true",
        help="Automatically answer 'yes' to all interactive prompts. Useful for\n"
        "non-interactive/batch execution.",
    )


def add_prep_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the prep subcommand."""
    parser.add_argument(
        "-i",
        "--input-genomes",
        nargs="+",
        required=False,
        default=[],
        help="Input genomes (FASTA or GenBank files). Required if no annotation directories are provided.",
    )
    parser.add_argument(
        "-a",
        "--annotation-dirs",
        nargs="+",
        help="Annotation directories (Prokka or Bakta output directories). Required if no input genomes are provided.",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory for prepared data."
    )
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        default=config.DEFAULT_THREADS,
        help=f"Number of threads to use [Default is {config.DEFAULT_THREADS}].",
    )
    parser.add_argument(
        "-gcm",
        "--gene-calling-method",
        choices=config.GENE_CALLING_METHODS,
        default="pyrodigal",
        help="Gene calling method [Default is pyrodigal].",
    )
    parser.add_argument(
        "-l",
        "--locus-tag-length",
        type=int,
        default=config.DEFAULT_LOCUS_TAG_LENGTH,
        help=f"Length of locus tags [Default is {config.DEFAULT_LOCUS_TAG_LENGTH}].",
    )
    parser.add_argument(
        "-m", "--meta-mode", action="store_true", help="Use meta mode for gene calling."
    )
    parser.add_argument(
        "-rlt",
        "--rename-locus-tags",
        action="store_true",
        help="Rename locus tags in GenBank files and annotation directories (Prokka/Bakta).",
    )
    parser.add_argument(
        "-mm",
        "--max-memory",
        type=int,
        default=config.DEFAULT_MAX_MEMORY,
        help=config.MAX_MEMORY_HELP,
    )
    parser.add_argument(
        "-ml",
        "--min-length",
        type=int,
        default=config.DEFAULT_MIN_LENGTH,
        help=f"Minimum length of domain/inter-domain unit to consider [Default is {config.DEFAULT_MIN_LENGTH}].",
    )
    parser.add_argument(
        "-sds",
        "--skip-domain-splitting",
        action="store_true",
        help="Skip domain-annotation based splitting of coding sequences.",
    )
    parser.add_argument(
        "-ntdo",
        "--no-trim-domain-overlaps",
        action="store_true",
        help="Do not trim the minor (<10%%) overlap between an accepted Pfam domain "
        "and higher-scoring domains already accepted for the same protein. By "
        "default overlaps are trimmed so each accepted domain yields exactly one, "
        "correctly labelled, chunk; this flag restores the behavior of v1.2.2 and "
        "earlier.",
    )
    parser.add_argument(
        "-rg",
        "--run-genomad",
        action="store_true",
        help="Whether to run geNomad - takes time but results can be used to achieve higher quality orthogroup predictions.",
    )
    parser.add_argument(
        "-gs",
        "--genomad-splits",
        type=int,
        default=config.DEFAULT_GENOMAD_SPLITS,
        help="Number of chunks geNomad splits its database into (geNomad --splits) - more "
        "splits lower memory usage at the cost of speed; only used with --run-genomad "
        f"[Default is {config.DEFAULT_GENOMAD_SPLITS}].",
    )
    parser.add_argument(
        "-x",
        "--ignore-upperbound-limit",
        action="store_true",
        help=f"Ignore the upper bound limit for the number of genomes to process [Default is {config.DEFAULT_MAX_GENOMES}].",
    )
    parser.add_argument(
        "-emg",
        "--extract-mge-genomes",
        action="store_true",
        help="Extract plasmids and phages identified by geNomad as separate genome files for downstream analysis.",
    )
    parser.add_argument(
        "-y",
        "--auto",
        action="store_true",
        help="Automatically answer 'yes' to all interactive prompts. Useful for\n"
        "non-interactive/batch execution.",
    )


def add_domscape_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the domscape subcommand."""
    parser.add_argument(
        "-i",
        "--input-genomes",
        nargs="+",
        default=[],
        help="Input genomes (FASTA or GenBank files).",
    )
    parser.add_argument(
        "-p",
        "--input-proteomes",
        nargs="+",
        default=[],
        help="Input proteomes (protein FASTA files). Sample names are taken from\n"
        "the file names without their extension.",
    )
    parser.add_argument(
        "-a",
        "--annotation-dirs",
        nargs="+",
        default=[],
        help="Annotation directories (Prokka or Bakta output directories).",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory."
    )
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        default=config.DEFAULT_THREADS,
        help=f"Number of threads to use [Default is {config.DEFAULT_THREADS}].",
    )
    parser.add_argument(
        "-gcm",
        "--gene-calling-method",
        choices=config.GENE_CALLING_METHODS,
        default="pyrodigal",
        help="Gene calling method for genomes provided as FASTA files\n"
        "[Default is pyrodigal].",
    )
    parser.add_argument(
        "-m", "--meta-mode", action="store_true", help="Use meta mode for gene calling."
    )
    parser.add_argument(
        "-ml",
        "--min-length",
        type=int,
        default=config.DEFAULT_MIN_LENGTH,
        help="Minimum length of proteins (in amino acids) to keep from genomes or\n"
        "annotation directories - not applied to input proteomes\n"
        f"[Default is {config.DEFAULT_MIN_LENGTH}].",
    )
    parser.add_argument(
        "-bs",
        "--batch-size",
        type=int,
        default=config.DEFAULT_DOMSCAPE_BATCH_SIZE,
        help="Number of proteomes to search against Pfam together in a single\n"
        "batch. Larger batches are more efficient but use more memory\n"
        f"[Default is {config.DEFAULT_DOMSCAPE_BATCH_SIZE}].",
    )
    parser.add_argument(
        "-mm",
        "--max-memory",
        type=int,
        default=config.DEFAULT_MAX_MEMORY,
        help=config.MAX_MEMORY_HELP,
    )
    parser.add_argument(
        "-y",
        "--auto",
        action="store_true",
        help="Automatically answer 'yes' to all interactive prompts. Useful for\n"
        "non-interactive/batch execution.",
    )


def add_setup_annotation_dbs_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the setup-dbs subcommand."""
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        help=f"Number of threads to use for parallel processing.\n"
        f"Default is {config.DEFAULT_THREADS} threads.",
        required=False,
        default=config.DEFAULT_THREADS,
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Force overwrite existing databases if they already exist.\n"
        "Use with caution as this will delete existing data.",
        required=False,
        default=False,
    )
    parser.add_argument(
        "-sg",
        "--skip-genomad",
        action="store_true",
        help="Skip downloading the geNomad database. The geNomad database is\n"
        "downloaded by default, but it is large and only needed for\n"
        "'bofasa prep --run-genomad'.",
        required=False,
        default=False,
    )
    parser.add_argument(
        "-pv",
        "--pfam-version",
        choices=["current"] + sorted(config.PFAM_RELEASES),
        help="Pfam release to download. Default is the current release. Use\n"
        "37.4 to reproduce the expected results of the test case.",
        required=False,
        default="current",
    )
    parser.add_argument(
        "-y",
        "--auto",
        action="store_true",
        help="Automatically answer 'yes' to all interactive prompts. Useful for\n"
        "non-interactive/batch execution.",
    )


def run_bofasa_analysis(args: argparse.Namespace) -> None:
    """Run the main bofasa analysis."""
    # Import here to avoid circular imports and heavy dependencies
    import os
    import shutil
    import subprocess
    from .utils import (
        create_logger_object, close_logger_object, memory_limit, 
        check_orthofinder_setup, setup_ready_directory, log_parameters_to_file
    )
    from .core import (
        determine_protein_orthogroups, resolve_orthogroups_using_phylogenetics,
        combine_orthofinder_results
    )
    from .analysis import (
        determine_ortholog_group_contexts,
        create_profile_hmms_and_consensus_seqs,
        create_final_report,
        create_final_visual,
    )
    from .alignment import create_protein_alignments, create_near_scc_resolved_domain_protein_alignments

    if args.mge_min_coverage is not None and not (0.0 <= args.mge_min_coverage <= 100.0):
        sys.stderr.write("Error: --mge-min-coverage must be a percentage between 0 and 100.\n")
        sys.exit(1)

    # Check OrthoFinder setup early in the workflow
    sys.stdout.write("Checking OrthoFinder setup...\n")
    orthofinder_ok, orthofinder_message = check_orthofinder_setup()
    if not orthofinder_ok:
        sys.stdout.write(orthofinder_message + '\n')
        sys.exit(1)
    sys.stdout.write("✓ OrthoFinder is properly set up\n")

    # Check if output directory already exists
    if os.path.exists(args.output_dir) and os.listdir(args.output_dir):
        sys.stdout.write(f"\nWARNING: The output directory '{args.output_dir}' already exists and contains files!\n")
        sys.stdout.write("This process may overwrite existing files.\n")

        if args.auto:
            sys.stdout.write("--auto specified: continuing with existing output directory...\n")
        else:
            while True:
                response = input("Do you want to continue and potentially overwrite existing files? (y/N): ").strip().lower()
                if response in ['y', 'yes']:
                    sys.stdout.write("Continuing with existing output directory...\n")
                    break
                elif response in ['n', 'no', '']:
                    sys.stdout.write("Exiting. Please choose a different output directory or remove existing files.\n")
                    sys.exit(1)
                else:
                    sys.stdout.write("Please enter 'y' for yes or 'n' for no.\n")

    # Debug: print(args.threads)  # Removed debug print

    # Set recursion limit
    max_recursion_depth = args.max_recursion_depth
    sys.setrecursionlimit(max_recursion_depth)

    # Set memory limit
    if args.max_memory:
        memory_limit(args.max_memory)

    # Create output directory first
    setup_ready_directory([args.output_dir], overwrite_mode="skip")

    # Create logger
    log_file = os.path.join(args.output_dir, "bofasa.log")
    logger = create_logger_object(log_file)
    
    # Check if logger was created successfully
    if logger is None:
        sys.stdout.write("Warning: Failed to create logger. Continuing without logging.\n")
        logger = None

    try:
        if logger:
            logger.info("Starting bofasa analysis")
            logger.info(f"Input directory: {args.bofasa_prep_dir}")
            logger.info(f"Output directory: {args.output_dir}")
            logger.info(f"Threads: {args.threads}")
            logger.info(f"Memory limit requested: {args.max_memory}GB" if args.max_memory else "Memory limit requested: none")
        else:
            sys.stdout.write("Starting bofasa analysis\n")
            sys.stdout.write(f"Input directory: {args.bofasa_prep_dir}\n")
            sys.stdout.write(f"Output directory: {args.output_dir}\n")
            sys.stdout.write(f"Threads: {args.threads}\n")
            sys.stdout.write(f"Memory limit requested: {args.max_memory}GB\n" if args.max_memory else "Memory limit requested: none\n")

        # Log parameters
        if logger:
            logger.info("Saving parameters for future records.")
            parameters_file = os.path.join(args.output_dir, "Parameter_Inputs.txt")
            parameter_values = [
                args.bofasa_prep_dir, args.output_dir, args.surrounding_bp,
                args.mcl_inflation, args.ultra_sens,
                args.fixation_index_cutoff, args.dog_jaccard,
                args.skip_phylo_refine, args.rooting_seeds,
                args.skip_merge_back, args.og_consensus, args.core_genome,
                args.near_scc_prop, args.threads, args.max_memory,
                args.mge_min_coverage, args.aligner
            ]
            parameter_names = [
                "bofasa_prep directory",
                "Output directory",
                "Surrounding BP for syntenic conservation assessment",
                "MCL inflation",
                "Use ultra-sensitivity mode for DIAMOND searching in OrthoFinder?",
                (
                    "Fixation index cutoff for domain ortholog group re-merging "
                    "following phylogenetic splitting"
                ),
                (
                    "Jaccard index cutoff for domain ortholog groups shared between "
                    "full protein pairs"
                ),
                "Skip phylogenetic refinement of domain ortholog groups?",
                "Maximum number of nodes to try for rooting",
                "Skip fixation index-based re-mergining of split ortholog partitions?",
                "Determine consensus sequences for ortholog groups?",
                "Create core genome alignment(s) for phylogenomics?",
                "Near-SCC conservation proportion?",
                "Number of threads/cores",
                "Maximum memory in GB",
                "Minimum coverage (%) for MGE-vs-MGE clustering hits",
                "Multiple sequence aligner",
            ]
            log_parameters_to_file(parameters_file, parameter_names, parameter_values)
            logger.info("Done saving parameters!")

        # Validate that the requested aligner is available
        import shutil
        if shutil.which(args.aligner) is None:
            msg = (f"Error: the requested aligner '{args.aligner}' was not found in PATH. "
                   f"Please install it or choose another via --aligner "
                   f"({', '.join(config.ALIGNMENT_METHODS)}).")
            if logger:
                logger.error(msg)
            sys.stderr.write(msg + '\n')
            sys.exit(1)
        if logger:
            logger.info(f"Multiple sequence aligner: {args.aligner}")
        else:
            sys.stdout.write(f"Multiple sequence aligner: {args.aligner}\n")

        # Validate input directory
        index_file = os.path.join(args.bofasa_prep_dir, "Info_on_Input_Genome_Files.txt")
        isfinder_file = os.path.join(args.bofasa_prep_dir, "Sample_IS_Element_Proteins.txt")
        phage_file = os.path.join(args.bofasa_prep_dir, "Sample_Phage_Proteins.txt")
        plasmid_file = os.path.join(args.bofasa_prep_dir, "Sample_Plasmid_Proteins.txt")
        
        try:
            assert os.path.isdir(args.bofasa_prep_dir)
            assert os.path.isfile(index_file) and os.path.isfile(isfinder_file)
            msg = 'geNomad was not run ...'
            if os.path.isfile(phage_file) and os.path.isfile(plasmid_file):
                msg = 'geNomad was successfully run on input genomes and will be used downstream ...'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
        except Exception as e:
            if logger:
                logger.error('Issue with validating directory with input genomes exists.')
            else:
                sys.stdout.write('Issue with validating directory with input genomes exists.\n')
            sys.exit(1)

        # Create checkpoint and final results directories
        checkdir = os.path.join(args.output_dir, "Checkpoint_Files/")
        findir = os.path.join(args.output_dir, "Final_Results/")
        setup_ready_directory([checkdir, findir], overwrite_mode="skip")

        # Step 1: Run OrthoFinder (on bacterial genomes only, excluding MGEs)
        step1_checkpoint_file = os.path.join(checkdir, "Step1.txt")
        msg = '\n--------------------\nStep 1\n--------------------\nRunning OrthoFinder on bacterial genomes (MGEs will be processed separately if bofasa prep was run with the -emg flag).'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        # Load MGE metadata to determine which samples are MGEs
        mge_metadata_file = os.path.join(args.bofasa_prep_dir, "Sample_MGE_Metadata.txt")
        sample_is_mge = set()
        if os.path.exists(mge_metadata_file):
            with open(mge_metadata_file) as mge_handle:
                for i, line in enumerate(mge_handle):
                    if i == 0:  # Skip header
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) == 2 and parts[1] == "True":
                        sample_is_mge.add(parts[0])
        
        # All files are in a common directory now
        all_fastas_dir = os.path.join(args.bofasa_prep_dir, "Domain_and_Interdomain_FASTAs/")
        orthofinder_results_dir = os.path.join(args.output_dir, "OrthoFinder_Results/")
        
        # Create temporary directory with copies of only bacterial genome files for OrthoFinder
        orthofinder_input_dir = os.path.join(args.output_dir, "OrthoFinder_Input_Bacterial_Genomes/")
        
        if not os.path.isfile(step1_checkpoint_file):
            # Create temporary directory with copies of bacterial genome files only
            setup_ready_directory([orthofinder_input_dir], overwrite_mode="overwrite")
            
            # Copy bacterial genomes (excluding MGEs)
            for fasta_file in sorted(os.listdir(all_fastas_dir)):
                if fasta_file.endswith('.ccds.faa'):
                    sample_name = fasta_file.replace('.ccds.faa', '')
                    # Only include if not an MGE
                    if sample_name not in sample_is_mge:
                        src = os.path.join(all_fastas_dir, fasta_file)
                        dst = os.path.join(orthofinder_input_dir, fasta_file)
                        shutil.copy2(src, dst)
            
            num_bacterial = len(os.listdir(orthofinder_input_dir))
            num_mge = len(sample_is_mge)
            msg = f"Running OrthoFinder on {num_bacterial} bacterial genomes ({num_mge} MGEs will be integrated separately)"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
            
            orthofinder_cmd = [
                'orthofinder', '-f', orthofinder_input_dir, '-o', orthofinder_results_dir,
                '-t', str(args.threads), '-og', '-I', str(args.mcl_inflation)
            ]
            if args.ultra_sens:
                orthofinder_cmd += ['-S', 'diamond_ultra_sens']

            # OrthoFinder refuses to write into an existing non-default output directory,
            # so remove leftovers from a previously interrupted run.
            if os.path.exists(orthofinder_results_dir):
                msg = f"Removing OrthoFinder results directory from a previous incomplete run: {orthofinder_results_dir}"
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')
                shutil.rmtree(orthofinder_results_dir)
            
            try:
                subprocess.run(orthofinder_cmd, check=True)
                if logger:
                    logger.info(f'Successfully ran: {" ".join(orthofinder_cmd)}')
            except subprocess.CalledProcessError as e:
                if logger:
                    logger.error(f'Had an issue running: {" ".join(orthofinder_cmd)}')
                    logger.error(str(e))
                else:
                    sys.stdout.write(f'Had an issue running: {" ".join(orthofinder_cmd)}\n')
                sys.exit(1)
            
            with open(step1_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 1 checkpoint found - skipping OrthoFinder execution"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Find OrthoFinder results files
        orthofinder_subdirs = []
        if os.path.isdir(orthofinder_results_dir):
            orthofinder_subdirs = [d for d in sorted(os.listdir(orthofinder_results_dir))
                                  if d.startswith('Results_') and os.path.isdir(os.path.join(orthofinder_results_dir, d))]
        if not orthofinder_subdirs:
            if logger:
                logger.error('Could not find OrthoFinder results directory.')
            else:
                sys.stdout.write('Could not find OrthoFinder results directory.\n')
            sys.exit(1)
        
        # If multiple Results_* directories exist, use the most recently modified one
        orthofinder_subdir = max(orthofinder_subdirs, key=lambda d: os.path.getmtime(os.path.join(orthofinder_results_dir, d)))
        if len(orthofinder_subdirs) > 1:
            msg = f'Warning: multiple OrthoFinder results directories found; using the most recent: {orthofinder_subdir}'
            if logger:
                logger.warning(msg)
            else:
                sys.stdout.write(msg + '\n')
        orthofinder_base = os.path.join(orthofinder_results_dir, orthofinder_subdir)
        
        orthofinder_tsv_file = os.path.join(orthofinder_base, "Orthogroups/Orthogroups.tsv")
        orthofinder_tsv_singletons_file = os.path.join(orthofinder_base, "Orthogroups/Orthogroups_UnassignedGenes.tsv")
        orthofinder_mod_tsv_file = os.path.join(orthofinder_base, "Orthogroups/Orthogroups_Modified.tsv")
        orthofinder_mod_tsv_singletons_file = os.path.join(orthofinder_base, "Orthogroups/Orthogroups_UnassignedGenes_Modified.tsv")
        orthofinder_graph_file = os.path.join(orthofinder_base, "WorkingDirectory/OrthoFinder_graph.txt")
        orthofinder_seqid_file = os.path.join(orthofinder_base, "WorkingDirectory/SequenceIDs.txt")
        orthofinder_fasta_dir = os.path.join(orthofinder_base, "Orthogroup_Sequences/")
        
        try:
            assert os.path.isfile(orthofinder_tsv_file) and os.path.isfile(orthofinder_tsv_singletons_file)
            assert os.path.isfile(orthofinder_graph_file) and os.path.isfile(orthofinder_seqid_file)
        except AssertionError:
            if logger:
                logger.error(
                    'Could not validate OrthoFinder results exist. '
                    'Perhaps investigate the OrthoFinder log files for further '
                    'logging information.'
                )
            else:
                sys.stdout.write(
                    'Could not validate Orthofinder results exist. '
                    'Perhaps investigate the OrthoFinder log files for further '
                    'logging information.\n'
                )
            sys.exit(1)

        # Step 1b: Integrate MGE proteins into ortholog groups
        step1b_checkpoint_file = os.path.join(checkdir, "Step1b.txt")
        mge_og_assignment_file = os.path.join(args.output_dir, "MGE_Orthogroup_Assignments.tsv")
        
        if not os.path.isfile(step1b_checkpoint_file):
            msg = '\n--------------------\nStep 1b\n--------------------\nIntegrating MGE (phage/plasmid) proteins into ortholog groups.'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
            
            # Check if there are any MGEs to integrate
            if sample_is_mge:
                # Import MGE integration function
                from bofasa.analysis import integrate_mge_proteins_into_orthogroups
                
                workspace_dir = os.path.join(args.output_dir, "MGE_Integration_Workspace/")
                setup_ready_directory([workspace_dir], overwrite_mode="overwrite")

                # Run MGE integration 
                integrate_mge_proteins_into_orthogroups(
                    input_dir=all_fastas_dir,
                    sample_is_mge=sample_is_mge,
                    orthofinder_tsv_file=orthofinder_tsv_file,
                    orthofinder_tsv_singletons_file=orthofinder_tsv_singletons_file,
                    orthofinder_mod_tsv_file=orthofinder_mod_tsv_file,
                    orthofinder_mod_tsv_singletons_file=orthofinder_mod_tsv_singletons_file,
                    workspace_dir=workspace_dir,
                    mge_og_assignment_file=mge_og_assignment_file,
                    orthogroup_sequences_dir=orthofinder_fasta_dir,
                    threads=args.threads,
                    ultra_sens=args.ultra_sens,
                    evalue_cutoff=1e-3,
                    min_coverage=args.mge_min_coverage,
                    log_object=logger
                )
                
                msg = "MGE protein integration completed"
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')
            else:
                msg = "No MGE proteins found or were requested to be extracted via the -emg flag in bofasa prep"
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')
            
            # Create Step 1b checkpoint
            with open(step1b_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 1b checkpoint found"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        if sample_is_mge:
            orthofinder_tsv_file = orthofinder_mod_tsv_file
            orthofinder_tsv_singletons_file = orthofinder_mod_tsv_singletons_file

        # Step 2: Process OrthoFinder results and determine more coarse domain ortholog groups
        step2_checkpoint_file = os.path.join(checkdir, "Step2.txt")
        msg = '\n--------------------\nStep 2\n--------------------\nDetermine refined domain-resolution ortholog groups.'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        resdog_dir = os.path.join(args.output_dir, "Resolve_Domain_Ortholog_Groups/")
        resulting_dogs_file = os.path.join(findir, "Domain_Ortholog_Groups.tsv")
        
        if not os.path.isfile(step2_checkpoint_file):
            if args.skip_phylo_refine:
                combine_orthofinder_results(
                    orthofinder_tsv_file, orthofinder_tsv_singletons_file, 
                    resulting_dogs_file, logger
                )
            else:
                setup_ready_directory([resdog_dir], overwrite_mode="overwrite")
                resolve_orthogroups_using_phylogenetics(
                    orthofinder_fasta_dir, orthofinder_tsv_file, orthofinder_tsv_singletons_file,
                    resdog_dir, resulting_dogs_file, logger,
                    skip_merge_back_flag=args.skip_merge_back,
                    rooting_seeds=args.rooting_seeds,
                    fixation_index_cutoff=args.fixation_index_cutoff,
                    threads=args.threads, 
                    aligner=args.aligner
                )
            
            with open(step2_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 2 checkpoint found - skipping domain ortholog group resolution"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 3: Determine protein-resolution ortholog groups from domain-resolution ortholog groups
        step3_checkpoint_file = os.path.join(checkdir, "Step3.txt")
        msg = '\n--------------------\nStep 3\n--------------------\nDetermine protein-resolution ortholog groups.'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        protein_cluster_dir = os.path.join(args.output_dir, "Protein_Clustering/")
        resulting_ogs_file = os.path.join(findir, "Protein_Ortholog_Groups.tsv")
        
        if not os.path.isfile(step3_checkpoint_file):
            setup_ready_directory([protein_cluster_dir], overwrite_mode="overwrite")
            determine_protein_orthogroups(
                resulting_dogs_file, protein_cluster_dir, resulting_ogs_file, 
                logger, dj=args.dog_jaccard, threads=args.threads
            )
            
            with open(step3_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 3 checkpoint found - skipping protein ortholog group determination"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 4: Determine surrounding contexts for each ortholog group
        step4_checkpoint_file = os.path.join(checkdir, "Step4.txt")
        msg = '\n--------------------\nStep 4\n--------------------\nAssessing neighborhoods for each ortholog group.'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        surround_info_dir = os.path.join(args.output_dir, "Surrounding_Context_Information/")
        og_context_info_file = os.path.join(surround_info_dir, "OrthoGroup_Contexts.tsv")
        
        if not os.path.isfile(step4_checkpoint_file):
            setup_ready_directory([surround_info_dir], overwrite_mode="overwrite")
            determine_ortholog_group_contexts(
                args.bofasa_prep_dir, og_context_info_file, surround_info_dir, 
                resulting_ogs_file, logger, 
                surrounding_bp=args.surrounding_bp, 
                threads=args.threads
            )
            
            with open(step4_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 4 checkpoint found - skipping ortholog group context determination"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 5: Create final summarizations
        step5_checkpoint_file = os.path.join(checkdir, "Step5.txt")
        msg = '\n--------------------\nStep 5\n--------------------\nCreating final consolidated report.'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        final_result_file = os.path.join(findir, "Orthogroup_Overview.xlsx")
        
        if not os.path.isfile(step5_checkpoint_file):
            create_final_report(
                args.bofasa_prep_dir, og_context_info_file, final_result_file, logger
            )
            
            with open(step5_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 5 checkpoint found - skipping final report creation"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 6: Create final visuals
        step6_checkpoint_file = os.path.join(checkdir, "Step6.txt")
        msg = (
            '\n--------------------\nStep 6\n--------------------\n'
            'Creating final summary visualization of conservation vs. '
            'neighborhood syntenic conservation.'
        )
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        tmp_result_file = os.path.join(surround_info_dir, "Simplified_Info_for_Plotting.tsv")
        final_result_plot = os.path.join(findir, "Orthogroup_Conservation_vs_ContextEntropy.html")
        
        if not os.path.isfile(step6_checkpoint_file):
            create_final_visual(
                args.bofasa_prep_dir, tmp_result_file, og_context_info_file, 
                final_result_plot, logger
            )
            
            with open(step6_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 6 checkpoint found - skipping final visualization creation"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 7: (Optional) Create profile-HMM database
        if args.og_consensus:
            step7_checkpoint_file = os.path.join(checkdir, "Step7.txt")
            msg = '\n--------------------\nStep 7\n--------------------\nConstructing ortholog group profile HMMs and consensus sequences.'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
            
            og_seqs_dir = os.path.join(args.output_dir, "Ortholog_Group_Sequences/")
            og_algn_dir = os.path.join(args.output_dir, "Ortholog_Group_Alignments/")
            og_hmms_dir = os.path.join(args.output_dir, "Ortholog_Group_HMMs/")
            og_cons_dir = os.path.join(args.output_dir, "Ortholog_Group_Consensus_Sequences/")
            concatenated_consensus_seqs_file = os.path.join(findir, "Orthogroup_Consensus_Sequences.faa")
            
            if not os.path.isfile(step7_checkpoint_file):
                setup_ready_directory([og_seqs_dir, og_algn_dir, og_hmms_dir, og_cons_dir], overwrite_mode="overwrite")
                create_protein_alignments(
                    args.bofasa_prep_dir, resulting_ogs_file, og_seqs_dir, og_algn_dir,
                    logger, threads=args.threads,
                    aligner=args.aligner
                )
                create_profile_hmms_and_consensus_seqs(
                    og_algn_dir, og_hmms_dir, og_cons_dir, concatenated_consensus_seqs_file, logger, threads=args.threads
                )
                
                with open(step7_checkpoint_file, 'w') as f:
                    f.write("DONE")
            else:
                msg = "Step 7 checkpoint found - skipping profile HMM and consensus sequence creation"
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')

        # Step 8: (Optional) Create core genome alignment
        if args.core_genome:
            step8_checkpoint_file = os.path.join(checkdir, "Step8.txt")
            msg = '\n--------------------\nStep 8\n--------------------\nCreating core genome alignment from (near) single-copy-core domain ortholog groups.'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
            
            dogs_seqs_dir = os.path.join(args.output_dir, "NearSCC_Coarse_Domain_Ortholog_Group_Sequences/")
            dogs_algn_dir = os.path.join(args.output_dir, "NearSCC_Coarse_Domain_Ortholog_Group_Alignments/")
            dogs_trim_dir = os.path.join(findir, "NearSCC_Coarse_Domain_Ortholog_Group_Trimmed_Alignments/")
            merged_core_genome_file = os.path.join(findir, "Core_Genome_Alignment.faa")
            
            if not os.path.isfile(step8_checkpoint_file):
                setup_ready_directory([dogs_seqs_dir, dogs_algn_dir, dogs_trim_dir], overwrite_mode="overwrite")
                create_near_scc_resolved_domain_protein_alignments(
                    args.bofasa_prep_dir, resulting_dogs_file, dogs_seqs_dir,
                    dogs_algn_dir, dogs_trim_dir, merged_core_genome_file,
                    logger, near_scc_prop=args.near_scc_prop, 
                    threads=args.threads,
                    aligner=args.aligner
                )
                
                with open(step8_checkpoint_file, 'w') as f:
                    f.write("DONE")
            else:
                msg = "Step 8 checkpoint found - skipping core genome alignment creation"
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')

        if logger:
            logger.info("bofasa analysis completed successfully")
        else:
            sys.stdout.write("bofasa analysis completed successfully\n")

    except Exception as e:
        if logger:
            logger.error(f"Error during bofasa analysis: {str(e)}")
        else:
            sys.stdout.write(f"Error during bofasa analysis: {str(e)}\n")
        raise
    finally:
        close_logger_object(logger)


def _extracted_mge_ids(sample_ids) -> set:
    """
    Return the IDs among sample_ids that are MGE genomes extracted by prep Step 2a.5,
    i.e. named '<parent>_plasmid_N' / '<parent>_phage_N' where <parent> is itself a sample.
    """
    import re
    mge_name_pattern = re.compile(r'^(.+)_(plasmid|phage)_\d+$')
    sample_ids = set(sample_ids)
    parent_ids = set(s for s in sample_ids if not mge_name_pattern.match(s))
    mge_ids = set()
    for s in sample_ids:
        match = mge_name_pattern.match(s)
        if match and match.group(1) in parent_ids:
            mge_ids.add(s)
    return mge_ids


def run_bofasa_prep(args: argparse.Namespace) -> None:
    """Run the bofasa preparation step."""
    # Import here to avoid circular imports and heavy dependencies
    import os
    import shutil
    from .utils import (
        create_logger_object, close_logger_object, memory_limit,
        check_bofasa_db_setup, check_genomad_setup, setup_ready_directory,
        log_parameters_to_file
    )
    from .processing import (
        run_gene_calling, process_genomes_as_genbanks, 
        process_annotation_directories
    )
    from .analysis import (
        determine_phages_and_plasmids, annotate_is_finder, 
        annotate_and_split_proteins_using_pfam
    )

    # Validate input arguments
    input_genomes = getattr(args, 'input_genomes', [])
    annotation_dirs = getattr(args, 'annotation_dirs', [])
    
    if not input_genomes and not annotation_dirs:
        sys.stdout.write("Error: Either input genomes (-i/--input-genomes) or annotation directories (-a/--annotation-dirs) must be provided.\n")
        sys.exit(1)
    
    if input_genomes and annotation_dirs:
        sys.stdout.write("Warning: Both input genomes and annotation directories provided. Processing both types of input.\n")

    if args.genomad_splits < 1:
        sys.stdout.write("Error: --genomad-splits must be at least 1.\n")
        sys.exit(1)
    
    # Validate extract-mge-genomes requires run-genomad
    if args.extract_mge_genomes and not args.run_genomad:
        sys.stdout.write("Error: --extract-mge-genomes requires --run-genomad to be enabled.\n")
        sys.stdout.write("Please add the --run-genomad flag to use --extract-mge-genomes.\n")
        sys.exit(1)

    sys.stdout.write("Checking bofasa database setup...\n")
    dbs_ok, dbs_message = check_bofasa_db_setup(require_genomad=args.run_genomad)
    if not dbs_ok:
        sys.stdout.write(f"Error: {dbs_message}\n")
        sys.exit(1)
    sys.stdout.write("✓ bofasa databases are properly set up\n")

    # Check if output directory already exists
    if os.path.exists(args.output_dir) and os.listdir(args.output_dir):
        sys.stdout.write(f"\nWARNING: The output directory '{args.output_dir}' already exists and contains files!\n")
        sys.stdout.write("This process may overwrite existing files.\n")

        if args.auto:
            sys.stdout.write("--auto specified: continuing with existing output directory...\n")
        else:
            while True:
                response = input("Do you want to continue and potentially overwrite existing files? (y/N): ").strip().lower()
                if response in ['y', 'yes']:
                    sys.stdout.write("Continuing with existing output directory...\n")
                    break
                elif response in ['n', 'no', '']:
                    sys.stdout.write("Exiting. Please choose a different output directory or remove existing files.\n")
                    sys.exit(1)
                else:
                    sys.stdout.write("Please enter 'y' for yes or 'n' for no.\n")

    # Check genomad setup early in the workflow (only needed with --run-genomad)
    genomad_ok, genomad_message = False, "not checked (--run-genomad not specified)"
    if args.run_genomad:
        sys.stdout.write("Checking genomad setup...\n")
        genomad_ok, genomad_message = check_genomad_setup()
        if not genomad_ok:
            sys.stdout.write(f"Error: {genomad_message}\n")
            sys.stdout.write("geNomad is required when --run-genomad is specified.\n")
            sys.stdout.write("Please run 'bofasa setup' (without --skip-genomad) first to configure geNomad, or drop --run-genomad.\n")
            sys.exit(1)
        else:
            sys.stdout.write("✓ genomad is properly set up\n")

    # Set memory limit
    if args.max_memory:
        memory_limit(args.max_memory)

    # Create output directory first
    setup_ready_directory([args.output_dir], overwrite_mode="skip")

    # Create checkpoint directory
    checkdir = os.path.join(args.output_dir, "Checkpoint_Files/")
    setup_ready_directory([checkdir], overwrite_mode="skip")

    # Create logger
    log_file = os.path.join(args.output_dir, "bofasa_prep.log")
    logger = create_logger_object(log_file)
    
    # Check if logger was created successfully
    if logger is None:
        sys.stdout.write("Warning: Failed to create logger. Continuing without logging.\n")
        logger = None

    try:
        if logger:
            logger.info("Starting bofasa preparation")
            if input_genomes:
                logger.info(f"Number of input genomes: {len(input_genomes)}")
            if annotation_dirs:
                logger.info(f"Number of annotation directories: {len(annotation_dirs)}")
            logger.info(f"Output directory: {args.output_dir}")
            logger.info(f"Threads: {args.threads}")
            logger.info(f"Gene calling method: {args.gene_calling_method}")
            logger.info(f"genomad setup status: {genomad_ok} - {genomad_message}")
        else:
            sys.stdout.write("Starting bofasa preparation\n")
            if input_genomes:
                sys.stdout.write(f"Input genomes: {input_genomes}\n")
            if annotation_dirs:
                sys.stdout.write(f"Annotation directories: {annotation_dirs}\n")
            sys.stdout.write(f"Output directory: {args.output_dir}\n")
            sys.stdout.write(f"Threads: {args.threads}\n")
            sys.stdout.write(f"Gene calling method: {args.gene_calling_method}\n")
            sys.stdout.write(f"genomad setup status: {genomad_ok} - {genomad_message}\n")

        # Log parameters
        if logger:
            logger.info("Saving parameters for future records.")
            parameters_file = os.path.join(args.output_dir, "Parameter_Inputs.txt")
            parameter_values = [
                input_genomes,
                annotation_dirs,
                args.output_dir,
                args.rename_locus_tags,
                args.meta_mode,
                args.run_genomad,
                args.genomad_splits,
                args.extract_mge_genomes,
                args.ignore_upperbound_limit,
                args.gene_calling_method,
                args.min_length,
                args.skip_domain_splitting,
                not args.no_trim_domain_overlaps,
                args.threads,
                args.max_memory
            ]
            parameter_names = [
                "Input Genomes Files",
                "Annotation Directories", 
                "Output Directory",
                "Rename Locus Tags?",
                "Draft Mode?",
                "Run geNomad?",
                "geNomad Splits",
                "Extract MGE Genomes?",
                "Ignore Upperbound Genome Limit?",
                "Gene Calling Method",
                "Minimum Length of Domain/Inter-Domain Unit",
                "Skip Domain Splitting?",
                "Trim Domain Overlaps?",
                "Number of Threads",
                "Maximum Memory in GB"
            ]
            log_parameters_to_file(parameters_file, parameter_names, parameter_values)
            logger.info("Done saving parameters!")

        # Step 1: Process genomes
        step1_checkpoint_file = os.path.join(checkdir, "Step1.txt")
        msg = '\n--------------------\nStep 1\n--------------------\nProcessing genomes / annotation folders provided as input.'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')
        
        gp_dir = os.path.join(args.output_dir, "Genome_Processing/")
        faa_dir = os.path.join(gp_dir, "Proteomes/")
        bed_dir = os.path.join(gp_dir, "BEDs/")
        wgs_dir = os.path.join(gp_dir, "Genomes/")

        # Initialize sample mappings
        sample_wgs = {}
        sample_proteomes = {}
        sample_beds = {}
        sample_is_mge = set()  # Track which samples are MGEs (plasmids/phages)

        if not os.path.isfile(step1_checkpoint_file):
            setup_ready_directory([gp_dir, faa_dir, bed_dir, wgs_dir], overwrite_mode="overwrite")

            # Process input genomes if provided
            if input_genomes:
                # Determine file types and process accordingly
                fasta_files = []
                genbank_files = []

                for genome_file in input_genomes:
                    if genome_file.split('.')[-1].lower() in ["gbk", "gb", "genbank", "gbff"]:
                        genbank_files.append(genome_file)
                    elif genome_file.split('.')[-1].lower() in ["fasta", "fa", "fna"]:
                        fasta_files.append(genome_file)
                    else:
                        if logger:
                            logger.warning(f"{genome_file} is not a valid genome file. Skipping...")
                        else:
                            sys.stdout.write(f"Warning: {genome_file} is not a valid genome file. Skipping...\n")

                # Process FASTA files with Prodigal
                if fasta_files:
                    if logger:
                        logger.info(f"Processing {len(fasta_files)} FASTA files with {args.gene_calling_method}...")
                    else:
                        sys.stdout.write(f"Processing {len(fasta_files)} FASTA files with {args.gene_calling_method}...\n")
                    
                    # Convert list of fasta files to dictionary with sample names
                    fasta_dict = {}
                    for fasta_file in fasta_files:
                        # Use the original filename without extension as sample name
                        sample_name = os.path.splitext(os.path.basename(fasta_file))[0]
                        fasta_dict[sample_name] = fasta_file
                        sample_wgs[sample_name] = fasta_file
                    
                    run_gene_calling(
                        fasta_dict,
                        gp_dir,
                        logger,
                        threads=args.threads,
                        locus_tag_length=args.locus_tag_length,
                        gene_calling_method=args.gene_calling_method,
                        meta_mode=args.meta_mode,
                        min_length=args.min_length
                    )

                    # Move files to proper locations
                    for sample in fasta_dict:
                        faa = os.path.join(faa_dir, f"{sample}.faa")
                        bed = os.path.join(bed_dir, f"{sample}.bed")
                        fna = os.path.join(wgs_dir, f"{sample}.fna")
                        try:
                            assert os.path.isfile(faa) and os.path.isfile(bed) and os.path.isfile(fna)
                            sample_wgs[sample] = "Genome_Processing/Genomes/" + f"{sample}.fna"
                            sample_proteomes[sample] = "Genome_Processing/Proteomes/" + f"{sample}.faa"
                            sample_beds[sample] = "Genome_Processing/BEDs/" + f"{sample}.bed"
                        except Exception as e:
                            msg = f'Unable to validate proper processing for sample {sample}, skipping it'
                            if logger:
                                logger.info(msg)
                            else:
                                sys.stdout.write(msg + '\n')
                            sample_wgs.pop(sample, None)
                            sample_proteomes.pop(sample, None)
                            sample_beds.pop(sample, None)

                # Process GenBank files
                if genbank_files:
                    if logger:
                        logger.info(f"Processing {len(genbank_files)} GenBank files...")
                    else:
                        sys.stdout.write(f"Processing {len(genbank_files)} GenBank files...\n")
                    
                    # Convert list of genbank files to dictionary with sample names
                    genbank_dict = {}
                    for genbank_file in genbank_files:
                        # Use the original filename without extension as sample name
                        sample_name = os.path.splitext(os.path.basename(genbank_file))[0]
                        genbank_dict[sample_name] = genbank_file
                    
                    process_genomes_as_genbanks(
                        genbank_dict,
                        gp_dir,
                        logger,
                        threads=args.threads,
                        locus_tag_length=args.locus_tag_length,
                        rename_locus_tags=args.rename_locus_tags,
                        min_length=args.min_length
                    )

                    # Move files to proper locations
                    for sample in genbank_dict:
                        faa = os.path.join(gp_dir, f"{sample}.faa")
                        bed = os.path.join(gp_dir, f"{sample}.bed")
                        fna = os.path.join(gp_dir, f"{sample}.fna")
                        try:
                            assert os.path.isfile(faa) and os.path.isfile(bed) and os.path.isfile(fna)
                            assert os.path.getsize(faa) > 0 and os.path.getsize(bed) > 0 and os.path.getsize(fna) > 0
                            renamed_faa = os.path.join(faa_dir, f"{sample}.faa")
                            renamed_bed = os.path.join(bed_dir, f"{sample}.bed")
                            renamed_fna = os.path.join(wgs_dir, f"{sample}.fna")
                            shutil.move(faa, renamed_faa)
                            shutil.move(bed, renamed_bed)
                            shutil.move(fna, renamed_fna)
                            sample_proteomes[sample] = "Genome_Processing/Proteomes/" + f"{sample}.faa"
                            sample_beds[sample] = "Genome_Processing/BEDs/" + f"{sample}.bed"
                            sample_wgs[sample] = "Genome_Processing/Genomes/" + f"{sample}.fna"
                        except Exception as e:
                            msg = f'Unable to validate proper processing for sample {sample}, skipping it'
                            if logger:
                                logger.info(msg)
                            else:
                                sys.stdout.write(msg + '\n')
                            sample_wgs.pop(sample, None)
                            sample_proteomes.pop(sample, None)
                            sample_beds.pop(sample, None)

            # Process annotation directories (Prokka/Bakta)
            if annotation_dirs:
                if logger:
                    logger.info(f"Processing {len(annotation_dirs)} annotation directories...")
                else:
                    sys.stdout.write(f"Processing {len(annotation_dirs)} annotation directories...\n")

                # Process annotation directories and get sample mappings
                annotation_results = process_annotation_directories(
                    annotation_dirs,
                    gp_dir,
                    logger,
                    locus_tag_length=args.locus_tag_length,
                    rename_locus_tags=args.rename_locus_tags,
                    threads=args.threads,
                    min_length=args.min_length
                )
                
                # Update sample mappings
                sample_wgs.update(annotation_results.get('sample_wgs', {}))
                sample_proteomes.update(annotation_results.get('sample_proteomes', {}))
                sample_beds.update(annotation_results.get('sample_beds', {}))

            # Create Step 1 checkpoint
            with open(step1_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 1 checkpoint found - skipping genome processing"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
            
            # Load existing data from previous run
            gp_dir = os.path.join(args.output_dir, "Genome_Processing/")
            faa_dir = os.path.join(gp_dir, "Proteomes/")
            bed_dir = os.path.join(gp_dir, "BEDs/")
            wgs_dir = os.path.join(gp_dir, "Genomes/")

            # Load from existing proteome files
            if os.path.exists(faa_dir):
                for faa_file in sorted(os.listdir(faa_dir)):
                    if faa_file.endswith('.faa'):
                        sample_name = os.path.splitext(faa_file)[0]
                        sample_proteomes[sample_name] = os.path.join('Genome_Processing/Proteomes', faa_file)
            
            # Load from existing bed files
            if os.path.exists(bed_dir):
                for bed_file in sorted(os.listdir(bed_dir)):
                    if bed_file.endswith('.bed'):
                        sample_name = os.path.splitext(bed_file)[0]
                        sample_beds[sample_name] = os.path.join('Genome_Processing/BEDs', bed_file)
            
            # Load genome files (check both .fna files and original input files)
            if os.path.exists(wgs_dir):
                for fna_file in sorted(os.listdir(wgs_dir)):
                    if fna_file.endswith('.fna'):
                        sample_name = os.path.splitext(fna_file)[0]
                        sample_wgs[sample_name] = os.path.join('Genome_Processing/Genomes', fna_file)

            # MGE genomes extracted in Step 2a.5 live in the same directories. They are
            # owned by Step 2a.5 (reloaded there if its checkpoint exists, regenerated
            # otherwise) and must not be picked up here as ordinary genomes - in
            # particular partial files left by an interrupted Step 2a.5.
            if args.extract_mge_genomes:
                all_ids = set(sample_wgs) | set(sample_proteomes) | set(sample_beds)
                for mge_id in _extracted_mge_ids(all_ids):
                    sample_wgs.pop(mge_id, None)
                    sample_proteomes.pop(mge_id, None)
                    sample_beds.pop(mge_id, None)

        msg = f'Found and successfully processed {len(sample_proteomes)} genomes, continuing ...'
        if logger:
            logger.info(msg)
        else:
            sys.stdout.write(msg + '\n')

        # Check if we have enough genomes for analysis
        if len(sample_proteomes) < 4:
            error_msg = f"Error: bofasa prep requires at least 4 genomes, but only {len(sample_proteomes)} were successfully processed."
            if logger:
                logger.error(error_msg)
            else:
                sys.stdout.write(error_msg + '\n')
            sys.stdout.write("Please provide at least 4 valid genome files or annotation directories.\n")
            sys.exit(1)

        if len(sample_proteomes) > config.DEFAULT_MAX_GENOMES and not args.ignore_upperbound_limit:
            error_msg = f"Error: Upper limit of genomes bofasa (default of {config.DEFAULT_MAX_GENOMES}) is designed to handle was exceeded."
            if logger:
                logger.error(error_msg)
            else:
                sys.stdout.write(error_msg + '\n')
            sys.stdout.write(f"Please provide at most {config.DEFAULT_MAX_GENOMES} valid genome files or annotation directories or alternatively issue the --ignore-upperbound-limit flag.\n")
            sys.exit(1)

        for sample in sample_proteomes:
            sample_proteomes[sample] = os.path.join(args.output_dir, sample_proteomes[sample])

        for sample in sample_beds:
            sample_beds[sample] = os.path.join(args.output_dir, sample_beds[sample])
        
        for sample in sample_wgs:
            sample_wgs[sample] = os.path.join(args.output_dir, sample_wgs[sample])

        # Step 2a: Running geNomad (if requested and genome files are available)
        step2a_checkpoint_file = os.path.join(checkdir, "Step2a.txt")
        if args.run_genomad and sample_wgs:
            if not os.path.isfile(step2a_checkpoint_file):
                msg = '\n--------------------\nStep 2a\n--------------------\nRunning geNomad for phage and plasmid identification.'
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')

                genomad_dir = os.path.join(args.output_dir, "geNomad_Annotations/")
                setup_ready_directory([genomad_dir], overwrite_mode="overwrite")

                phage_protein_listing_file = os.path.join(args.output_dir, "Sample_Phage_Proteins.txt")
                plasmid_protein_listing_file = os.path.join(args.output_dir, "Sample_Plasmid_Proteins.txt")
                
                determine_phages_and_plasmids(
                    sample_wgs, sample_beds, genomad_dir, 
                    phage_protein_listing_file, plasmid_protein_listing_file, 
                    logger, threads=args.threads,
                    genome_splits=args.genomad_splits
                )
                
                # Create Step 2a checkpoint
                with open(step2a_checkpoint_file, 'w') as f:
                    f.write("DONE")
            else:
                msg = "Step 2a checkpoint found - skipping geNomad execution"
                if logger:
                    logger.info(msg)
                else:
                    sys.stdout.write(msg + '\n')
            
            # Step 2a.5: Extract MGE genomes if requested
            if args.extract_mge_genomes:
                step2a5_checkpoint_file = os.path.join(checkdir, "Step2a5.txt")
                if not os.path.isfile(step2a5_checkpoint_file):
                    msg = '\n--------------------\nStep 2a.5\n--------------------\nExtracting plasmid and phage genomes as separate entities.'
                    if logger:
                        logger.info(msg)
                    else:
                        sys.stdout.write(msg + '\n')
                    
                    genomad_dir = os.path.join(args.output_dir, "geNomad_Annotations/")

                    # Remove MGE files left by an interrupted earlier Step 2a.5 so that a
                    # partial set cannot linger alongside the freshly extracted one.
                    gp_dir_mge = os.path.join(args.output_dir, "Genome_Processing")
                    leftover_ids = set()
                    for sub, ext in [("Genomes", ".fna"), ("Proteomes", ".faa"), ("BEDs", ".bed")]:
                        sub_dir = os.path.join(gp_dir_mge, sub)
                        if os.path.isdir(sub_dir):
                            leftover_ids.update(f[:-len(ext)] for f in sorted(os.listdir(sub_dir)) if f.endswith(ext))
                    leftover_mges = _extracted_mge_ids(leftover_ids | set(sample_wgs))
                    for mge_id in leftover_mges:
                        for sub, ext in [("Genomes", ".fna"), ("Proteomes", ".faa"), ("BEDs", ".bed")]:
                            leftover_file = os.path.join(gp_dir_mge, sub, mge_id + ext)
                            if os.path.isfile(leftover_file):
                                os.remove(leftover_file)
                        sample_wgs.pop(mge_id, None)
                        sample_proteomes.pop(mge_id, None)
                        sample_beds.pop(mge_id, None)
                    if leftover_mges:
                        msg = f'Removed files of {len(leftover_mges)} MGE genomes left by an interrupted previous run'
                        if logger:
                            logger.info(msg)
                        else:
                            sys.stdout.write(msg + '\n')

                    # Extract MGE genomes with annotations from parent genomes
                    # This avoids re-running gene calling on small MGE sequences
                    from .analysis import extract_mge_annotations_from_parent_genomes, update_mge_protein_listings_from_proteomes
                    mge_data = extract_mge_annotations_from_parent_genomes(
                        genomad_dir, sample_wgs, sample_proteomes, sample_beds, 
                        args.output_dir, logger
                    )
                    
                    # Add MGE genomes to sample tracking (skip empty proteome files)
                    num_mges = len(mge_data['mge_wgs'])
                    num_skipped = 0
                    if num_mges > 0:
                        for mge_id in mge_data['mge_wgs']:
                            proteome_file = mge_data['mge_proteomes'][mge_id]
                            # Only add if proteome file is not empty
                            if os.path.isfile(proteome_file) and os.path.getsize(proteome_file) > 0:
                                sample_wgs[mge_id] = mge_data['mge_wgs'][mge_id]
                                sample_proteomes[mge_id] = proteome_file
                                sample_beds[mge_id] = mge_data['mge_beds'][mge_id]
                                sample_is_mge.add(mge_id)  # Mark this sample as an MGE
                            else:
                                num_skipped += 1
                        
                        num_valid = num_mges - num_skipped
                        msg = f'Successfully extracted annotations for {num_valid} MGE genomes from parent genomes'
                        if num_skipped > 0:
                            msg += f' ({num_skipped} MGEs skipped due to no extractable genes)'
                        if logger:
                            logger.info(msg)
                        else:
                            sys.stdout.write(msg + '\n')
                        
                        # Update MGE protein listing files based on extracted proteomes
                        phage_protein_listing_file = os.path.join(args.output_dir, "Sample_Phage_Proteins.txt")
                        plasmid_protein_listing_file = os.path.join(args.output_dir, "Sample_Plasmid_Proteins.txt")
                        update_mge_protein_listings_from_proteomes(
                            args.output_dir, phage_protein_listing_file, 
                            plasmid_protein_listing_file, logger
                        )
                    
                    # Create Step 2a.5 checkpoint
                    with open(step2a5_checkpoint_file, 'w') as f:
                        f.write("DONE")
                else:
                    msg = "Step 2a.5 checkpoint found - skipping MGE genome extraction"
                    if logger:
                        logger.info(msg)
                    else:
                        sys.stdout.write(msg + '\n')
                    
                    # Need to reload MGE genomes from previous run. MGE genome/proteome/BED files are
                    # written into Genome_Processing/{Genomes,Proteomes,BEDs}/ and named
                    # '<parent>_plasmid_N' / '<parent>_phage_N' (the Extracted_*_Genomes/ dirs only hold
                    # copies of the sequences). The Step 1 reload above may already have picked them up.
                    import re
                    mge_name_pattern = re.compile(r'^(.+)_(plasmid|phage)_\d+$')
                    mge_metadata = {}
                    mge_metadata_file = os.path.join(args.output_dir, "Sample_MGE_Metadata.txt")
                    if os.path.exists(mge_metadata_file):
                        with open(mge_metadata_file) as mge_handle:
                            for i, line in enumerate(mge_handle):
                                if i == 0:  # Skip header
                                    continue
                                parts = line.strip().split('\t')
                                if len(parts) == 2:
                                    mge_metadata[parts[0]] = (parts[1] == "True")

                    candidate_ids = set()
                    for cand_dir, cand_ext in [(wgs_dir, '.fna'), (faa_dir, '.faa'), (bed_dir, '.bed')]:
                        if os.path.exists(cand_dir):
                            for cand_file in sorted(os.listdir(cand_dir)):
                                if cand_file.endswith(cand_ext):
                                    candidate_ids.add(cand_file[:-len(cand_ext)])
                    parent_ids = set(c for c in candidate_ids if not mge_name_pattern.match(c))

                    num_reloaded = 0
                    for mge_id in sorted(candidate_ids):
                        mge_match = mge_name_pattern.match(mge_id)
                        if not mge_match or mge_match.group(1) not in parent_ids:
                            continue
                        if mge_metadata and mge_id in mge_metadata and not mge_metadata[mge_id]:
                            continue
                        mge_fna = os.path.join(wgs_dir, f"{mge_id}.fna")
                        mge_faa = os.path.join(faa_dir, f"{mge_id}.faa")
                        mge_bed = os.path.join(bed_dir, f"{mge_id}.bed")
                        if (os.path.isfile(mge_fna) and os.path.isfile(mge_faa) and os.path.isfile(mge_bed)
                                and os.path.getsize(mge_faa) > 0):
                            sample_wgs[mge_id] = mge_fna
                            sample_proteomes[mge_id] = mge_faa
                            sample_beds[mge_id] = mge_bed
                            sample_is_mge.add(mge_id)  # Mark this sample as an MGE
                            num_reloaded += 1
                        else:
                            # Mirror the non-resume path: MGEs without extractable genes are skipped
                            sample_wgs.pop(mge_id, None)
                            sample_proteomes.pop(mge_id, None)
                            sample_beds.pop(mge_id, None)
                    msg = f'Reloaded {num_reloaded} MGE genomes from previous run'
                    if logger:
                        logger.info(msg)
                    else:
                        sys.stdout.write(msg + '\n')

        elif args.run_genomad and not sample_wgs:
            msg = (
                'Warning: geNomad requested but no genome files available '
                '(only annotation directories provided). Skipping geNomad step.'
            )
            if logger:
                logger.warning(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 2b: Annotate proteins with ISfinder databases
        step2b_checkpoint_file = os.path.join(checkdir, "Step2b.txt")
        if not os.path.isfile(step2b_checkpoint_file):
            msg = '\n--------------------\nStep 2b\n--------------------\nAnnotating IS elements using the ISFinder database.'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

            annot_dir = os.path.join(args.output_dir, "ISFinder_Annotations/")
            setup_ready_directory([annot_dir], overwrite_mode="overwrite")

            isfinder_protein_listing_file = os.path.join(args.output_dir, "Sample_IS_Element_Proteins.txt")
            annotate_is_finder(
                sample_proteomes, annot_dir, isfinder_protein_listing_file, 
                logger, threads=args.threads
            )
            
            # Create Step 2b checkpoint
            with open(step2b_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 2b checkpoint found - skipping ISfinder annotation"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        # Step 3: Annotate proteins with Pfam and split into domain/inter-domain units
        step3_checkpoint_file = os.path.join(checkdir, "Step3.txt")
        if not os.path.isfile(step3_checkpoint_file):
            msg = '\n--------------------\nStep 3\n--------------------\nAnnotating proteins with Pfam HMMs to split them into domain/inter-domain "chunks" or units.'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

            split_proteins_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_FASTAs/")
            domain_coords_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_Coordinates/")
            domain_coord_info_file = os.path.join(args.output_dir, "Sample_Domain_and_InterDomain_Information.txt")
            setup_ready_directory([split_proteins_dir, domain_coords_dir], overwrite_mode="overwrite")
            
            sample_ccds_proteomes = annotate_and_split_proteins_using_pfam(
                sample_proteomes, split_proteins_dir, domain_coords_dir, 
                domain_coord_info_file, logger, 
                sample_is_mge=sample_is_mge,
                minimal_length=args.min_length, 
                threads=args.threads, 
                skip_domain_splitting=args.skip_domain_splitting,
                trim_domain_overlaps=not args.no_trim_domain_overlaps
            )

            if not sample_ccds_proteomes:
                msg = ("Pfam annotation and domain splitting produced no domain-chopped proteomes - "
                       "see the log for details. Not writing the Step 3 checkpoint.")
                if logger:
                    logger.error(msg)
                sys.stderr.write(msg + '\n')
                sys.exit(1)

            # Save sample_is_mge metadata for downstream use
            mge_metadata_file = os.path.join(args.output_dir, "Sample_MGE_Metadata.txt")
            with open(mge_metadata_file, 'w') as mge_handle:
                mge_handle.write("Sample\tIs_MGE\n")
                for sample in sample_proteomes:
                    is_mge = "True" if sample in sample_is_mge else "False"
                    mge_handle.write(f"{sample}\t{is_mge}\n")
            
            # Create Step 3 checkpoint
            with open(step3_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 3 checkpoint found - skipping Pfam annotation and domain splitting"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')
            
            # Load existing data from previous run
            split_proteins_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_FASTAs/")
            domain_coords_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_Coordinates/")
            domain_coord_info_file = os.path.join(args.output_dir, "Sample_Domain_and_InterDomain_Information.txt")
            
            # Reconstruct sample_ccds_proteomes from existing files (single directory)
            sample_ccds_proteomes = {}
            if os.path.exists(split_proteins_dir):
                for ccds_file in sorted(os.listdir(split_proteins_dir)):
                    if ccds_file.endswith('.ccds.faa'):
                        sample_name = ccds_file.replace('.ccds.faa', '')
                        sample_ccds_proteomes[sample_name] = os.path.join(split_proteins_dir, ccds_file)
            
            # Reload sample_is_mge metadata
            mge_metadata_file = os.path.join(args.output_dir, "Sample_MGE_Metadata.txt")
            if os.path.exists(mge_metadata_file):
                with open(mge_metadata_file) as mge_handle:
                    for i, line in enumerate(mge_handle):
                        if i == 0:  # Skip header
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) == 2 and parts[1] == "True":
                            sample_is_mge.add(parts[0])

        # Step 4: Finalize setting up for final bofasa analysis
        step4_checkpoint_file = os.path.join(checkdir, "Step4.txt")
        if not os.path.isfile(step4_checkpoint_file):
            msg = '\n--------------------\nStep 4\n--------------------\nNearly there, putting together final files.'
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

            # Create info file
            genome_info_listing_file = os.path.join(args.output_dir, "Info_on_Input_Genome_Files.txt")
            with open(genome_info_listing_file, "w") as f:
                header = ['sample', 'ccds_faa', 'proteome', 'bed_coord', 'genome']
                f.write('\t'.join(header) + '\n')
                
                for sample in sample_wgs:
                    if sample not in sample_ccds_proteomes or sample not in sample_proteomes or sample not in sample_beds:
                        msg = f'Warning: sample {sample} has no domain-chopped proteome (or proteome/BED) - excluding it from {os.path.basename(genome_info_listing_file)}'
                        if logger:
                            logger.warning(msg)
                        else:
                            sys.stdout.write(msg + '\n')
                        continue
                    fna = 'Genome_Processing/' + sample_wgs[sample].split('/Genome_Processing/')[-1]
                    faa = 'Genome_Processing/' + sample_proteomes[sample].split('/Genome_Processing/')[-1]
                    ccds_faa = 'Domain_and_Interdomain_FASTAs/' + sample_ccds_proteomes[sample].split('/Domain_and_Interdomain_FASTAs/')[-1]
                    bed = 'Genome_Processing/' + sample_beds[sample].split('Genome_Processing/')[-1]
                    sample_info = [str(x) for x in [sample, ccds_faa, faa, bed, fna]]
                    f.write('\t'.join(sample_info) + '\n')
            
            # Create Step 4 checkpoint
            with open(step4_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 4 checkpoint found - skipping final file creation"
            if logger:
                logger.info(msg)
            else:
                sys.stdout.write(msg + '\n')

        if logger:
            logger.info("bofasa preparation completed successfully")
        else:
            sys.stdout.write("bofasa preparation completed successfully\n")

    except Exception as e:
        if logger:
            logger.error(f"Error during bofasa preparation: {str(e)}")
        else:
            sys.stdout.write(f"Error during bofasa preparation: {str(e)}\n")
        raise
    finally:
        close_logger_object(logger)


def run_bofasa_domscape(args: argparse.Namespace) -> None:
    """Run the bofasa domscape program - Pfam domain counts across genomes."""
    # Import here to avoid circular imports and heavy dependencies
    import os
    import shutil
    from .utils import (
        create_logger_object, close_logger_object, memory_limit,
        check_bofasa_db_setup, setup_ready_directory, log_parameters_to_file
    )
    from .processing import (
        run_gene_calling, process_genomes_as_genbanks,
        process_annotation_directories
    )
    from .domscape import annotate_pfam_domains, create_pfam_count_matrix

    input_genomes = args.input_genomes or []
    input_proteomes = args.input_proteomes or []
    annotation_dirs = args.annotation_dirs or []

    if not input_genomes and not input_proteomes and not annotation_dirs:
        sys.stdout.write("Error: At least one of input genomes (-i), input proteomes (-p), or annotation directories (-a) must be provided.\n")
        sys.exit(1)

    if args.batch_size < 1:
        sys.stdout.write("Error: --batch-size must be at least 1.\n")
        sys.exit(1)

    # Sort inputs by type and make sure sample names are unique across all inputs
    fasta_files = {}
    genbank_files = {}
    proteome_files = {}
    sample_sources = {}
    duplicate_samples = set()

    def register_sample(sample_name: str, source: str) -> bool:
        if sample_name in sample_sources:
            duplicate_samples.add(sample_name)
            return False
        sample_sources[sample_name] = source
        return True

    for genome_file in input_genomes:
        sample_name = os.path.splitext(os.path.basename(genome_file))[0]
        suffix = genome_file.split('.')[-1].lower()
        if suffix in ["gbk", "gb", "genbank", "gbff"]:
            if register_sample(sample_name, genome_file):
                genbank_files[sample_name] = genome_file
        elif suffix in ["fasta", "fa", "fna"]:
            if register_sample(sample_name, genome_file):
                fasta_files[sample_name] = genome_file
        else:
            sys.stdout.write(f"Warning: {genome_file} is not a valid genome file. Skipping...\n")

    for proteome_file in input_proteomes:
        sample_name = os.path.splitext(os.path.basename(proteome_file))[0]
        if not os.path.isfile(proteome_file):
            sys.stdout.write(f"Warning: {proteome_file} does not exist. Skipping...\n")
            continue
        if register_sample(sample_name, proteome_file):
            proteome_files[sample_name] = os.path.abspath(proteome_file)

    for annotation_dir in annotation_dirs:
        register_sample(os.path.basename(annotation_dir.rstrip('/')), annotation_dir)

    if duplicate_samples:
        sys.stdout.write("Error: The following sample names occur more than once across the inputs (sample names are derived from file/directory names):\n")
        for sample_name in sorted(duplicate_samples):
            sys.stdout.write(f"  {sample_name}\n")
        sys.exit(1)

    sys.stdout.write("Checking bofasa database setup...\n")
    dbs_ok, dbs_message = check_bofasa_db_setup(require_isfinder=False)
    if not dbs_ok:
        sys.stdout.write(f"Error: {dbs_message}\n")
        sys.exit(1)
    sys.stdout.write("✓ bofasa databases are properly set up\n")

    # Check if output directory already exists
    if os.path.exists(args.output_dir) and os.listdir(args.output_dir):
        sys.stdout.write(f"\nWARNING: The output directory '{args.output_dir}' already exists and contains files!\n")
        sys.stdout.write("Completed steps will be resumed from checkpoints and other files may be overwritten.\n")

        if args.auto:
            sys.stdout.write("--auto specified: continuing with existing output directory...\n")
        else:
            while True:
                response = input("Do you want to continue and potentially overwrite existing files? (y/N): ").strip().lower()
                if response in ['y', 'yes']:
                    sys.stdout.write("Continuing with existing output directory...\n")
                    break
                elif response in ['n', 'no', '']:
                    sys.stdout.write("Exiting. Please choose a different output directory or remove existing files.\n")
                    sys.exit(1)
                else:
                    sys.stdout.write("Please enter 'y' for yes or 'n' for no.\n")

    # Set memory limit
    if args.max_memory:
        memory_limit(args.max_memory)

    setup_ready_directory([args.output_dir], overwrite_mode="skip")
    checkdir = os.path.join(args.output_dir, "Checkpoint_Files/")
    setup_ready_directory([checkdir], overwrite_mode="skip")

    log_file = os.path.join(args.output_dir, "bofasa_domscape.log")
    logger = create_logger_object(log_file)
    if logger is None:
        sys.stdout.write("Error: Failed to create logger.\n")
        sys.exit(1)

    try:
        logger.info("Starting bofasa domscape")
        logger.info(f"Number of input genomes: {len(fasta_files) + len(genbank_files)}")
        logger.info(f"Number of input proteomes: {len(proteome_files)}")
        logger.info(f"Number of annotation directories: {len(annotation_dirs)}")

        parameters_file = os.path.join(args.output_dir, "Parameter_Inputs.txt")
        parameter_values = [
            input_genomes, input_proteomes, annotation_dirs, args.output_dir,
            args.gene_calling_method, args.meta_mode, args.min_length,
            args.batch_size, args.threads, args.max_memory
        ]
        parameter_names = [
            "Input Genomes Files", "Input Proteome Files", "Annotation Directories",
            "Output Directory", "Gene Calling Method", "Meta Mode?",
            "Minimum Protein Length", "Batch Size", "Number of Threads",
            "Maximum Memory in GB"
        ]
        log_parameters_to_file(parameters_file, parameter_names, parameter_values)

        # Step 1: Process genomes / annotation directories into proteomes
        step1_checkpoint_file = os.path.join(checkdir, "Step1.txt")
        sample_listing_file = os.path.join(args.output_dir, "Sample_Proteomes.txt")
        logger.info('\n--------------------\nStep 1\n--------------------\nProcessing genomes / annotation folders / proteomes provided as input.')

        gp_dir = os.path.join(args.output_dir, "Genome_Processing/")
        faa_dir = os.path.join(gp_dir, "Proteomes/")
        bed_dir = os.path.join(gp_dir, "BEDs/")
        wgs_dir = os.path.join(gp_dir, "Genomes/")

        sample_proteomes = {}
        if not os.path.isfile(step1_checkpoint_file):
            sample_proteomes.update(proteome_files)

            if fasta_files or genbank_files or annotation_dirs:
                setup_ready_directory([gp_dir, faa_dir, bed_dir, wgs_dir], overwrite_mode="overwrite")

            if fasta_files:
                logger.info(f"Processing {len(fasta_files)} FASTA files with {args.gene_calling_method}...")
                # Locus tag prefixes must be unique per genome, so grow their length with the number of genomes
                locus_tag_length = config.DEFAULT_LOCUS_TAG_LENGTH
                while 26 ** locus_tag_length < len(fasta_files):
                    locus_tag_length += 1
                run_gene_calling(
                    fasta_files,
                    gp_dir,
                    logger,
                    threads=args.threads,
                    locus_tag_length=locus_tag_length,
                    gene_calling_method=args.gene_calling_method,
                    meta_mode=args.meta_mode,
                    min_length=args.min_length
                )
                for sample in fasta_files:
                    faa = os.path.join(faa_dir, f"{sample}.faa")
                    if os.path.isfile(faa) and os.path.getsize(faa) > 0:
                        sample_proteomes[sample] = os.path.abspath(faa)
                    else:
                        logger.warning(f"Unable to validate proper processing for sample {sample}, skipping it")

            if genbank_files:
                logger.info(f"Processing {len(genbank_files)} GenBank files...")
                process_genomes_as_genbanks(
                    genbank_files,
                    gp_dir,
                    logger,
                    threads=args.threads,
                    rename_locus_tags=False,
                    min_length=args.min_length
                )
                for sample in genbank_files:
                    try:
                        for suffix, dest_dir in [('faa', faa_dir), ('bed', bed_dir), ('fna', wgs_dir)]:
                            src = os.path.join(gp_dir, f"{sample}.{suffix}")
                            assert os.path.isfile(src) and os.path.getsize(src) > 0
                            shutil.move(src, os.path.join(dest_dir, f"{sample}.{suffix}"))
                        sample_proteomes[sample] = os.path.abspath(os.path.join(faa_dir, f"{sample}.faa"))
                    except Exception:
                        logger.warning(f"Unable to validate proper processing for sample {sample}, skipping it")

            if annotation_dirs:
                logger.info(f"Processing {len(annotation_dirs)} annotation directories...")
                annotation_results = process_annotation_directories(
                    annotation_dirs,
                    gp_dir,
                    logger,
                    rename_locus_tags=False,
                    threads=args.threads,
                    min_length=args.min_length
                )
                for sample, faa in annotation_results.get('sample_proteomes', {}).items():
                    sample_proteomes[sample] = os.path.abspath(os.path.join(args.output_dir, faa))

            with open(sample_listing_file, 'w') as sl_handle:
                sl_handle.write("Sample\tProteome\n")
                for sample in sorted(sample_proteomes):
                    sl_handle.write(f"{sample}\t{sample_proteomes[sample]}\n")

            with open(step1_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            logger.info("Step 1 checkpoint found - loading processed proteomes")
            with open(sample_listing_file) as sl_handle:
                next(sl_handle)
                for line in sl_handle:
                    sample, faa = line.rstrip('\n').split('\t')
                    sample_proteomes[sample] = faa

        logger.info(f"Found and successfully processed {len(sample_proteomes)} samples, continuing ...")
        if not sample_proteomes:
            logger.error("No samples were successfully processed - exiting.")
            sys.stderr.write("Error: No samples were successfully processed. See the log for details.\n")
            sys.exit(1)

        # Step 2: Annotate proteomes with Pfam using gathering thresholds
        step2_checkpoint_file = os.path.join(checkdir, "Step2.txt")
        hits_dir = os.path.join(args.output_dir, "Pfam_Hits/")
        logger.info('\n--------------------\nStep 2\n--------------------\nAnnotating proteomes with Pfam HMMs using gathering thresholds.')
        setup_ready_directory([hits_dir], overwrite_mode="skip")
        if not os.path.isfile(step2_checkpoint_file):
            sample_hits_files = annotate_pfam_domains(
                sample_proteomes, hits_dir, logger,
                threads=args.threads, batch_size=args.batch_size
            )
            if len(sample_hits_files) == len(sample_proteomes):
                with open(step2_checkpoint_file, 'w') as f:
                    f.write("DONE")
        else:
            logger.info("Step 2 checkpoint found - skipping Pfam annotation")
            sample_hits_files = {
                sample: os.path.join(hits_dir, f"{sample}.pfam_hits.tsv")
                for sample in sample_proteomes
            }

        if not sample_hits_files:
            logger.error("Pfam annotation failed for all samples - exiting.")
            sys.stderr.write("Error: Pfam annotation failed for all samples. See the log for details.\n")
            sys.exit(1)

        missing_samples = sorted(set(sample_proteomes) - set(sample_hits_files))
        if missing_samples:
            msg = (f"WARNING: Pfam annotation failed for {len(missing_samples)} of {len(sample_proteomes)} samples; "
                   f"they will be missing from the count matrix (rerun with the same output directory to retry them): "
                   + ', '.join(missing_samples))
            logger.warning(msg)
            sys.stderr.write(msg + '\n')

        # Step 3: Build the Pfam by genome count matrix
        logger.info('\n--------------------\nStep 3\n--------------------\nCreating Pfam domain protein count matrix.')
        matrix_file = os.path.join(args.output_dir, "Pfam_Protein_Count_Matrix.tsv")
        pfam_info_file = os.path.join(args.output_dir, "Pfam_Domain_Info.tsv")
        protein_listing_file = os.path.join(args.output_dir, "Pfam_Domain_Proteins.tsv")
        create_pfam_count_matrix(sample_hits_files, matrix_file, pfam_info_file, protein_listing_file, logger)

        sys.stdout.write(f"Pfam protein count matrix: {matrix_file}\n")
        if missing_samples:
            msg = (f"bofasa domscape finished INCOMPLETE: {len(missing_samples)} samples lacked Pfam annotations "
                   f"and are absent from the matrix - exiting with non-zero status.")
            logger.error(msg)
            sys.stderr.write(msg + '\n')
            sys.exit(1)
        logger.info("bofasa domscape completed successfully")

    except Exception as e:
        logger.error(f"Error during bofasa domscape: {str(e)}")
        raise
    finally:
        close_logger_object(logger)


def run_setup_annotation_dbs(args: argparse.Namespace) -> None:
    """Run the setup annotation databases step."""
    # Import here to avoid circular imports
    import os
    from .setup_dbs import setup_annotation_databases

    # Get output directory from environment variable
    output_dir = os.getenv("BOFASA_DB_PATH")
    if not output_dir:
        sys.stdout.write("Error: BOFASA_DB_PATH environment variable is not set.\n")
        sys.stdout.write(
            "Please set BOFASA_DB_PATH to the directory where you want to store annotation databases.\n"
        )
        sys.exit(1)

    output_dir = output_dir.strip()
    if not output_dir:
        sys.stdout.write("Error: BOFASA_DB_PATH environment variable is empty.\n")
        sys.exit(1)

    try:
        setup_ok = setup_annotation_databases(
            output_dir, args.threads, args.force, args.auto, args.skip_genomad,
            pfam_version=args.pfam_version,
        )
    except Exception as e:
        sys.stdout.write(f"Error during setup annotation databases: {str(e)}\n")
        sys.exit(1)
    if not setup_ok:
        sys.stderr.write("Error: bofasa setup did not complete successfully - see messages above.\n")
        sys.exit(1)


def main() -> None:
    """Main entry point for the bofasa command."""
    # Create parser and parse arguments
    parser = create_main_parser()

    # Handle version flag
    if len(sys.argv) > 1 and ("-v" in set(sys.argv) or "--version" in set(sys.argv)):
        sys.stdout.write(config.get_version() + '\n')
        sys.exit(0)

    # Let argparse handle help automatically through our custom parser
    args = parser.parse_args()

    # Handle case where no command is provided
    if not args.command:
        print_bofasa_help()
        sys.exit(1)

    # Dispatch to appropriate function
    if args.command == "run":
        run_bofasa_analysis(args)
    elif args.command == "prep":
        run_bofasa_prep(args)
    elif args.command == "setup":
        run_setup_annotation_dbs(args)
    elif args.command == "domscape":
        run_bofasa_domscape(args)
    else:
        sys.stdout.write(f"Unknown command: {args.command}\n")
        print_bofasa_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
