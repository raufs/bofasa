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
        help_type = getattr(self, 'help_type', None)
        
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
        description="bacterial ortholog finding and syntenic analysis (bofasa)",
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
        "  bofasa setup --threads 8 --force"
    )

    # Add test subcommand
    test_parser = subparsers.add_parser(
        "test", 
        help="Run bofasa prep and run on test genomes",
        description="Run bofasa prep and run on test genomes to verify installation and functionality.",
        formatter_class=BofasaHelpFormatter,
    )
    add_test_arguments(test_parser)
    test_parser.epilog = (
        "Examples:\n"
        "  bofasa test -o test_results\n"
        "  bofasa test -o test_results --threads 8 --cleanup"
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
        "will be important for parsing out sample names.",
        formatter_class=BofasaHelpFormatter,
    )
    add_prep_arguments(prep_parser)
    prep_parser.epilog = (
        "Examples:\n"
        "  bofasa prep -i genome1.fasta genome2.gbk -o prepared_data\n"
        "  bofasa prep -i *.fasta -o prepared_data -c 8 --gene-calling-method prodigal"
    )

    # Add run subcommand
    run_parser = subparsers.add_parser(
        "run",
        help="Execute the main bofasa orthology analysis",
        description="Execute the main bofasa orthology analysis.",
        formatter_class=BofasaHelpFormatter,
    )
    add_run_arguments(run_parser)
    run_parser.epilog = (
        "Examples:\n"
        "  bofasa run -i prepared_data -o analysis_results\n"
        "  bofasa run -i prepared_data -o analysis_results --threads 8 --og-consensus --core-genome\n"
        "  bofasa run -i prepared_data -o analysis_results --surrounding-bp 5000 --dog-jaccard 0.3"
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
        default=10000,
        help="The number of basepairs to look up and downstream for syntenic\n"
        "analysis [Default is 10000].",
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
        default=0.25,
        help="The Jaccard index threshold for domain ortholog group overlap between\n"
        "two pairs of proteins needed to consider them as sharing an edge. Will\n"
        "be automatically lowered for cases where single-copy domain ortholog\n"
        "groups are observed [Default is 0.25].",
    )
    parser.add_argument(
        "-fic",
        "--fixation-index-cutoff",
        type=float,
        default=0.25,
        help="Fixation index cutoff for domain ortholog group re-merging following\n"
        "phylogenetic splitting [Default is 0.25].",
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
        default=1,
        help="The maximum number of nodes to try for rooting. A random sampling is\n"
        "performed [Default is 1; uses midpoint rooting].",
    )
    parser.add_argument(
        "-qa",
        "--quality-alignments",
        action="store_true",
        help="Prioritize quality over speed for constructing domain-resolution\n"
        "ortholog group protein alignments. Uses MUSCLE align mode instead of super5 mode.",
    )
    parser.add_argument(
        "-us",
        "--ultra-sens",
        action="store_true",
        help="Use ultra-sensitivity parameters for running OrthoFinder DIAMOND\n"
        "searches.",
    )
    parser.add_argument(
        "-mi",
        "--mcl-inflation",
        type=float,
        default=1.2,
        help="MCL inflation parameter for determining coarse domain-resolution\n"
        "ortholog groups via OrthoFinder [Default is 1.2].",
    )
    parser.add_argument(
        "-ns",
        "--near-scc-prop",
        type=float,
        default=0.95,
        help="Proportion of genomes which single-copy coarse domain ortholog\n"
        "groups need to be found in for constructing core genome alignment\n"
        "for phylogenomics [Default is 0.95].",
    )
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        default=4,
        help="Total number of cores/threads to use [Default is 4].",
    )
    parser.add_argument(
        "-mrd",
        "--max-recursion-depth",
        type=int,
        default=5000,
        help="The maximum recursion depth.",
    )
    parser.add_argument(
        "-mm",
        "--max-memory",
        type=int,
        default=32,
        help="Uses resource module to set soft memory limit. Provide in Giga-bytes\n"
        "[Default is 32].",
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
        help=f"Memory limit in GB [Default is {config.DEFAULT_MAX_MEMORY}].",
    )
    parser.add_argument(
        "-ml",
        "--min-length",
        type=int,
        default=20,
        help="Minimum length of domain/inter-domain unit to consider [Default is 20].",
    )
    parser.add_argument(
        "-sds",
        "--skip-domain-splitting",
        action="store_true",
        help="Skip domain-annotation based splitting of coding sequences.",
    )
    parser.add_argument(
        "-rg",
        "--run-genomad",
        action="store_true",
        help="Whether to run geNomad - takes time but results can be used to achieve higher quality orthogroup predictions.",
    )
    parser.add_argument(
        "-x",
        "--ignore-upperbound-limit",
        action="store_true",
        help="Ignore the upper bound limit for the number of genomes to process [Default is 200].",
    )


def add_setup_annotation_dbs_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the setup-dbs subcommand."""
    parser.add_argument(
        "-t",
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


def add_test_arguments(parser):
    """Add arguments for the test subcommand."""
    parser.add_argument(
        "-o", "--output-dir", required=True, help="Output directory for test results."
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use [Default is 4].",
    )
    parser.add_argument(
        "-gcm",
        "--gene-calling-method",
        choices=["pyrodigal", "prodigal"],
        default="pyrodigal",
        help="Gene calling method [Default is pyrodigal].",
    )
    parser.add_argument(
        "-l",
        "--locus-tag-length",
        type=int,
        default=3,
        help="Length of locus tags [Default is 3].",
    )
    parser.add_argument(
        "-m", "--meta-mode", action="store_true", help="Use meta mode for gene calling."
    )
    parser.add_argument(
        "-rlt",
        "--rename-locus-tags",
        action="store_true",
        help="Rename locus tags in GenBank files.",
    )
    parser.add_argument(
        "-mm",
        "--max-memory",
        type=int,
        default=32,
        help="Memory limit in GB [Default is 32].",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        help="Clean up temporary files after test completion.",
    )


def run_bofasa_analysis(args: argparse.Namespace) -> None:
    """Run the main bofasa analysis."""
    # Import here to avoid circular imports and heavy dependencies
    import os
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
        create_protein_alignments,
        create_profile_hmms_and_consensus_seqs,
        concatenate_consensus_alignment,
        create_near_scc_resolved_domain_protein_alignments,
        create_final_report,
        create_final_visual,
    )

    # Check OrthoFinder setup early in the workflow
    print("Checking OrthoFinder setup...")
    orthofinder_ok, orthofinder_message = check_orthofinder_setup()
    if not orthofinder_ok:
        print(f"Error: {orthofinder_message}")
        print("OrthoFinder is essential for the run command and must be properly configured.")
        sys.exit(1)
    print("✓ OrthoFinder is properly set up")

    # Check if output directory already exists
    if os.path.exists(args.output_dir) and os.listdir(args.output_dir):
        print(f"\nWARNING: The output directory '{args.output_dir}' already exists and contains files!")
        print("This process may overwrite existing files.")
        
        while True:
            response = input("Do you want to continue and potentially overwrite existing files? (y/N): ").strip().lower()
            if response in ['y', 'yes']:
                print("Continuing with existing output directory...")
                break
            elif response in ['n', 'no', '']:
                print("Exiting. Please choose a different output directory or remove existing files.")
                sys.exit(1)
            else:
                print("Please enter 'y' for yes or 'n' for no.")

    print(args.threads)

    # Set recursion limit
    max_recursion_depth = getattr(args, 'max_recursion_depth', 5000)
    sys.setrecursionlimit(max_recursion_depth)

    # Set memory limit
    if args.max_memory:
        memory_limit(args.max_memory)

    # Create output directory first
    os.makedirs(args.output_dir, exist_ok=True)

    # Create logger
    log_file = os.path.join(args.output_dir, "bofasa.log")
    logger = create_logger_object(log_file)
    
    # Check if logger was created successfully
    if logger is None:
        print("Warning: Failed to create logger. Continuing without logging.")
        logger = None

    try:
        if logger:
            logger.info("Starting bofasa analysis")
            logger.info(f"Input directory: {args.bofasa_prep_dir}")
            logger.info(f"Output directory: {args.output_dir}")
            logger.info(f"Threads: {args.threads}")
            logger.info(f"Memory limit: {args.max_memory}GB")
        else:
            print("Starting bofasa analysis")
            print(f"Input directory: {args.bofasa_prep_dir}")
            print(f"Output directory: {args.output_dir}")
            print(f"Threads: {args.threads}")
            print(f"Memory limit: {args.max_memory}GB")

        # Log parameters
        if logger:
            logger.info("Saving parameters for future records.")
            parameters_file = os.path.join(args.output_dir, "Parameter_Inputs.txt")
            parameter_values = [
                args.bofasa_prep_dir, args.output_dir, getattr(args, 'surrounding_bp', 10000),
                getattr(args, 'mcl_inflation', 1.2), getattr(args, 'ultra_sens', False),
                getattr(args, 'fixation_index_cutoff', 0.25), getattr(args, 'dog_jaccard', 0.25),
                getattr(args, 'skip_phylo_refine', False), getattr(args, 'rooting_seeds', 1),
                getattr(args, 'skip_merge_back', False), getattr(args, 'quality_alignments', False),
                getattr(args, 'og_consensus', False), getattr(args, 'core_genome', False),
                getattr(args, 'near_scc_prop', 0.95), args.threads, args.max_memory
            ]
            parameter_names = [
                "bofasa_prep directory", "Output directory", "Surrounding BP for syntenic conservation assessment",
                "MCL inflation", "Use ultra-sensitivity mode for DIAMOND searching in OrthoFinder?",
                "Fixation index cutoff for domain ortholog group re-merging following phylogenetic splitting",
                "Jaccard index cutoff for domain ortholog groups shared between full protein pairs",
                "Skip phylogenetic refinement of domain ortholog groups?", "Maximum number of nodes to try for rooting",
                "Skip fixation index-based re-mergining of split ortholog partitions?",
                "High-quality alignment method?", "Determine consensus sequences for ortholog groups?",
                "Create core genome alignment(s) for phylogenomics?", "Near-SCC conservation proportion?",
                "Number of threads/cores", "Maximum memory in GB"
            ]
            log_parameters_to_file(parameters_file, parameter_names, parameter_values)
            logger.info("Done saving parameters!")

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
                print(msg)
        except Exception as e:
            if logger:
                logger.error('Issue with validating directory with input genomes exists.')
            else:
                print('Issue with validating directory with input genomes exists.')
            sys.exit(1)

        # Create checkpoint and final results directories
        checkdir = os.path.join(args.output_dir, "Checkpoint_Files/")
        findir = os.path.join(args.output_dir, "Final_Results/")
        setup_ready_directory([checkdir, findir])

        # Step 1: Run OrthoFinder
        step1_checkpoint_file = os.path.join(checkdir, "Step1.txt")
        msg = '\n--------------------\nStep 1\n--------------------\nRunning OrthoFinder.'
        if logger:
            logger.info(msg)
        else:
            print(msg)
        
        orthofinder_input_dir = os.path.join(args.bofasa_prep_dir, "Domain_and_Interdomain_FASTAs/")
        orthofinder_results_dir = os.path.join(args.output_dir, "OrthoFinder_Results/")
        
        if not os.path.isfile(step1_checkpoint_file):
            orthofinder_cmd = [
                'orthofinder', '-f', orthofinder_input_dir, '-o', orthofinder_results_dir,
                '-t', str(args.threads), '-og', '-I', str(getattr(args, 'mcl_inflation', 1.2))
            ]
            if getattr(args, 'ultra_sens', False):
                orthofinder_cmd += ['-S', 'diamond_ultra_sens']
            
            try:
                subprocess.run(orthofinder_cmd, check=True)
                if logger:
                    logger.info(f'Successfully ran: {" ".join(orthofinder_cmd)}')
            except subprocess.CalledProcessError as e:
                if logger:
                    logger.error(f'Had an issue running: {" ".join(orthofinder_cmd)}')
                    logger.error(str(e))
                else:
                    print(f'Had an issue running: {" ".join(orthofinder_cmd)}')
                sys.exit(1)
            
            with open(step1_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 1 checkpoint found - skipping OrthoFinder execution"
            if logger:
                logger.info(msg)
            else:
                print(msg)

        # Find OrthoFinder results files
        orthofinder_subdirs = [d for d in os.listdir(orthofinder_results_dir) 
                              if os.path.isdir(os.path.join(orthofinder_results_dir, d))]
        if not orthofinder_subdirs:
            if logger:
                logger.error('Could not find OrthoFinder results directory.')
            else:
                print('Could not find OrthoFinder results directory.')
            sys.exit(1)
        
        orthofinder_subdir = orthofinder_subdirs[0]
        orthofinder_base = os.path.join(orthofinder_results_dir, orthofinder_subdir)
        
        orthofinder_tsv_file = os.path.join(orthofinder_base, "Orthogroups/Orthogroups.tsv")
        orthofinder_tsv_singletons_file = os.path.join(orthofinder_base, "Orthogroups/Orthogroups_UnassignedGenes.tsv")
        orthofinder_graph_file = os.path.join(orthofinder_base, "WorkingDirectory/OrthoFinder_graph.txt")
        orthofinder_seqid_file = os.path.join(orthofinder_base, "WorkingDirectory/SequenceIDs.txt")
        orthofinder_fasta_dir = os.path.join(orthofinder_base, "Orthogroup_Sequences/")
        
        try:
            assert os.path.isfile(orthofinder_tsv_file) and os.path.isfile(orthofinder_tsv_singletons_file)
            assert os.path.isfile(orthofinder_graph_file) and os.path.isfile(orthofinder_seqid_file)
        except AssertionError:
            if logger:
                logger.error('Could not validate OrthoFinder results exist. Perhaps investigate the OrthoFinder log files for further logging information.')
            else:
                print('Could not validate Orthofinder results exist. Perhaps investigate the OrthoFinder log files for further logging information.')
            sys.exit(1)

        # Step 2: Process OrthoFinder results and determine more coarse domain ortholog groups
        step2_checkpoint_file = os.path.join(checkdir, "Step2.txt")
        msg = '\n--------------------\nStep 2\n--------------------\nDetermine refined domain-resolution ortholog groups.'
        if logger:
            logger.info(msg)
        else:
            print(msg)
        
        resdog_dir = os.path.join(args.output_dir, "Resolve_Domain_Ortholog_Groups/")
        resulting_dogs_file = os.path.join(findir, "Domain_Ortholog_Groups.tsv")
        
        if not os.path.isfile(step2_checkpoint_file):
            if getattr(args, 'skip_phylo_refine', False):
                combine_orthofinder_results(
                    orthofinder_tsv_file, orthofinder_tsv_singletons_file, 
                    resulting_dogs_file, logger
                )
            else:
                setup_ready_directory([resdog_dir])
                resolve_orthogroups_using_phylogenetics(
                    orthofinder_fasta_dir, orthofinder_tsv_file, orthofinder_tsv_singletons_file,
                    resdog_dir, resulting_dogs_file, logger,
                    skip_merge_back_flag=getattr(args, 'skip_merge_back', False),
                    rooting_seeds=getattr(args, 'rooting_seeds', 1),
                    fixation_index_cutoff=getattr(args, 'fixation_index_cutoff', 0.25),
                    threads=args.threads
                )
            
            with open(step2_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 2 checkpoint found - skipping domain ortholog group resolution"
            if logger:
                logger.info(msg)
            else:
                print(msg)

        # Step 3: Determine protein-resolution ortholog groups from domain-resolution ortholog groups
        step3_checkpoint_file = os.path.join(checkdir, "Step3.txt")
        msg = '\n--------------------\nStep 3\n--------------------\nDetermine protein-resolution ortholog groups.'
        if logger:
            logger.info(msg)
        else:
            print(msg)
        
        protein_cluster_dir = os.path.join(args.output_dir, "Protein_Clustering/")
        resulting_ogs_file = os.path.join(findir, "Protein_Ortholog_Groups.tsv")
        
        if not os.path.isfile(step3_checkpoint_file):
            setup_ready_directory([protein_cluster_dir])
            determine_protein_orthogroups(
                resulting_dogs_file, protein_cluster_dir, resulting_ogs_file, 
                logger, dj=getattr(args, 'dog_jaccard', 0.25), threads=args.threads
            )
            
            with open(step3_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 3 checkpoint found - skipping protein ortholog group determination"
            if logger:
                logger.info(msg)
            else:
                print(msg)

        # Step 4: Determine surrounding contexts for each ortholog group
        step4_checkpoint_file = os.path.join(checkdir, "Step4.txt")
        msg = '\n--------------------\nStep 4\n--------------------\nAssessing neighborhoods for each ortholog group.'
        if logger:
            logger.info(msg)
        else:
            print(msg)
        
        surround_info_dir = os.path.join(args.output_dir, "Surrounding_Context_Information/")
        og_context_info_file = os.path.join(surround_info_dir, "OrthoGroup_Contexts.tsv")
        
        if not os.path.isfile(step4_checkpoint_file):
            setup_ready_directory([surround_info_dir])
            determine_ortholog_group_contexts(
                args.bofasa_prep_dir, og_context_info_file, surround_info_dir, 
                resulting_ogs_file, logger, 
                surrounding_bp=getattr(args, 'surrounding_bp', 10000), 
                threads=args.threads
            )
            
            with open(step4_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 4 checkpoint found - skipping ortholog group context determination"
            if logger:
                logger.info(msg)
            else:
                print(msg)

        # Step 5: Create final summarizations
        step5_checkpoint_file = os.path.join(checkdir, "Step5.txt")
        msg = '\n--------------------\nStep 5\n--------------------\nCreating final consolidated report.'
        if logger:
            logger.info(msg)
        else:
            print(msg)
        
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
                print(msg)

        # Step 6: Create final visuals
        step6_checkpoint_file = os.path.join(checkdir, "Step6.txt")
        msg = '\n--------------------\nStep 6\n--------------------\nCreating final summary visualization of conservation vs. neighborhood syntenic conservation.'
        if logger:
            logger.info(msg)
        else:
            print(msg)
        
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
                print(msg)

        # Step 7: (Optional) Create profile-HMM database
        if getattr(args, 'og_consensus', False):
            step7_checkpoint_file = os.path.join(checkdir, "Step7.txt")
            msg = '\n--------------------\nStep 7\n--------------------\nConstructing ortholog group profile HMMs and consensus sequences.'
            if logger:
                logger.info(msg)
            else:
                print(msg)
            
            og_seqs_dir = os.path.join(args.output_dir, "Ortholog_Group_Sequences/")
            og_algn_dir = os.path.join(args.output_dir, "Ortholog_Group_Alignments/")
            og_hmms_dir = os.path.join(args.output_dir, "Ortholog_Group_HMMs/")
            og_cons_dir = os.path.join(args.output_dir, "Ortholog_Group_Consensus_Sequences/")
            concatenated_consensus_seqs_file = os.path.join(findir, "Orthogroup_Consensus_Sequences.faa")
            
            if not os.path.isfile(step7_checkpoint_file):
                setup_ready_directory([og_seqs_dir, og_algn_dir, og_hmms_dir, og_cons_dir])
                create_protein_alignments(
                    args.bofasa_prep_dir, resulting_ogs_file, og_seqs_dir, og_algn_dir,
                    logger, threads=args.threads
                )
                create_profile_hmms_and_consensus_seqs(
                    og_algn_dir, og_hmms_dir, og_cons_dir, logger, threads=args.threads
                )
                concatenate_consensus_alignment(
                    og_cons_dir, concatenated_consensus_seqs_file, logger
                )
                
                with open(step7_checkpoint_file, 'w') as f:
                    f.write("DONE")
            else:
                msg = "Step 7 checkpoint found - skipping profile HMM and consensus sequence creation"
                if logger:
                    logger.info(msg)
                else:
                    print(msg)

        # Step 8: (Optional) Create core genome alignment
        if getattr(args, 'core_genome', False):
            step8_checkpoint_file = os.path.join(checkdir, "Step8.txt")
            msg = '\n--------------------\nStep 8\n--------------------\nCreating core genome alignment from (near) single-copy-core domain ortholog groups.'
            if logger:
                logger.info(msg)
            else:
                print(msg)
            
            dogs_seqs_dir = os.path.join(args.output_dir, "NearSCC_Coarse_Domain_Ortholog_Group_Sequences/")
            dogs_algn_dir = os.path.join(args.output_dir, "NearSCC_Coarse_Domain_Ortholog_Group_Alignments/")
            dogs_trim_dir = os.path.join(findir, "NearSCC_Coarse_Domain_Ortholog_Group_Trimmed_Alignments/")
            merged_core_genome_file = os.path.join(findir, "Core_Genome_Alignment.faa")
            
            if not os.path.isfile(step8_checkpoint_file):
                setup_ready_directory([dogs_seqs_dir, dogs_algn_dir, dogs_trim_dir])
                create_near_scc_resolved_domain_protein_alignments(
                    args.bofasa_prep_dir, resulting_dogs_file, dogs_seqs_dir,
                    dogs_algn_dir, dogs_trim_dir, merged_core_genome_file,
                    logger, near_scc_prop=getattr(args, 'near_scc_prop', 0.95), 
                    threads=args.threads
                )
                
                with open(step8_checkpoint_file, 'w') as f:
                    f.write("DONE")
            else:
                msg = "Step 8 checkpoint found - skipping core genome alignment creation"
                if logger:
                    logger.info(msg)
                else:
                    print(msg)

        if logger:
            logger.info("bofasa analysis completed successfully")
        else:
            print("bofasa analysis completed successfully")

    except Exception as e:
        if logger:
            logger.error(f"Error during bofasa analysis: {str(e)}")
        else:
            print(f"Error during bofasa analysis: {str(e)}")
        raise
    finally:
        close_logger_object(logger)


def run_bofasa_prep(args: argparse.Namespace) -> None:
    """Run the bofasa preparation step."""
    # Import here to avoid circular imports and heavy dependencies
    import os
    import shutil
    from .utils import (
        create_logger_object, close_logger_object, memory_limit, 
        check_genomad_setup, setup_ready_directory, create_locus_tag_options,
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
        print("Error: Either input genomes (-i/--input-genomes) or annotation directories (-a/--annotation-dirs) must be provided.")
        sys.exit(1)
    
    if input_genomes and annotation_dirs:
        print("Warning: Both input genomes and annotation directories provided. Processing both types of input.")

    # Check if output directory already exists
    if os.path.exists(args.output_dir) and os.listdir(args.output_dir):
        print(f"\nWARNING: The output directory '{args.output_dir}' already exists and contains files!")
        print("This process may overwrite existing files.")
        
        while True:
            response = input("Do you want to continue and potentially overwrite existing files? (y/N): ").strip().lower()
            if response in ['y', 'yes']:
                print("Continuing with existing output directory...")
                break
            elif response in ['n', 'no', '']:
                print("Exiting. Please choose a different output directory or remove existing files.")
                sys.exit(1)
            else:
                print("Please enter 'y' for yes or 'n' for no.")

    # Check genomad setup early in the workflow
    print("Checking genomad setup...")
    genomad_ok, genomad_message = check_genomad_setup()
    if not genomad_ok:
        print(f"Warning: {genomad_message}")
        print("genomad is required for phage/plasmid annotation in the analysis step.")
        print("You can continue with prep, but the analysis step may fail.")
        response = input("Do you want to continue with prep anyway? (y/N): ")
        if response.lower() not in ['y', 'yes']:
            print("Exiting. Please run 'bofasa setup' first to configure genomad.")
            sys.exit(1)
        print("Continuing with prep despite genomad issues...")
    else:
        print("✓ genomad is properly set up")

    # Set memory limit
    if args.max_memory:
        memory_limit(args.max_memory)

    # Create output directory first
    os.makedirs(args.output_dir, exist_ok=True)

    # Create checkpoint directory
    checkdir = os.path.join(args.output_dir, "Checkpoint_Files/")
    setup_ready_directory([checkdir])

    # Create logger
    log_file = os.path.join(args.output_dir, "bofasa_prep.log")
    logger = create_logger_object(log_file)
    
    # Check if logger was created successfully
    if logger is None:
        print("Warning: Failed to create logger. Continuing without logging.")
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
            print("Starting bofasa preparation")
            if input_genomes:
                print(f"Input genomes: {input_genomes}")
            if annotation_dirs:
                print(f"Annotation directories: {annotation_dirs}")
            print(f"Output directory: {args.output_dir}")
            print(f"Threads: {args.threads}")
            print(f"Gene calling method: {args.gene_calling_method}")
            print(f"genomad setup status: {genomad_ok} - {genomad_message}")

        # Log parameters
        if logger:
            logger.info("Saving parameters for future records.")
            parameters_file = os.path.join(args.output_dir, "Parameter_Inputs.txt")
            parameter_values = [
                input_genomes,
                annotation_dirs,
                args.output_dir,
                getattr(args, 'rename_locus_tags', False),
                getattr(args, 'meta_mode', False),
                getattr(args, 'run_genomad', False),
                getattr(args, 'ignore_upperbound_limit', False),
                args.gene_calling_method,
                getattr(args, 'min_length', 20),
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
                "Ignore Upperbound Genome Limit?",
                "Gene Calling Method",
                "Minimum Length of Domain/Inter-Domain Unit",
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
            print(msg)
        
        if not os.path.isfile(step1_checkpoint_file):
            gp_dir = os.path.join(args.output_dir, "Genome_Processing/")
            faa_dir = os.path.join(gp_dir, "Proteomes/")
            bed_dir = os.path.join(gp_dir, "BEDs/")

            setup_ready_directory([gp_dir, faa_dir, bed_dir])

        # Initialize sample mappings
        sample_wgs = {}
        sample_proteomes = {}
        sample_beds = {}

        # Process input genomes if provided
        if input_genomes:
            # Determine file types and process accordingly
            fasta_files = []
            genbank_files = []

            for genome_file in input_genomes:
                if os.path.splitext(genome_file)[1].lower() in [".gbk", ".gb", ".genbank", ".gbff"]:
                    genbank_files.append(genome_file)
                elif os.path.splitext(genome_file)[1].lower() in [".fasta", ".fa", ".fna"]:
                    fasta_files.append(genome_file)
                else:
                    if logger:
                        logger.warning(f"{genome_file} is not a valid genome file. Skipping...")
                    else:
                        print(f"Warning: {genome_file} is not a valid genome file. Skipping...")

            # Process FASTA files with Prodigal
            if fasta_files:
                if logger:
                    logger.info(f"Processing {len(fasta_files)} FASTA files with {args.gene_calling_method}...")
                else:
                    print(f"Processing {len(fasta_files)} FASTA files with {args.gene_calling_method}...")
                
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
                    meta_mode=getattr(args, 'meta_mode', False),
                )

                # Move files to proper locations
                for sample in fasta_dict:
                    faa = os.path.join(gp_dir, f"{sample}.faa")
                    bed = os.path.join(gp_dir, f"{sample}.coords.bed")
                    try:
                        assert os.path.isfile(faa) and os.path.isfile(bed)
                        renamed_faa = os.path.join(faa_dir, f"{sample}.faa")
                        renamed_bed = os.path.join(bed_dir, f"{sample}.bed")
                        shutil.move(faa, renamed_faa)
                        shutil.move(bed, renamed_bed)
                        sample_proteomes[sample] = renamed_faa
                        sample_beds[sample] = renamed_bed
                    except Exception as e:
                        msg = f'Unable to validate proper processing for sample {sample}, skipping it'
                        if logger:
                            logger.info(msg)
                        else:
                            print(msg)
                        if sample in sample_wgs:
                            del sample_wgs[sample]

            # Process GenBank files
            if genbank_files:
                if logger:
                    logger.info(f"Processing {len(genbank_files)} GenBank files...")
                else:
                    print(f"Processing {len(genbank_files)} GenBank files...")
                
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
                    rename_locus_tags=getattr(args, 'rename_locus_tags', False),
                )

                # Move files to proper locations
                for sample in genbank_dict:
                    faa = os.path.join(gp_dir, f"{sample}.faa")
                    bed = os.path.join(gp_dir, f"{sample}.coords.bed")
                    fna = os.path.join(gp_dir, f"{sample}.fna")
                    try:
                        assert os.path.isfile(faa) and os.path.isfile(bed) and os.path.isfile(fna)
                        assert os.path.getsize(faa) > 0 and os.path.getsize(bed) > 0 and os.path.getsize(fna) > 0
                        renamed_faa = os.path.join(faa_dir, f"{sample}.faa")
                        renamed_bed = os.path.join(bed_dir, f"{sample}.bed")
                        shutil.move(faa, renamed_faa)
                        shutil.move(bed, renamed_bed)
                        sample_proteomes[sample] = renamed_faa
                        sample_beds[sample] = renamed_bed
                        sample_wgs[sample] = fna
                    except Exception as e:
                        msg = f'Unable to validate proper processing for sample {sample}, skipping it'
                        if logger:
                            logger.info(msg)
                        else:
                            print(msg)
                        if sample in sample_wgs:
                            del sample_wgs[sample]

        # Process annotation directories (Prokka/Bakta)
        if annotation_dirs:
            if logger:
                logger.info(f"Processing {len(annotation_dirs)} annotation directories...")
            else:
                print(f"Processing {len(annotation_dirs)} annotation directories...")
            
            # Process annotation directories and get sample mappings
            annotation_results = process_annotation_directories(
                annotation_dirs,
                gp_dir,
                logger,
                locus_tag_length=args.locus_tag_length,
                rename_locus_tags=getattr(args, 'rename_locus_tags', False),
                threads=args.threads,
            )
            
            # Update sample mappings
            sample_wgs.update(annotation_results.get('sample_wgs', {}))
            sample_proteomes.update(annotation_results.get('sample_proteomes', {}))
            sample_beds.update(annotation_results.get('sample_beds', {}))

        # Create Step 1 checkpoint only after all processing is complete
        if not os.path.isfile(step1_checkpoint_file):
            # Create Step 1 checkpoint
            with open(step1_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 1 checkpoint found - skipping genome processing"
            if logger:
                logger.info(msg)
            else:
                print(msg)
            
            # Load existing data from previous run
            gp_dir = os.path.join(args.output_dir, "Genome_Processing/")
            faa_dir = os.path.join(gp_dir, "Proteomes/")
            bed_dir = os.path.join(gp_dir, "BEDs/")
            
            # Reconstruct sample mappings from existing files
            sample_wgs = {}
            sample_proteomes = {}
            sample_beds = {}
            
            # Load from existing proteome files
            if os.path.exists(faa_dir):
                for faa_file in os.listdir(faa_dir):
                    if faa_file.endswith('.faa'):
                        sample_name = os.path.splitext(faa_file)[0]
                        sample_proteomes[sample_name] = os.path.join(faa_dir, faa_file)
            
            # Load from existing bed files
            if os.path.exists(bed_dir):
                for bed_file in os.listdir(bed_dir):
                    if bed_file.endswith('.bed'):
                        sample_name = os.path.splitext(bed_file)[0]
                        sample_beds[sample_name] = os.path.join(bed_dir, bed_file)
            
            # Load genome files (check both .fna files and original input files)
            for sample in sample_proteomes:
                fna_file = os.path.join(gp_dir, f"{sample}.fna")
                if os.path.exists(fna_file):
                    sample_wgs[sample] = fna_file
                else:
                    # If no .fna file, try to find original input file
                    for input_file in input_genomes:
                        if os.path.splitext(os.path.basename(input_file))[0] == sample:
                            sample_wgs[sample] = input_file
                            break

        msg = f'Found and successfully processed {len(sample_proteomes)} genomes, continuing ...'
        if logger:
            logger.info(msg)
        else:
            print(msg)

        # Check if we have enough genomes for analysis
        if len(sample_proteomes) < 4:
            error_msg = f"Error: bofasa prep requires at least 4 genomes, but only {len(sample_proteomes)} were successfully processed."
            if logger:
                logger.error(error_msg)
            else:
                print(error_msg)
            print("Please provide at least 4 valid genome files or annotation directories.")
            sys.exit(1)

        # Step 2a: Running geNomad (if requested and genome files are available)
        step2a_checkpoint_file = os.path.join(checkdir, "Step2a.txt")
        if getattr(args, 'run_genomad', False) and sample_wgs:
            if not os.path.isfile(step2a_checkpoint_file):
                msg = '\n--------------------\nStep 2a\n--------------------\nRunning geNomad for phage and plasmid identification.'
                if logger:
                    logger.info(msg)
                else:
                    print(msg)

                genomad_dir = os.path.join(args.output_dir, "geNomad_Annotations/")
                setup_ready_directory([genomad_dir])

                phage_protein_listing_file = os.path.join(args.output_dir, "Sample_Phage_Proteins.txt")
                plasmid_protein_listing_file = os.path.join(args.output_dir, "Sample_Plasmid_Proteins.txt")
                
                determine_phages_and_plasmids(
                    sample_wgs, sample_beds, genomad_dir, 
                    phage_protein_listing_file, plasmid_protein_listing_file, 
                    logger, threads=args.threads
                )
                
                # Create Step 2a checkpoint
                with open(step2a_checkpoint_file, 'w') as f:
                    f.write("DONE")
            else:
                msg = "Step 2a checkpoint found - skipping geNomad execution"
                if logger:
                    logger.info(msg)
                else:
                    print(msg)
        elif getattr(args, 'run_genomad', False) and not sample_wgs:
            msg = 'Warning: geNomad requested but no genome files available (only annotation directories provided). Skipping geNomad step.'
            if logger:
                logger.warning(msg)
            else:
                print(msg)

        # Step 2b: Annotate proteins with ISfinder databases
        step2b_checkpoint_file = os.path.join(checkdir, "Step2b.txt")
        if not os.path.isfile(step2b_checkpoint_file):
            msg = '\n--------------------\nStep 2b\n--------------------\nAnnotating IS elements using the ISFinder database.'
            if logger:
                logger.info(msg)
            else:
                print(msg)

            annot_dir = os.path.join(args.output_dir, "ISFinder_Annotations/")
            setup_ready_directory([annot_dir])

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
                print(msg)

        # Step 3: Annotate proteins with Pfam and split into domain/inter-domain units
        step3_checkpoint_file = os.path.join(checkdir, "Step3.txt")
        if not os.path.isfile(step3_checkpoint_file):
            msg = '\n--------------------\nStep 3\n--------------------\nAnnotating proteins with Pfam HMMs to split them into domain/inter-domain "chunks" or units.'
            if logger:
                logger.info(msg)
            else:
                print(msg)

            split_proteins_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_FASTAs/")
            domain_coords_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_Coordinates/")
            domain_coord_info_file = os.path.join(args.output_dir, "Sample_Domain_and_InterDomain_Information.txt")
            setup_ready_directory([split_proteins_dir, domain_coords_dir])
            
            sample_ccds_proteomes = annotate_and_split_proteins_using_pfam(
                sample_proteomes, split_proteins_dir, domain_coords_dir, 
                domain_coord_info_file, logger, 
                minimal_length=getattr(args, 'min_length', 20), 
                threads=args.threads, 
                skip_domain_splitting=getattr(args, 'skip_domain_splitting', False)
            )
            
            # Create Step 3 checkpoint
            with open(step3_checkpoint_file, 'w') as f:
                f.write("DONE")
        else:
            msg = "Step 3 checkpoint found - skipping Pfam annotation and domain splitting"
            if logger:
                logger.info(msg)
            else:
                print(msg)
            
            # Load existing data from previous run
            split_proteins_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_FASTAs/")
            domain_coords_dir = os.path.join(args.output_dir, "Domain_and_Interdomain_Coordinates/")
            domain_coord_info_file = os.path.join(args.output_dir, "Sample_Domain_and_InterDomain_Information.txt")
            
            # Reconstruct sample_ccds_proteomes from existing files
            sample_ccds_proteomes = {}
            if os.path.exists(split_proteins_dir):
                for ccds_file in os.listdir(split_proteins_dir):
                    if ccds_file.endswith('.faa'):
                        sample_name = os.path.splitext(ccds_file)[0]
                        sample_ccds_proteomes[sample_name] = os.path.join(split_proteins_dir, ccds_file)

        # Step 4: Finalize setting up for final bofasa analysis
        step4_checkpoint_file = os.path.join(checkdir, "Step4.txt")
        if not os.path.isfile(step4_checkpoint_file):
            msg = '\n--------------------\nStep 4\n--------------------\nNearly there, putting together final files.'
            if logger:
                logger.info(msg)
            else:
                print(msg)

            # Create info file
            genome_info_listing_file = os.path.join(args.output_dir, "Info_on_Input_Genome_Files.txt")
            with open(genome_info_listing_file, "w") as f:
                header = ['sample', 'ccds_faa', 'proteome', 'bed_coord', 'genome']
                f.write('\t'.join(header) + '\n')
                
                for sample in sample_wgs:
                    fna = sample_wgs[sample]
                    faa = sample_proteomes[sample]
                    ccds_faa = sample_ccds_proteomes[sample]
                    bed = sample_beds[sample]
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
                print(msg)

        if logger:
            logger.info("bofasa preparation completed successfully")
        else:
            print("bofasa preparation completed successfully")

    except Exception as e:
        if logger:
            logger.error(f"Error during bofasa preparation: {str(e)}")
        else:
            print(f"Error during bofasa preparation: {str(e)}")
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
        print("Error: BOFASA_DB_PATH environment variable is not set.")
        print(
            "Please set BOFASA_DB_PATH to the directory where you want to store annotation databases."
        )
        sys.exit(1)

    output_dir = output_dir.strip()
    if not output_dir:
        print("Error: BOFASA_DB_PATH environment variable is empty.")
        sys.exit(1)

    try:
        setup_annotation_databases(output_dir, args.threads, args.force)
    except Exception as e:
        print(f"Error during setup annotation databases: {str(e)}")
        raise


def run_bofasa_test(args):
    """Run the bofasa test workflow."""
    # Import here to avoid circular imports and heavy dependencies
    import os
    import tempfile
    import shutil
    from .utils import create_logger_object, close_logger_object, memory_limit
    
    # Set memory limit
    if args.max_memory:
        memory_limit(args.max_memory)

    # Create output directory first
    os.makedirs(args.output_dir, exist_ok=True)

    # Create logger
    log_file = os.path.join(args.output_dir, "bofasa_test.log")
    logger = create_logger_object(log_file)
    
    # Check if logger was created successfully
    if logger is None:
        print("Warning: Failed to create logger. Continuing without logging.")
        logger = None

    try:
        if logger:
            logger.info("Starting bofasa test workflow")
            logger.info(f"Output directory: {args.output_dir}")
            logger.info(f"Threads: {args.threads}")
            logger.info(f"Gene calling method: {args.gene_calling_method}")
        else:
            print("Starting bofasa test workflow")
            print(f"Output directory: {args.output_dir}")
            print(f"Threads: {args.threads}")
            print(f"Gene calling method: {args.gene_calling_method}")

        # Create temporary directory for test data
        with tempfile.TemporaryDirectory() as temp_dir:
            if logger:
                logger.info(f"Created temporary directory: {temp_dir}")
            else:
                print(f"Created temporary directory: {temp_dir}")
            
            # Generate test genomes
            test_genomes = create_test_genomes(temp_dir, logger)
            if logger:
                logger.info(f"Created {len(test_genomes)} test genomes")
            else:
                print(f"Created {len(test_genomes)} test genomes")
            
            # Step 1: Run bofasa prep
            prep_dir = os.path.join(args.output_dir, "prep_output")
            if logger:
                logger.info("Running bofasa prep...")
            else:
                print("Running bofasa prep...")
            
            # Create prep args
            prep_args = type('Args', (), {
                'input_genomes': test_genomes,
                'output_dir': prep_dir,
                'threads': args.threads,
                'gene_calling_method': args.gene_calling_method,
                'locus_tag_length': args.locus_tag_length,
                'meta_mode': args.meta_mode,
                'rename_locus_tags': args.rename_locus_tags,
                'max_memory': args.max_memory
            })()
            
            run_bofasa_prep(prep_args)
            if logger:
                logger.info("bofasa prep completed successfully")
            else:
                print("bofasa prep completed successfully")
            
            # Step 2: Run bofasa run
            run_dir = os.path.join(args.output_dir, "run_output")
            if logger:
                logger.info("Running bofasa run...")
            else:
                print("Running bofasa run...")
            
            # Create run args
            run_args = type('Args', (), {
                'bofasa_prep_dir': prep_dir,
                'output_dir': run_dir,
                'surrounding_bp': 10000,
                'og_consensus': True,
                'core_genome': True,
                'dog_jaccard': 0.25,
                'fixation_index_cutoff': 0.25,
                'skip_merge_back': False,
                'skip_phylo_refine': False,
                'rooting_seeds': 1,
                'refine': False,
                'n_refinements': 100,
                'ultra_sens': False,
                'mcl_inflation': 1.2,
                'near_scc_prop': 0.80,
                'threads': args.threads,
                'max_recursion_depth': 5000,
                'max_memory': args.max_memory
            })()
            
            run_bofasa_analysis(run_args)
            logger.info("bofasa run completed successfully")
            
            # Cleanup if requested
            if args.cleanup:
                logger.info("Cleaning up temporary files...")
                if os.path.exists(prep_dir):
                    shutil.rmtree(prep_dir)
                logger.info("Cleanup completed")
            
            logger.info("bofasa test workflow completed successfully")
            logger.info(f"Results available in: {args.output_dir}")

    except Exception as e:
        logger.error(f"Error during bofasa test workflow: {str(e)}")
        raise
    finally:
        close_logger_object(logger)


def create_test_genomes(temp_dir, logger):
    """Create test genome files for testing."""
    import random
    import os
    
    test_genomes = []
    
    # Create 3 test genomes with different characteristics
    for i in range(3):
        genome_file = os.path.join(temp_dir, f"test_genome_{i+1}.fasta")
        
        # Generate random DNA sequence (simplified for testing)
        # In a real implementation, you might want more realistic sequences
        dna_length = random.randint(100000, 200000)  # 100-200 kb for faster testing
        dna_sequence = ''.join(random.choices(['A', 'T', 'G', 'C'], k=dna_length))
        
        # Write FASTA file directly
        with open(genome_file, 'w') as f:
            f.write(f">test_genome_{i+1} Test genome {i+1} for bofasa testing\n")
            # Write sequence in chunks of 80 characters (standard FASTA format)
            for j in range(0, len(dna_sequence), 80):
                f.write(dna_sequence[j:j+80] + '\n')
        
        test_genomes.append(genome_file)
        logger.info(f"Created test genome: {genome_file}")
    
    return test_genomes


def main() -> None:
    """Main entry point for the bofasa command."""
    # Create parser and parse arguments
    parser = create_main_parser()

    # Handle version flag
    if len(sys.argv) > 1 and ("-v" in set(sys.argv) or "--version" in set(sys.argv)):
        print(config.get_version())
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
    elif args.command == "test":
        run_bofasa_test(args)
    else:
        print(f"Unknown command: {args.command}")
        print_bofasa_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
