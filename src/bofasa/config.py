"""
Configuration settings for bofasa.

This module contains global configuration variables and settings used throughout
the bofasa package.
"""

import random
from typing import Set, Dict, Any, Optional
from .utils import get_version

# Set random seed for reproducibility
random.seed(12345)

# Version information - dynamically retrieved from package metadata
VERSION: str = get_version()

# Global variables for orthology analysis
single_copy_dogs: Set[str] = set([])
largely_idr_dogs: Set[str] = set([])
protein_dogs: Dict[str, Dict[str, int]] = {}
dog_conservation: Dict[str, Any] = {}
tree_obj: Optional[Any] = None

# Default parameters
DEFAULT_DOG_JACCARD: float = 0.5
DEFAULT_FIXATION_INDEX_CUTOFF: float = 0.1
DEFAULT_ROOTING_SEEDS: int = 1
DEFAULT_MCL_INFLATION: float = 1.2
DEFAULT_NEAR_SCC_PROP: float = 0.95
DEFAULT_THREADS: int = 4
DEFAULT_MAX_MEMORY: int = 16
DEFAULT_SURROUNDING_BP: int = 10000
DEFAULT_TRIMAL_OPTIONS: str = "-strict -keepseqs"
DEFAULT_MAX_RECURSION_DEPTH: int = 5000
DEFAULT_MAX_GENOMES: int = 200

# File extensions and patterns
FASTA_EXTENSIONS: list[str] = [".fasta", ".fa", ".fna", ".fas"]
GENBANK_EXTENSIONS: list[str] = [".gbk", ".gb", ".genbank"]
PROTEIN_EXTENSIONS: list[str] = [".faa", ".fasta", ".fa"]

# Supported software for comparison
SUPPORTED_SOFTWARE: Set[str] = {
    "bofasa",
    "panx",
    "sonicparanoid",
    "pirate",
    "panaroo",
    "scarap",
    "orthofinder_hog",
    "orthofinder_og",
}

# Gene calling methods
GENE_CALLING_METHODS: list[str] = ["pyrodigal", "prodigal"]

# Alignment methods
ALIGNMENT_METHODS: list[str] = ["muscle"]

# Default locus tag length
DEFAULT_LOCUS_TAG_LENGTH: int = 3

# Default minimum length of domain/inter-domain unit to consider
DEFAULT_MIN_LENGTH: int = 20

# Memory limits (in GB)
MIN_MEMORY_LIMIT: int = 4
MAX_MEMORY_LIMIT: int = 256

# File size limits (in MB)
MAX_FILE_SIZE: int = 1024 * 1024 * 1024  # 1GB

# Logging configuration
LOG_FORMAT: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_LEVEL: str = "INFO"

# Subcommand colors for help text
SUBCMD_SETUP_COLOR: str =  "\033[93m\033[1m"  # Yellow
SUBCMD_PREP_COLOR: str = "\033[95m\033[1m"    # Magenta
SUBCMD_RUN_COLOR: str = "\033[95m\033[1m"    # Magenta
SUBCMD_END_COLOR: str = "\033[0m"            # Reset

def reset_global_variables() -> None:
    """Reset global variables to their initial state."""
    global single_copy_dogs, largely_idr_dogs, protein_dogs, dog_conservation, tree_obj

    single_copy_dogs = set([])
    largely_idr_dogs = set([])
    protein_dogs = {}
    dog_conservation = {}
    tree_obj = None
