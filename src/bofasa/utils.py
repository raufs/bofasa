"""
Utility functions for bofasa.

This module contains general utility functions for logging, command execution,
file operations, and other helper functions.
"""

import os
import sys
import subprocess
import resource
import logging
import traceback
import shutil
import itertools
from typing import Any, List, Optional, Set, Tuple
import pandas as pd
from . import config
import importlib.metadata

try:
    package_name = "bofasa"
    package_version = str(importlib.metadata.version(package_name))
except importlib.metadata.PackageNotFoundError:
    package_version = "NA"

def load_table_in_pandas_dataframe(
    input_file: str, 
    numeric_columns: List[str], 
    cut_last_columns: Optional[int] = None
) -> pd.DataFrame:
    """
    Load a table into a pandas DataFrame with proper data types.

    Parameters:
    -----------
    input_file : str
        Path to input file
    numeric_columns : List[str]
        List of column names that should be numeric
    cut_last_columns : Optional[int], default=None
        Number of columns to cut from the end

    Returns:
    --------
    pd.DataFrame
        Loaded DataFrame with proper data types
    """
    try:
        df: pd.DataFrame = pd.read_csv(input_file, sep='\t')
        
        # Convert numeric columns
        for col in numeric_columns:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Cut last columns if specified
        if cut_last_columns is not None and cut_last_columns > 0:
            df = df.iloc[:, :-cut_last_columns]
        
        return df
        
    except Exception as e:
        raise ValueError(f"Error loading table from {input_file}: {str(e)}")


def _iter_progress(iterable, total=None, description: Optional[str] = None):
    """
    Iterate with a nice progress display using Rich's track when available,
    falling back to tqdm otherwise.
    """
    try:
        # Lightweight import so Rich remains optional at runtime
        from rich.progress import track as _rich_track

        return _rich_track(iterable, total=total, description=description)
    except Exception:
        # Fallback to tqdm if Rich isn't available
        import tqdm  # type: ignore

        # tqdm uses desc instead of description
        return tqdm.tqdm(iterable, total=total, desc=description)


class Colors:
    """ANSI color codes for terminal output."""
    BOLD = '\033[1m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    END = '\033[0m'
    UNDERLINE = '\033[4m'


def print_colored_help(help_text: str, command_name: Optional[str] = None) -> None:
    """
    Print colored help text with bofasa branding.

    This function formats and displays help text with ANSI color codes
    and bofasa branding for terminal output.

    Parameters:
    -----------
    help_text : str
        The help text to display
    command_name : str, optional
        Optional command name for specific styling

    Returns:
    --------
    None
        Prints colored help text to stdout
    """
    # Check dependencies for help commands
    if command_name:
        if command_name == "prep":
            # For prep help, check genomad but only warn
            genomad_ok: bool
            genomad_message: str
            genomad_ok, genomad_message = check_genomad_setup()
            if not genomad_ok:
                sys.stdout.write(f"\nWarning: {genomad_message}\n")
                sys.stdout.write("genomad is required for phage/plasmid annotation in the analysis step.\n")
        elif command_name == "run":
            # For run help, check OrthoFinder and error if missing
            orthofinder_ok: bool
            orthofinder_message: str
            orthofinder_ok, orthofinder_message = check_orthofinder_setup()
            if not orthofinder_ok:
                sys.stdout.write(orthofinder_message + '\n')
                sys.exit(1)
    
    # Check if we're in a terminal that supports colors
    if not sys.stdout.isatty():
        sys.stdout.write(help_text + '\n')
        return
    
    # Define color scheme
    title_color: str = Colors.CYAN + Colors.BOLD
    subtitle_color: str = Colors.YELLOW + Colors.BOLD
    section_color: str = Colors.GREEN + Colors.BOLD
    option_color: str = Colors.MAGENTA + Colors.BOLD
    description_color: str = Colors.WHITE
    arg_description_color: str = Colors.CYAN
    end_color: str = Colors.END
    
    # Add bofasa logo
    logo: str = f"""{title_color}
██████╗  ██████╗ ███████╗ █████╗ ███████╗ █████╗
██╔══██╗██╔═══██╗██╔════╝██╔══██╗██╔════╝██╔══██╗
██████╔╝██║   ██║█████╗  ███████║███████╗███████║
██╔══██╗██║   ██║██╔══╝  ██╔══██║╚════██║██╔══██║
██████╔╝╚██████╔╝██║     ██║  ██║███████║██║  ██║
╚═════╝  ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝{end_color}"""
    
    if command_name:
        command_header: str = f"\n{subtitle_color}Program: bofasa {command_name}{end_color}"
        sys.stdout.write(logo + command_header + '\n')
    else:
        sys.stdout.write(logo + '\n')
    
    # Process help text with colors
    lines: List[str] = help_text.split('\n')
    in_options: bool = False
    
    for line in lines:
        if line.startswith('usage:') or line.startswith('positional arguments:') or line.startswith('optional arguments:'):
            sys.stdout.write(f"{title_color}{line}{end_color}\n")
        elif line.startswith('  -') or line.startswith('  --'):
            parts: List[str] = line.split('  ', 2)
            if len(parts) >= 3:
                option_part: str = parts[1]
                description_part: str = parts[2]
                sys.stdout.write(f"  {option_color}{option_part}{end_color}  {description_color}{description_part}{end_color}\n")
            else:
                sys.stdout.write(f"{option_color}{line}{end_color}\n")
        elif line.strip() and not line.startswith('  '):
            sys.stdout.write(f"{section_color}{line}{end_color}\n")
        else:
            sys.stdout.write(f"{description_color}{line}{end_color}\n")


def print_bofasa_help() -> None:
    """
    Print the main bofasa help message.

    Returns:
    --------
    None
        Prints help text to stdout
    """
    help_text: str = f"""
BOFASA: Bacterial Ortholog Finder and Synteny Analysis

A comprehensive tool for identifying ortholog groups and analyzing 
syntenic conservation in bacterial genomes representing multiple
species.

Commands:
  {config.SUBCMD_SETUP_COLOR}setup{config.SUBCMD_END_COLOR}    Set up annotation databases
  {config.SUBCMD_PREP_COLOR}prep{config.SUBCMD_END_COLOR}     Prepare genomic data for analysis
  {config.SUBCMD_RUN_COLOR}run{config.SUBCMD_END_COLOR}      Run the main BOFASA analysis pipeline

For detailed help on any command, use: bofasa <command> --help

Examples:
  bofasa prep -i genome1.fna genome2.fna -o prep_output/
  bofasa run -i prep_output/ -o analysis_results/ -c 8

For more information, visit: https://github.com/raufs/bofasa
"""
    print_colored_help(help_text)


def run_cmd(
    cmd: List[str],
    logObject: Any,
    check_files: List[str] = [],
    check_directories: List[str] = [],
    stdout: Any = subprocess.DEVNULL,
    stderr: Any = subprocess.DEVNULL,
) -> None:
    """
    Run a command and handle logging.

    Parameters:
    -----------
    cmd : List[str]
        Command to run as a list of strings
    logObject : Any
        Logger object for recording output
    check_files : List[str], default=[]
        List of files to check for existence after command execution
    check_directories : List[str], default=[]
        List of directories to check for existence after command execution
    stdout : Any, default=subprocess.DEVNULL
        Standard output destination
    stderr : Any, default=subprocess.DEVNULL
        Standard error destination

    Returns:
    --------
    None
        Executes the command and logs results
    """
    try:
        cmd_str: str = ' '.join(cmd)
        if logObject:
            logObject.info(f"Running command: {cmd_str}")
        
        result: subprocess.CompletedProcess = subprocess.run(
            cmd, stdout=stdout, stderr=stderr, check=True, text=True
        )
        
        if logObject:
            logObject.info(f"Command completed successfully: {cmd_str}")
        
        # Check for required files and directories
        for file_path in check_files:
            if not os.path.isfile(file_path):
                raise FileNotFoundError(f"Required file not found: {file_path}")
        
        for dir_path in check_directories:
            if not os.path.isdir(dir_path):
                raise FileNotFoundError(f"Required directory not found: {dir_path}")
                
    except subprocess.CalledProcessError as e:
        error_msg: str = f"Command failed with return code {e.returncode}: {cmd_str}"
        if logObject:
            logObject.error(error_msg)
        else:
            sys.stderr.write(error_msg + '\n')
        raise
    except Exception as e:
        error_msg: str = f"Error running command {cmd_str}: {str(e)}"
        if logObject:
            logObject.error(error_msg)
        else:
            sys.stderr.write(error_msg + '\n')
        raise


def multi_process(input_data: List[Any]) -> None:
    """
    Execute a command in a multiprocessing context.

    Parameters:
    -----------
    input_data : List[Any]
        List containing command and logger object

    Returns:
    --------
    None
        Executes the command
    """
    try:
        cmd: List[str] = input_data[:-1]
        log_object: Any = input_data[-1]
        run_cmd(cmd, log_object)
    except Exception as e:
        if len(input_data) > 0 and hasattr(input_data[-1], 'error'):
            input_data[-1].error(f"Error in multi_process: {str(e)}")
        else:
            sys.stderr.write(f"Error in multi_process: {str(e)}\n")


def get_version() -> str:
    """
    Get the current version of bofasa package.

    This function retrieves the version information from the package metadata.

    Parameters:
    -----------
    None

    Returns:
    --------
    str
        Version string from package metadata, or "unknown" if not available
    """
    return package_version 


def single_linkage_cluster(
    pairs: List[Tuple[str, str, float]], all_lts: Set[str], paired_lts: Set[str]
) -> List[Set[str]]:
    """
    Perform single-linkage clustering on pairs.

    Solution adapted from a common approach to merge lists with shared elements.

    Parameters:
    - pairs: list of 2-tuples (or lists) representing links to cluster
    - all_lts: set of all labels
    - paired_lts: set of labels that appear in any pair

    Returns:
    - List of clustered groups (each group is a list of merged labels)
    """
    try:
        merged_groups: List[List[str]] = [list(p[:2]) for p in pairs]
        flat_items: Set[str] = set(itertools.chain.from_iterable(merged_groups))

        for item in flat_items:
            components = [grp for grp in merged_groups if item in grp]
            for comp in components:
                merged_groups.remove(comp)
            merged_groups.append(list(set(itertools.chain.from_iterable(components))))

        # Add singleton groups for any unpaired labels
        for label in all_lts:
            if label not in paired_lts:
                merged_groups.append([label])

        # Convert inner lists to sets for compatibility with callers
        return [set(group) for group in merged_groups]
    except Exception:
        msg = 'Issue running single linkage clustering!'
        sys.stderr.write(msg + '\n')
        sys.stderr.write(traceback.format_exc() + '\n')
        raise

def create_logger_object(log_file: str) -> Optional[logging.Logger]:
    """
    Create and configure a logger object for application logging.

    This function sets up a logger with both file and console handlers
    for comprehensive logging capabilities.

    Parameters:
    -----------
    log_file : str
        Path to the log file where messages will be written

    Returns:
    --------
    logging.Logger
        Configured logger object with file and console handlers

    Raises:
    --------
    Exception
        If logger creation fails
    """
    try:
        # Create logger
        logger: logging.Logger = logging.getLogger("bofasa")
        logger.setLevel(logging.INFO)

        # Create file handler
        file_handler: logging.FileHandler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)

        # Create console handler
        console_handler: logging.StreamHandler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)

        # Create formatter
        formatter: logging.Formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add handlers to logger
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        return logger

    except Exception as e:
        sys.stderr.write(f"Error creating logger: {str(e)}\n")
        sys.stderr.write(traceback.format_exc() + "\n")
        return None


def close_logger_object(logObject: Optional[logging.Logger]) -> None:
    """
    Close and cleanup logger object handlers.

    This function properly closes and removes all handlers from a logger
    object to prevent resource leaks.

    Parameters:
    -----------
    logObject : logging.Logger
        Logger object to close and cleanup

    Returns:
    --------
    None
        Removes and closes all handlers
    """
    try:
        if logObject:
            for handler in logObject.handlers[:]:
                logObject.removeHandler(handler)
                handler.close()
    except Exception as e:
        # Note: This function doesn't have access to log_object, so we keep stderr for this case
        sys.stderr.write(f"Error closing logger: {str(e)}\n")


def log_parameters_to_file(parameter_file: str, parameter_names: List[str], parameter_values: List[Any]) -> None:
    """
    Log parameters to a file for record keeping.

    Parameters:
    -----------
    parameter_file : str
        Path to the parameter log file
    parameter_names : List[str]
        List of parameter names
    parameter_values : List[Any]
        List of parameter values

    Returns:
    --------
    None
        Writes parameters to the specified file
    """
    try:
        with open(parameter_file, 'w') as f:
            f.write("Parameter\tValue\n")
            for name, value in zip(parameter_names, parameter_values):
                f.write(f"{name}\t{value}\n")
    except Exception as e:
        sys.stderr.write(f"Error logging parameters: {str(e)}\n")


def memory_limit(mem: int) -> None:
    """
    Set memory limit for the process.

    Parameters:
    -----------
    mem : int
        Memory limit in gigabytes

    Returns:
    --------
    None
        Sets the memory limit for the current process
    """
    try:
        # Import config here to avoid circular imports
        from .config import MIN_MEMORY_LIMIT, MAX_MEMORY_LIMIT
        
        # Validate memory limit is within reasonable bounds
        if mem < MIN_MEMORY_LIMIT:
            sys.stdout.write(f"Warning: Requested memory limit ({mem}GB) is below minimum ({MIN_MEMORY_LIMIT}GB)\n")
            sys.stdout.write(f"Using minimum limit of {MIN_MEMORY_LIMIT}GB\n")
            mem = MIN_MEMORY_LIMIT
        elif mem > MAX_MEMORY_LIMIT:
            sys.stdout.write(f"Warning: Requested memory limit ({mem}GB) is above maximum ({MAX_MEMORY_LIMIT}GB)\n")
            sys.stdout.write(f"Using maximum limit of {MAX_MEMORY_LIMIT}GB\n")
            mem = MAX_MEMORY_LIMIT
        
        max_virtual_memory: int = mem * 1000000000
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        
        # Handle RLIM_INFINITY (unlimited) values
        soft_gb = "unlimited" if soft == resource.RLIM_INFINITY else f"{soft // 1000000000}GB"
        hard_gb = "unlimited" if hard == resource.RLIM_INFINITY else f"{hard // 1000000000}GB"
        
        # Check if the requested limit exceeds the system's hard limit
        if hard != resource.RLIM_INFINITY and max_virtual_memory > hard:
            sys.stdout.write(f"Warning: Requested memory limit ({mem}GB) exceeds system maximum ({hard_gb})\n")
            sys.stdout.write(f"Using system maximum instead\n")
            max_virtual_memory = hard
        
        # Check if the requested limit is lower than the current soft limit
        if soft != resource.RLIM_INFINITY and max_virtual_memory < soft:
            sys.stdout.write(f"Reducing memory limit from {soft_gb} to {mem}GB\n")
        elif soft == resource.RLIM_INFINITY:
            sys.stdout.write(f"Setting memory limit from unlimited to {mem}GB\n")
        
        # When setting a lower limit, we need to set both soft and hard limits
        # to the same value to ensure the limit is enforced
        try:
            if soft == resource.RLIM_INFINITY or max_virtual_memory < soft:
                resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, max_virtual_memory))
            else:
                resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, hard))
            
            sys.stdout.write(f"Memory limit set to: {mem}GB\n")
        except ValueError as e:
            if "current limit exceeds maximum limit" in str(e):
                sys.stdout.write(f"Warning: Unable to set memory limit to {mem}GB on this system\n")
                sys.stdout.write("This may be due to system restrictions (common on macOS)\n")
                sys.stdout.write("Memory usage will not be limited by this process\n")
            else:
                raise
    except Exception as e:
        sys.stderr.write(f"Error setting memory limit: {str(e)}\n")


def setup_ready_directory(directories: List[str], overwrite_mode: str = "skip") -> None:
    """
    Create directories with optional overwrite behavior.
    
    Parameters:
    -----------
    directories : List[str]
        List of directory paths to create
    overwrite_mode : str, optional
        How to handle existing directories:
        - "skip": Skip existing directories (default)
        - "overwrite": Delete and recreate existing directories
        - "ask": Ask user for each existing directory

    Returns:
    --------
    None
        Creates all specified directories
    """
    try:
        for directory in directories:
            if os.path.exists(directory):
                if overwrite_mode == "skip":
                    continue
                elif overwrite_mode == "overwrite":
                    shutil.rmtree(directory)
                    os.makedirs(directory, exist_ok=True)
                elif overwrite_mode == "ask":
                    response = input(f"Directory '{directory}' already exists. Delete it? (y/n): ").strip().lower()
                    if response in ['y', 'yes']:
                        shutil.rmtree(directory)
                        os.makedirs(directory, exist_ok=True)
                    else:
                        sys.stdout.write(f"Skipping directory '{directory}'\n")
                else:
                    raise ValueError(f"Invalid overwrite_mode: {overwrite_mode}. Must be 'skip', 'overwrite', or 'ask'")
            else:
                os.makedirs(directory, exist_ok=True)
    except Exception as e:
        sys.stderr.write(f"Error creating directories: {str(e)}\n")
        raise


def create_locus_tag_options(locus_tag_length: int) -> List[str]:
    """
    Create locus tag options for different lengths.

    Parameters:
    -----------
    locus_tag_length : int
        Length of locus tags to generate

    Returns:
    --------
    List[str]
        List of locus tag options
    """
    try:
        alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        possible_locustags = sorted(list(set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))])))
        return possible_locustags
        
    except Exception as e:
        sys.stderr.write(f"Error creating locus tag options: {str(e)}\n")
        return []


def assess_job_intensity(faa_file: str) -> bool:
    """
    Assess the intensity of a job based on FASTA file characteristics.

    Parameters:
    -----------
    faa_file : str
        Path to FASTA file to assess

    Returns:
    --------
    bool
        True if job is considered heavy/intense, False otherwise
    """
    try:
        sequence_count: int = 0
        total_length: int = 0
        
        with open(faa_file, 'r') as handle:
            for line in handle:
                if line.startswith('>'):
                    sequence_count += 1
                else:
                    total_length += len(line.strip())
        
        # Consider job intense if:
        # >200 sequences or the average sequence length >2,000 sequences
        avg_length: float = total_length / sequence_count if sequence_count > 0 else 0
        
        return (sequence_count > 200 or avg_length > 2000)
                
    except Exception:
        # Default to simple job if assessment fails
        return False


def check_genomad_setup() -> Tuple[bool, str]:
    """
    Check if genomad is properly set up and accessible.

    Returns:
    --------
    Tuple[bool, str]
        (is_available, message) - True if genomad is available, False otherwise
    """
    try:
        # Check if genomad command is available
        result: subprocess.CompletedProcess = subprocess.run(
            ['genomad', '--version'], 
            capture_output=True, 
            text=True, 
            timeout=10
        )
        
        if result.returncode != 0:
            return False, "genomad command failed"
            
    except subprocess.TimeoutExpired:
        return False, "genomad command timed out"
    except FileNotFoundError:
        return False, "genomad command not found in PATH"
    except Exception as e:
        return False, f"Error checking genomad: {str(e)}"
    
    # Check if genomad database is set up using database_location_paths.txt
    bofasa_db_dir = os.getenv("BOFASA_DB_PATH", "").strip()
    if not bofasa_db_dir:
        return False, "BOFASA_DB_PATH environment variable not set"
    
    db_locations = os.path.join(bofasa_db_dir, 'database_location_paths.txt')
    if not os.path.isfile(db_locations):
        return False, f"Database locations file not found: {db_locations}"
    
    # Find genomad database path from the file
    genomad_db_path = None
    try:
        with open(db_locations) as odb:
            for line in odb:
                line = line.strip()
                ls = line.split('\t')
                if ls[0] == 'genomad':
                    genomad_db_path = ls[2]
                    break
    except Exception as e:
        return False, f"Error reading database locations file: {str(e)}"
    
    if not genomad_db_path:
        return False, "geNomad database path not found in database locations file"
    
    # Check if the database directory exists and contains required files
    if not os.path.exists(genomad_db_path) or not os.path.isdir(genomad_db_path):
        return False, f"geNomad database not found: {genomad_db_path}"
    
    # Check if the database directory contains the expected files
    expected_files = ["genomad_db", "version.txt"]
    has_required_files = any(
        os.path.exists(os.path.join(genomad_db_path, f)) for f in expected_files
    )
    
    if not has_required_files:
        return False, f"geNomad database incomplete: {genomad_db_path}"
    
    return True, "genomad is properly set up"


def check_orthofinder_setup() -> Tuple[bool, str]:
    """
    Check if OrthoFinder is properly set up and accessible.
    
    This function performs a comprehensive check of the OrthoFinder installation
    by running the help command and examining both the return code and output
    for any error messages. Provides detailed troubleshooting guidance if issues
    are detected.

    Returns:
    --------
    Tuple[bool, str]
        (is_available, message) - True if OrthoFinder is available, False otherwise
        The message contains detailed information about the status or error.
    """
    try:
        # Check if orthofinder command is available by trying to run it with help
        result: subprocess.CompletedProcess = subprocess.run(
            ['orthofinder', '-h'], 
            capture_output=True, 
            text=True, 
            timeout=10
        )
        
        # Check for error messages in the output (similar to lsaBGC-Pan approach)
        combined_output = (result.stdout + result.stderr).lower()
        
        if 'error' in combined_output:
            error_msg = """
[ERROR] OrthoFinder was not installed (properly) or set up correctly!

> Set the bin/ folder in your conda environment to have the highest priority for path lookup:
$ export PATH=$CONDA_PREFIX/bin:$PATH

> Then, test if you are able to see the help message for OrthoFinder without errors:
$ orthofinder -h

Finally, please try re-running bofasa again!
"""
            return False, error_msg
        
        # Verify successful execution with proper output
        if result.returncode == 0 and ("OrthoFinder" in result.stdout or "orthofinder" in combined_output):
            return True, "OrthoFinder is properly set up"
        else:
            error_msg = """
[ERROR] OrthoFinder command failed to execute properly!

> Set the bin/ folder in your conda environment to have the highest priority for path lookup:
$ export PATH=$CONDA_PREFIX/bin:$PATH

> Then, test if you are able to see the help message for OrthoFinder without errors:
$ orthofinder -h

Finally, please try re-running bofasa again!
"""
            return False, error_msg
            
    except subprocess.TimeoutExpired:
        error_msg = """
[ERROR] OrthoFinder command timed out!

This may indicate an installation or configuration issue.

> Try setting the bin/ folder in your conda environment to have the highest priority:
$ export PATH=$CONDA_PREFIX/bin:$PATH

> Test the OrthoFinder installation:
$ orthofinder -h

Finally, please try re-running bofasa again!
"""
        return False, error_msg
        
    except FileNotFoundError:
        error_msg = """
[ERROR] OrthoFinder was not installed (properly) or set up correctly!

> Set the bin/ folder in your conda environment to have the highest priority for path lookup:
$ export PATH=$CONDA_PREFIX/bin:$PATH

> Then, test if you are able to see the help message for OrthoFinder without errors:
$ orthofinder -h

Finally, please try re-running bofasa again!
"""
        return False, error_msg
        
    except Exception as e:
        error_msg = f"""
[ERROR] Unexpected error while checking OrthoFinder: {str(e)}

> Try setting the bin/ folder in your conda environment to have the highest priority:
$ export PATH=$CONDA_PREFIX/bin:$PATH

> Test the OrthoFinder installation:
$ orthofinder -h

Finally, please try re-running bofasa again!
"""
        return False, error_msg
