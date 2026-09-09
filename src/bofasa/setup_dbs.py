"""
Setup functions for bofasa.

This module contains functions for setting up annotation databases
required for bofasa analysis.
"""

import os
import sys
import subprocess
import gzip
import shutil
from . import config
from Bio import SeqIO
from .utils import setup_ready_directory


def setup_annotation_databases(
    output_dir: str, threads: int = config.DEFAULT_THREADS, force: bool = False,
    auto: bool = False
) -> None:
    """
    Setup annotation databases.

    Args:
        output_dir: Output directory for databases
        threads: Number of threads to use
        force: Whether to force overwrite existing databases
        auto: Whether to automatically answer 'yes' to interactive prompts

    Returns:
        None: Creates databases in output directory
    """
    # Check if output directory already exists and has content
    if os.path.exists(output_dir) and os.listdir(output_dir):
        sys.stdout.write(
            f"\nWARNING: The directory {output_dir} already exists and contains files!\n"
        )
        sys.stdout.write("This setup process may overwrite existing database files.\n")

        if auto and not force:
            sys.stdout.write("--auto specified: continuing with setup...\n")
        elif not force:
            while True:
                response = (
                    input(
                        "Do you want to continue and potentially overwrite existing files? (y/N): "
                    )
                    .strip()
                    .lower()
                )
                if response in ['y', 'yes']:
                    break
                elif response in ['n', 'no', '']:
                    sys.stdout.write(
                        "Setup cancelled. Please choose a different directory or use --force to skip this prompt.\n"
                    )
                    return
                else:
                    sys.stdout.write("Please enter 'y' for yes or 'n' for no.\n")

        sys.stdout.write("Proceeding with setup...\n\n")

    # Create output directory
    setup_ready_directory([output_dir], overwrite_mode="skip")

    sys.stdout.write(f"Setting up annotation databases in {output_dir}\n")
    sys.stdout.write(f"Using {threads} threads\n")

    # Database paths file
    db_paths_file = os.path.join(output_dir, "database_location_paths.txt")

    if os.path.exists(db_paths_file) and not force:
        sys.stdout.write("Database paths file already exists. Use --force to overwrite.\n")
        return

    # Setup geNomad database
    sys.stdout.write("Setting up geNomad database...\n")
    genomad_db_dir = os.path.join(output_dir, "genomad_db")
    
    # If force is True and directory exists, remove it
    if force and os.path.exists(genomad_db_dir):
        sys.stdout.write(f"Removing existing geNomad database directory: {genomad_db_dir}\n")
        shutil.rmtree(genomad_db_dir)
    
    setup_ready_directory([genomad_db_dir], overwrite_mode="overwrite")

    # Download geNomad database
    genomad_cmd = ["genomad", "download-database", genomad_db_dir]

    try:
        result = subprocess.run(genomad_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stdout.write(
                f"Error setting up geNomad database: Command returned exit status {result.returncode}\n"
            )
            sys.stdout.write(f"Command: {' '.join(genomad_cmd)}\n")
            sys.stdout.write(f"stdout: {result.stdout}\n")
            sys.stdout.write(f"stderr: {result.stderr}\n")
            sys.stdout.write("\nPlease ensure geNomad is installed and accessible.\n")
            sys.stdout.write("You can install geNomad with: conda install -c bioconda genomad\n")
            return
        sys.stdout.write("geNomad database setup completed\n")
    except FileNotFoundError:
        sys.stdout.write("Error: 'genomad' command not found.\n")
        sys.stdout.write("Please install geNomad first:\n")
        sys.stdout.write("  conda install -c bioconda genomad\n")
        sys.stdout.write("  or visit: https://github.com/apcamargo/genomad\n")
        return
    except Exception as e:
        sys.stdout.write(f"Error setting up geNomad database: {e}\n")
        sys.stdout.write("Please ensure geNomad is installed and accessible\n")
        return

    # Setup ISfinder database
    sys.stdout.write("Setting up ISfinder database...\n")
    isfinder_dir = os.path.join(output_dir, "isfinder_db")
    
    # If force is True and directory exists, remove it
    if force and os.path.exists(isfinder_dir):
        sys.stdout.write(f"Removing existing ISfinder database directory: {isfinder_dir}\n")
        shutil.rmtree(isfinder_dir)
    
    setup_ready_directory([isfinder_dir], overwrite_mode="overwrite")

    # Download ISfinder database and create DIAMOND database
    isfinder_fasta = os.path.join(isfinder_dir, "ISfinder.faa")
    isfinder_proc_fasta = os.path.join(isfinder_dir, "ISfinder_proc.faa")
    isfinder_dmnd = os.path.join(isfinder_dir, "ISfinder.dmnd")

    # Download ISfinder database from GitHub repository
    isfinder_url = "https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/refs/heads/master/IS.faa"
    
    try:
        sys.stdout.write(f"Downloading ISfinder database from {isfinder_url}...\n")
        
        # If force is True and file exists, remove it
        if force and os.path.exists(isfinder_fasta):
            sys.stdout.write(f"Removing existing ISfinder FASTA file: {isfinder_fasta}\n")
            os.remove(isfinder_fasta)
        
        # Use curl to download the file
        curl_cmd = ["curl", "-L", "-o", isfinder_fasta, isfinder_url]
        result = subprocess.run(curl_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stdout.write(
                f"Error downloading ISfinder database: Command returned exit status {result.returncode}\n"
            )
            sys.stdout.write(f"Command: {' '.join(curl_cmd)}\n")
            sys.stdout.write(f"stdout: {result.stdout}\n")
            sys.stdout.write(f"stderr: {result.stderr}\n")
            raise Exception("Failed to download ISfinder database")
        
        sys.stdout.write("ISfinder database downloaded successfully\n")
        
    except Exception as e:
        sys.stdout.write(f"Error downloading ISfinder database: {e}\n")
        sys.stdout.write("Please ensure you have internet connectivity and curl is available\n")
        return

    if os.path.exists(isfinder_fasta):
        # If force is True and DIAMOND database exists, remove it
        if force and os.path.exists(isfinder_dmnd):
            sys.stdout.write(f"Removing existing ISfinder DIAMOND database: {isfinder_dmnd}\n")
            os.remove(isfinder_dmnd)
        
        oipf = open(isfinder_proc_fasta, 'w')
        with open(isfinder_fasta, 'r') as f:
            for rec in SeqIO.parse(f, 'fasta'):
                if '~~~Passenger Gene~~~' in rec.description:
                    continue
                oipf.write(f">{rec.description}\n{rec.seq}\n")
        oipf.close()

        try:
            # Create DIAMOND database
            diamond_cmd = [
                "diamond",
                "makedb",
                "--in",
                isfinder_proc_fasta,
                "--db",
                isfinder_dmnd,
                "--threads",
                str(threads),
            ]
            result = subprocess.run(diamond_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                sys.stdout.write(
                    f"Error creating ISfinder DIAMOND database: Command returned exit status {result.returncode}\n"
                )
                sys.stdout.write(f"Command: {' '.join(diamond_cmd)}\n")
                sys.stdout.write(f"stdout: {result.stdout}\n")
                sys.stdout.write(f"stderr: {result.stderr}\n")
            else:
                sys.stdout.write("ISfinder DIAMOND database created\n")
        except Exception as e:
            sys.stdout.write(f"Error creating ISfinder DIAMOND database: {e}\n")
    else:
        sys.stdout.write("ISfinder FASTA file not found. Please download manually.\n")

    # Setup Pfam database
    sys.stdout.write("Setting up Pfam database...\n")
    pfam_dir = os.path.join(output_dir, "pfam_db")
    
    # If force is True and directory exists, remove it
    if force and os.path.exists(pfam_dir):
        sys.stdout.write(f"Removing existing Pfam database directory: {pfam_dir}\n")
        shutil.rmtree(pfam_dir)
    
    setup_ready_directory([pfam_dir], overwrite_mode="overwrite")

    pfam_hmm = os.path.join(pfam_dir, "Pfam-A.hmm")
    pfam_hmm_gz = os.path.join(pfam_dir, "Pfam-A.hmm.gz")

    # Download Pfam database
    pfam_url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

    try:
        sys.stdout.write(f"Downloading Pfam database from {pfam_url}...\n")

        # If force is True and files exist, remove them
        if force:
            if os.path.exists(pfam_hmm_gz):
                sys.stdout.write(f"Removing existing Pfam gzipped file: {pfam_hmm_gz}\n")
                os.remove(pfam_hmm_gz)
            if os.path.exists(pfam_hmm):
                sys.stdout.write(f"Removing existing Pfam HMM file: {pfam_hmm}\n")
                os.remove(pfam_hmm)

        # Use curl to download the file
        curl_cmd = ["curl", "-L", "-o", pfam_hmm_gz, pfam_url]
        result = subprocess.run(curl_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            sys.stdout.write(
                f"Error downloading Pfam database: Command returned exit status {result.returncode}\n"
            )
            sys.stdout.write(f"Command: {' '.join(curl_cmd)}\n")
            sys.stdout.write(f"stdout: {result.stdout}\n")
            sys.stdout.write(f"stderr: {result.stderr}\n")
            raise Exception("Failed to download Pfam database")

        # Extract gzipped file
        with gzip.open(pfam_hmm_gz, 'rb') as f_in:
            with open(pfam_hmm, 'wb') as f_out:
                f_out.write(f_in.read())

        # Remove gzipped file
        os.remove(pfam_hmm_gz)

        # Count records for Z parameter
        pfam_z = 0
        with open(pfam_hmm, 'r') as f:
            for line in f:
                if line.startswith('NAME'):
                    pfam_z += 1

        sys.stdout.write(f"Pfam database setup completed with {pfam_z} records\n")

    except Exception as e:
        sys.stdout.write(f"Error setting up Pfam database: {e}\n")
        sys.stdout.write(f"Please download manually from: {pfam_url}\n")

    # Write database paths file
    with open(db_paths_file, 'w') as f:
        f.write("Database\tType\tPath\tAdditional_Info\n")
        # geNomad database path should include the genomad_db subdirectory
        genomad_db_path = os.path.join(genomad_db_dir, "genomad_db")
        f.write(f"genomad\tgenomad\t{genomad_db_path}\t\n")
        if os.path.exists(isfinder_dmnd):
            f.write(f"isfinder\tdiamond\t{isfinder_dmnd}\t\n")
        if os.path.exists(pfam_hmm):
            f.write(f"pfam\thmm\t{pfam_hmm}\t{pfam_z}\n")

    sys.stdout.write("Annotation databases setup completed successfully\n")
    sys.stdout.write(f"Database paths saved to: {db_paths_file}\n")
