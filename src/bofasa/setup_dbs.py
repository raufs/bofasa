"""
Setup functions for bofasa.

This module contains functions for setting up annotation databases
required for bofasa analysis.
"""

import os
import subprocess
import gzip
import shutil
from . import config


def setup_annotation_databases(
    output_dir: str, threads: int = config.DEFAULT_THREADS, force: bool = False
) -> None:
    """
    Setup annotation databases.

    Args:
        output_dir: Output directory for databases
        threads: Number of threads to use
        force: Whether to force overwrite existing databases

    Returns:
        None: Creates databases in output directory
    """
    # Check if output directory already exists and has content
    if os.path.exists(output_dir) and os.listdir(output_dir):
        print(
            f"\nWARNING: The directory {output_dir} already exists and contains files!"
        )
        print("This setup process may overwrite existing database files.")

        if not force:
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
                    print(
                        "Setup cancelled. Please choose a different directory or use --force to skip this prompt."
                    )
                    return
                else:
                    print("Please enter 'y' for yes or 'n' for no.")

        print("Proceeding with setup...\n")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    print(f"Setting up annotation databases in {output_dir}")
    print(f"Using {threads} threads")

    # Database paths file
    db_paths_file = os.path.join(output_dir, "database_location_paths.txt")

    if os.path.exists(db_paths_file) and not force:
        print("Database paths file already exists. Use --force to overwrite.")
        return

    # Setup geNomad database
    print("Setting up geNomad database...")
    genomad_db_dir = os.path.join(output_dir, "genomad_db")
    
    # If force is True and directory exists, remove it
    if force and os.path.exists(genomad_db_dir):
        print(f"Removing existing geNomad database directory: {genomad_db_dir}")
        shutil.rmtree(genomad_db_dir)
    
    os.makedirs(genomad_db_dir, exist_ok=True)

    # Download geNomad database
    genomad_cmd = ["genomad", "download-database", genomad_db_dir]

    try:
        result = subprocess.run(genomad_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(
                f"Error setting up geNomad database: Command returned exit status {result.returncode}"
            )
            print(f"Command: {' '.join(genomad_cmd)}")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            print("\nPlease ensure geNomad is installed and accessible.")
            print("You can install geNomad with: conda install -c bioconda genomad")
            return
        print("geNomad database setup completed")
    except FileNotFoundError:
        print("Error: 'genomad' command not found.")
        print("Please install geNomad first:")
        print("  conda install -c bioconda genomad")
        print("  or visit: https://github.com/apcamargo/genomad")
        return
    except Exception as e:
        print(f"Error setting up geNomad database: {e}")
        print("Please ensure geNomad is installed and accessible")
        return

    # Setup ISfinder database
    print("Setting up ISfinder database...")
    isfinder_dir = os.path.join(output_dir, "isfinder_db")
    
    # If force is True and directory exists, remove it
    if force and os.path.exists(isfinder_dir):
        print(f"Removing existing ISfinder database directory: {isfinder_dir}")
        shutil.rmtree(isfinder_dir)
    
    os.makedirs(isfinder_dir, exist_ok=True)

    # Download ISfinder database and create DIAMOND database
    isfinder_fasta = os.path.join(isfinder_dir, "ISfinder.faa")
    isfinder_dmnd = os.path.join(isfinder_dir, "ISfinder.dmnd")

    # Download ISfinder database from GitHub repository
    isfinder_url = "https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/refs/heads/master/IS.faa"
    
    try:
        print(f"Downloading ISfinder database from {isfinder_url}...")
        
        # If force is True and file exists, remove it
        if force and os.path.exists(isfinder_fasta):
            print(f"Removing existing ISfinder FASTA file: {isfinder_fasta}")
            os.remove(isfinder_fasta)
        
        # Use curl to download the file
        curl_cmd = ["curl", "-L", "-o", isfinder_fasta, isfinder_url]
        result = subprocess.run(curl_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(
                f"Error downloading ISfinder database: Command returned exit status {result.returncode}"
            )
            print(f"Command: {' '.join(curl_cmd)}")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
            raise Exception("Failed to download ISfinder database")
        
        print("ISfinder database downloaded successfully")
        
    except Exception as e:
        print(f"Error downloading ISfinder database: {e}")
        print("Please ensure you have internet connectivity and curl is available")
        return

    if os.path.exists(isfinder_fasta):
        print("ISfinder FASTA file already exists - skipping download (checkpoint)")
        # If force is True and DIAMOND database exists, remove it
        if force and os.path.exists(isfinder_dmnd):
            print(f"Removing existing ISfinder DIAMOND database: {isfinder_dmnd}")
            os.remove(isfinder_dmnd)
        
        try:
            # Create DIAMOND database
            diamond_cmd = [
                "diamond",
                "makedb",
                "--in",
                isfinder_fasta,
                "--db",
                isfinder_dmnd,
                "--threads",
                str(threads),
            ]
            result = subprocess.run(diamond_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print(
                    f"Error creating ISfinder DIAMOND database: Command returned exit status {result.returncode}"
                )
                print(f"Command: {' '.join(diamond_cmd)}")
                print(f"stdout: {result.stdout}")
                print(f"stderr: {result.stderr}")
            else:
                print("ISfinder DIAMOND database created")
        except FileNotFoundError:
            print("Error: 'diamond' command not found.")
            print("Please install DIAMOND first:")
            print("  conda install -c bioconda diamond")
        except Exception as e:
            print(f"Error creating ISfinder DIAMOND database: {e}")
    else:
        print("ISfinder FASTA file not found. Please download manually.")

    # Setup Pfam database
    print("Setting up Pfam database...")
    pfam_dir = os.path.join(output_dir, "pfam_db")
    
    # If force is True and directory exists, remove it
    if force and os.path.exists(pfam_dir):
        print(f"Removing existing Pfam database directory: {pfam_dir}")
        shutil.rmtree(pfam_dir)
    
    os.makedirs(pfam_dir, exist_ok=True)

    pfam_hmm = os.path.join(pfam_dir, "Pfam-A.hmm")
    pfam_hmm_gz = os.path.join(pfam_dir, "Pfam-A.hmm.gz")

    # Download Pfam database
    pfam_url = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"

    try:
        print(f"Downloading Pfam database from {pfam_url}...")

        # If force is True and files exist, remove them
        if force:
            if os.path.exists(pfam_hmm_gz):
                print(f"Removing existing Pfam gzipped file: {pfam_hmm_gz}")
                os.remove(pfam_hmm_gz)
            if os.path.exists(pfam_hmm):
                print(f"Removing existing Pfam HMM file: {pfam_hmm}")
                os.remove(pfam_hmm)

        # Use curl to download the file
        curl_cmd = ["curl", "-L", "-o", pfam_hmm_gz, pfam_url]
        result = subprocess.run(curl_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(
                f"Error downloading Pfam database: Command returned exit status {result.returncode}"
            )
            print(f"Command: {' '.join(curl_cmd)}")
            print(f"stdout: {result.stdout}")
            print(f"stderr: {result.stderr}")
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

        print(f"Pfam database setup completed with {pfam_z} records")

    except Exception as e:
        print(f"Error setting up Pfam database: {e}")
        print("Please download manually from:", pfam_url)

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

    print("Annotation databases setup completed successfully")
    print(f"Database paths saved to: {db_paths_file}")
