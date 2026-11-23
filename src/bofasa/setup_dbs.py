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
import email.utils
import threading
import concurrent.futures
import time
from typing import List, Optional, Tuple
from . import config
from Bio import SeqIO
from .utils import setup_ready_directory


_CURL_BASE: List[str] = [
    "curl", "-f", "-s", "-S", "-L", "--retry", "5", "--retry-delay", "5", "--connect-timeout", "60"
]


def _probe_remote(url: str) -> Tuple[Optional[int], Optional[float], bool]:
    """
    Return (size in bytes or None, Last-Modified as epoch seconds or None,
    whether byte-range requests are honored) for a remote file.
    """
    size, mtime = None, None
    try:
        res = subprocess.run(
            ["curl", "-s", "-I", "-L", "--max-time", "60", url], capture_output=True, text=True
        )
        # With redirects, curl prints a header block per hop - only the final one is relevant
        for line in res.stdout.splitlines():
            if line.upper().startswith("HTTP/"):
                size, mtime = None, None
                continue
            key, _, val = line.partition(":")
            key, val = key.strip().lower(), val.strip()
            if key == "content-length" and val.isdigit():
                size = int(val) if int(val) > 0 else None
            elif key == "last-modified":
                try:
                    mtime = email.utils.parsedate_to_datetime(val).timestamp()
                except Exception:
                    pass
    except Exception:
        pass

    supports_ranges = False
    if size is not None and size >= config.DOWNLOAD_SEGMENT_MIN_SIZE:
        # Servers do not always advertise range support, so test with an actual 1-byte
        # range request (reading at most 2 bytes so a server ignoring the range does not
        # start a full download)
        try:
            proc = subprocess.Popen(
                ["curl", "-s", "-L", "--max-time", "60", "-r", "0-0", url], stdout=subprocess.PIPE
            )
            data = proc.stdout.read(2)
            proc.kill()
            proc.wait()
            supports_ranges = len(data) == 1
        except Exception:
            supports_ranges = False
    return size, mtime, supports_ranges


def _curl_ok(cmd: List[str]) -> bool:
    try:
        return subprocess.run(cmd).returncode == 0
    except Exception as e:
        sys.stderr.write(f"Error: {e}\n")
        return False


def _download_file(url: str, output_file: str, segments: int = config.DEFAULT_DOWNLOAD_SEGMENTS) -> bool:
    """
    Download url to output_file - in parallel byte-range segments when the file is
    large and the server allows it - printing periodic progress.

    Args:
        url: URL to download
        output_file: Path to write the download to
        segments: Number of parallel connections used for large files (1 disables segmenting)

    Returns:
        bool: True if the download succeeded (and matches the size reported by the server).
    """
    segments = max(1, segments)
    size, mtime, supports_ranges = _probe_remote(url)
    if os.path.isfile(output_file):
        os.remove(output_file)

    use_segments = segments > 1 and supports_ranges and size is not None
    if use_segments:
        chunk = -(-size // segments)  # ceiling division
        bounds = [(start, min(start + chunk, size) - 1) for start in range(0, size, chunk)]
        parts = [f"{output_file}.part{i}" for i in range(len(bounds))]
    else:
        bounds, parts = [], [output_file]

    # Report progress every 20 seconds while downloading
    done = threading.Event()

    def report() -> None:
        start = time.time()
        while not done.wait(20):
            got = sum(os.path.getsize(f) for f in parts if os.path.isfile(f))
            total = f" of {size / 1e6:,.0f}" if size else ""
            sys.stdout.write(
                f"  ... downloaded {got / 1e6:,.0f}{total} MB ({(time.time() - start) / 60:.1f} min elapsed)\n"
            )
            sys.stdout.flush()

    reporter = threading.Thread(target=report, daemon=True)
    reporter.start()

    try:
        if use_segments:
            sys.stdout.write(f"  Downloading in {len(bounds)} parallel segments...\n")

            def fetch(i: int) -> bool:
                start, end = bounds[i]
                for _ in range(3):
                    _curl_ok(_CURL_BASE + ["-r", f"{start}-{end}", "-o", parts[i], url])
                    if os.path.isfile(parts[i]) and os.path.getsize(parts[i]) == end - start + 1:
                        return True
                return False

            with concurrent.futures.ThreadPoolExecutor(max_workers=len(bounds)) as ex:
                success = all(ex.map(fetch, range(len(bounds))))
            if success:
                with open(output_file, "wb") as tgt:
                    for part in parts:
                        with open(part, "rb") as src:
                            shutil.copyfileobj(src, tgt)
            for part in parts:
                if os.path.isfile(part):
                    os.remove(part)
            success = success and os.path.isfile(output_file) and os.path.getsize(output_file) == size
        else:
            success = _curl_ok(_CURL_BASE + ["-o", output_file, url]) and os.path.isfile(output_file)
            if success and size is not None and os.path.getsize(output_file) != size:
                sys.stdout.write(
                    f"Warning: size of {output_file} differs from the size reported by the server for {url}\n"
                )
    finally:
        done.set()

    if success and mtime is not None:
        os.utime(output_file, (mtime, mtime))  # keep the server's timestamp (like curl -R)
    return success


def setup_annotation_databases(
    output_dir: str, threads: int = config.DEFAULT_THREADS, force: bool = False,
    auto: bool = False, skip_genomad: bool = False, pfam_version: str = "current",
    download_segments: int = config.DEFAULT_DOWNLOAD_SEGMENTS
) -> bool:
    """
    Setup annotation databases.

    Args:
        output_dir: Output directory for databases
        threads: Number of threads to use
        force: Whether to force overwrite existing databases
        auto: Whether to automatically answer 'yes' to interactive prompts
        skip_genomad: Whether to skip downloading the geNomad database
        pfam_version: Pfam release to download - "current" for the current release,
            or a specific release listed in config.PFAM_RELEASES (e.g. "37.4")
        download_segments: Number of parallel connections used to download each large
            file in byte-range segments, if the server supports it (1 disables this)

    Returns:
        bool: True if setup completed successfully (or an existing setup was kept),
        False if any database failed to set up or the user cancelled.
    """
    if pfam_version != "current" and pfam_version not in config.PFAM_RELEASES:
        sys.stderr.write(
            f"Error: unsupported Pfam version '{pfam_version}'. Options are: "
            f"current, {', '.join(sorted(config.PFAM_RELEASES))}\n"
        )
        return False

    # Use an absolute path so all recorded database locations are absolute
    output_dir = os.path.abspath(output_dir)

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
                    return False
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
        # Don't report success for an existing but incomplete/broken setup (e.g. left
        # by an earlier failed run): validate what the manifest points to.
        from .utils import check_bofasa_db_setup
        dbs_ok, dbs_message = check_bofasa_db_setup(require_genomad=not skip_genomad)
        expected_z = config.PFAM_RELEASES.get(pfam_version)
        if dbs_ok and expected_z is not None:
            existing_z = None
            with open(db_paths_file) as f:
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) >= 4 and fields[0] == 'pfam':
                        existing_z = fields[3]
            if existing_z != str(expected_z):
                sys.stderr.write(
                    f"Error: the existing Pfam database in {output_dir} has {existing_z} HMMs, "
                    f"which does not match Pfam {pfam_version} ({expected_z} HMMs).\n"
                    "Re-run with --force to replace the databases.\n"
                )
                return False
        if dbs_ok:
            sys.stdout.write(
                "Databases are already set up and validated. Use --force to re-download/overwrite.\n"
            )
            return True
        sys.stderr.write(
            f"Error: an existing database setup was found in {output_dir} but it is incomplete or invalid:\n"
            f"  {dbs_message}\n"
            "Re-run with --force to rebuild the databases.\n"
        )
        return False

    # Setup geNomad database (downloaded by default, skipped with --skip-genomad)
    # The geNomad database path should include the genomad_db subdirectory created
    # by 'genomad download-database'.
    genomad_db_dir = os.path.join(output_dir, "genomad_db")
    genomad_db_path = os.path.join(genomad_db_dir, "genomad_db")

    if skip_genomad:
        sys.stdout.write(
            "Skipping geNomad database setup (--skip-genomad specified).\n"
            "Note: 'bofasa prep --run-genomad' will not be usable until the geNomad\n"
            "database is set up.\n"
        )
        # Keep recording a previously downloaded geNomad database (e.g. with --force --skip-genomad)
        if os.path.isdir(genomad_db_path) and os.listdir(genomad_db_path):
            sys.stdout.write(f"Keeping existing geNomad database: {genomad_db_path}\n")
        else:
            genomad_db_path = None
    else:
        sys.stdout.write("Setting up geNomad database...\n")

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
                sys.stdout.write(
                    "Alternatively, rerun with --skip-genomad to set up the other databases only.\n"
                )
                return False
            sys.stdout.write("geNomad database setup completed\n")
        except FileNotFoundError:
            sys.stdout.write("Error: 'genomad' command not found.\n")
            sys.stdout.write("Please install geNomad first:\n")
            sys.stdout.write("  conda install -c bioconda genomad\n")
            sys.stdout.write("  or visit: https://github.com/apcamargo/genomad\n")
            sys.stdout.write(
                "Alternatively, rerun with --skip-genomad to set up the other databases only.\n"
            )
            return False
        except Exception as e:
            sys.stdout.write(f"Error setting up geNomad database: {e}\n")
            sys.stdout.write("Please ensure geNomad is installed and accessible\n")
            sys.stdout.write(
                "Alternatively, rerun with --skip-genomad to set up the other databases only.\n"
            )
            return False

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

    # Download ISfinder database from GitHub repository (pinned to the most recent
    # commit of the repository, which has not been updated since October 2020)
    isfinder_url = config.ISFINDER_URL
    
    try:
        sys.stdout.write(f"Downloading ISfinder database from {isfinder_url}...\n")
        
        # If force is True and file exists, remove it
        if force and os.path.exists(isfinder_fasta):
            sys.stdout.write(f"Removing existing ISfinder FASTA file: {isfinder_fasta}\n")
            os.remove(isfinder_fasta)
        
        if not _download_file(isfinder_url, isfinder_fasta, segments=download_segments):
            raise Exception("Failed to download ISfinder database")

        # Validate that the download looks like a FASTA file (not e.g. an HTML error page)
        first_line = ''
        with open(isfinder_fasta, 'r', errors='replace') as f:
            for line in f:
                if line.strip():
                    first_line = line.strip()
                    break
        if not first_line.startswith('>'):
            raise Exception(f"Downloaded file {isfinder_fasta} does not look like a FASTA file")
        
        sys.stdout.write("ISfinder database downloaded successfully\n")
        
    except Exception as e:
        sys.stdout.write(f"Error downloading ISfinder database: {e}\n")
        sys.stdout.write("Please ensure you have internet connectivity and curl is available\n")
        return False

    setup_ok = True

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
                setup_ok = False
            else:
                sys.stdout.write("ISfinder DIAMOND database created\n")
        except Exception as e:
            sys.stdout.write(f"Error creating ISfinder DIAMOND database: {e}\n")
            setup_ok = False
    else:
        sys.stdout.write("ISfinder FASTA file not found. Please download manually.\n")
        setup_ok = False

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
    if pfam_version == "current":
        pfam_url = config.PFAM_CURRENT_RELEASE_URL
    else:
        pfam_url = config.PFAM_RELEASE_URL_TEMPLATE.format(version=pfam_version)
    pfam_z = 0

    try:
        sys.stdout.write(f"Downloading Pfam database (release: {pfam_version}) from {pfam_url}...\n")

        # If force is True and files exist, remove them
        if force:
            if os.path.exists(pfam_hmm_gz):
                sys.stdout.write(f"Removing existing Pfam gzipped file: {pfam_hmm_gz}\n")
                os.remove(pfam_hmm_gz)
            if os.path.exists(pfam_hmm):
                sys.stdout.write(f"Removing existing Pfam HMM file: {pfam_hmm}\n")
                os.remove(pfam_hmm)

        if not _download_file(pfam_url, pfam_hmm_gz, segments=download_segments):
            raise Exception("Failed to download Pfam database")

        # Extract gzipped file (streamed to a temporary file, then moved into place)
        pfam_hmm_tmp = pfam_hmm + ".tmp"
        with gzip.open(pfam_hmm_gz, 'rb') as f_in:
            with open(pfam_hmm_tmp, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.replace(pfam_hmm_tmp, pfam_hmm)

        # Remove gzipped file
        os.remove(pfam_hmm_gz)

        # Count records for Z parameter
        with open(pfam_hmm, 'r') as f:
            for line in f:
                if line.startswith('NAME'):
                    pfam_z += 1

        if pfam_z <= 0:
            raise Exception(f"No HMM records found in {pfam_hmm}")

        expected_z = config.PFAM_RELEASES.get(pfam_version)
        if expected_z is not None and pfam_z != expected_z:
            raise Exception(
                f"Expected {expected_z} HMM records for Pfam {pfam_version} but found {pfam_z}"
            )

        sys.stdout.write(f"Pfam database setup completed with {pfam_z} records\n")

    except Exception as e:
        sys.stdout.write(f"Error setting up Pfam database: {e}\n")
        sys.stdout.write(f"Please download manually from: {pfam_url}\n")
        setup_ok = False

    # Write database paths file, recording absolute paths for every database
    with open(db_paths_file, 'w') as f:
        f.write("Database\tType\tPath\tAdditional_Info\n")
        if genomad_db_path is not None:
            f.write(f"genomad\tgenomad\t{os.path.abspath(genomad_db_path)}\t\n")
        if os.path.exists(isfinder_dmnd):
            f.write(f"isfinder\tdiamond\t{os.path.abspath(isfinder_dmnd)}\t\n")
        if os.path.exists(pfam_hmm) and pfam_z > 0:
            f.write(f"pfam\thmm\t{os.path.abspath(pfam_hmm)}\t{pfam_z}\n")

    if not setup_ok:
        sys.stdout.write("Annotation database setup finished with errors - see messages above.\n")
        sys.stdout.write(f"Database paths (for successfully set up databases) saved to: {db_paths_file}\n")
        return False

    sys.stdout.write("Annotation databases setup completed successfully\n")
    sys.stdout.write(f"Database paths saved to: {db_paths_file}\n")
    return True
