#!python

### Program: extract_genomad_proteomes.py
### Author: Rauf Salamzade, Aamuktha Kottapalli, Lindsay R. Kalan
### Affiliation: University of Wisconsin - Madison, McMaster University

# BSD 3-Clause License
#
# Copyright (c) 2024, Rauf Salamzade
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import argparse
import os
from Bio import SeqIO
import sys
import shutil
from pathlib import Path

from rich_argparse import RawTextRichHelpFormatter


def create_parser():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="""
    Program: extract_genomad_proteomes.py
    Author: Rauf Salamzade
    Affiliation: University of Wisconsin - Madison, McMaster University

    Extract virus/prophage proteome files from multiple geNomad results directories.

    This script searches for virus protein FASTA files (*_virus_proteins.faa) in the
    summary folders of geNomad result directories and consolidates them into a single
    output directory or file.
    """,
        formatter_class=RawTextRichHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--genomad_dirs",
        nargs="+",
        required=True,
        help="Path(s) to geNomad results directories (space-separated for multiple).",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        required=True,
        help="Path to output directory where proteome files will be written.",
    )
    parser.add_argument(
        "-c",
        "--consolidate",
        action="store_true",
        help="For each sample (geNomad directory) just output a single FASTA file - i.e.\n"
             "just use the \"_virus_proteins.faa\" file in the summary directory of geNomad\n"
             "results (default: separate files for each prophage from each sample).",
    )

    args = parser.parse_args()
    return args


def find_virus_proteome_files(genomad_dir):
    """
    Find virus proteome files in a geNomad results directory.
    
    Args:
        genomad_dir (str): Path to geNomad results directory.
        
    Returns:
        list: List of paths to virus proteome .faa files.
    """
    proteome_files = []
    genomad_path = Path(genomad_dir)
    
    if not genomad_path.exists():
        sys.stderr.write(f"Warning: Directory does not exist: {genomad_dir}\n")
        return proteome_files
    
    # Look for summary directories
    for item in genomad_path.iterdir():
        if item.is_dir() and item.name.endswith("_summary"):
            # Look for virus protein files in the summary directory
            for file in item.iterdir():
                if file.is_file() and file.name.endswith("_virus_proteins.faa"):
                    proteome_files.append(file)
    
    return proteome_files


def extract_sample_name(genomad_dir):
    """
    Extract a sample name from the geNomad directory path.
    
    Args:
        genomad_dir (str): Path to geNomad results directory
        
    Returns:
        str: Sample name
    """
    # Use the directory name as the sample name
    return Path(genomad_dir).name

def split_genomad_virus_file(genomad_virus_protein_file, output_dir, sample_name):
    """
    Split a geNomad virus protein file into multiple files - one for each virus/phage.
    
    Args:
        genomad_virus_protein_file (str): Path to geNomad virus protein file
        output_dir (Path): Path to output directory
        sample_name (str): Sample name
        
    Returns:
        None
    """

    prophage_name_to_id = {}
    prophage_id = 1
    with open(genomad_virus_protein_file) as ogvpf:
        for rec in SeqIO.parse(ogvpf, 'fasta'):
            prophage_name = '_'.join(rec.id.split('_')[:-1])
            if not prophage_name in prophage_name_to_id:
                prophage_name_to_id[prophage_name] = str(prophage_id)
                prophage_id += 1

    with open(genomad_virus_protein_file) as ogvpf:
        for rec in SeqIO.parse(ogvpf, 'fasta'):
            prophage_name = '_'.join(rec.id.split('_')[:-1])
            prophage_id = prophage_name_to_id[prophage_name]
            output_faa = output_dir / (sample_name + '_' + prophage_id + '.faa')
            of_handle = open(output_faa, 'a+')
            of_handle.write('>' + rec.description + '\n' + str(rec.seq) + '\n')
            of_handle.close()

def main():
    """
    Main workflow for extracting virus proteomes from geNomad results.
    """
    # Parse arguments
    myargs = create_parser()
    
    genomad_dirs = myargs.genomad_dirs
    output_dir = Path(myargs.output_dir)
    consolidate_flag = myargs.consolidate
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    sys.stdout.write("=" * 80 + "\n")
    sys.stdout.write("extract_genomad_proteomes.py\n")
    sys.stdout.write("=" * 80 + "\n")
    sys.stdout.write(f"Processing {len(genomad_dirs)} geNomad results directories...\n")
    
    # Collect all proteome files
    all_proteome_files = []
    samples_processed = 0
    
    for genomad_dir in genomad_dirs:
        sample_name = extract_sample_name(genomad_dir)
        proteome_files = find_virus_proteome_files(genomad_dir)
        
        if not proteome_files:
            sys.stderr.write(f"Warning: No virus proteome files found in {genomad_dir}\n")
            continue
        
        for proteome_file in proteome_files:
            all_proteome_files.append((sample_name, proteome_file))
        
        samples_processed += 1
    
    if not all_proteome_files:
        sys.stderr.write("Error: No virus proteome files found in any of the provided directories.\n")
        sys.exit(1)
    
    sys.stdout.write(f"Found {len(all_proteome_files)} virus proteome file(s) from {samples_processed} sample(s)\n")
    
    # Process files based on consolidate option
    if consolidate_flag:
        sys.stdout.write(f"Copying proteome files to: {output_dir}\n")
        for sample_name, proteome_file in all_proteome_files:
            output_file = output_dir / f"{sample_name}.faa"
            shutil.copy2(proteome_file, output_file)
        sys.stdout.write(f"Successfully copied {len(all_proteome_files)} proteome file(s) to {output_dir}\n")

    else:        
        sys.stdout.write(f"Splitting genomad virus protein files into individual files for each virus/phage\n")        
        for sample_name, proteome_file in all_proteome_files:
            split_genomad_virus_file(proteome_file, output_dir, sample_name)
        sys.stdout.write(f"Successfully split {len(all_proteome_files)} genomad virus protein files into individual files for each virus/phage\n")
    
    sys.stdout.write("=" * 80 + "\n")
    sys.stdout.write("Extraction complete!\n")


if __name__ == "__main__":
    main()

