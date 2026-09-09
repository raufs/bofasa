#!python

### Program: extract_focal_protein_dogs
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
import sys
from collections import defaultdict

from rich_argparse import RawTextRichHelpFormatter

from bofasa import utils


def create_parser():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="""
Program: extract_focal_protein_dogs
Authors: Rauf Salamzade, Aamuktha Kottapalli, Lindsay R. Kalan 
Affiliation: University of Wisconsin - Madison, McMaster University

extract_focal_protein_dogs: Extract domain ortholog groups for a focal 
protein from BOFASA results. Creates a mini Domain Ortholog Groups file 
containing only DOGs found in the focal protein and their distribution 
across domain chunks (proteins) that have at least one of those DOGs.

Note: Domain chunks in BOFASA have the format sample|protein|annotation|iterator.
You can provide either the base protein name (e.g., "PROTEIN_001") or the 
full domain chunk name (e.g., "sample1|PROTEIN_001|pfam_domain|1").
""",
        formatter_class=RawTextRichHelpFormatter,
    )

    parser.add_argument(
        '-i',
        '--results-dir',
        help='BOFASA results directory (contains Final_Results subdirectory).',
        required=True,
    )
    parser.add_argument(
        '-p',
        '--protein-id',
        help='Identifier of the focal protein to extract DOGs for. Can be either '
             'the base protein name (e.g., "protein123") or a full domain chunk '
             'name (e.g., "sample|protein123|pfam_domain|1"). If a base name is '
             'provided, all domain chunks from that protein will be matched.',
        required=True,
    )
    parser.add_argument(
        '-o',
        '--output-file',
        help='Output file - default is standard output.',
        required=False,
        default=None,
    )
    parser.add_argument(
        '-v',
        '--version',
        action='store_true',
        help="Get version and exit.",
        required=False,
        default=False,
    )

    args = parser.parse_args()
    return args


def extract_focal_protein_dogs():
    """
    Void function which runs primary workflow for program.
    """

    # get version
    version = utils.get_version()

    if len(sys.argv) > 1 and ('-v' in set(sys.argv) or '--version' in set(sys.argv)):
        sys.stdout.write(version + '\n')
        sys.exit(0)

    # PARSE ARGUMENTS
    myargs = create_parser()

    results_dir = os.path.abspath(myargs.results_dir) + '/'
    protein_id = myargs.protein_id
    outfile = myargs.output_file

    # Validate results directory
    final_results_dir = os.path.join(results_dir, 'Final_Results')
    if not os.path.isdir(final_results_dir):
        sys.stderr.write(
            f"Error: Final_Results directory not found in {results_dir}\n"
        )
        sys.exit(1)

    # Path to Domain Ortholog Groups file
    dogs_file = os.path.join(final_results_dir, 'Domain_Ortholog_Groups.tsv')
    if not os.path.isfile(dogs_file):
        sys.stderr.write(
            f"Error: Domain_Ortholog_Groups.tsv not found in {final_results_dir}\n"
        )
        sys.exit(1)

    # START WORKFLOW
    sys.stderr.write(f"Reading Domain Ortholog Groups from: {dogs_file}\n")
    sys.stderr.write(f"Focal protein: {protein_id}\n")

    # Determine if protein_id is a base name or full domain chunk name
    # Domain chunks have format: sample|protein|annotation|iterator
    is_full_domain_chunk = protein_id.count('|') >= 3
    
    def matches_focal_protein(domain_chunk):
        """Check if a domain chunk matches the focal protein identifier."""
        domain_chunk = domain_chunk.strip()
        if not domain_chunk:
            return False
        
        if is_full_domain_chunk:
            # Exact match for full domain chunk name
            return domain_chunk == protein_id
        else:
            # Match base protein name (second field after splitting by |)
            parts = domain_chunk.split('|')
            if len(parts) >= 2:
                return parts[1] == protein_id
            return False

    # Step 1: Read the Domain Ortholog Groups file and find DOGs containing the focal protein
    focal_dogs = set()  # DOGs that contain the focal protein
    focal_domain_chunks = set()  # All domain chunks matching the focal protein
    dog_data = {}  # Store all DOG data: dog_id -> [sample1_proteins, sample2_proteins, ...]
    samples = []  # Sample names in order

    try:
        with open(dogs_file) as inf:
            for i, line in enumerate(inf):
                line = line.rstrip('\n')
                ls = line.split('\t')

                if i == 0:  # Header line
                    samples = ls[1:]  # First column is "OG/Sample"
                    continue

                dog_id = ls[0]
                dog_data[dog_id] = ls[1:]  # Store sample columns

                # Check if focal protein is in this DOG
                for sample_proteins in ls[1:]:
                    domain_chunks = [
                        p.strip() for p in sample_proteins.split(',') if p.strip()
                    ]
                    for dc in domain_chunks:
                        if matches_focal_protein(dc):
                            focal_dogs.add(dog_id)
                            focal_domain_chunks.add(dc)

    except Exception as e:
        sys.stderr.write(f"Error reading Domain_Ortholog_Groups.tsv: {str(e)}\n")
        sys.exit(1)

    if not focal_dogs:
        sys.stderr.write(
            f"Warning: Focal protein '{protein_id}' not found in any Domain Ortholog Groups\n"
        )
        sys.stderr.write("Creating empty output file.\n")
    else:
        sys.stderr.write(
            f"Found {len(focal_domain_chunks)} domain chunk(s) matching focal protein\n"
        )
        sys.stderr.write(
            f"Found {len(focal_dogs)} domain ortholog groups containing the focal protein\n"
        )

    # Step 2: Find all base proteins (sample|protein) that have at least one of the focal DOGs
    proteins_with_focal_dogs = set()  # Store as (sample, protein) tuples
    for dog_id in focal_dogs:
        sample_proteins_list = dog_data[dog_id]
        for sample_proteins in sample_proteins_list:
            domain_chunks = [p.strip() for p in sample_proteins.split(',') if p.strip()]
            for dc in domain_chunks:
                dc_parts = dc.split('|')
                if len(dc_parts) >= 2:
                    # Store the sample and protein name as a tuple
                    proteins_with_focal_dogs.add((dc_parts[0], dc_parts[1]))

    sys.stderr.write(
        f"Found {len(proteins_with_focal_dogs)} unique proteins with at least one focal DOG\n"
    )

    # Step 3: Filter DOGs to only include domain chunks from proteins with focal DOGs
    filtered_dogs = {}

    for dog_id in focal_dogs:
        sample_proteins_list = dog_data[dog_id]
        filtered_sample_data = []

        for i, sample_proteins in enumerate(sample_proteins_list):
            sample_name = samples[i]
            # Remove any suffix from sample name (e.g., .ccds) to get base sample name
            base_sample_name = sample_name.split('.')[0]
            
            # Get domain chunks in this sample for this DOG
            domain_chunks = [p.strip() for p in sample_proteins.split(',') if p.strip()]

            # Keep only domain chunks whose base protein has at least one focal DOG
            filtered_chunks = []
            for dc in domain_chunks:
                dc_parts = dc.split('|')
                if len(dc_parts) >= 2:
                    dc_sample = dc_parts[0]
                    dc_protein = dc_parts[1]
                    
                    # Only keep if this protein has at least one focal DOG
                    if (dc_sample, dc_protein) in proteins_with_focal_dogs:
                        filtered_chunks.append(dc)

            # Join back with comma separation
            filtered_sample_data.append(', '.join(filtered_chunks))

        filtered_dogs[dog_id] = filtered_sample_data

    # Step 4: Write output
    outf_handle = sys.stdout
    if outfile is not None:
        outf_handle = open(outfile, 'w')

    try:
        # Write header
        outf_handle.write('OG/Sample\t' + '\t'.join(samples) + '\n')

        # Write filtered DOGs in sorted order
        for dog_id in sorted(filtered_dogs.keys()):
            filtered_sample_data = filtered_dogs[dog_id]
            outf_handle.write(dog_id + '\t' + '\t'.join(filtered_sample_data) + '\n')

        sys.stderr.write(
            f"Successfully wrote {len(filtered_dogs)} domain ortholog groups to output\n"
        )

    finally:
        if outfile is not None:
            outf_handle.close()


if __name__ == '__main__':
    extract_focal_protein_dogs()
