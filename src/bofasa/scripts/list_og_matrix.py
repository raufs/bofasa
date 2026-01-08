#!/usr/bin/env python3
"""
Script to list OrthoFinder matrices in simple table form.

Takes an OrthoFinder Orthogroups.tsv file and outputs a two-column
tab-separated table with OG ID and protein ID.
"""

import argparse
import sys
import os


def list_og_matrix(input_file, output_file=None):
    """
    Parse OrthoFinder Orthogroups.tsv and output simple two-column table.

    Args:
        input_file: Path to Orthogroups.tsv file
        output_file: Optional output file (default: stdout)
    """
    if not os.path.isfile(input_file):
        sys.stderr.write(f"Error: File not found: {input_file}\n")
        sys.exit(1)

    # Open output
    if output_file:
        outf = open(output_file, "w")
    else:
        outf = sys.stdout

    try:
        # Write header
        outf.write("og_id\tprotein_id\n")

        # Parse input file
        with open(input_file) as inf:
            for i, line in enumerate(inf):
                # Skip header row
                if i == 0:
                    continue

                line = line.rstrip("\n")
                ls = line.split("\t")

                if len(ls) < 1:
                    continue

                og_id = ls[0]

                # Parse all columns after the first (which contain protein IDs)
                for col in ls[1:]:
                    # Split by comma (proteins are comma-separated within each column)
                    for protein_id in col.split(","):
                        protein_id = protein_id.strip()
                        if protein_id != "":
                            outf.write(f"{og_id}\t{protein_id}\n")

    finally:
        if output_file:
            outf.close()


def main():
    parser = argparse.ArgumentParser(
        description="List OrthoFinder matrices in simple two-column table format.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Output to stdout
  list_og_matrix.py -i Orthogroups.tsv
  
  # Output to file
  list_og_matrix.py -i Orthogroups.tsv -o og_table.tsv
  
  # Using with OrthoFinder results directory
  list_og_matrix.py -i OrthoFinder_Results/Results_*/Orthogroups/Orthogroups.tsv
""",
    )

    parser.add_argument(
        "-i", "--input", required=True, help="Path to OrthoFinder Orthogroups.tsv\nor Orthogroups_Unassigned.tsv file"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=False,
        default=None,
        help="Output file (default: stdout)",
    )

    args = parser.parse_args()

    list_og_matrix(args.input, args.output)


if __name__ == "__main__":
    main()
