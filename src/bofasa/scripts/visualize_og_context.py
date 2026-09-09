#!/usr/bin/env python3

### Program: visualize_og_context.py
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
import subprocess
import sys
import pickle
import re
import warnings
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional

from rich_argparse import RawTextRichHelpFormatter

from bofasa import utils

try:
    from Bio import Phylo, SeqIO
    from Bio.Phylo.BaseTree import Tree, Clade
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

def create_parser():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="""
Program: visualize_og_context.py

visualize_og_context.py: Visualize the genomic context surrounding instances 
of a focal protein ortholog group across samples. Similar in functionality 
to lsaBGC-See, this tool creates visualizations showing genes across a 
species tree or gene tree.

The tool will:
1. Identify all instances of the focal protein OG across samples
2. Extract surrounding genomic context information
3. Assign colors to neighboring OGs for visualization
4. Create output files for visualization in R (using gggenes & ggtree) and iTOL
5. Handle paralogs by duplicating tree leaves for multi-instance samples
6. Optionally build a phylogeny from the focal OG if no tree is provided

Output Files:
=============
- OG_Context_Visualization.gggenes.txt: Input for gggenes R package
- OG_Context_Visualization.iTol.txt: iTOL annotation track
- species_tree.modified.nwk: Modified tree with duplicate leaves (if paralogs)
- OG_info.txt: Summary information about the focal OG
- Color_Scheme.txt: Color assignments for OGs in the context

Still in testing/development, use with caution!
""",
        formatter_class=RawTextRichHelpFormatter,
    )

    parser.add_argument(
        '-pd',
        '--prep-dir',
        help='BOFASA prep directory (output from "bofasa prep").',
        required=True,
    )
    parser.add_argument(
        '-rd',
        '--run-dir',
        help='BOFASA run directory (output from "bofasa run").',
        required=True,
    )
    parser.add_argument(
        '-og',
        '--focal-og',
        help='Focal protein ortholog group identifier to visualize context around.',
        required=True,
    )
    parser.add_argument(
        '-o',
        '--output-dir',
        help='Output directory for visualization files.',
        required=True,
    )
    parser.add_argument(
        '-s',
        '--species-tree',
        help='Species phylogeny in Newick format. If not provided, a gene tree\n'
             'will be constructed from the focal OG sequences.',
        required=False,
        default=None,
    )
    parser.add_argument(
        '-w',
        '--context-window',
        type=int,
        default=5,
        help='Number of genes to include upstream and downstream of focal OG\n'
             '[Default is 5].',
        required=False,
    )
    parser.add_argument(
        '-bp',
        '--context-bp',
        type=int,
        default=None,
        help='Alternative to --context-window: specify context as number of\n'
             'base pairs instead of gene count. Overrides --context-window.',
        required=False,
    )
    parser.add_argument(
        '-sc',
        '--single-copy-only',
        action='store_true',
        help='Only include samples where the focal OG is single-copy.',
        required=False,
        default=False,
    )
    parser.add_argument(
        '-ml',
        '--min-length',
        type=int,
        default=None,
        help='Minimum protein length (amino acids) for focal OG instances to\n'
             'include. Useful for filtering out truncated copies.',
        required=False,
    )
    parser.add_argument(
        '-sl',
        '--sample-list',
        help='File containing sample IDs to include (one per line). If not\n'
             'provided, all samples will be included.',
        required=False,
        default=None,
    )
    parser.add_argument(
        '-c',
        '--cpus',
        type=int,
        default=1,
        help='Number of CPUs to use for alignment and tree building\n'
             '[Default is 1].',
        required=False,
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


def load_prep_results(prep_dir: str) -> Dict:
    """
    Load results from bofasa prep directory.
    
    Args:
        prep_dir (str): Path to bofasa prep directory
        
    Returns:
        dict: Dictionary containing prep results
    """
    prep_data = {
        'sample_genbanks': {},
        'sample_proteins': {},
        'sample_metadata': {},
    }
    
    # Check for required directories
    genbank_dir = os.path.join(prep_dir, 'Sample_Annotation_Files')
    protein_dir = os.path.join(prep_dir, 'Sample_Protein_Files')
    bed_dir = os.path.join(prep_dir, 'Genome_Processing', 'BEDs')
    
    if not os.path.isdir(genbank_dir):
        sys.stderr.write(
            f"Warning: Sample_Annotation_Files directory not found in prep directory\n"
        )
    
    if not os.path.isdir(protein_dir):
        sys.stderr.write(
            f"Warning: Sample_Protein_Files directory not found in prep directory\n"
        )
    
    # Load GenBank files
    if os.path.isdir(genbank_dir):
        for filename in os.listdir(genbank_dir):
            if filename.endswith('.gbk') or filename.endswith('.genbank') or filename.endswith('.gbff'):
                sample_id = filename.rsplit('.', 1)[0]
                prep_data['sample_genbanks'][sample_id] = os.path.join(genbank_dir, filename)
    
    # If no GenBank files found, try BED files (for -emg mode)
    if not prep_data['sample_genbanks'] and os.path.isdir(bed_dir):
        sys.stderr.write(
            f"No GenBank files found. Using BED files from Genome_Processing/BEDs/\n"
        )
        for filename in os.listdir(bed_dir):
            if filename.endswith('.bed'):
                sample_id = filename.rsplit('.', 1)[0]
                prep_data['sample_genbanks'][sample_id] = os.path.join(bed_dir, filename)
    
    # Load protein FASTA files
    if os.path.isdir(protein_dir):
        for filename in os.listdir(protein_dir):
            if filename.endswith('.faa') or filename.endswith('.fasta'):
                sample_id = filename.rsplit('.', 1)[0]
                prep_data['sample_proteins'][sample_id] = os.path.join(protein_dir, filename)
    
    return prep_data


def load_run_results(run_dir: str) -> Dict:
    """
    Load results from bofasa run directory.
    
    Args:
        run_dir (str): Path to bofasa run directory
        
    Returns:
        dict: Dictionary containing run results
    """
    run_data = {
        'og_matrix': None,
        'dog_matrix': None,
        'og_data': {},
        'dog_data': {},
    }
    
    final_results_dir = os.path.join(run_dir, 'Final_Results')
    
    if not os.path.isdir(final_results_dir):
        sys.stderr.write(
            f"Error: Final_Results directory not found in {run_dir}\n"
        )
        sys.exit(1)
    
    # Load Protein Ortholog Groups
    og_file = os.path.join(final_results_dir, 'Protein_Ortholog_Groups.tsv')
    if os.path.isfile(og_file):
        run_data['og_matrix'] = og_file
    else:
        sys.stderr.write(
            f"Error: Protein_Ortholog_Groups.tsv not found in {final_results_dir}\n"
        )
        sys.exit(1)
    
    # Load Domain Ortholog Groups
    dog_file = os.path.join(final_results_dir, 'Domain_Ortholog_Groups.tsv')
    if os.path.isfile(dog_file):
        run_data['dog_matrix'] = dog_file
    
    return run_data


def parse_og_matrix(og_file: str, focal_og: str) -> Tuple[Dict, List]:
    """
    Parse the Protein Ortholog Groups matrix and extract focal OG information.
    
    Args:
        og_file (str): Path to Protein_Ortholog_Groups.tsv
        focal_og (str): Focal OG identifier
        
    Returns:
        tuple: (focal_og_data, sample_names)
            - focal_og_data: dict mapping sample -> list of protein IDs
            - sample_names: list of sample names in order
    """
    focal_og_data = {}
    sample_names = []
    og_found = False
    
    try:
        with open(og_file, 'r') as f:
            for i, line in enumerate(f):
                line = line.rstrip('\n')
                parts = line.split('\t')
                
                if i == 0:
                    # Header line: OG/Sample, sample1, sample2, ...
                    sample_names = parts[1:]
                    continue
                
                og_id = parts[0]
                
                if og_id == focal_og:
                    og_found = True
                    # Parse protein IDs for each sample
                    for j, sample_proteins in enumerate(parts[1:]):
                        sample_id = sample_names[j]
                        if sample_proteins.strip():
                            # Multiple proteins are comma-separated
                            proteins = [p.strip() for p in sample_proteins.split(',')]
                            focal_og_data[sample_id] = proteins
                    break
    
    except Exception as e:
        sys.stderr.write(f"Error parsing OG matrix: {str(e)}\n")
        sys.exit(1)
    
    if not og_found:
        sys.stderr.write(f"Error: Focal OG '{focal_og}' not found in OG matrix\n")
        sys.exit(1)
    
    return focal_og_data, sample_names


def parse_full_og_matrix(og_file: str) -> Tuple[Dict, List]:
    """
    Parse the complete Protein Ortholog Groups matrix.
    
    Args:
        og_file (str): Path to Protein_Ortholog_Groups.tsv
        
    Returns:
        tuple: (all_og_data, sample_names)
            - all_og_data: dict mapping OG -> dict(sample -> list of protein IDs)
            - sample_names: list of sample names in order
    """
    all_og_data = {}
    sample_names = []
    
    try:
        with open(og_file, 'r') as f:
            for i, line in enumerate(f):
                line = line.rstrip('\n')
                parts = line.split('\t')
                
                if i == 0:
                    sample_names = parts[1:]
                    continue
                
                og_id = parts[0]
                all_og_data[og_id] = {}
                
                for j, sample_proteins in enumerate(parts[1:]):
                    sample_id = sample_names[j]
                    if sample_proteins.strip():
                        proteins = [p.strip() for p in sample_proteins.split(',')]
                        all_og_data[og_id][sample_id] = proteins
    
    except Exception as e:
        sys.stderr.write(f"Error parsing full OG matrix: {str(e)}\n")
        sys.exit(1)
    
    return all_og_data, sample_names


def load_bed_annotations(bed_file: str) -> Dict:
    """
    Load gene annotations from a BED file.
    
    Args:
        bed_file (str): Path to BED file
        
    Returns:
        dict: Dictionary mapping protein ID -> gene information
    """
    gene_info = {}
    
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                
                contig = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                protein_id = parts[3]
                strand_symbol = parts[5]
                
                # Convert strand symbol to BioPython convention
                strand = 1 if strand_symbol == '+' else -1
                
                gene_info[protein_id] = {
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'contig': contig,
                    'product': 'Unknown',
                    'gene': '',
                }
    
    except Exception as e:
        sys.stderr.write(f"Warning: Error parsing BED file {bed_file}: {str(e)}\n")
    
    return gene_info


def load_genbank_annotations(genbank_file: str) -> Dict:
    """
    Load gene annotations from a GenBank file.
    
    Args:
        genbank_file (str): Path to GenBank file
        
    Returns:
        dict: Dictionary mapping protein ID -> gene information
    """
    if not BIOPYTHON_AVAILABLE:
        sys.stderr.write(
            "Warning: BioPython not available. Cannot parse GenBank files.\n"
        )
        return {}
    
    gene_info = {}
    
    try:
        for record in SeqIO.parse(genbank_file, 'genbank'):
            for feature in record.features:
                if feature.type == 'CDS':
                    # Extract protein ID
                    protein_id = None
                    if 'protein_id' in feature.qualifiers:
                        protein_id = feature.qualifiers['protein_id'][0]
                    elif 'locus_tag' in feature.qualifiers:
                        protein_id = feature.qualifiers['locus_tag'][0]
                    
                    if protein_id:
                        gene_info[protein_id] = {
                            'start': int(feature.location.start),
                            'end': int(feature.location.end),
                            'strand': feature.location.strand,
                            'contig': record.id,
                            'product': feature.qualifiers.get('product', ['Unknown'])[0],
                            'gene': feature.qualifiers.get('gene', [''])[0],
                        }
    
    except Exception as e:
        sys.stderr.write(f"Warning: Error parsing GenBank file {genbank_file}: {str(e)}\n")
    
    return gene_info


def load_annotation_file(annotation_file: str) -> Dict:
    """
    Load gene annotations from either a GenBank or BED file.
    
    Args:
        annotation_file (str): Path to annotation file (.gbk or .bed)
        
    Returns:
        dict: Dictionary mapping protein ID -> gene information
    """
    if annotation_file.endswith('.bed'):
        return load_bed_annotations(annotation_file)
    else:
        return load_genbank_annotations(annotation_file)


def extract_genomic_context(
    sample_genbanks: Dict[str, str],
    focal_og_data: Dict[str, List[str]],
    all_og_data: Dict[str, Dict[str, List[str]]],
    context_window: int,
    context_bp: Optional[int] = None,
) -> Dict:
    """
    Extract genomic context around focal OG instances.
    
    Args:
        sample_genbanks: Dict mapping sample ID -> GenBank file path
        focal_og_data: Dict mapping sample -> list of focal protein IDs
        all_og_data: Complete OG matrix data
        context_window: Number of genes to include up/downstream
        context_bp: Alternative: number of base pairs for context
        
    Returns:
        dict: Context information for visualization
    """
    context_data = {}
    
    # Build reverse mapping: protein_id -> OG
    protein_to_og = {}
    for og_id, og_samples in all_og_data.items():
        for sample_id, proteins in og_samples.items():
            for protein_id in proteins:
                protein_to_og[protein_id] = og_id
    
    # Process each sample
    for sample_id, focal_proteins in focal_og_data.items():
        if sample_id not in sample_genbanks:
            sys.stderr.write(f"Warning: No annotation file found for sample {sample_id}\n")
            continue
        
        # Load gene annotations (GenBank or BED)
        gene_info = load_annotation_file(sample_genbanks[sample_id])
        
        if not gene_info:
            sys.stderr.write(f"Warning: No gene annotations loaded for {sample_id}\n")
            continue
        
        # Organize genes by contig and position
        contig_genes = defaultdict(list)
        for protein_id, info in gene_info.items():
            contig_genes[info['contig']].append((protein_id, info))
        
        # Sort genes by position on each contig
        for contig in contig_genes:
            contig_genes[contig].sort(key=lambda x: x[1]['start'])
        
        # Extract context for each focal protein
        for focal_protein in focal_proteins:
            if focal_protein not in gene_info:
                continue
            
            focal_info = gene_info[focal_protein]
            contig = focal_info['contig']
            genes_on_contig = contig_genes[contig]
            
            # Find focal protein index
            focal_idx = None
            for idx, (pid, _) in enumerate(genes_on_contig):
                if pid == focal_protein:
                    focal_idx = idx
                    break
            
            if focal_idx is None:
                continue
            
            # Extract context genes
            context_genes = []
            
            if context_bp is not None:
                # Extract by base pair window
                focal_start = focal_info['start']
                focal_end = focal_info['end']
                focal_center = (focal_start + focal_end) // 2
                
                for protein_id, info in genes_on_contig:
                    gene_center = (info['start'] + info['end']) // 2
                    if abs(gene_center - focal_center) <= context_bp:
                        og_id = protein_to_og.get(protein_id, 'Unknown')
                        context_genes.append({
                            'protein_id': protein_id,
                            'og_id': og_id,
                            'start': info['start'],
                            'end': info['end'],
                            'strand': info['strand'],
                            'product': info['product'],
                            'is_focal': protein_id == focal_protein,
                        })
            else:
                # Extract by gene count window
                start_idx = max(0, focal_idx - context_window)
                end_idx = min(len(genes_on_contig), focal_idx + context_window + 1)
                
                for idx in range(start_idx, end_idx):
                    protein_id, info = genes_on_contig[idx]
                    og_id = protein_to_og.get(protein_id, 'Unknown')
                    context_genes.append({
                        'protein_id': protein_id,
                        'og_id': og_id,
                        'start': info['start'],
                        'end': info['end'],
                        'strand': info['strand'],
                        'product': info['product'],
                        'is_focal': protein_id == focal_protein,
                    })
            
            # Store context data - use unique key for each instance
            context_key = f"{sample_id}_{focal_protein}"
            context_data[context_key] = {
                'sample_id': sample_id,
                'focal_protein': focal_protein,
                'contig': contig,
                'genes': context_genes,
            }
    
    return context_data


def _generate_n_colors(n: int) -> List[str]:
    """
    Generate n visually distinct hex colours.

    Uses the curated tab20/tab20b/tab20c palettes for up to 60 OGs, then
    falls back to golden-ratio HSV stepping so that any number of OGs
    receive a unique colour (no wrapping/repeating).
    """
    import colorsys

    # 60 hand-picked visually distinct colours (tab20 + tab20b + tab20c)
    base = [
        '#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd',
        '#8c564b','#e377c2','#7f7f7f','#bcbd22','#17becf',
        '#aec7e8','#ffbb78','#98df8a','#ff9896','#c5b0d5',
        '#c49c94','#f7b6d2','#c7c7c7','#dbdb8d','#9edae5',
        '#393b79','#5254a3','#6b6ecf','#9c9ede','#637939',
        '#8ca252','#b5cf6b','#cedb9c','#8c6d31','#bd9e39',
        '#e7ba52','#e7cb94','#843c39','#ad494a','#d6616b',
        '#e7969c','#7b4173','#a55194','#ce6dbd','#de9ed6',
        '#3182bd','#6baed6','#9ecae1','#c6dbef','#e6550d',
        '#fd8d3c','#fdae6b','#fdd0a2','#31a354','#74c476',
        '#a1d99b','#c7e9c0','#756bb1','#9e9ac8','#bcbddc',
        '#dadaeb','#636363','#969696','#bdbdbd','#d9d9d9',
    ]

    if n <= len(base):
        # Spread selection across the palette for maximum contrast
        step = len(base) / n
        return [base[int(i * step)] for i in range(n)]

    # More OGs than palette entries – golden-ratio HSV wheel (no repeats)
    golden = 0.618033988749895
    h = 0.0
    extra = []
    for _ in range(n - len(base)):
        r, g, b = colorsys.hsv_to_rgb(h % 1.0, 0.65, 0.82)
        extra.append('#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255)))
        h += golden
    return base + extra


def assign_og_colors(context_data: Dict, color_scheme: str = 'tab20') -> Dict:
    """
    Assign colors to OGs for visualization.
    
    Args:
        context_data: Context information from extract_genomic_context
        color_scheme: Color scheme name (currently ignored; kept for API compat)
        
    Returns:
        dict: Mapping of OG ID -> color
    """
    all_ogs = set()
    focal_og = None

    for context in context_data.values():
        for gene in context['genes']:
            all_ogs.add(gene['og_id'])
            if gene['is_focal']:
                focal_og = gene['og_id']

    og_colors = {}

    if focal_og:
        og_colors[focal_og] = '#FF0000'
        all_ogs.discard(focal_og)

    other_ogs = sorted(all_ogs)
    palette   = _generate_n_colors(len(other_ogs))
    for og_id, color in zip(other_ogs, palette):
        og_colors[og_id] = color

    return og_colors


def get_parent_sample_name(sample_id: str) -> str:
    """
    Extract parent sample name from MGE-specific sample IDs.
    
    When bofasa prep is run with -emg flag, phages and plasmids get separate
    names like 'Staphylococcus_aureus_phage_3' or 'Staphylococcus_epidermidis_plasmid_1'.
    This function extracts the parent sample name ('Staphylococcus_aureus').
    
    Args:
        sample_id: Sample identifier (may include _phage_X or _plasmid_X suffix)
        
    Returns:
        str: Parent sample name
    """
    pattern = r'_(phage|plasmid)_\d+$'
    parent = re.sub(pattern, '', sample_id)
    return parent


def get_unique_sample_ids(context_data: Dict) -> Dict[str, str]:
    """
    Generate unique sample IDs for cases where samples have multiple instances.
    
    For samples with multiple instances of the focal OG, this creates
    unique identifiers like "sample1_instance1", "sample1_instance2".
    Handles MGE-specific sample names by mapping to parent samples.
    
    Args:
        context_data: Context information from extract_genomic_context
        
    Returns:
        dict: Mapping of context_key -> unique_sample_id
    """
    # Count instances per parent sample
    sample_counts = defaultdict(int)
    
    for context_key in sorted(context_data.keys()):
        sample_id = context_data[context_key]['sample_id']
        parent_sample = get_parent_sample_name(sample_id)
        sample_counts[parent_sample] += 1
    
    # Assign unique IDs
    unique_ids = {}
    instance_counters = defaultdict(int)
    
    for context_key in sorted(context_data.keys()):
        sample_id = context_data[context_key]['sample_id']
        parent_sample = get_parent_sample_name(sample_id)
        
        if sample_counts[parent_sample] > 1:
            # Multiple instances - add suffix
            instance_counters[parent_sample] += 1
            unique_ids[context_key] = f"{parent_sample}_instance{instance_counters[parent_sample]}"
        else:
            # Single instance - keep parent sample name
            unique_ids[context_key] = parent_sample
    
    return unique_ids


def modify_tree_for_paralogs(tree_file: str, context_data: Dict, output_file: str) -> Dict[str, str]:
    """
    Modify a phylogenetic tree to add duplicate leaves for samples with multiple instances.
    
    Similar to lsaBGC-See's modifyPhylogenyForSamplesWithMultipleBGCs functionality.
    
    Args:
        tree_file: Input Newick tree file
        context_data: Context information
        output_file: Output modified tree file
        
    Returns:
        dict: Mapping of context_key -> unique_sample_id_in_tree
    """
    if not BIOPYTHON_AVAILABLE:
        sys.stderr.write(
            "Warning: BioPython not available. Cannot modify tree for paralogs.\n"
        )
        return get_unique_sample_ids(context_data)
    
    try:
        # Read the tree
        tree = Phylo.read(tree_file, 'newick')
        
        # Get unique sample IDs
        unique_ids = get_unique_sample_ids(context_data)
        
        # Find which parent samples have multiple instances
        parent_sample_counts = defaultdict(int)
        for context_key in context_data.keys():
            sample_id = context_data[context_key]['sample_id']
            parent_sample = get_parent_sample_name(sample_id)
            parent_sample_counts[parent_sample] += 1
        
        samples_to_duplicate = {s for s, c in parent_sample_counts.items() if c > 1}
        
        # Modify tree by duplicating terminals for multi-instance samples
        for terminal in list(tree.get_terminals()):
            parent_sample = terminal.name
            
            if parent_sample in samples_to_duplicate:
                # Find parent of this terminal
                parent_clade = None
                for clade in tree.find_clades():
                    if terminal in clade.clades:
                        parent_clade = clade
                        break
                
                if parent_clade is None:
                    continue
                
                # Get branch length
                branch_length = terminal.branch_length if terminal.branch_length else 0.0
                
                # Remove original terminal
                parent_clade.clades.remove(terminal)
                
                # Create new terminals for each instance
                num_instances = parent_sample_counts[parent_sample]
                
                # Add duplicated terminals with slightly modified branch lengths
                for i in range(1, num_instances + 1):
                    new_terminal = Clade(
                        branch_length=branch_length + (i * 0.00001),
                        name=f"{parent_sample}_instance{i}"
                    )
                    parent_clade.clades.append(new_terminal)
        
        # Write modified tree, then midpoint-root it via ete3
        Phylo.write(tree, output_file, 'newick')
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', SyntaxWarning)
                from ete3 import Tree as EteTree
            _et = EteTree(output_file)
            _et.set_outgroup(_et.get_midpoint_outgroup())
            _et.write(outfile=output_file, format=1)
        except Exception:
            pass  # Non-fatal; tree is still usable unrooted

        sys.stderr.write(
            f"  Modified tree written with duplicated leaves for {len(samples_to_duplicate)} samples\n"
        )

        return unique_ids
    
    except Exception as e:
        sys.stderr.write(f"Warning: Error modifying tree: {str(e)}\n")
        sys.stderr.write("Proceeding with unique sample IDs but unmodified tree.\n")
        return get_unique_sample_ids(context_data)


def write_context_heatmap(context_data: Dict, og_colors: Dict, focal_og: str,
                          output_file: str, unique_sample_ids: Optional[Dict[str, str]] = None,
                          tree_file: Optional[str] = None):
    """
    Write a heatmap PNG showing OG presence per sample, optionally aligned to a tree.

    Columns = unique OGs ordered by their median relative position to the focal gene.
    Cells are coloured by OG colour when present, light grey when absent.
    The focal OG column is outlined in black for emphasis.

    Args:
        context_data:       Context information from extract_genomic_context.
        og_colors:          OG -> hex colour mapping.
        focal_og:           The focal OG ID to centre columns around.
        output_file:        Output PNG path.
        unique_sample_ids:  Optional context_key -> unique sample ID mapping.
        tree_file:          Optional Newick tree path; samples are pruned to
                            those present in context_data.
    """
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', SyntaxWarning)
            from ete3 import Tree as EteTree, TreeStyle, RectFace, TextFace, NodeStyle, faces
    except ImportError:
        sys.stderr.write("Warning: ete3 not available. Skipping heatmap PNG.\n")
        return

    if unique_sample_ids is None:
        unique_sample_ids = get_unique_sample_ids(context_data)

    # ── build per-sample OG presence + relative-position index ─────────────
    # Each context_key must already have a unique molecule name (guaranteed by
    # get_unique_sample_ids / modify_tree_for_paralogs upstream).  Plain dict
    # assignment is intentional: a collision here means the ID pipeline failed.
    sample_og_presence: Dict[str, set] = {}
    og_rel_positions: Dict[str, List[int]] = defaultdict(list)

    for context_key, context in context_data.items():
        sample_id = unique_sample_ids.get(context_key, context['sample_id'])
        genes = context['genes']

        focal_idx = next(
            (i for i, g in enumerate(genes) if g['og_id'] == focal_og), None
        )

        # Strand sign: if the focal gene is on the reverse strand, flip the
        # direction so that negative offsets always mean "upstream" and positive
        # offsets always mean "downstream" in the biological sense.
        if focal_idx is not None:
            focal_strand = genes[focal_idx].get('strand', 1)
            strand_sign  = 1 if focal_strand == 1 else -1
        else:
            strand_sign = 1

        og_set: set = set()
        for i, gene in enumerate(genes):
            og_id = gene['og_id']
            og_set.add(og_id)
            if focal_idx is not None:
                og_rel_positions[og_id].append(strand_sign * (i - focal_idx))

        if sample_id in sample_og_presence:
            sys.stderr.write(
                f"Warning: duplicate molecule name '{sample_id}' in heatmap data — "
                f"unique ID assignment may have failed.\n"
            )
        sample_og_presence[sample_id] = og_set

    # ── column order: by median relative position ───────────────────────────
    def _median(vals):
        s = sorted(vals)
        return s[len(s) // 2]

    all_ogs = sorted(
        {og for s in sample_og_presence.values() for og in s},
        key=lambda og: _median(og_rel_positions.get(og, [0]))
    )

    # ── build / prune tree ──────────────────────────────────────────────────
    samples_in_data = set(sample_og_presence.keys())

    if tree_file and os.path.isfile(tree_file):
        try:
            etree = EteTree(tree_file)
            tree_tips = {leaf.name for leaf in etree.get_leaves()}
            keep = tree_tips & samples_in_data
            if not keep:
                raise ValueError("No tree tips match any sample ID.")
            etree.prune(list(keep), preserve_branch_length=True)
            etree.set_outgroup(etree.get_midpoint_outgroup())
        except Exception as e:
            sys.stderr.write(f"Warning: could not use tree for heatmap ({e}). "
                             f"Using star topology.\n")
            etree = None
    else:
        etree = None

    if etree is None:
        etree = EteTree()
        etree.name = "root"
        for sid in sorted(samples_in_data):
            etree.add_child(name=sid)
    else:
        # Safety net: any sample present in the data but absent from the
        # (possibly pruned) tree gets appended as an extra leaf so every
        # row always appears in the heatmap.
        tree_leaf_names = {leaf.name for leaf in etree.get_leaves()}
        for sid in sorted(samples_in_data - tree_leaf_names):
            sys.stderr.write(f"  Note: '{sid}' not in tree, appended as extra row.\n")
            etree.add_child(name=sid)

    # ── attach heatmap faces to each leaf ───────────────────────────────────
    CELL_W, CELL_H = 16, 16
    ABSENT_COLOR   = '#E8E8E8'
    FOCAL_BORDER   = '#000000'

    for leaf in etree.get_leaves():
        nstyle = NodeStyle()
        nstyle['size'] = 0
        leaf.set_style(nstyle)

        og_set = sample_og_presence.get(leaf.name, set())

        for col_idx, og_id in enumerate(all_ogs):
            is_focal = (og_id == focal_og)
            if og_id in og_set:
                bg = og_colors.get(og_id, '#CCCCCC')
                fg = FOCAL_BORDER if is_focal else bg
            else:
                bg = ABSENT_COLOR
                fg = FOCAL_BORDER if is_focal else ABSENT_COLOR

            rect = RectFace(CELL_W, CELL_H, bgcolor=bg, fgcolor=fg)
            rect.margin_left  = 1
            rect.margin_right = 1
            leaf.add_face(rect, column=col_idx, position='aligned')

        name_face = TextFace(f'  {leaf.name}', fsize=8)
        leaf.add_face(name_face, column=len(all_ogs), position='aligned')

    # ── tree style + column headers ─────────────────────────────────────────
    ts = TreeStyle()
    ts.show_leaf_name          = False
    ts.mode                    = 'r'
    ts.branch_vertical_margin  = 2
    ts.show_scale              = False

    # Auto-scale: fit the tree into ~180 px regardless of branch length units.
    leaves = etree.get_leaves()
    max_dist = max((etree.get_distance(lf) for lf in leaves), default=0)
    ts.scale = int(180 / max_dist) if max_dist > 0 else 60

    # Column headers: rotated OG ID label only (no coloured square row).
    # rotation=270 renders text bottom-to-top so labels don't overlap.
    for col_idx, og_id in enumerate(all_ogs):
        lbl = TextFace(og_id, fsize=7)
        lbl.rotation     = 270
        lbl.margin_left  = 1
        lbl.margin_right = 1
        lbl.hz_align     = 1
        ts.aligned_header.add_face(lbl, column=col_idx)

    # Right margin: wide enough for the longest sample name label.
    max_label_len = max((len(sid) for sid in sample_og_presence), default=20)
    ts.margin_right = max_label_len * 6 + 20

    # Top margin: rotated OG labels are rendered horizontally by Qt then
    # rotated, so their effective height equals their text width (~5 px per
    # character at fsize=7). Give the header room to show them in full.
    max_og_name_len = max((len(og) for og in all_ogs), default=10)
    ts.margin_top = max_og_name_len * 5 + 20

    try:
        # SVG is vector-based so text and lines are always sharp regardless of
        # zoom or print size.  The w parameter still controls the layout width
        # that ete3 uses internally when placing faces.
        img_w = max(800, len(all_ogs) * (CELL_W + 2) + 500)
        etree.render(output_file, tree_style=ts, w=img_w, units='px')
        sys.stderr.write(f"  Wrote heatmap: {output_file}\n")

        # Also render PNG for smaller file size
        png_file = output_file.replace('.svg', '.png')
        etree.render(png_file, tree_style=ts, w=img_w, units='px')
        sys.stderr.write(f"  Wrote heatmap PNG: {png_file}\n")
    except Exception as e:
        sys.stderr.write(f"Warning: ete3 render failed ({e}). "
                         f"Is a display / Qt available?\n")


def write_neighborhood_genbanks(context_data: Dict, focal_og: str,
                                 output_dir: str, prep_dir: str,
                                 unique_sample_ids: Optional[Dict[str, str]] = None):
    """
    Write one GenBank file per focal-gene neighborhood.

    Each record contains:
    - source feature spanning the full extracted region
    - CDS features with protein translations, OG annotation in /note,
      and /focal_gene qualifier on the focal CDS
    - The nucleotide sequence of the extracted contig region as the record
      sequence (origin section in the GenBank file)

    Files are written to <output_dir>/Neighborhood_GenBanks/.
    """
    if not BIOPYTHON_AVAILABLE:
        sys.stderr.write("Warning: BioPython not available — skipping neighborhood GenBanks.\n")
        return

    from Bio.SeqRecord import SeqRecord as BioSeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Seq import Seq

    if unique_sample_ids is None:
        unique_sample_ids = get_unique_sample_ids(context_data)

    genomes_dir   = os.path.join(prep_dir, 'Genome_Processing', 'Genomes')
    proteomes_dir = _proteomes_dir(prep_dir)
    gbk_dir       = os.path.join(output_dir, 'Neighborhood_GenBanks')
    os.makedirs(gbk_dir, exist_ok=True)

    # Lazy genome / proteome caches keyed by original sample_id
    genome_cache:  Dict[str, Dict[str, str]] = {}
    protein_cache: Dict[str, Dict[str, str]] = {}

    def _load_fasta(path: str) -> Dict[str, str]:
        seqs: Dict[str, str] = {}
        cur_id = None
        buf: List[str] = []
        with open(path) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    if cur_id:
                        seqs[cur_id] = ''.join(buf)
                    cur_id = line[1:].split()[0]
                    buf = []
                elif cur_id:
                    buf.append(line.strip())
            if cur_id:
                seqs[cur_id] = ''.join(buf)
        return seqs

    def _ensure_genome(sample_id: str) -> Dict[str, str]:
        if sample_id not in genome_cache:
            fna = os.path.join(genomes_dir, f'{sample_id}.fna')
            genome_cache[sample_id] = _load_fasta(fna) if os.path.isfile(fna) else {}
        return genome_cache[sample_id]

    def _ensure_proteins(sample_id: str) -> Dict[str, str]:
        if sample_id not in protein_cache:
            faa = os.path.join(proteomes_dir, f'{sample_id}.faa')
            seqs = _load_fasta(faa) if os.path.isfile(faa) else {}
            # Strip trailing stop codon symbol
            protein_cache[sample_id] = {k: v.rstrip('*') for k, v in seqs.items()}
        return protein_cache[sample_id]

    written = 0
    for context_key, context in sorted(context_data.items()):
        sample_id      = unique_sample_ids.get(context_key, context['sample_id'])
        orig_sample_id = context['sample_id']
        genes          = context['genes']
        contig         = context['contig']

        if not genes:
            continue

        region_start = min(g['start'] for g in genes)
        region_end   = max(g['end']   for g in genes)

        # Nucleotide sequence for the region
        contig_seqs = _ensure_genome(orig_sample_id)
        contig_seq  = contig_seqs.get(contig, '')
        nt_seq      = contig_seq[region_start:region_end] if contig_seq else 'N' * (region_end - region_start)

        prot_seqs = _ensure_proteins(orig_sample_id)

        # Build SeqRecord
        rec_id = re.sub(r'[^\w\-]', '_', f'{sample_id}_{focal_og}')
        record = BioSeqRecord(
            Seq(nt_seq),
            id=rec_id,
            name=rec_id[:16],
            description=(
                f'Genomic neighborhood of {focal_og} in {sample_id} '
                f'| contig={contig} region={region_start+1}..{region_end}'
            ),
        )
        record.annotations['molecule_type'] = 'DNA'

        # source feature
        record.features.append(SeqFeature(
            FeatureLocation(0, len(nt_seq), strand=1),
            type='source',
            qualifiers={
                'organism':   [sample_id],
                'mol_type':   ['genomic DNA'],
                'note':       [f'contig={contig} abs_start={region_start+1} abs_end={region_end}'],
            },
        ))

        # CDS features
        for gene in genes:
            feat_start = gene['start'] - region_start
            feat_end   = gene['end']   - region_start

            qualifiers: Dict[str, List[str]] = {
                'locus_tag':  [gene['protein_id']],
                'protein_id': [gene['protein_id']],
                'product':    [gene['product']],
                'note':       [f'OG={gene["og_id"]}'],
            }
            if gene.get('is_focal'):
                qualifiers['note'].append(f'focal_og={focal_og}')

            prot_seq = prot_seqs.get(gene['protein_id'], '')
            if prot_seq:
                qualifiers['translation'] = [prot_seq]

            record.features.append(SeqFeature(
                FeatureLocation(feat_start, feat_end, strand=gene['strand']),
                type='CDS',
                qualifiers=qualifiers,
            ))

        gbk_file = os.path.join(gbk_dir, f'{rec_id}.gbk')
        with open(gbk_file, 'w') as out_fh:
            SeqIO.write(record, out_fh, 'genbank')
        written += 1

    sys.stderr.write(f"  Wrote {written} neighborhood GenBank files to: {gbk_dir}\n")


def write_gggenes_output(context_data: Dict, og_colors: Dict, output_file: str, 
                        unique_sample_ids: Optional[Dict[str, str]] = None):
    """
    Write output for gggenes R package visualization.
    
    Args:
        context_data: Context information
        og_colors: OG color assignments
        output_file: Output file path
        unique_sample_ids: Optional mapping of context_key -> unique sample ID
    """
    if unique_sample_ids is None:
        unique_sample_ids = get_unique_sample_ids(context_data)
    
    try:
        with open(output_file, 'w') as f:
            # Header
            f.write('molecule\tgene\tstart\tend\tstrand\tog\tcolor\tproduct\n')
            
            # Write each gene
            for context_key, context in sorted(context_data.items()):
                # Use unique sample ID for cases with multiple instances
                sample_id = unique_sample_ids.get(context_key, context['sample_id'])
                
                for gene in context['genes']:
                    og_id = gene['og_id']
                    color = og_colors.get(og_id, '#CCCCCC')
                    strand_str = '1' if gene['strand'] == 1 else '0'
                    product = gene['product'].replace('\t', ' ').replace('\n', ' ').replace('\r', ' ')
                    
                    f.write(
                        f"{sample_id}\t{gene['protein_id']}\t{gene['start']}\t"
                        f"{gene['end']}\t{strand_str}\t{og_id}\t{color}\t"
                        f"{product}\n"
                    )
    
    except Exception as e:
        sys.stderr.write(f"Error writing gggenes output: {str(e)}\n")
        sys.exit(1)


def write_itol_output(context_data: Dict, og_colors: Dict, output_file: str,
                     unique_sample_ids: Optional[Dict[str, str]] = None):
    """
    Write iTOL annotation track file.
    
    Args:
        context_data: Context information
        og_colors: OG color assignments
        output_file: Output file path
        unique_sample_ids: Optional mapping of context_key -> unique sample ID
    """
    if unique_sample_ids is None:
        unique_sample_ids = get_unique_sample_ids(context_data)
    
    try:
        with open(output_file, 'w') as f:
            f.write('DATASET_DOMAINS\n')
            f.write('SEPARATOR TAB\n')
            f.write('DATASET_LABEL\tOG_Context\n')
            f.write('COLOR\t#ff0000\n')
            f.write('DATA\n')
            
            # Write domain annotations for each sample
            for context_key, context in sorted(context_data.items()):
                # Use unique sample ID for cases with multiple instances
                sample_id = unique_sample_ids.get(context_key, context['sample_id'])
                
                # Calculate domain positions relative to focal gene
                domains = []
                for gene in context['genes']:
                    og_id = gene['og_id']
                    color = og_colors.get(og_id, '#CCCCCC')
                    domains.append(f"{gene['start']}-{gene['end']}|{color}|{og_id}")
                
                f.write(f"{sample_id}\t{'|'.join(domains)}\n")
    
    except Exception as e:
        sys.stderr.write(f"Error writing iTOL output: {str(e)}\n")
        sys.exit(1)


def _proteomes_dir(prep_dir: str) -> str:
    """Return the path to per-sample proteome FASTAs inside a prep directory."""
    return os.path.join(prep_dir, 'Genome_Processing', 'Proteomes')


def _aln_length(fasta_file: str) -> int:
    """Return the aligned length (columns) of the first sequence in a FASTA."""
    with open(fasta_file) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith('>'):
                return len(line.replace('-', '').replace('*', ''))
    return 0


def _extract_og_sequences(og_id: str, og_matrix_file: str,
                           proteomes_dir: str, out_fasta: str,
                           label_as_sample: bool = False) -> int:
    """
    Write a FASTA of all proteins belonging to og_id to out_fasta.

    Args:
        og_id:           OG identifier to extract.
        og_matrix_file:  Path to Protein_Ortholog_Groups.tsv.
        proteomes_dir:   Directory containing per-sample .faa files.
        out_fasta:       Output FASTA path.
        label_as_sample: If True, sequence headers are just >{sample_id}
                         (useful for species-tree concatenation); otherwise
                         headers are >{sample_id}|{protein_id}.

    Returns:
        Number of sequences written.
    """
    # Parse OG row
    sample_to_proteins: Dict[str, List[str]] = {}
    with open(og_matrix_file) as fh:
        samples: List[str] = []
        for i, line in enumerate(fh):
            parts = line.rstrip('\n').split('\t')
            if i == 0:
                samples = parts[1:]
                continue
            if parts[0] != og_id:
                continue
            for j, cell in enumerate(parts[1:]):
                cell = cell.strip()
                if cell:
                    sample_to_proteins[samples[j]] = [
                        p.strip() for p in cell.split(',') if p.strip()
                    ]
            break

    if not sample_to_proteins:
        return 0

    # Build set of wanted protein IDs for fast lookup
    wanted: Dict[str, str] = {}  # protein_id -> sample_id
    for sample_id, prots in sample_to_proteins.items():
        for p in prots:
            wanted[p] = sample_id

    written = 0
    with open(out_fasta, 'w') as out_fh:
        for faa_name in sorted(os.listdir(proteomes_dir)):
            if not faa_name.endswith('.faa'):
                continue
            sample_id = faa_name[:-4]
            header = None
            buf: List[str] = []
            with open(os.path.join(proteomes_dir, faa_name)) as pfh:
                for line in pfh:
                    if line.startswith('>'):
                        if header is not None and buf:
                            out_fh.write(header + '\n')
                            out_fh.writelines(buf)
                            buf = []
                            header = None
                        prot_id = line[1:].split()[0]
                        if prot_id in wanted:
                            if label_as_sample:
                                header = f'>{wanted[prot_id]}'
                            else:
                                header = f'>{wanted[prot_id]}|{prot_id}'
                            written += 1
                    elif header is not None:
                        buf.append(line)
                if header is not None and buf:
                    out_fh.write(header + '\n')
                    out_fh.writelines(buf)
    return written


def _run_muscle_trimal(in_faa: str, aligned_faa: str, trimmed_faa: str,
                        trim_log: str, cpus: int = 1) -> Tuple[str, int]:
    """
    Align with muscle super5, trim with trimAl -automated1.

    Returns (path_to_use, alignment_length_aa).
    If trimmed length < 10 aa, reverts to the untrimmed file and logs a note.
    """
    subprocess.run(
        ['muscle', '-super5', in_faa, '-output', aligned_faa,
         '-threads', str(cpus)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

    try:
        subprocess.run(
            ['trimal', '-in', aligned_faa, '-out', trimmed_faa, '-automated1'],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        aln_len = _aln_length(trimmed_faa)
    except (subprocess.CalledProcessError, FileNotFoundError):
        aln_len = 0

    if aln_len >= 10:
        return trimmed_faa, aln_len

    msg = (f"Trimmed alignment is only {aln_len} aa (<10); "
           f"using untrimmed alignment instead.\n")
    sys.stderr.write(f"  Note: {msg}")
    with open(trim_log, 'w') as lf:
        lf.write(msg)
    return aligned_faa, _aln_length(aligned_faa)


def build_gene_tree(focal_og: str, run_dir: str, prep_dir: str,
                    output_dir: str, cpus: int = 1) -> Optional[str]:
    """
    Build a gene tree for focal_og using muscle super5 + trimAl + FastTree2.

    Sequences are collected from the prep proteomes directory using the
    membership listed in Protein_Ortholog_Groups.tsv.

    Returns the path to the Newick tree, or None on failure.
    """
    og_matrix = os.path.join(run_dir, 'Final_Results', 'Protein_Ortholog_Groups.tsv')
    proteomes  = _proteomes_dir(prep_dir)
    tree_dir   = os.path.join(output_dir, 'gene_tree')
    os.makedirs(tree_dir, exist_ok=True)

    sys.stderr.write(f"\nBuilding gene tree for {focal_og}...\n")

    raw_faa     = os.path.join(tree_dir, f'{focal_og}.faa')
    aligned_faa = os.path.join(tree_dir, f'{focal_og}_aligned.faa')
    trimmed_faa = os.path.join(tree_dir, f'{focal_og}_trimmed.faa')
    trim_log    = os.path.join(tree_dir, f'{focal_og}_trimming.log')
    tree_nwk    = os.path.join(tree_dir, f'{focal_og}_gene_tree.nwk')

    n_seqs = _extract_og_sequences(focal_og, og_matrix, proteomes, raw_faa)
    if n_seqs < 3:
        sys.stderr.write(f"  Only {n_seqs} sequences found — skipping gene tree.\n")
        return None
    sys.stderr.write(f"  Extracted {n_seqs} sequences.\n")

    try:
        use_faa, aln_len = _run_muscle_trimal(
            raw_faa, aligned_faa, trimmed_faa, trim_log, cpus
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.stderr.write(f"  Alignment failed: {e}\n")
        return None

    sys.stderr.write(f"  Alignment length used: {aln_len} aa.\n")

    try:
        with open(tree_nwk, 'w') as tree_out:
            subprocess.run(
                ['FastTree', '-lg', '-quiet', use_faa],
                check=True, stdout=tree_out, stderr=subprocess.DEVNULL,
            )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.stderr.write(f"  FastTree failed: {e}\n")
        return None

    sys.stderr.write(f"  Gene tree written: {tree_nwk}\n")
    return tree_nwk


def _subsample_and_fasttree(concat_faa: str, tree_nwk: str,
                             max_sites: int = 10000) -> Optional[str]:
    """
    Optionally subsample columns of a concatenated protein FASTA to max_sites,
    then run FastTree2.  Returns the tree path on success, None on failure.
    """
    import random

    # Read sequences
    seqs: Dict[str, str] = {}
    cur = None
    with open(concat_faa) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                cur = line[1:].strip()
                seqs[cur] = ''
            elif cur:
                seqs[cur] += line.replace('*', '-')

    if not seqs:
        return None

    total_sites = len(next(iter(seqs.values())))
    sys.stderr.write(f"  Concatenated alignment: {total_sites} sites.\n")

    if total_sites > max_sites:
        rng_cols = sorted(random.sample(range(total_sites), max_sites))
        seqs = {s: ''.join(seq[c] for c in rng_cols) for s, seq in seqs.items()}
        sys.stderr.write(f"  Randomly selected {max_sites} sites.\n")
        # Overwrite with subsampled version
        with open(concat_faa, 'w') as out_fh:
            for s, seq in seqs.items():
                out_fh.write(f'>{s}\n{seq}\n')

    try:
        with open(tree_nwk, 'w') as tree_out:
            subprocess.run(
                ['FastTree', '-lg', '-quiet', concat_faa],
                check=True, stdout=tree_out, stderr=subprocess.DEVNULL,
            )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        sys.stderr.write(f"  FastTree failed: {e}\n")
        return None

    sys.stderr.write(f"  Species tree written: {tree_nwk}\n")
    return tree_nwk


def build_species_tree_from_scc(run_dir: str, prep_dir: str,
                                 output_dir: str, cpus: int = 1,
                                 max_sites: int = 10000) -> Optional[str]:
    """
    Build a species tree from a core-genome alignment.

    Priority:
      1. If bofasa run was invoked with --core_genome, the pre-built
         concatenated alignment is at Final_Results/Core_Genome_Alignment.faa.
         Use that directly (apply max_sites subsampling then FastTree2).
      2. Otherwise, build the alignment on the fly from single-copy core (SCC)
         OGs in Protein_Ortholog_Groups.tsv:
           a. Identify SCC OGs (one protein per sample, no blanks/commas).
           b. For each: extract → muscle super5 → trimAl.
           c. Concatenate horizontally; subsample to max_sites if needed.
           d. FastTree2.

    Returns the path to the Newick tree, or None on failure.
    """
    import random

    og_matrix  = os.path.join(run_dir, 'Final_Results', 'Protein_Ortholog_Groups.tsv')
    proteomes  = _proteomes_dir(prep_dir)
    tree_dir   = os.path.join(output_dir, 'species_tree_scc')
    os.makedirs(tree_dir, exist_ok=True)

    tree_nwk = os.path.join(tree_dir, 'species_tree_scc.nwk')

    # ── Priority 1: use pre-built core genome alignment if available ────────
    prebuilt_aln = os.path.join(run_dir, 'Final_Results', 'Core_Genome_Alignment.faa')
    if os.path.isfile(prebuilt_aln):
        sys.stderr.write(
            "\nBuilding species tree from pre-built core genome alignment "
            f"({prebuilt_aln})...\n"
        )
        import shutil
        working_copy = os.path.join(tree_dir, 'Core_Genome_Alignment.faa')
        shutil.copy2(prebuilt_aln, working_copy)
        return _subsample_and_fasttree(working_copy, tree_nwk, max_sites)

    sys.stderr.write("\nBuilding species tree from single-copy core OGs...\n")

    # ── 1. Parse OG matrix ──────────────────────────────────────────────────
    all_samples: List[str] = []
    scc_ogs: List[str] = []
    og_rows: Dict[str, Dict[str, str]] = {}  # og_id -> {sample: protein}

    with open(og_matrix) as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip('\n').split('\t')
            if i == 0:
                all_samples = parts[1:]
                continue
            og_id = parts[0]
            cells = parts[1:]
            # SCC: every sample has exactly one protein
            if len(cells) < len(all_samples):
                cells += [''] * (len(all_samples) - len(cells))
            if all(c.strip() and ',' not in c for c in cells):
                scc_ogs.append(og_id)
                og_rows[og_id] = {
                    all_samples[j]: cells[j].strip()
                    for j in range(len(all_samples))
                }

    sys.stderr.write(f"  Found {len(scc_ogs)} single-copy core OGs "
                     f"across {len(all_samples)} samples.\n")

    if len(scc_ogs) < 1:
        sys.stderr.write("  No SCC OGs found — cannot build species tree.\n")
        return None

    # ── 2–3. Align each SCC OG and concatenate ──────────────────────────────
    concat: Dict[str, List[str]] = {s: [] for s in all_samples}
    used_ogs = 0

    for og_id in scc_ogs:
        raw_faa     = os.path.join(tree_dir, f'{og_id}.faa')
        aligned_faa = os.path.join(tree_dir, f'{og_id}_aligned.faa')
        trimmed_faa = os.path.join(tree_dir, f'{og_id}_trimmed.faa')
        trim_log    = os.path.join(tree_dir, f'{og_id}_trimming.log')

        n = _extract_og_sequences(og_id, og_matrix, proteomes, raw_faa,
                                   label_as_sample=True)
        if n < len(all_samples):
            continue

        try:
            use_faa, aln_len = _run_muscle_trimal(
                raw_faa, aligned_faa, trimmed_faa, trim_log, cpus
            )
        except (subprocess.CalledProcessError, FileNotFoundError):
            continue

        if aln_len < 1:
            continue

        # Read aligned sequences keyed by sample name
        seqs: Dict[str, str] = {}
        cur_label = None
        with open(use_faa) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('>'):
                    cur_label = line[1:].strip()
                    seqs[cur_label] = ''
                elif cur_label:
                    seqs[cur_label] += line.replace('*', '-')

        # Only keep OG if every sample has a sequence in the alignment
        if not all(s in seqs for s in all_samples):
            continue

        for s in all_samples:
            concat[s].append(seqs[s])
        used_ogs += 1

    sys.stderr.write(f"  Used {used_ogs} SCC OGs in concatenated alignment.\n")

    if used_ogs == 0:
        sys.stderr.write("  No usable SCC alignments — cannot build species tree.\n")
        return None

    # ── 4–5. Write concatenated FASTA then subsample + FastTree ─────────────
    concat_faa = os.path.join(tree_dir, 'scc_concat.faa')
    with open(concat_faa, 'w') as out_fh:
        for s in all_samples:
            out_fh.write(f'>{s}\n{"".join(concat[s])}\n')

    return _subsample_and_fasttree(concat_faa, tree_nwk, max_sites)


def visualize_og_context():
    """
    Main function to visualize OG context.
    """
    # Get version
    version = utils.get_version()

    if len(sys.argv) > 1 and ('-v' in set(sys.argv) or '--version' in set(sys.argv)):
        sys.stdout.write(version + '\n')
        sys.exit(0)

    # Parse arguments
    myargs = create_parser()

    prep_dir = os.path.abspath(myargs.prep_dir)
    run_dir = os.path.abspath(myargs.run_dir)
    focal_og = myargs.focal_og
    output_dir = os.path.abspath(myargs.output_dir)
    species_tree = myargs.species_tree
    context_window = myargs.context_window
    context_bp = myargs.context_bp
    single_copy_only = myargs.single_copy_only
    min_length = myargs.min_length
    sample_list = myargs.sample_list
    cpus = myargs.cpus

    # Validate inputs
    if not os.path.isdir(prep_dir):
        sys.stderr.write(f"Error: Prep directory not found: {prep_dir}\n")
        sys.exit(1)
    
    if not os.path.isdir(run_dir):
        sys.stderr.write(f"Error: Run directory not found: {run_dir}\n")
        sys.exit(1)
    
    # Create output directory
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
            sys.stderr.write(f"Created output directory: {output_dir}\n")
        except Exception as e:
            sys.stderr.write(f"Error creating output directory: {str(e)}\n")
            sys.exit(1)

    # Load prep and run results
    sys.stderr.write("\nLoading bofasa results...\n")
    prep_data = load_prep_results(prep_dir)
    run_data = load_run_results(run_dir)
    
    sys.stderr.write(f"  Found {len(prep_data['sample_genbanks'])} GenBank files\n")
    sys.stderr.write(f"  Found {len(prep_data['sample_proteins'])} protein files\n")

    # Parse OG matrices
    sys.stderr.write(f"\nParsing ortholog group data for focal OG: {focal_og}\n")
    focal_og_data, sample_names = parse_og_matrix(run_data['og_matrix'], focal_og)
    
    sys.stderr.write(f"  Focal OG found in {len(focal_og_data)} samples\n")
    
    # Filter for single-copy if requested
    if single_copy_only:
        focal_og_data = {
            sample: proteins 
            for sample, proteins in focal_og_data.items() 
            if len(proteins) == 1
        }
        sys.stderr.write(f"  After single-copy filter: {len(focal_og_data)} samples\n")
    
    # Load full OG matrix for context
    sys.stderr.write("\nLoading complete OG matrix for context analysis...\n")
    all_og_data, _ = parse_full_og_matrix(run_data['og_matrix'])
    sys.stderr.write(f"  Loaded {len(all_og_data)} ortholog groups\n")

    # Extract genomic context
    sys.stderr.write("\nExtracting genomic context...\n")
    context_data = extract_genomic_context(
        prep_data['sample_genbanks'],
        focal_og_data,
        all_og_data,
        context_window,
        context_bp,
    )
    
    sys.stderr.write(f"  Extracted context for {len(context_data)} instances\n")
    
    if not context_data:
        sys.stderr.write("Error: No context data extracted. Check input files.\n")
        sys.exit(1)

    # ── Always build a gene tree for the focal OG ───────────────────────────
    gene_tree_nwk = build_gene_tree(focal_og, run_dir, prep_dir, output_dir, cpus)

    # ── If no species tree provided, build one from SCC OGs ─────────────────
    if not species_tree:
        scc_tree_nwk = build_species_tree_from_scc(run_dir, prep_dir, output_dir, cpus)
        if scc_tree_nwk and os.path.isfile(scc_tree_nwk):
            sys.stderr.write(
                f"\nNo species tree provided — using SCC-based tree: {scc_tree_nwk}\n"
            )
            species_tree = scc_tree_nwk
        elif gene_tree_nwk and os.path.isfile(gene_tree_nwk):
            sys.stderr.write(
                f"\nSCC tree unavailable — falling back to gene tree for ordering.\n"
            )
            species_tree = gene_tree_nwk

    # Assign colors to OGs
    sys.stderr.write("\nAssigning colors to ortholog groups...\n")
    og_colors = assign_og_colors(context_data)
    sys.stderr.write(f"  Assigned colors to {len(og_colors)} OGs\n")

    # Handle samples with multiple instances (paralogs)
    sys.stderr.write("\nHandling samples with multiple instances...\n")

    # Count instances per PARENT sample (strips _phage_N/_plasmid_N suffixes),
    # matching the same logic used by get_unique_sample_ids.  This ensures that
    # two different MGE sample_ids from the same genome (e.g. Staph_phage_1 and
    # Staph_plasmid_2) are treated as multiple instances and the tree is modified
    # to give them separate leaves.
    parent_instance_counts = defaultdict(int)
    for context in context_data.values():
        parent = get_parent_sample_name(context['sample_id'])
        parent_instance_counts[parent] += 1

    samples_with_paralogs = {p: c for p, c in parent_instance_counts.items() if c > 1}

    if samples_with_paralogs:
        sys.stderr.write(f"  Found {len(samples_with_paralogs)} parent samples with multiple instances:\n")
        for sample_id, count in sorted(samples_with_paralogs.items()):
            sys.stderr.write(f"    - {sample_id}: {count} instances\n")
    else:
        sys.stderr.write("  All samples have single instances (no paralogs)\n")
    
    # Generate unique sample IDs and modify tree if needed
    unique_sample_ids = None
    if species_tree and os.path.isfile(species_tree):
        # Validate that tree tip labels match expected parent sample names.
        # The species tree must use genome-level labels (e.g. "Staphylococcus_aureus"),
        # NOT MGE-level labels (e.g. "Staphylococcus_aureus_phage_1"), because
        # get_parent_sample_name strips the _phage_N/_plasmid_N suffix when building
        # molecule IDs. Warn early if there is a mismatch.
        try:
            tree_tips = set()
            with open(species_tree) as _tf:
                for _line in _tf:
                    tree_tips.update(re.findall(r'[^(),;:\[\]]+(?=[:,)\[]|$)', _line.strip()))
            tree_tips = {t.strip() for t in tree_tips if t.strip()}
            expected_parents = {
                get_parent_sample_name(ctx['sample_id'])
                for ctx in context_data.values()
            }
            unmatched = expected_parents - tree_tips
            if unmatched:
                sys.stderr.write(
                    f"\nWARNING: {len(unmatched)} sample parent name(s) not found as tree tips.\n"
                    f"  This usually means the tree has MGE-level tip labels instead of\n"
                    f"  genome-level labels. Heatmap tree alignment will be incorrect.\n"
                    f"  Unmatched: {', '.join(sorted(unmatched)[:5])}"
                    f"{'...' if len(unmatched) > 5 else ''}\n"
                )
        except Exception:
            pass  # Don't let validation block the main workflow

        if samples_with_paralogs:
            sys.stderr.write("\nModifying species tree for samples with multiple instances...\n")
            modified_tree_file = os.path.join(output_dir, 'species_tree.modified.nwk')
            unique_sample_ids = modify_tree_for_paralogs(
                species_tree, 
                context_data, 
                modified_tree_file
            )
            sys.stderr.write(f"  Modified tree saved to: {modified_tree_file}\n")
            sys.stderr.write("  NOTE: Use the modified tree for visualization!\n")
        else:
            sys.stderr.write("\nNo tree modification needed (no paralogs detected)\n")
            unique_sample_ids = get_unique_sample_ids(context_data)
    else:
        # No tree provided, just generate unique IDs
        unique_sample_ids = get_unique_sample_ids(context_data)

    # Determine which tree file to use for heatmap (modified if paralogs exist)
    heatmap_tree_file = None
    if species_tree and os.path.isfile(species_tree):
        if samples_with_paralogs:
            heatmap_tree_file = os.path.join(output_dir, 'species_tree.modified.nwk')
        else:
            heatmap_tree_file = os.path.abspath(species_tree)

    # Write output files
    sys.stderr.write("\nWriting output files...\n")

    heatmap_file = os.path.join(output_dir, 'OG_Context_Heatmap.svg')
    write_context_heatmap(context_data, og_colors, focal_og, heatmap_file,
                          unique_sample_ids, tree_file=heatmap_tree_file)

    gggenes_file = os.path.join(output_dir, 'OG_Context_Visualization.gggenes.txt')
    write_gggenes_output(context_data, og_colors, gggenes_file, unique_sample_ids)
    sys.stderr.write(f"  Wrote gggenes file: {gggenes_file}\n")
    
    itol_file = os.path.join(output_dir, 'OG_Context_Visualization.iTol.txt')
    write_itol_output(context_data, og_colors, itol_file, unique_sample_ids)
    sys.stderr.write(f"  Wrote iTOL file: {itol_file}\n")

    write_neighborhood_genbanks(context_data, focal_og, output_dir, prep_dir,
                                 unique_sample_ids)

    # Write color scheme
    color_file = os.path.join(output_dir, 'Color_Scheme.txt')
    with open(color_file, 'w') as f:
        f.write('OG\tColor\n')
        for og_id, color in sorted(og_colors.items()):
            f.write(f"{og_id}\t{color}\n")
    sys.stderr.write(f"  Wrote color scheme: {color_file}\n")
    
    # Write summary info
    info_file = os.path.join(output_dir, 'OG_info.txt')
    with open(info_file, 'w') as f:
        f.write(f"Focal OG: {focal_og}\n")
        f.write(f"Samples with focal OG: {len(focal_og_data)}\n")
        f.write(f"Total instances visualized: {len(context_data)}\n")
        f.write(f"Context window: {context_window} genes\n")
        if context_bp:
            f.write(f"Context window (bp): {context_bp} bp\n")
        f.write(f"Unique OGs in context: {len(og_colors)}\n")
        f.write(f"\nParalog Information:\n")
        f.write(f"Samples with multiple instances: {len(samples_with_paralogs)}\n")
        if samples_with_paralogs:
            f.write(f"\nDetails:\n")
            for sample_id, count in sorted(samples_with_paralogs.items()):
                f.write(f"  {sample_id}: {count} instances\n")
        if species_tree and samples_with_paralogs:
            f.write(f"\nNote: Species tree was modified to include duplicate leaves.\n")
            f.write(f"Use species_tree.modified.nwk for visualization.\n")
    sys.stderr.write(f"  Wrote summary info: {info_file}\n")

    sys.stderr.write("\n" + "="*60 + "\n")
    sys.stderr.write("Visualization complete!\n")
    sys.stderr.write("\nOutput files:\n")
    sys.stderr.write(f"  - {heatmap_file}  (primary visualization - SVG)\n")
    sys.stderr.write(f"  - {heatmap_file.replace('.svg', '.png')}  (raster format - smaller file size)\n")
    sys.stderr.write(f"  - {os.path.join(output_dir, 'Neighborhood_GenBanks')}/ (per-instance GenBank files)\n")
    if gene_tree_nwk and os.path.isfile(gene_tree_nwk):
        sys.stderr.write(f"  - {gene_tree_nwk}  (focal OG gene tree)\n")
    sys.stderr.write(f"  - {gggenes_file}  (tab-separated; import into R/gggenes if needed)\n")
    sys.stderr.write(f"  - {itol_file}\n")
    sys.stderr.write(f"  - {color_file}\n")
    sys.stderr.write(f"  - {info_file}\n")
    if species_tree and samples_with_paralogs:
        modified_tree_file = os.path.join(output_dir, 'species_tree.modified.nwk')
        sys.stderr.write(f"  - {modified_tree_file} (modified tree with duplicate leaves)\n")

    if samples_with_paralogs:
        sys.stderr.write("\n" + "!"*60 + "\n")
        sys.stderr.write("IMPORTANT: Samples with multiple instances detected!\n")
        sys.stderr.write("Sample IDs in output files have been modified:\n")
        for sample_id, count in sorted(samples_with_paralogs.items()):
            instances = ", ".join([f"{sample_id}_instance{i}" for i in range(1, count+1)])
            sys.stderr.write(f"  {sample_id} -> {instances}\n")
        if species_tree:
            sys.stderr.write("\nThe species tree has been modified with duplicate leaves.\n")
        sys.stderr.write("!"*60 + "\n")

    sys.stderr.write("="*60 + "\n")


if __name__ == '__main__':
    visualize_og_context()
