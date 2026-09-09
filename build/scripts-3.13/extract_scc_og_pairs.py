#!python

"""
Program: extract_scc_og_pairs.py
Authors: Rauf Salamzade, Aamuktha Kottapalli, Lindsay R. Kalan 
Affiliation: University of Wisconsin - Madison, McMaster University
"""

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
Program: extract_scc_og_pairs
Authors: Rauf Salamzade, Aamuktha Kottapalli, Lindsay R. Kalan 
Affiliation: University of Wisconsin - Madison, McMaster University

extract_scc_og_pairs: a simple script to list pairs of proteins in the 
same ortholog group for method comparison. Only ortholog groups 
which are scc are regarded.

Options for software supported include:
===================================================================
- bofasa 
- panx
- sonicparanoid
- panaroo
- scarap
- pirate
- orthofinder_hog (OrthoFinder v2.5.4/v2.5.5 results - using 
                   phylogenetically refined hierarchical 
                   ortholog groups from the N0.tsv file)
- orthofinder_og (OrthoFinder v2.5.4/v2.5.5 results - using 
                  coarse ortholog groups from the 
                  Orthogroups.tsv file)
""",
        formatter_class=RawTextRichHelpFormatter,
    )

    parser.add_argument(
        '-i',
        '--results-dir',
        help='Resulting directory produced by orthology prediction software.',
        required=True,
    )
    parser.add_argument(
        '-m',
        '--method',
        help='The name of the orthology prediction software - see help function '
             'for supported methods.',
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
        '-pi',
        '--pirate-input-dir',
        help='Input directory with GFFs for PIRATE - needed because PIRATE '
             'renames genes.',
        required=False,
        default=None,
    )
    parser.add_argument(
        '-pgi',
        '--pirate-gff-identifier',
        help='The identifier for CDS feautres in GFFs in the input directory for PIRATE. [Default is "ID"]',
        required=False,
        default="ID",
    )
    parser.add_argument(
        '-emg',
        '--exclude-mge-genomes',
        action='store_true',
        help='Exclude MGE (mobile genetic element) samples when determining single-copy core.\n'
             'Filters out samples ending with "_plasmid", "_phage", or ".ccds" suffix.\n'
             'Useful when bofasa prep was run with -emg flag.',
        required=False,
        default=False,
    )
    parser.add_argument(
        '-mm',
        '--max_memory',
        type=int,
        help='Uses resource module to set soft memory limit. Provide in Giga-bytes.\nGenerally memory shouldn\'t be a major concern unless working\nwith hundreds of large metagenomes. [currently\nexperimental; default is None].',
        default=None,
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


valid_software = set(
    [
        'bofasa',
        'panx',
        'sonicparanoid',
        'pirate',
        'panaroo',
        'scarap',
        'orthofinder_hog',
        'orthofinder_og',
    ]
)


def extract_og_pairs():
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

    input_dir = os.path.abspath(myargs.results_dir) + '/'
    outfile = myargs.output_file
    method = myargs.method
    max_memory = myargs.max_memory
    pirate_input_dir = myargs.pirate_input_dir
    pirate_gff_identifier = myargs.pirate_gff_identifier
    exclude_mge_genomes = myargs.exclude_mge_genomes

    # START WORKFLOW

    # set max memory limit
    if max_memory is not None:
        try:
            utils.memory_limit(max_memory)
        except Exception as e:
            sys.stderr.write(f"Error setting memory limit: {str(e)}\n")

    try:
        assert method in valid_software
    except Exception:
        msg = 'Method %s not in set of accepted methods.' % method
        sys.stderr.write(msg)
        sys.exit(1)

    outf_handle = sys.stdout
    if outfile is not None:
        outf_handle = open(outfile, 'w')

    if method == 'bofasa':
        result_file = input_dir + '/Final_Results/Protein_Ortholog_Groups.tsv'
        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'Error: Issue validating the existence of file: %s' % result_file
            sys.stderr.write(msg)
            sys.exit(1)

        # Determine which sample columns to include
        sample_indices_to_include = []
        with open(result_file) as orf:
            for i, line in enumerate(orf):
                line = line.rstrip('\n')
                ls = line.split('\t')
                
                if i == 0:
                    # Header line - determine which samples to include
                    for j, sample_name in enumerate(ls[1:], start=1):
                        sample_name = '_'.join(sample_name.strip().split('_')[:-1])
                        is_mge = (sample_name.endswith('_phage') or sample_name.endswith('_plasmid'))
                        
                        if exclude_mge_genomes and is_mge:
                            continue
                        else:
                            sample_indices_to_include.append(j)
                    
                    if exclude_mge_genomes:
                        num_excluded = len(ls) - 1 - len(sample_indices_to_include)
                        if num_excluded > 0:
                            sys.stderr.write(f"Excluding {num_excluded} MGE sample(s) from analysis\n")
                    continue
                
                # Data lines - check for single-copy core
                og_lts = set([])
                scc_flag = True
                for j in sample_indices_to_include:
                    lts = ls[j].strip()
                    if lts == '' or ',' in lts:
                        scc_flag = False
                    for lt in lts.split(','):
                        if lt.strip() == '':
                            continue
                        og_lts.add(lt.strip())
                
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    elif method == 'panx':
        result_file = input_dir + '/allclusters_final.tsv'
        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'Error: Issue validating the existence of file: %s' % result_file
            sys.stderr.write(msg)
            sys.exit(1)

        gbk_dir = input_dir + '/input_GenBank/'

        all_genomes = set([])
        for f in os.listdir(gbk_dir):
            genome = '.'.join(f.split('.')[:-1])
            all_genomes.add(genome)

        with open(result_file) as orf:
            for i, line in enumerate(orf):
                line = line.rstrip('\n')
                ls = line.split('\t')
                genome_counts = defaultdict(int)
                og_lts = set([])
                for lts in ls:
                    lts = lts.strip()
                    genome, lt = lts.split('|')
                    genome_counts[genome] += 1
                    og_lts.add(lt)

                scc_flag = True
                for g in all_genomes:
                    if genome_counts[g] != 1:
                        scc_flag = False
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    elif method == 'sonicparanoid':
        info_file = input_dir + 'last_run_info.txt'
        try:
            assert os.path.isfile(info_file)
        except Exception:
            msg = 'Error: Issue validating the existence of file: %s' % info_file
            sys.stderr.write(msg)
            sys.exit(1)

        result_file = None
        with open(info_file) as oif:
            for line in oif:
                line = line.rstrip('\n')
                ls = line.split('\t')
                if ls[0] == 'Directory with ortholog groups:':
                    result_file = ls[1] + 'flat.ortholog_groups.tsv'

        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'Error: Issue validating the existence of file: %s' % result_file
            sys.stderr.write(msg)
            sys.exit(1)

        with open(result_file) as orf:
            for i, line in enumerate(orf):
                if i == 0:
                    continue
                line = line.rstrip('\n')
                ls = line.split('\t')
                og_lts = set([])
                scc_flag = True
                for lts in ls[1:]:
                    lts = lts.strip()
                    if lts == '*' or ',' in lts:
                        scc_flag = False
                    if lts == '*':
                        continue
                    for lt in lts.split(','):
                        if lt.strip() == '':
                            continue
                        og_lts.add(lt.strip())
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    elif method == 'panaroo':
        result_file = input_dir + 'gene_presence_absence.csv'

        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'Error: Issue validating the existence of file: %s' % result_file
            sys.stderr.write(msg)
            sys.exit(1)

        with open(result_file) as orf:
            for i, line in enumerate(orf):
                if i == 0:
                    continue
                line = line.rstrip('\n')
                ls = line.split(',')
                og_lts = set([])
                scc_flag = True
                for lts in ls[3:]:
                    lts = lts.strip()
                    if lts == '' or ';' in lts:
                        scc_flag = False
                    for lt in lts.split(';'):
                        if lt.strip() == '':
                            continue
                        og_lts.add(lt)
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    elif method == 'pirate':
        try:
            assert os.path.isdir(pirate_input_dir)
        except Exception:
            msg = 'PIRATE input directory with GFFs was not found!'
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        sample_ids_in_order = defaultdict(list)
        for f in os.listdir(pirate_input_dir):
            s = '_'.join(f.split('.')[:-1]).replace('-', '_')
            with open(pirate_input_dir + f) as ogff:
                for line in ogff:
                    if line.startswith('#'):
                        continue
                    line = line.strip()
                    ls = line.split('\t')
                    if len(ls) < 6:
                        continue
                    if ls[2] != 'CDS':
                        continue
                    gid = ls[-1].split(pirate_gff_identifier + '=')[1].split(';')[0]
                    sample_ids_in_order[s].append(gid)

        g2l_file = input_dir + 'genome2loci.tab'
        try:
            assert os.path.isfile(g2l_file)
        except Exception:
            msg = 'File not found %s' % g2l_file
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        gid_to_gname = {}
        with open(g2l_file) as og2f:
            for line in og2f:
                line = line.strip()
                ls = line.split()
                s = ls[1]
                gid = int(ls[0].split('_')[-1]) - 1
                gname = sample_ids_in_order[s][gid]
                gid_to_gname[ls[0]] = gname

        result_file = input_dir + 'PIRATE.gene_families.ordered.tsv'
        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'File not found %s' % result_file
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        with open(result_file) as orf:
            for i, line in enumerate(orf):
                if i == 0:
                    continue
                line = line.rstrip('\n')
                ls = line.split('\t')
                og_lts = set([])
                scc_flag = True
                for lts in ls[22:]:
                    lts = lts.strip()
                    if lts == '' or ';' in lts:
                        scc_flag = False
                    for lt in lts.split(';'):
                        if lt.strip() == '':
                            continue
                        if lt.startswith('('):
                            lt = lt[1:-1]
                            for slt in lt.split(':'):
                                og_lts.add(gid_to_gname[slt.strip()])
                        else:
                            og_lts.add(gid_to_gname[lt.strip()])
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    elif method == 'scarap':
        result_file = input_dir + 'pangenome.tsv'
        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'File not found %s' % result_file
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        og_lts = defaultdict(set)
        all_samples = set([])
        og_sample_counts = defaultdict(lambda: defaultdict(int))
        with open(result_file) as orf:
            for i, line in enumerate(orf):
                line = line.rstrip('\n')
                lt, g, og = line.split('\t')
                og_lts[og].add(lt)
                og_sample_counts[og][g] += 1
                all_samples.add(g)

        scc_ogs = set([])
        for og in og_sample_counts:
            scc_flag = True
            for s in all_samples:
                if og_sample_counts[og][s] != 1:
                    scc_flag = False
            if scc_flag:
                scc_ogs.add(og)

        for og in og_lts:
            if og in scc_ogs:
                outf_handle.write('\t'.join(sorted(og_lts[og])) + '\n')

    elif method == 'orthofinder_og':
        subdirs = []
        for sd in os.listdir(input_dir):
            subdirs.append(sd)

        try:
            assert len(subdirs) == 1
        except Exception:
            msg = 'None or multiple results found within OrthoFinder directory!'
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        result_file = input_dir + subdirs[0] + '/Orthogroups/Orthogroups.tsv'

        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'File not found %s' % result_file
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        with open(result_file) as orf:
            for i, line in enumerate(orf):
                if i == 0:
                    continue
                line = line.rstrip('\n')
                ls = line.split('\t')
                og_lts = set([])
                scc_flag = True
                for lts in ls[1:]:
                    lts = lts.strip()
                    if lts == '' or ',' in lts:
                        scc_flag = False
                    for lt in lts.split(','):
                        if lt.strip() == '':
                            continue
                        og_lts.add(lt.strip())
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    elif method == 'orthofinder_hog':

        subdirs = []
        for sd in os.listdir(input_dir):
            subdirs.append(sd)

        try:
            assert len(subdirs) == 1
        except Exception:
            msg = 'None or multiple results found within OrthoFinder directory!'
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        result_file = (
            input_dir + subdirs[0] + '/Phylogenetic_Hierarchical_Orthogroups/N0.tsv'
        )

        try:
            assert os.path.isfile(result_file)
        except Exception:
            msg = 'File not found %s' % result_file
            sys.stderr.write(msg + '\n')
            sys.exit(1)

        with open(result_file) as orf:
            for i, line in enumerate(orf):
                if i == 0:
                    continue
                line = line.rstrip('\n')
                ls = line.split('\t')
                og_lts = set([])
                scc_flag = True
                for lts in ls[3:]:
                    lts = lts.strip()
                    if lts == '' or ',' in lts:
                        scc_flag = False
                    for lt in lts.split(','):
                        if lt.strip() == '':
                            continue
                        og_lts.add(lt.strip())
                if scc_flag:
                    outf_handle.write('\t'.join(sorted(og_lts)) + '\n')

    if outfile is not None:
        outf_handle.close()


if __name__ == '__main__':
    extract_og_pairs()
