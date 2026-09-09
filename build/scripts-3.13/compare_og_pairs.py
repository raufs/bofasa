#!python

"""
Program: compare_og_pairs.py
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
from typing import Set, Tuple

from rich_argparse import RawTextRichHelpFormatter


def create_parser() -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="""
Program: compare_og_pairs
Authors: Rauf Salamzade, Aamuktha Kottapalli, Lindsay R. Kalan 
Affiliation: University of Wisconsin - Madison, McMaster University

compare_og_pairs: a simple script to compare results from two 
different orthology prediction methods.
""",
        formatter_class=RawTextRichHelpFormatter,
    )

    parser.add_argument(
        '-a',
        '--method-a',
        help='File listing orthology results for method (a) as pair of proteins '
             'co-found in ortholog groups per line.',
        required=True,
    )
    parser.add_argument(
        '-b',
        '--method-b',
        help='File listing orthology results for method (b) as pair of proteins '
             'co-found in ortholog groups per line.',
        required=True,
    )
    parser.add_argument(
        '-an',
        '--method-a-name',
        help='Name of method (a).',
        required=False,
        default='Method A',
    )
    parser.add_argument(
        '-bn',
        '--method-b-name',
        help='Name of method (b).',
        required=False,
        default='Method B',
    )

    args = parser.parse_args()
    return args


def compare_og_pairs() -> None:
    """Void function which runs primary workflow for program."""
    # PARSE ARGUMENTS
    myargs = create_parser()

    method_a_file = myargs.method_a
    method_b_file = myargs.method_b

    method_a_name = myargs.method_a_name
    method_b_name = myargs.method_b_name

    # START WORKFLOW
    method_a_set: Set[Tuple[str, str]] = set([])
    method_b_set: Set[Tuple[str, str]] = set([])

    with open(method_a_file) as omaf:
        for line in omaf:
            line = line.strip()
            ls = line.split('\t')
            method_a_set.add(tuple(sorted(ls)))

    with open(method_b_file) as ombf:
        for line in ombf:
            line = line.strip()
            ls = line.split('\t')
            method_b_set.add(tuple(sorted(ls)))

    intersection = len(method_a_set.intersection(method_b_set))
    union = len(method_a_set.union(method_b_set))
    unique_a = len(method_a_set.difference(method_b_set))
    unique_b = len(method_b_set.difference(method_a_set))
    jaccard_index = intersection / union

    sys.stdout.write(
        '\t'.join(
            [
                'Method_A',
                'Method_B',
                'Jaccard_Index',
                'Method_A_Unique',
                'Method_B_Unique',
                'Intersection',
                'Union',
            ]
        ) + '\n'
    )
    sys.stdout.write(
        '\t'.join(
            [
                str(x)
                for x in [
                    method_a_name,
                    method_b_name,
                    jaccard_index,
                    unique_a,
                    unique_b,
                    intersection,
                    union,
                ]
            ]
        ) + '\n'
    )


if __name__ == '__main__':
    compare_og_pairs()
