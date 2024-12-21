#!/usr/bin/env python3

### Program: processNCBIGenBank.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison, Department of Medical Microbiology and Immunology

# BSD 3-Clause License
#
# Copyright (c) 2021, Kalan-Lab
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

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
from bofasa import util

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: processNCBIGenBank.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Process NCBI Genbanks to create proteome + coords BED with specific locus tag.
	""", formatter_class=argparse.RawTextHelpFormatter)

	parser.add_argument('-i', '--input_ncbi_genbank', help='Path to genomic assembly in GenBank format.', required=True)
	parser.add_argument('-o', '--outdir', help='Path to output directory where files should be written. Should already be created!', required=True)
	parser.add_argument('-s', '--sample_name', help='Sample name', default='Sample', required=False)
	parser.add_argument('-l', '--locus_tag', help='Locus tag', default=None, required=False)

	args = parser.parse_args()
	return args


def processAndReformatNCBIGenbanks():
	"""
	Void function which runs primary workflow for program.
	"""

	"""
	PARSE REQUIRED INPUTS
	"""
	myargs = create_parser()

	input_genbank_file = os.path.abspath(myargs.input_ncbi_genbank)
	outdir = os.path.abspath(myargs.outdir) + '/'

	try:
		assert (util.is_genbank(input_genbank_file))
	except:
		raise RuntimeError('Issue with input Genbank file from NCBI.')

	outdirs = [outdir]
	for outdir in outdirs:
		if not os.path.isdir(outdir):
			sys.stderr.write("Output directory %s does not exist! Please create and retry program." % outdir)

	"""
	PARSE OPTIONAL INPUTS
	"""

	sample_name = myargs.sample_name
	locus_tag = myargs.locus_tag

	"""
	START WORKFLOW
	"""

	# Step 1: Process NCBI Genbank and (re)create genbank/proteome with updated locus tags.
	try:
		bed_outfile = outdir + sample_name + '.coords.bed'
		pro_outfile = outdir + sample_name + '.faa'
		fna_outfile = outdir + sample_name + '.fna'
		map_outfile = outdir + sample_name + '.name_map.txt'

		bed_outfile_handle = open(bed_outfile, 'w')
		pro_outfile_handle = open(pro_outfile, 'w')
		map_outfile_handle = open(map_outfile, 'w')
		fna_outfile_handle = open(fna_outfile, 'w')


		locus_tag_iterator = 1
		oigf = None
		if input_genbank_file.endswith('.gz'):
			oigf = gzip.open(input_genbank_file, 'rt')
		else:
			oigf = open(input_genbank_file)
		for rec in SeqIO.parse(oigf, 'genbank'):
			scaffold = rec.id
			scaffold_length = len(str(rec.seq))
			fna_outfile_handle.write('>' + scaffold + '\n' + str(rec.seq) + '\n')
			for feature in rec.features:
				if feature.type == "CDS":
					start = min([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')]) + 1
					end = max([int(x.strip('>').strip('<')) for x in str(feature.location)[1:].split(']')[0].split(':')])
					direction = str(feature.location).split('(')[1].split(')')[0]

					old_locus_tag = None
					prot_seq = None
					try:
						old_locus_tag = feature.qualifiers.get('locus_tag')[0]
					except:
						pass
					try:
						prot_seq = str(feature.qualifiers.get('translation')[0]).replace('*', '')
					except:
						msg = "Currently only full Genbanks with translations available for each CDS is accepted."
						sys.stderr.write(msg + '\n')

					if locus_tag != None:
						new_locus_tag = locus_tag + '_'
						if locus_tag_iterator < 10:
							new_locus_tag += '00000' + str(locus_tag_iterator)
						elif locus_tag_iterator < 100:
							new_locus_tag += '0000' + str(locus_tag_iterator)
						elif locus_tag_iterator < 1000:
							new_locus_tag += '000' + str(locus_tag_iterator)
						elif locus_tag_iterator < 10000:
							new_locus_tag += '00' + str(locus_tag_iterator)
						elif locus_tag_iterator < 100000:
							new_locus_tag += '0' + str(locus_tag_iterator)
						else:
							new_locus_tag += str(locus_tag_iterator)
						locus_tag_iterator += 1
						feature.qualifiers['locus_tag'] = new_locus_tag
					else:
						new_locus_tag = old_locus_tag
						assert(new_locus_tag != None)

					prot_score = '1'
					if (scaffold_length-end) < 1000 or start < 1000:
						prot_score = '0'
					bed_outfile_handle.write('\t'.join([str(x) for x in [scaffold, start, end, new_locus_tag, prot_score, direction]]) + '\n')
					pro_outfile_handle.write('>' + str(new_locus_tag) + ' ' + rec.id + ' ' + str(start) + ' ' + str(end) + ' ' + str(direction) + '\n' + prot_seq + '\n')
					map_outfile_handle.write(str(old_locus_tag) + '\t' + str(new_locus_tag) + '\n')
				
		oigf.close()
		bed_outfile_handle.close()
		pro_outfile_handle.close()
		map_outfile_handle.close()
		fna_outfile_handle.close()
	except:
		raise RuntimeError("Issue processing GenBank file.")

if __name__ == '__main__':
	processAndReformatNCBIGenbanks()