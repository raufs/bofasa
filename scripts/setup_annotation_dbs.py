#!/usr/bin/env python3

### Program: setup_annotation_dbs.py
### Author: Rauf Salamzade
### Kalan Lab
### UW Madison & McMaster University
### This code is largely adapted from the prepTG program in zol.

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

import os
import sys
import argparse
from Bio import SeqIO
import shutil
import traceback

def create_parser():
	""" Parse arguments """
	parser = argparse.ArgumentParser(description="""
	Program: setup_annotation_dbs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology
		
	Downloads annotation databases for ISfinder, Pfam, and geNomad.

	Location of where to download databases is controlled by setting the environmental
	variable BOFASA_DB_PATH e.g.:
						
	$ export BOFASA_DB_PATH=/path/to/database_location/
	""", formatter_class=argparse.RawTextHelpFormatter)
	args = parser.parse_args()
	return args

def setup_annot_dbs():
	myargs = create_parser()

	download_path = None
	if str(os.getenv("BOFASA_DB_PATH")) != 'None':
		download_path = str(os.getenv("BOFASA_DB_PATH"))
		
	if download_path == None or not os.path.isdir(download_path):
		sys.stderr.write('Issues validing database download directory exists.\n')
		sys.exit(1)

	try:
		assert(os.path.isdir(download_path))
	except:
		sys.stderr.write('Error: Provided directory for downloading annotation files does not exist!\n')

	try:
		shutil.rmtree(download_path)
		os.mkdir(download_path)
		readme_outf = open(download_path + 'README.txt', 'w')
		readme_outf.write('Default space for downloading annotation databases.\n')
		readme_outf.close()
	except Exception as e:
		sys.stderr.write('Issues clearing contents of db/ to re-try downloads.\n')
		sys.stderr.write(str(e) + '\n')
		sys.exit(1)

	listing_file = download_path + 'database_location_paths.txt'
	issues_file = download_path + 'issues.txt'
	listing_handle = open(listing_file, 'w')
	issues_handle = open(issues_file, 'w')

	
	# Final annotation files
	pfam_phmm_file = download_path + 'Pfam-A.hmm'
	is_faa_file = download_path + 'isfinder.dmnd'
	genomad_db_tar = download_path + 'genomad_db_v1.7.tar.gz'
	download_links = ['https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz',
					  'https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.faa',
					  'https://zenodo.org/records/10594875/files/genomad_db_v1.7.tar.gz?download=1']
	download_files = ['Pfam-A.hmm.gz', 'IS.faa', 'genomad_db_v1.7.tar.gz']

	# Download
	print('Starting download of files!')
	os.chdir(download_path)
	try:
		for i, dl in enumerate(download_links):
			df = download_files[i]
			curl_download_dbs_cmd = ['curl', dl, '-o', df]
			os.system(' '.join(curl_download_dbs_cmd))
	except Exception as e:
		sys.stderr.write('Error occurred during downloading!\n')
		issues_handle.write('Error occurred during downloading with axel.\n')
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(str(e) + '\n')

	try:
		print('Setting up Pfam database!')
		os.system(' '.join(['gunzip', 'Pfam-A.hmm.gz']))
		assert(os.path.isfile(pfam_phmm_file))
		name = None
		desc = None
		pfam_descriptions_file = download_path + 'pfam_descriptions.txt'
		pdf_handle = open(pfam_descriptions_file, 'w')
		with open(pfam_phmm_file) as oppf:
			for line in oppf:
				line = line.strip()
				ls = line.split()
				if ls[0].strip() == 'NAME':
					name = ' '.join(ls[1:]).strip()
				elif ls[0].strip() == 'DESC':
					desc = ' '.join(ls[1:]).strip()
					pdf_handle.write(name + '\t' + desc + '\n')
		pdf_handle.close()
		os.system(' '.join(['hmmpress', pfam_phmm_file]))
		z = 0
		with open(pfam_phmm_file) as oppf:
			for line in oppf:
				if line.startswith('NAME'): z += 1
		listing_handle.write('pfam\t' + pfam_descriptions_file + '\t' + pfam_phmm_file + '\t' + str(z) + '\n')
	except Exception as e:
		sys.stderr.write('Issues setting up Pfam database.\n')
		issues_handle.write('Issues setting up Pfam database.\n')
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(str(e) + '\n')

	try:
		print('Setting up ISFinder database!')
		is_faa_path = download_path + 'IS.faa'
		os.system(' '.join(['diamond', 'makedb', '--in', is_faa_path, '-d', is_faa_file]))
		is_descriptions_file = download_path + 'is_descriptions.txt'
		idf_handle = open(is_descriptions_file, 'w')
		with open(is_faa_path) as oif:
			for rec in SeqIO.parse(oif, 'fasta'):
				idf_handle.write(rec.id + '\t' + rec.description + '\n')
		idf_handle.close()
		assert(os.path.isfile(is_faa_path))
		assert(os.path.isfile(is_descriptions_file))
		listing_handle.write('isfinder\t' + is_descriptions_file + '\t' + is_faa_file + '\tNA\n')
		os.system(' '.join(['rm', '-rf', is_faa_path]))
	except Exception as e:
		sys.stderr.write('Issues setting up ISFinder database.\n')
		issues_handle.write('Issues setting up ISFinder database.\n')
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(str(e) + '\n')

	try:
		print('Setting up geNomad databases!')
		assert(os.path.isfile(genomad_db_tar))
		os.system(' '.join(['tar', '-zxvf', genomad_db_tar]))
		genomad_db_dir = download_path + 'genomad_db/'
		assert(os.path.isdir(genomad_db_dir))
		listing_handle.write('genomad\tNA\t' + genomad_db_dir + '\tNA\n')
	except Exception as e:
		sys.stderr.write('Issues setting up geNomad databases.\n')
		issues_handle.write('Issues setting up geNomad databases.\n')
		sys.stderr.write(traceback.format_exc())
		sys.stderr.write(str(e) + '\n')

	listing_handle.close()
	issues_handle.close()

	print("Done setting up annotation databases!")
	print("Information on final files used by bofasa can be found at:\n%s" % listing_file)
	sys.exit(0)
setup_annot_dbs()