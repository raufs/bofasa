import os
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
import logging
import subprocess
from operator import itemgetter
from collections import defaultdict
import traceback
import numpy as np
import gzip
import copy
import itertools
import multiprocessing
import pickle
import resource
import pkg_resources  # part of setuptools
import shutil
import statistics
from scipy import stats

version = pkg_resources.require("bofasa")[0].version


def multiProcess(input):
	"""
	Description:
	This is a generalizable function to be used with multiprocessing to parallelize list of commands. Inputs should
	correspond to space separated command (as list), with last item in list corresponding to a logging object handle for
	logging progress.
	********************************************************************************************************************
	Parameters:
	- input: A list corresponding to a command to run with the last item in the list corresponding to a logging object
	         for the function.
	********************************************************************************************************************
	"""
	input_cmd = input[:-1]
	logObject = input[-1]
	logObject.info('Running the following command: %s' % ' '.join(input_cmd))
	try:
		subprocess.call(' '.join(input_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
						executable='/bin/bash')
		logObject.info('Successfully ran: %s' % ' '.join(input_cmd))
	except Exception as e:
		logObject.error('Had an issue running: %s' % ' '.join(input_cmd))
		sys.stderr.write('Had an issue running: %s' % ' '.join(input_cmd))
		logObject.error(e)
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)


def getVersion():
	"""
	Description:
	Parses the version of the zol suite from the setup.py program.
	********************************************************************************************************************
	"""
	return(str(version))

def createLoggerObject(log_file):
	"""
	Description:
	This function creates a logging object.
	********************************************************************************************************************
	Parameters:
	- log_file: Path to file to which to write logging.
	********************************************************************************************************************
	Returns:
	- logger: A logging object.
	********************************************************************************************************************
	"""

	logger = logging.getLogger('task_logger')
	logger.setLevel(logging.DEBUG)
	# create file handler which logs even debug messages
	fh = logging.FileHandler(log_file)
	fh.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', "%Y-%m-%d %H:%M")
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	return logger


def closeLoggerObject(logObject):
	"""
	Description:
	This function closes a logging object.
	********************************************************************************************************************
	Parameters:
	- logObject: A logging object.
	********************************************************************************************************************
	"""

	handlers = logObject.handlers[:]
	for handler in handlers:
		handler.close()
		logObject.removeHandler(handler)


def logParametersToFile(parameter_file, parameter_names, parameter_values):
	"""
	Description:
	This function serves to create a parameters input file for major programs, e.g. fai and zol.
	********************************************************************************************************************
	Parameters:
	- parameter_file: The path to the file where to write parameter information to. Will overwrite each time.
	- parameter_names: A list containing the parameter names.
	- parameter_values: A list in the same order as parameter_names which contains the respective arguments provided.
	********************************************************************************************************************
	"""
	parameter_handle = open(parameter_file, 'w')
	for i, pv in enumerate(parameter_values):
		pn = parameter_names[i]
		parameter_handle.write(pn + ': ' + str(pv) + '\n')
	parameter_handle.close()


def memory_limit(mem):
	"""
	Description:
	Experimental function to limit memory.
	********************************************************************************************************************
	Parameters:
	- mem: The memory limit in GB.
	********************************************************************************************************************
	"""
	max_virtual_memory = mem*1000000000
	soft, hard = resource.getrlimit(resource.RLIMIT_AS)
	resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, hard))
	print(resource.getrlimit(resource.RLIMIT_AS))


def setupReadyDirectory(directories):
	"""
	Description:
	This is a generalizable function to create directories.
	********************************************************************************************************************
	Parameters:
	- dictionaries: A list of paths to directories to create or recreate (after removing).
	********************************************************************************************************************
	"""
	try:
		assert (type(directories) is list)
		for d in directories:
			if os.path.isdir(d):
				os.system('rm -rf %s' % d)
			os.system('mkdir %s' % d)
	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def processGenomesUsingProdigal(sample_genomes, results_directory, logObject, cpus=1, locus_tag_length=3, gene_calling_method="pyrodigal", meta_mode=False, avoid_locus_tags=set([])):
	"""
	Description:
	This function oversees processing of input genomes to create proteome and GenBank files using p(y)rodigal.
	********************************************************************************************************************
	Parameters:
	- sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in FASTA format.
	- prodigal_outdir: Workspace where gene-calling results will be written to.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	- locus_tag_length: The length of the locus tags to generate.
	- use_prodigal: Whether to use prodigal instead of pyrodigal.
	- avoid_locus_tags: Whether to avoid using certain locus tags.
	********************************************************************************************************************
	"""
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]).difference(
				avoid_locus_tags)))

		prodigal_cmds = []
		for i, sample in enumerate(sample_genomes):
			sample_assembly = sample_genomes[sample]
			sample_locus_tag = ''.join(list(possible_locustags[i]))

			prodigal_cmd = ['runProdigalAndMakeInputsForBofasa.py', '-i', sample_assembly, '-s', sample,
							'-l', sample_locus_tag, '-o', results_directory]
			if gene_calling_method == 'prodigal':
				prodigal_cmd += ['-p']
			if meta_mode:
				prodigal_cmd += ['-m']
			prodigal_cmds.append(prodigal_cmd + [logObject])

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, prodigal_cmds)
		p.close()

		for sample in sample_genomes:
			try:
				assert (os.path.isfile(results_directory + sample + '.faa') and os.path.isfile(results_directory + sample + '.coords.bed'))
			except:
				sys.stderr.write("Unable to validate successful .bed & .faa creation for sample %s\n" % sample)
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
	except Exception as e:
		logObject.error(
			"Problem with creating commands for running prodigal via script runProdigalAndMakeInputsForBofasa.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def processGenomesAsGenbanks(sample_genomes, results_directory, logObject, cpus=1, locus_tag_length=3, avoid_locus_tags=set([]),
							 rename_locus_tags=False):
	"""
	Description:
	This function oversees processing of input genomes as GenBanks with CDS features already available.
	********************************************************************************************************************
	Parameters:
	- sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in GenBank format with CDS
	                  features available.
	- proteomes_directory: Directory where final proteome files (in FASTA format) for target genomes will be saved.
	- genbanks_directory: Directory where final GenBank files for target genomes will be saved.
	- gene_name_mapping_outdir: Directory where mapping files for original locus tags to new locus tags will be saved.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	- locus_tag_length: The length of the locus tags to generate.
	- avoid_locus_tags: Whether to avoid using certain locus tags.
	- rename_locus_tags: Whether to rename locus tags.
	********************************************************************************************************************
	Returns:
	- sample_genomes_updated: Dictionary mapping sample names to paths of final/processed sample GenBanks.
	********************************************************************************************************************
	"""

	sample_genomes_updated = {}
	process_cmds = []
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))]).difference(
				avoid_locus_tags)))
		lacking_cds_gbks = set([])

		for i, sample in enumerate(sample_genomes):
			sample_locus_tag = possible_locustags[i]
			sample_genbank = sample_genomes[sample]
			process_cmd = ['processNCBIGenBankAndCreateInputs.py', '-i', sample_genbank, '-s', sample,
						   '-o', results_directory]
			if rename_locus_tags:
				process_cmd += ['-l', sample_locus_tag]
			process_cmds.append(process_cmd + [logObject])

		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, process_cmds)
		p.close()

		for sample in sample_genomes:
			if sample in lacking_cds_gbks:
				continue
			try:
				assert (os.path.isfile(results_directory + sample + '.faa') and
						os.path.isfile(results_directory + sample + '.coords.bed') and
						os.path.isfile(results_directory + sample + '.txt'))
			except:
				sys.stderr.write("Unable to validate successful genbank/predicted-proteome creation for sample %s" % sample)
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
	except Exception as e:
		logObject.error("Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return sample_genomes_updated


def determineGenomeFormat(inputs):
	"""
	Description:
	This function determines whether a target sample genome is provided in GenBank or FASTA format.
	********************************************************************************************************************
	Parameters:
	- input a list which can be expanded to the following:
		- sample: The sample identifier / name.
		- genome_File: The path to the genome file.
		- format_assess_dir: The directory where the genome type information for the sample will be written.
		- logObject: A logging object.
	********************************************************************************************************************
	"""

	sample, genome_file, format_assess_dir, logObject = inputs
	try:
		gtype = 'unknown'
		if is_fasta(genome_file):
			gtype = 'fasta'
		if is_genbank(genome_file, check_for_cds=True):
			if gtype == 'fasta':
				gtype = 'unknown'
			else:
				gtype = 'genbank'
		sample_res_handle = open(format_assess_dir + sample + '.txt', 'w')
		sample_res_handle.write(sample + '\t' + str(gtype) + '\n')
		sample_res_handle.close()
	except:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def parseSampleGenomes(genome_listing_file, format_assess_dir, format_predictions_file, logObject, cpus=1):
	"""
	Description:
	This function parses the input sample target genomes and determines whether they are all provided in the same format
	and whether everything aligns with expectations.
	********************************************************************************************************************
	Parameters:
	- genome_listing_file: A tab separated file with two columns: (1) sample name, (2) path to genome file.
	- format_assess_dir: The directory/workspace where genome format information will be saved.
	- format_predictions_file: The file where to concatenate genome format information.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	********************************************************************************************************************
	Returns:
	- sample_genomes: A dictionary which maps sample names to genome file paths (note, unknown format files will be
	                  dropped).
	- format_prediction: The format prediction for genome files.
	********************************************************************************************************************
	"""
	try:
		sample_genomes = {}
		assess_inputs = []
		with open(genome_listing_file) as oglf:
			for line in oglf:
				line = line.strip()
				ls = line.split('\t')
				sample, genome_file = ls
				assess_inputs.append([sample, genome_file, format_assess_dir, logObject])
				try:
					assert (os.path.isfile(genome_file))
				except:
					logObject.warning(
						"Problem with finding genome file %s for sample %s, skipping" % (genome_file, sample))
					continue
				if sample in sample_genomes:
					logObject.warning('Skipping genome %s for sample %s because a genome file was already provided for this sample' % (genome_file, sample))
					continue
				sample_genomes[sample] = genome_file

		p = multiprocessing.Pool(cpus)
		p.map(determineGenomeFormat, assess_inputs)
		p.close()

		os.system('find %s -maxdepth 1 -type f | xargs cat >> %s' % (format_assess_dir, format_predictions_file))

		format_prediction = 'mixed'
		gtypes = set([])
		with open(format_predictions_file) as ofpf:
			for line in ofpf:
				line = line.strip()
				sample, gtype = line.split('\t')
				if gtype == 'unknown':
					sys.stderr.write('unsure about format for genome %s for sample %s, skipping inclusion...\n' % (sample_genomes[sample], sample))
					logObject.warning('unsure about format for genome %s for sample %s, skipping inclusion...' % (sample_genomes[sample], sample))
					del sample_genomes[sample]
				else:
					gtypes.add(gtype)

		if len(gtypes) == 1:
			format_prediction = list(gtypes)[0]

		return ([sample_genomes, format_prediction])

	except Exception as e:
		logObject.error("Problem with creating commands for running Prodigal. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def extractGeneContexts(inputs):
	sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, logObject = inputs
	try:
		assert(os.path.isfile(coords_file))

		scaffold_features = defaultdict(list)
		with open(coords_file) as ocf:
			for line in ocf:
				line = line.strip()
				scaffold, start, end, name, score, strand = line.split('\t')
				start = int(start); end = int(end)
				scaffold_features[scaffold].append([start, end, name, score, strand])

		output_handle = open(output_file, 'w')
		output_handle.write('\t'.join(['Protein', 'OG', 'Near scaffold edge?', 'Length', 'Upstream genes', 'Downstream genes', 'Upstream OGs', 'Downstream OGs']))
		for scaffold in scaffold_features:
			scaffold_features_sorted = sorted(scaffold_features[scaffold], key=itemgetter(0))
			max_features = len(scaffold_features[scaffold])
			for cds_index, cds in enumerate(scaffold_features_sorted):
				cds_start = cds[0]
				cds_end = cds[1]
				cds_name = cds[2]
				cds_og = gene_to_og[cds_name]
				left_boundary = cds_start-surrounding_bp
				right_boundary = cds_end+surrounding_bp
				left_side_genes_and_ogs = set([])
				right_side_genes_and_ogs = set([])
				near_scaffold_edge = False
				if cds[3] == '0':
					near_scaffold_edge = True

				limit_reached = False
				cds_iter_index = cds_index-1
				while not limit_reached:
					try:
						cds_iter = scaffold_features_sorted[cds_iter_index]
						if cds_iter_index < 0:
							near_scaffold_edge = True
							limit_reached = True
						elif cds_iter[1] < left_boundary:
							if cds_iter[3] == '0':
								near_scaffold_edge = True
							limit_reached = True
						elif cds_iter[0] < left_boundary:
							if cds_iter[3] == '0':
								near_scaffold_edge = True
							left_side_genes_and_ogs.add(tuple([cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]))
							limit_reached = True
						else:
							left_side_genes_and_ogs.add(tuple([cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]))
					except:
						near_scaffold_edge = True
						limit_reached = True
					cds_iter_index -= 1

				limit_reached = False
				cds_iter_index = cds_index+1
				while not limit_reached:
					try:
						cds_iter = scaffold_features_sorted[cds_iter_index]
						if cds_iter_index > max_features-1:
							near_scaffold_edge = True
							limit_reached = True
						elif cds_iter[0] > right_boundary:
							if cds_iter[3] == '0':
								near_scaffold_edge = True
							limit_reached = True
						elif cds_iter[1] > right_boundary:
							if cds_iter[3] == '0':
								near_scaffold_edge = True
							right_side_genes_and_ogs.add(tuple([cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]))
							limit_reached = True
						else:
							right_side_genes_and_ogs.add(tuple([cds_iter[2], gene_to_og[cds_iter[2]], cds_iter[0]]))
					except:
						near_scaffold_edge = True
						limit_reached = True
					cds_iter_index += 1

				upstream_genes = []
				downstream_genes = []
				upstream_ogs = []
				downstream_ogs = []
				if cds[2] == '1':
					for go in sorted(left_side_genes_and_ogs, key=itemgetter(2)):
						upstream_genes.append(go[0])
						upstream_ogs.append(go[1])
					for go in sorted(right_side_genes_and_ogs,key=itemgetter(2)):
						downstream_genes.append(go[0])
						downstream_ogs.append(go[1])
				else:
					for go in sorted(right_side_genes_and_ogs, key=itemgetter(2)):
						upstream_genes.append(go[0])
						upstream_ogs.append(go[1])
					for go in sorted(left_side_genes_and_ogs,key=itemgetter(2)):
						downstream_genes.append(go[0])
						downstream_ogs.append(go[1])
				output_handle.write('\t'.join([cds_name, cds_og, str(near_scaffold_edge), str(abs(cds_end-cds_start+1)), 
								   ', '.join(upstream_genes), ', '.join(downstream_genes), 
								   ', '.join(upstream_ogs), ', '.join(downstream_ogs)]) + '\n')
		output_handle.close()

	except Exception as e:
		logObject.error("Problem with determining context of genes for sample %s" % sample)
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determineOrthologGroupContexts(prepare_input_dir, og_context_info_file, surround_info_dir, of_tsv, of_singletons_tsv, logObject, surrounding_bp=10000, cpus=1):
	"""
	Description:
	This function loads coordinate information of CDSs for each input genome into a dictionary and also 
	********************************************************************************************************************
	Parameters:
	- prepare_input_dir: The results directory from a prior run of prepare_input for preparing genomes for bofasa.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	********************************************************************************************************************
	"""
	try:
		og_genes = defaultdict(set)
		og_samples = defaultdict(set)
		gene_to_og = {}
		if os.path.isfile(of_tsv):
			samples = []
			with open(of_tsv) as oot:
				for i, line in enumerate(oot):
					line = line.strip('\n')
					ls = line.split('\t')
					if i == 0:
						samples = ls[1:]
					else:
						og = ls[0]
						for j, gs in enumerate(ls[1:]):
							sample = samples[j]
							for g in gs.split(', '):
								g = g.strip()
								if g != '':
									og_genes[og].add(g)
									og_samples[og].add(sample)
									gene_to_og[g] = og

		if os.path.isfile(of_singletons_tsv):
			samples = []
			with open(of_singletons_tsv) as oost:
				for i, line in enumerate(oost):
					line = line.strip('\n')
					ls = line.split('\t')
					if i == 0: 
						samples = ls[1:]
					else:
						og = ls[0]
						for j, gs in enumerate(ls[1:]):
							for g in gs.split(', '):
								if g != '':
									og_genes[og].add(g)
									og_samples[og].add(sample)
									gene_to_og[g] = og

		try:
			assert(len(og_genes) > 0)
		except Exception as e:
			logObject.error("Difficulties parsing OrthoFinder results! Probably it didn't run successfully :(")
			logObject.error(traceback.format_exc())
			sys.stderr.write("Difficulties parsing OrthoFinder results! Probably it didn't run successfully :(\n")
			sys.stderr.write(traceback.format_exc() + '\n')

		listing_file = prepare_input_dir + 'Genome_Proteomes_Files.txt'
		assert(os.path.isfile(listing_file))
		genome_params = []
		with open(listing_file) as olf:
			for line in olf:
				line = line.strip()
				sample, proteome_file, coords_file = line.split('\t')
				coords_file = prepare_input_dir + coords_file
				output_file = surround_info_dir + sample + '.tsv'
				genome_params.append([sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, logObject])

		p = multiprocessing.Pool(cpus)
		p.map(extractGeneContexts, genome_params)
		p.close()

		og_surrounding_nogs = defaultdict(lambda: defaultdict(int))
		og_completed_surrounding_nogs = defaultdict(lambda: defaultdict(int))
		og_contexts = defaultdict(list)
		og_contexts_nses = defaultdict(int)
		og_gene_lengths = defaultdict(list)
		context_nogs = defaultdict(list)
		context_nogs_complete = defaultdict(list)
		for genome_info in genome_params:
			sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, logObject = genome_info
			with open(output_file) as oof:
				for i, line in enumerate(oof):
					if i == 0: continue
					line = line.strip('\n')
					protein, og, nse, length, upstream_genes, downstream_genes, upstream_ogs, downstream_ogs = line.split('\t')
					og_gene_lengths[og].append(float(length))
					complete_status = 'C:'
					if nse == 'True':
						og_contexts_nses[og] += 1
						complete_status = 'I:'
					context = complete_status
					if len(upstream_ogs) > 0:
						context += upstream_ogs + ', '
					context += og
					if len(downstream_ogs) > 0:
						context += ', ' + downstream_ogs
					
					og_contexts[og].append(context)
					context_ogs = set([])
					complete_context_ogs = set([])
					for ogc in upstream_ogs.split(', '):
						if nse == 'False':
							complete_context_ogs.add(ogc)
						if not ogc in context_ogs:
							og_surrounding_nogs[og][ogc] += 1
							if nse == 'False':
								og_completed_surrounding_nogs[og][ogc] += 1
						context_ogs.add(ogc)
					for ogc in downstream_ogs.split(', '):
						if nse == 'False':
							complete_context_ogs.add(ogc)
						if not ogc in context_ogs:
							og_surrounding_nogs[og][ogc] += 1
							if nse == 'False':
								og_completed_surrounding_nogs[og][ogc] += 1
						context_ogs.add(ogc)
					context_nogs[og].append(len(context_ogs))
					context_nogs_complete[og].append(len(context_ogs))

		og_context_info_handle = open(og_context_info_file, 'w')
		og_context_info_handle.write('\t'.join(['OG', 'Median OG length (bp)', 'Proportion contexts near scaffold edge', 'Number of genomes with OG', 
										        'Number of instances', 'Context variability score', 'Context variability score - complete contexts', 
												'Context entropy score', 'Context entropy score - complete contexts', 'Number of distinct neighbor OGs', 
												'Number of distinct OGs from complete contexts', 'Avg. number of distinct neighbor OGs', 
												'Avg. number of distinct neighbor OGs from complete contexts', 'Contexts']) + '\n')
		for og in sorted(og_contexts):
			median_length = og_gene_lengths[og][0]
			if len(og_gene_lengths) > 1:
				median_length = statistics.median(og_gene_lengths[og])
			num_samples = len(og_samples[og])
			num_contexts = len(og_contexts[og])
			prop_nse = round(og_contexts_nses[og]/float(num_contexts),2)
			nog_freqs = []
			total_nog = 0
			nog_freqs_complete = []
			total_nog_complete = 0
			for nog in og_surrounding_nogs[og]:
				nog_freqs.append(og_surrounding_nogs[og][nog])
				total_nog += 1
				if nog in og_completed_surrounding_nogs[og]:
					nog_freqs_complete.append(og_completed_surrounding_nogs[og][nog])
					total_nog_complete += 1
			avg_nog = round(statistics.mean(context_nogs[og]),2)
			avg_nog_complete = round(statistics.mean(context_nogs_complete[og]),2)
			context_entropy = 'NA'
			context_var_score = 'NA'
			context_entropy_complete = 'NA'
			context_var_score_complete = 'NA'
			if total_nog > 0:
				context_var_score = round(avg_nog/total_nog, 2)
				if total_nog > 1:
					context_entropy = round(stats.entropy([x/total_nog for x in nog_freqs]), 2)
			if total_nog_complete > 0:
				context_var_score_complete = round(avg_nog_complete/total_nog_complete, 2)
				if total_nog_complete > 1:
					context_entropy_complete = round(stats.entropy([x/total_nog_complete for x in nog_freqs_complete]), 2)

			og_context_info_handle.write('\t'.join([str(x) for x in [og, round(median_length,2), prop_nse, num_samples, num_contexts, context_var_score, 
								         context_var_score_complete, context_entropy, context_entropy_complete, total_nog, total_nog_complete,
										 avg_nog, avg_nog_complete, '; '.join(og_contexts[og])]]) + '\n')
		og_context_info_handle.close()
			
	except Exception as e:
		logObject.error("Problem determining contexts of ortholog groups. Exiting now...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)



def is_fasta(fasta):
	"""
	Description:
	Function to validate if a file is actually a FASTA file.
	********************************************************************************************************************
	Parameters:
	- fasta: A file that should be in FASTA format.
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether file is in FASTA format.
	********************************************************************************************************************
	"""
	try:
		recs = 0
		if fasta.endswith('.gz'):
			with gzip.open(fasta, 'rt') as ogf:
				for rec in SeqIO.parse(ogf, 'fasta'):
					recs += 1
					break
		else:
			with open(fasta) as of:
				for rec in SeqIO.parse(of, 'fasta'):
					recs += 1
					break
		if recs > 0:
			return True
		else:
			return False
	except:
		return False

def is_genbank(gbk, check_for_cds=False):
	"""
	Description:
	Function to validate if a file is actually a GenBank file.
	********************************************************************************************************************
	Parameters:
	- gbk: A file that should be in GenBank format.
	- check_for_cds: Whether to also check that the GenBank contains CDS features.
	********************************************************************************************************************
	Returns:
	- True or False statement depending on whether file is in GenBank format.
	********************************************************************************************************************
	"""
	try:
		recs = 0
		cds_flag = False
		assert (gbk.endswith('.gbk') or gbk.endswith('.gbff') or gbk.endswith('.gbk.gz') or gbk.endswith('.gbff.gz'))
		if gbk.endswith('.gz'):
			with gzip.open(gbk, 'rt') as ogf:
				for rec in SeqIO.parse(ogf, 'genbank'):
					if check_for_cds:
						for feature in rec.features:
							if feature.type == 'CDS':
								cds_flag = True
					if not check_for_cds or cds_flag:
						recs += 1
						break
		else:
			with open(gbk) as ogf:
				for rec in SeqIO.parse(ogf, 'genbank'):
					if check_for_cds:
						for feature in rec.features:
							if feature.type == 'CDS':
								cds_flag = True
					if not check_for_cds or cds_flag:
						recs += 1
						break
		if recs > 0:
			return True
		else:
			return False
	except:
		return False

def createProteinAlignments(prot_dir, prot_algn_dir, logObject, use_super5=True, cpus=1):
	"""
	Description:
	This function creates protein alignments from a directory of protein sequences.

	Function originally developed for zol.
	*******************************************************************************************************************
	Parameters:
	- prot_dir: A directory containing protein sequences.
	- prot_algn_dir: A directory to write protein alignments.
	- logObject: A logging object.
	- use_super5: Whether to use the SUPER5 algorithm for MUSCLE alignment.
	- cpus: The number of CPUs to use for alignment.
	"""
	try:
		for pf in os.listdir(prot_dir):
			prefix = '.fa'.join(pf.split('.fa')[:-1])
			prot_file = prot_dir + pf
			prot_algn_file = prot_algn_dir + prefix + '.msa.faa'
			align_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(cpus)]
			if use_super5:
				align_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(cpus)]
			try:
				subprocess.call(' '.join(align_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
								executable='/bin/bash')
				assert(os.path.isfile(prot_algn_file))
				logObject.info('Successfully ran: %s' % ' '.join(align_cmd))
			except Exception as e:
				logObject.error('Had an issue running MUSCLE: %s' % ' '.join(align_cmd))
				sys.stderr.write('Had an issue running MUSCLE: %s\n' % ' '.join(align_cmd))
				logObject.error(e)
				sys.exit(1)
	except Exception as e:
		sys.stderr.write('Issues with creating protein alignments.\n')
		logObject.error('Issues with creating protein alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createProfileHMMsAndConsensusSeqs(prot_algn_dir, phmm_dir, cons_dir, logObject, cpus=1):
	"""
	Description:
	This function creates profile HMMs and emits consensus sequences based on protein MSAs using HMMER.
	
	Function originally developed for zol.
	*******************************************************************************************************************
	Parameters:
	- prot_algn_dir: The directory containing protein alignments.
	- phmm_dir: The directory where profile HMMs in HMMER3 HMM format will be saved.
	- cons_dir: The directory where consensus sequences in FASTA format will be saved.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	*******************************************************************************************************************
	"""
	try:
		hmmbuild_cmds = []
		hmmemit_cmds = []
		for paf in os.listdir(prot_algn_dir):
			prefix = '.msa.faa'.join(paf.split('.msa.faa')[:-1])
			prot_algn_file = prot_algn_dir + paf
			prot_hmm_file = phmm_dir + prefix + '.hmm'
			prot_cons_file = cons_dir + prefix + '.cons.faa'
			hmmbuild_cmds.append(['hmmbuild', '--amino', '--cpu', '2', '-n', prefix, prot_hmm_file, prot_algn_file, logObject])
			hmmemit_cmds.append(['hmmemit', '-c', '-o', prot_cons_file, prot_hmm_file, logObject])
		p = multiprocessing.Pool(cpus)
		p.map(multiProcess, hmmbuild_cmds)
		p.map(multiProcess, hmmemit_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating profile HMMs and consensus sequences.\n')
		logObject.error('Issues with creating profile HMMs and consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def concatenateConsensusAlignment(og_cons_dir, orthofinder_tsv_singletons_file, prepare_input_dir, concatenated_consensus_seqs_file, logObject):
	"""
	"""
	try:
		concatenated_consensus_seqs_handle = open(concatenated_consensus_seqs_file, 'w')
		for f in os.listdir(og_cons_dir):
			cons_seq_file = og_cons_dir + f
			og = f.split('.cons.faa')[0]
			with open(cons_seq_file) as ocsf:
				for rec in SeqIO.parse(ocsf, 'fasta'):
					concatenated_consensus_seqs_handle.write('>' + og + '\n' + str(rec.seq) + '\n')

		# TODO: Consider making look up by sample first to make faster - this is just safer in case orthofinder changes names
		gene_to_og = {}
		with open(orthofinder_tsv_singletons_file) as otsf:
			for i, line in enumerate(otsf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i > 0:
					og = ls[0]
					for gs in ls[1:]:
						for g in gs.split(', '):
							g = g.strip()
							if g != '':
								gene_to_og[g] = og
		
		orthofinder_input_dir = prepare_input_dir + 'Formatted_Proteomes/'
		for f in os.listdir(orthofinder_input_dir):
			sample_proteome = orthofinder_input_dir + f
			if f.endswith('.faa'):
				with open(sample_proteome) as osp:
					for rec in SeqIO.parse(osp, 'fasta'):
						if rec.id in gene_to_og:
							concatenated_consensus_seqs_handle.write('>' + gene_to_og[rec.id] + '\n' + str(rec.seq).replace('*', '') + '\n')
		concatenated_consensus_seqs_handle.close()

	except Exception as e:
		sys.stderr.write('Issues with concatenating consensus sequences of multi-sequence ortholog groups and singleton protein sequences.\n')
		logObject.error('Issues with concatenating consensus sequences of multi-sequence ortholog groups and singleton protein sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)


def trimAlignments(prot_algn_dir, codo_algn_dir, prot_algn_trim_dir, codo_algn_trim_dir, logObject, cpus=1):
	"""
	Description:
	This function trims protein and codon alignments using TrimAl.
	*******************************************************************************************************************
	Parameters:
	- prot_algn_dir: The directory containing the protein alignments.
	- codo_algn_dir: The directory containing the codon alignments.
	- prot_algn_trim_dir: The directory where the trimmed protein alignments will be saved.
	- codo_algn_trim_dir: The directory where the trimmed codon alignments will be saved.
	- logObject: A logging object.
	- cpus: The number of CPUs to use for trimming the alignments.
	*******************************************************************************************************************
	"""
	try:
		trim_cmds = []
		for paf in os.listdir(prot_algn_dir):
			prefix = '.msa.faa'.join(paf.split('.msa.faa')[:-1])
			prot_algn_file = prot_algn_dir + paf
			prot_algn_trim_file = prot_algn_trim_dir + paf
			codo_algn_file = codo_algn_dir + prefix + '.msa.fna'
			codo_algn_trim_file = codo_algn_trim_dir + prefix + '.msa.fna'
			trim_cmds.append(['trimal', '-in', prot_algn_file, '-out', prot_algn_trim_file, '-keepseqs', '-gt', '0.9', logObject])
			trim_cmds.append(['trimal', '-in', codo_algn_file, '-out', codo_algn_trim_file, '-keepseqs', '-gt', '0.9', logObject])
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, trim_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with trimming protein/codon alignments.\n')
		logObject.error('Issues with trimming protein/codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createGeneTrees(codo_algn_trim_dir, tree_dir, logObject, cpus=1):
	"""
	Description:
	This function creates gene trees from trimmed codon alignments using FastTree2 for ortholog groups.
	*******************************************************************************************************************
	Parameters:
	- codo_algn_trim_dir: The directory containing trimmed codon alignments.
	- tree_dir: The directory where trees in Newick format will be saved.
	- logObject: A logging object.
	- cpus: The number of CPUs to use.
	*******************************************************************************************************************
	"""
	try:
		fasttree_cmds = []
		for catf in os.listdir(codo_algn_trim_dir):
			if not catf.endswith('.msa.fna'): continue
			prefix = '.msa.faa'.join(catf.split('.msa.fna')[:-1])
			codo_algn_trim_file = codo_algn_trim_dir + catf
			tree_file = tree_dir + prefix + '.tre'
			fasttree_cmds.append(['fasttree', '-nt', codo_algn_trim_file, '>', tree_file, logObject])
		p = multiprocessing.Pool(cpus)
		p.map(util.multiProcess, fasttree_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating gene-trees.\n')
		logObject.error('Issues with creating gene-trees.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)