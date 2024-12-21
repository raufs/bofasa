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
import tqdm
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
import decimal 
import pyhmmer

version = pkg_resources.require("bofasa")[0].version

def runCmd(cmd, logObject, check_files=[], check_directories=[], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL):
	if logObject != None:
		logObject.info('Running %s' % ' '.join(cmd))
	try:
		subprocess.call(' '.join(cmd), shell=True, stdout=stdout, stderr=stderr,
						executable='/bin/bash')
		for cf in check_files:
			assert (os.path.isfile(cf))
		for cd in check_directories:
			assert (os.path.isdir(cd))
		if logObject != None:
			logObject.info('Successfully ran: %s' % ' '.join(cmd))
	except:
		if logObject != None:
			logObject.error('Had an issue running: %s' % ' '.join(cmd))
			logObject.error(traceback.format_exc())
		raise RuntimeError('Had an issue running: %s' % ' '.join(cmd))

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


def createLocusTagOptions(locus_tag_length):
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=3))])))
		return(possible_locustags)
	except:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def processGenomesUsingProdigal(sample_genomes, prodigal_outdir, prodigal_proteomes, prodigal_genbanks, logObject,
								threads=1, locus_tag_length=3, gene_calling_method="pyrodigal", meta_mode=False):
	"""
	Description:
	This function oversees processing of input genomes to create proteome and GenBank files using p(y)rodigal.
	********************************************************************************************************************
	Parameters:
	- sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in FASTA format.
	- prodigal_outdir: Workspace where prodigal (intermediate) results should be written to directly.
	- prodigal_proteomes: Directory where final proteome files (in FASTA format) for target genomes will be saved.
	- prodigal_genbanks_directory: Directory where final GenBank files for target genomes will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- locus_tag_length: The length of the locus tags to generate.
	- gene_calling_method: Whether to use pyrodigal (default), prodigal, or prodigal-gv.
	- meta_mode: Whether to run pyrodigal/prodigal in metagenomics mode.
	- avoid_locus_tags: Whether to avoid using certain locus tags.
	********************************************************************************************************************
	"""
	try:
		possible_locustags = createLocusTagOptions(locus_tag_length)

		prodigal_cmds = []
		for i, sample in enumerate(sorted(sample_genomes)):
			sample_assembly = sample_genomes[sample]
			sample_locus_tag = ''.join(list(possible_locustags[i]))

			prodigal_cmd = ['runProdigalAndMakeProperGenbank.py', '-i', sample_assembly, '-s', sample, '-gcm', gene_calling_method,
							'-l', sample_locus_tag, '-o', prodigal_outdir]
			if meta_mode:
				prodigal_cmd += ['-m']
			prodigal_cmds.append(prodigal_cmd + [logObject])

		msg = "Running %s for %d genomes" % (gene_calling_method, len(prodigal_cmds)) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, prodigal_cmds), total=len(prodigal_cmds)):
			pass
		p.close()

		for sample in sample_genomes:
			try:
				assert (os.path.isfile(prodigal_outdir + sample + '.faa') and os.path.isfile(
					prodigal_outdir + sample + '.gbk'))
				os.system('mv %s %s' % (prodigal_outdir + sample + '.gbk', prodigal_genbanks))
				os.system('mv %s %s' % (prodigal_outdir + sample + '.faa', prodigal_proteomes))
			except:
				sys.stderr.write("Unable to validate successful genbank/predicted-proteome creation for sample %s\n" % sample)
				sys.stderr.write(traceback.format_exc())
				sys.exit(1)
	except Exception as e:
		logObject.error(
			"Problem with creating commands for running prodigal via script runProdigalAndMakeProperGenbank.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
		
def processGenomesAsGenbanks(sample_genomes, proteomes_directory, genbanks_directory, gene_name_mapping_outdir,
							 logObject, threads=1, locus_tag_length=3, avoid_locus_tags=set([]),
							 rename_locus_tags=False, rename_problem_gbks=False):
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
	- threads: The number of threads to use.
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

		for i, sample in enumerate(sorted(sample_genomes)):
			sample_locus_tag = possible_locustags[i]
			sample_genbank = sample_genomes[sample]
			process_cmd = ['processNCBIGenBank.py', '-i', sample_genbank, '-s', sample, 
						   '-g', genbanks_directory, '-p', proteomes_directory, '-n', 
						   gene_name_mapping_outdir]
			if rename_problem_gbks:
				process_cmd += ['-r', '-l', sample_locus_tag]
			elif rename_locus_tags:
				process_cmd += ['-l', sample_locus_tag]
			process_cmds.append(process_cmd + [logObject])
				
		msg = "Attempting to process/re-format %d genomes provided as GenBank files" % len(process_cmds) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, process_cmds), total=len(process_cmds)):
			pass
		p.close()

		successfully_processed = 0
		for sample in sample_genomes:
			faa_file = proteomes_directory + sample + '.faa'
			gbk_file = genbanks_directory + sample + '.gbk'
			map_file = gene_name_mapping_outdir + sample + '.txt'
			try:
				assert (os.path.isfile(faa_file) and os.path.isfile(gbk_file) and os.path.isfile(map_file))
				assert (os.path.getsize(faa_file) > 0 and os.path.getsize(gbk_file) > 0 and os.path.getsize(map_file) > 0)
				sample_genomes_updated[sample] = genbanks_directory + sample + '.gbk'
				successfully_processed += 1
			except AssertionError:
				if os.path.isfile(faa_file):
					os.system('rm -f ' + faa_file)
				if os.path.isfile(gbk_file):
					os.system('rm -f ' + gbk_file)
				if os.path.isfile(map_file):
					os.system('rm -f ' + map_file)
				sys.stderr.write("Unable to validate successful genbank reformatting/predicted-proteome creation for sample %s\n" % sample)
				pass

		sys.stdout.write('Successfully processed %s genomes!\n' % successfully_processed)
	
	except Exception as e:
		logObject.error("Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return sample_genomes_updated

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

def determineOrthologGroupContexts(prepare_input_dir, og_context_info_file, surround_info_dir, of_tsv, of_singletons_tsv, logObject, surrounding_bp=10000, threads=1):
	"""
	Description:
	This function loads coordinate information of CDSs for each input genome into a dictionary and also 
	********************************************************************************************************************
	Parameters:
	- prepare_input_dir: The results directory from a prior run of prepare_input for preparing genomes for bofasa.
	- logObject: A logging object.
	- threads: The number of threads to use.
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

		p = multiprocessing.Pool(threads)
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

def createProteinAlignments(prot_dir, prot_algn_dir, logObject, use_super5=True, threads=1):
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
	- threads: The number of threads to use for alignment.
	"""
	try:
		for pf in os.listdir(prot_dir):
			prefix = '.fa'.join(pf.split('.fa')[:-1])
			prot_file = prot_dir + pf
			prot_algn_file = prot_algn_dir + prefix + '.msa.faa'
			align_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(threads)]
			if use_super5:
				align_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-amino', '-threads', str(threads)]
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

def createProfileHMMsAndConsensusSeqs(prot_algn_dir, phmm_dir, cons_dir, logObject, threads=1):
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
	- threads: The number of threads to use.
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
		p = multiprocessing.Pool(threads)
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


def trimAlignments(prot_algn_dir, codo_algn_dir, prot_algn_trim_dir, codo_algn_trim_dir, logObject, threads=1):
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
	- threads: The number of threads to use for trimming the alignments.
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
		p = multiprocessing.Pool(threads)
		p.map(util.multiProcess, trim_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with trimming protein/codon alignments.\n')
		logObject.error('Issues with trimming protein/codon alignments.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createGeneTrees(codo_algn_trim_dir, tree_dir, logObject, threads=1):
	"""
	Description:
	This function creates gene trees from trimmed codon alignments using FastTree2 for ortholog groups.
	*******************************************************************************************************************
	Parameters:
	- codo_algn_trim_dir: The directory containing trimmed codon alignments.
	- tree_dir: The directory where trees in Newick format will be saved.
	- logObject: A logging object.
	- threads: The number of threads to use.
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
		p = multiprocessing.Pool(threads)
		p.map(util.multiProcess, fasttree_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating gene-trees.\n')
		logObject.error('Issues with creating gene-trees.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determinePhagesAndPlasmids(sample_wgs, sample_beds, genomad_dir, phage_protein_listing_file, plasmid_protein_listing_file, logObject, threads=1, genome_splits=8):
	"""
	Description:
	This function runs annotation of sample genomes for phages and plasmis using geNomad to determine proteins 
	belonging to such elements. 
	*******************************************************************************************************************
	Parameters:
	- sample_wgs: Dictionary mapping sample names (keys) to genome file paths (values).
	- sample_beds: Dictionary mapping sample names (keys) to gene-coordinate describing bed file paths (values).
	- genomad_dir: Directory where to store geNomad results.
	- phage_protein_listing_file: Final resulting file with information of which proteins are predicted to exist on 
								  (pro)phages.
	- plasmid_protein_listing_file: Final resulting file with information of which proteins are predicted to exist on 
								  plasmids.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- genome_splits: The number of splits for geNomad to limit memory usage.
	*******************************************************************************************************************
	"""
	bofasa_db_dir = str(os.getenv("BOFASA_DB_PATH")).strip()
	db_locations = None
	try:
		bofasa_db_dir = os.path.abspath(bofasa_db_dir) + '/'
		db_locations = bofasa_db_dir + 'database_location_paths.txt'
		assert(os.path.isfile(db_locations))
	except:
		pass
	if db_locations == None or not os.path.isfile(db_locations):
		msg = 'Databases do not appear to be setup or setup properly! Please run setup_annotation_dbs.py prior to run bofasa_prep, exiting ...'
		sys.stderr.write('Error: ' + msg + '\n')
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())		

	try:
		genomad_db_dir = None
		with open(db_locations) as odb:
			for line in odb:
				line = line.strip()
				ls = line.split('\t')
				if ls[0] == 'genomad':
					genomad_db_dir = ls[2]
		assert(os.path.isdir(genomad_db_dir))

		phpf_handle = open(phage_protein_listing_file, 'w')
		plpf_handle = open(plasmid_protein_listing_file, 'w')

		for sample in sample_wgs:
			input_genome = sample_wgs[sample]
			genomad_results = genomad_dir + sample + '/'
			genomad_cmd = ['genomad', 'end-to-end', '--cleanup', '--threads', str(threads), '--splits', 
				 	       str(genome_splits), input_genome, genomad_results, genomad_db_dir]
			prophage_coords_tsv = None
			plasmid_coords_tsv = None
			try:
				runCmd(genomad_cmd, None, check_directories=[genomad_results])

				for subdir, dirs, files in os.walk(genomad_results):
					for file in files:
						filepath = subdir + os.sep + file
						if filepath.endswith("_plasmid_summary.tsv"):
							plasmid_coords_tsv = filepath
						elif filepath.endswith("_virus_summary.tsv"):
							prophage_coords_tsv = filepath
				assert(prophage_coords_tsv != None and os.path.isfile(prophage_coords_tsv))
				assert(plasmid_coords_tsv != None and os.path.isfile(plasmid_coords_tsv))
			except:
				msg = 'Issue with running the command: %s.' % ' '.join(genomad_cmd)
				sys.stderr.write(msg + '\n')
				logObject.error(msg)
				sys.stderr.write(traceback.format_exc() + '\n')
				sys.exit(1)			

			full_plasmid_scaffs = set([])
			full_phage_scaffs = set([])
			phage_coords = defaultdict(set)
			with open(prophage_coords_tsv) as opaf:
				for i, line in enumerate(opaf):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					if ls[3] == 'NA':
						scaffold = ls[0]
						full_phage_scaffs.add(scaffold)
					else:
						scaffold = '|'.join(ls[0].split('|')[:-1])
						start = int(ls[3].split('-')[0])
						end = int(ls[3].split('-')[1])
						for pos in range(start, end+1):
							phage_coords[scaffold].add(pos)

			with open(plasmid_coords_tsv) as opaf:
				for i, line in enumerate(opaf):
					if i == 0: continue
					line = line.strip()
					ls = line.split('\t')
					scaffold = ls[0]
					full_plasmid_scaffs.add(scaffold)

			bed_file = sample_beds[sample]
			with open(bed_file) as obf:
				for line in obf:
					line = line.strip()
					scaffold, start, end, final_lt, prot_score, direction = line.split('\t')
					prot_coords = set(range(start, end+1))

					if (scaffold in full_phage_scaffs) or (scaffold in phage_coords and len(prot_coords.intersection(phage_coords[scaffold])) > 0):
						phpf_handle.write(sample + '\t' + final_lt + '\n')
					if scaffold in full_plasmid_scaffs:
						plpf_handle.write(sample + '\t' + final_lt + '\n')
		phpf_handle.close()
		plpf_handle.close()
	except:
		msg = 'Issue with determining phages/plasmids in assemblies using geNomad!'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(traceback.format_exc())
		sys.exit(1)		


def annotateIsFinder(sample_proteomes, annot_dir, isfinder_protein_listing_file, logObject, threads=1, max_annotation_evalue=1e-3):
	"""
	Description:
	This function runs annotation of sample proteomes using ISFinder IS elements via DIAMOND blastp.
	*******************************************************************************************************************
	Parameters:
	- sample_proteomes: Dictionary mapping sample names (keys) to proteome file paths (values).
	- annot_dir: Directory where to store DIAMOND BLASTp results.
	- isfinder_protein_listing_file: Final resulting file with information of which proteins exhibit homology to IS 
	                                 elements.
	- logObject: A logging object.
	- threads: The number of threads to use.
	- max_annotation_evalue: The maximum e-value to consider proteins as homologous to transposons.
	*******************************************************************************************************************
	"""
	bofasa_db_dir = str(os.getenv("BOFASA_DB_PATH")).strip()
	db_locations = None
	try:
		bofasa_db_dir = os.path.abspath(bofasa_db_dir) + '/'
		db_locations = bofasa_db_dir + 'database_location_paths.txt'
		assert(os.path.isfile(db_locations))
	except:
		pass
	if db_locations == None or not os.path.isfile(db_locations):
		msg = 'Databases do not appear to be setup or setup properly! Please run setup_annotation_dbs.py prior to run bofasa_prep, exiting ...'
		sys.stderr.write('Error: ' + msg + '\n')
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())		

	try:
		isfinder_dmnd_path = None
		with open(db_locations) as odb:
			for line in odb:
				line = line.strip()
				ls = line.split('\t')
				if ls[0] == 'isfinder':
					isfinder_dmnd_path = ls[2]
		assert(os.path.isfile(isfinder_dmnd_path))
		
		dmnd_search_cmds = []
		for sample in sample_proteomes:
			faa_file = sample_proteomes[sample]
			annotation_result_file = annot_dir + sample + '.isfinder_diamond_blastp.txt'
			search_cmd = ['diamond', 'blastp', '--ignore-warnings', '-p', str(1), '-d', isfinder_dmnd_path,
						  '-q', faa_file, '-o', annotation_result_file, logObject]
			dmnd_search_cmds.append(search_cmd)

		msg = "Running %d DIAMOND blastp jobs for IS element annotation" % len(dmnd_search_cmds)
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, dmnd_search_cmds), total=len(dmnd_search_cmds)):
			pass
		p.close()

		outf = open(isfinder_protein_listing_file, 'w')
		for rf in os.listdir(annot_dir):
			sample = rf.split('.isfinder_diamond_blastp.txt')[0]
			
			best_hits_by_bitscore = defaultdict(lambda: [[], [], 0.0])
			# parse DIAMOND BLASTp based results
			with open(annot_dir + rf) as oarf:
				for line in oarf:
					line = line.strip()
					ls = line.split('\t')
					query = ls[0]
					hit = ls[1]
					bitscore = float(ls[11])
					evalue = decimal.Decimal(ls[10])
					if evalue > max_annotation_evalue: continue
					if bitscore > best_hits_by_bitscore[query][2]:
						best_hits_by_bitscore[query] = [[hit], [evalue], bitscore]
					elif bitscore == best_hits_by_bitscore[query][2]:
						best_hits_by_bitscore[query][0].append(hit)
						best_hits_by_bitscore[query][1].append(evalue)

			for p in best_hits_by_bitscore:
				outf.write(sample + '\t' + p + '\t' + ', '.join(best_hits_by_bitscore[p][0]) + str(statistics.mean(best_hits_by_bitscore[p][1])))
		outf.close()

	except Exception as e:
		msg = 'Issue with annotating IS elements!'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.exit(1)

def split_by_idx(S, list_of_indices):
        """
        Function taken from https://stackoverflow.com/questions/10851445/splitting-a-string-by-list-of-indices
        """
        left, right = 0, list_of_indices[0]
        yield S[left:right]
        left = right
        for right in list_of_indices[1:]:
                yield S[left:right]
                left = right
        yield S[left:]

def createChoppedProteomes(inputs):
	"""
	Description:
	Create a chopped CDS GenBank file from a regular GenBank file - core function for batchCreateChoppedGenbanks().
	********************************************************************************************************************
	Parameters:
	- inputs:
		- prot_file: The original proteome file.
		- ccds_prot_file: The chopped up proteome FASTA file to create.
		- pfam_db_file: The Pfam HMM DB file.
		- pfam_z: The Pfam record count - for accurate E-value estimation.
		- minimal_length: The minimum length in amino acids for a domain matching or intra-domain region to be kept and
	                      tagged as a chopped CDS feature [Default is 20].
		- logObject: A logging object.
	- threads: The number of threads to use [Default is 1].
	********************************************************************************************************************
	"""
	prot_file, ccds_prot_file, pfam_db_file, pfam_z, minimal_length, logObject = inputs
	try:
		sample = '.'.join(prot_file.split('.')[:-1])

		# align Pfam domains and remove overlap similar to BiG-SCAPE
		alphabet = pyhmmer.easel.Alphabet.amino()
		sequences = []
		with pyhmmer.easel.SequenceFile(prot_file, digital=True, alphabet=alphabet) as seq_file:
			sequences = list(seq_file)

		target_dom_hits = defaultdict(list)
		with pyhmmer.plan7.HMMFile(pfam_db_file) as hmm_file:
			for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="trusted", Z=int(pfam_z), cpus=1):
				for hit in hits:
					for domain in hit.domains.included:
						target_dom_hits[hit.name.decode()].append([hits.query_name.decode(), domain.alignment.target_from, domain.alignment.target_to, domain.score, domain.i_evalue])

		# chop up FASTA based on mostly non-overlapping domains, 10% leaway is given
		breakpoints = defaultdict(list)
		dom_start_names = defaultdict(lambda: 'NA')
		for tg in target_dom_hits:
			tg_dom_name_iter = defaultdict(int)
			accounted_coords = set([])
			for dom_align_info in sorted(target_dom_hits[tg], key=itemgetter(3), reverse=True):
				dom_name, start, end, score, i_evalue = dom_align_info
				overlap_coords = accounted_coords.intersection(set(range(start, end+1)))
				if len(overlap_coords)/float(end-start+1) >= 0.1: continue
				accounted_coords = accounted_coords.union(set(range(start, end+1)))
				breakpoints[tg].append(start)
				breakpoints[tg].append(end+1)
				dom_start_names[tg + '|' + str(start)] = tg + '|' + dom_name + '|' + str(tg_dom_name_iter[dom_name]+1) 
				tg_dom_name_iter[dom_name] += 1

		cpf_handle = open(ccds_prot_file, 'w')
		with open(prot_file) as ocf:
			for rec in SeqIO.parse(ocf, 'fasta'):
				tg = rec.id
				tg_seq = str(rec.seq)
				prev_end_coord = 1
				tg_interdomain_index = 1
				if not tg in breakpoints and len(tg_seq) >= minimal_length:
					dn = sample + '|' + tg + '|full_protein|1'
					cpf_handle.write('>' + dn + '\n' + str(tg_seq) + '\n')
				else:
					for tg_seq_chunk in split_by_idx(tg_seq, ([0] + sorted(breakpoints[tg]))):
						if tg_seq_chunk.strip() == '': continue
						end_coord = prev_end_coord + len(tg_seq_chunk) - 1
						if len(tg_seq_chunk) >= minimal_length:
							dn = sample + '|' + dom_start_names[tg + '|' + str(prev_end_coord-1)]
							if dn == 'NA':
								dn = sample + '|' + tg + '|inter-domain_region|' + str(tg_interdomain_index)
								tg_interdomain_index += 1
							cpf_handle.write('>' + dn + '\n' + str(tg_seq) + '\n')
						prev_end_coord = end_coord + 1
		cpf_handle.close()
	except:
		msg = 'An issue occurred with creating chopped up version of proteome file %s.' % prot_file
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.exit(1)

def annotateAndSplitProteinsUsingPfam(sample_proteomes, split_proteins_dir, domain_coord_info_file, logObject, minimal_length=20, threads=1)
	"""	
	Description:
	Create chopped CDS GenBank files from regular GenBank input with CDS features.
	********************************************************************************************************************
	Parameters:
	- sample_proteomes: Dictionary mapping sample names (keys) to proteome file paths (values).
	- split_proteins_dir: Directory where to write resulting cCDS proteomes.
	- domain_coord_info_file: Resulting domain/inter-domain information file.
	- logObject: A logging object.
	- minimal_length: The minimum length in amino acids for a domain matching or intra-domain region to be kept and
	                  tagged as a chopped CDS feature
	- threads: The number of threads to use [Default is 1].
	********************************************************************************************************************
	Returns:
	- sample_ccds_proteomes: A dictionary mapping sample names (keys) to chopped proteome file paths (values).
	********************************************************************************************************************
	"""

	bofasa_db_dir = str(os.getenv("BOFASA_DB_PATH")).strip()
	db_locations = None
	try:
		bofasa_db_dir = os.path.abspath(bofasa_db_dir) + '/'
		db_locations = bofasa_db_dir + 'database_location_paths.txt'
		assert(os.path.isfile(db_locations))
	except:
		pass
	if db_locations == None or not os.path.isfile(db_locations):
		msg = 'Databases do not appear to be setup or setup properly! Please run setup_annotation_dbs.py prior to run bofasa_prep, exiting ...'
		sys.stderr.write('Error: ' + msg + '\n')
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())		

	try:
		pfam_hmm_path = None
		pfam_z = None
		with open(db_locations) as odb:
			for line in odb:
				line = line.strip()
				ls = line.split('\t')
				if ls[0] == 'pfam':
					pfam_hmm_path = ls[2]
					pfam_z = int(ls[3])
		assert(os.path.isfile(pfam_hmm_path) and pfam_z != None)

		prot_mod_inputs = []
		for sample in sample_proteomes:
			ccds_prot_file = split_proteins_dir + sample + '.ccds.faa'
			prot_file = sample_proteomes[sample]
			prot_mod_inputs.append([prot_file, ccds_prot_file, pfam_hmm_path, pfam_z, minimal_length, logObject])

		msg = "Creating domain-chopped up version of GenBank files for %d gene clusters" % len(prot_mod_inputs) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(createChoppedProteomes, prot_mod_inputs), total=len(prot_mod_inputs)):
			pass
		p.close()

		dci_handle = open(domain_coord_info_file, 'w')
		sample_ccds_proteomes = {}
		for f in os.listdir(split_proteins_dir):
			sample = f.split('.ccds.faa')[0]
			ccds_prot_file = split_proteins_dir + f
			sample_ccds_proteomes[sample] = ccds_prot_file
			with open(ccds_prot_file) as ogf:
				for rec in SeqIO.parse(ogf, 'fasta'):
					name = rec.id
					sample, prot, dom, index = name.split('|')
					dci_handle.write('\t'.join([sample, prot, dom, index, ccds_prot_file]) + '\n')
		dci_handle.close()
		return(sample_ccds_proteomes)
	except:
		msg = 'An issue occurred with creating chopped up proteomes.' 
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
