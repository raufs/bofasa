import os
import sys
from Bio import SeqIO
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
import resource
import pkg_resources  # part of setuptools
import statistics
from scipy import stats
import decimal 
import pyhmmer
from ete3 import Tree
import pandas as pd
import plotly.express as px
from scipy import stats
from scipy.spatial import distance
import random

random.seed(12345)

single_copy_dogs = set([])
largely_idr_dogs = set([])
protein_dogs = defaultdict(lambda: defaultdict(int))
dog_conservation = {}
tree_obj = None

version = pkg_resources.require("bofasa")[0].version

def castToNumeric(x):
	"""
	Description:
	This function attempts to cast a variable into a float. A special exception is whether "< 3 segregating sites!" is
	the value of the variable, which will simply be retained as a string.
	********************************************************************************************************************
	Parameters:
	- x: Input variable.
	********************************************************************************************************************
	Returns:
	- A float casting of the variable's value if numeric or "nan" if not.
	********************************************************************************************************************
	"""
	try:
		if x == '< 3 segregating sites!':
			return(x)
		else:
			x = float(x)
			return (x)
	except:
		return float('nan')

def generate_og_name(i):
	try:
		pid = None
		if (i+1) < 10:
			pid = 'OG00000'+str(i+1)
		elif (i+1) < 100:
			pid = 'OG0000'+str(i+1)
		elif (i+1) < 1000:
			pid = 'OG000'+str(i+1)
		elif (i+1) < 10000:
			pid = 'OG00'+str(i+1)
		elif (i+1) < 100000:
			pid = 'OG0'+str(i+1)
		else:
			pid = 'OG' + str(i+1)
		assert(pid != None)
		return(pid)
	except:
		msg = 'Issue generating OG identifier for number: %s' % str(i)
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')

def runSetOfProteinComparisons(inputs):
	try:
		input_listing_file, result_file, dj = inputs
		
		outf_handle = open(result_file, 'w')

		with open(input_listing_file) as oilf:
			for line in oilf:
				line = line.strip()
				p1, p2 = line.split('\t')

				p1dogs = protein_dogs[p1]
				p2dogs = protein_dogs[p2]
				intersect_dogs = (set(p1dogs.keys())).intersection(set(p2dogs.keys()))
				union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))
				sc_dogs = single_copy_dogs.intersection(intersect_dogs)

				threshold = dj
				if len(sc_dogs) >= 1:
					for sd in sc_dogs:
						sd_conservation = dog_conservation[sd]
						updated_threshold = dj - (dj*sd_conservation)
						if updated_threshold < threshold:
							threshold = updated_threshold

				union_count = 0
				intersect_count = 0
				for d in union_dogs:
					if d in largely_idr_dogs: continue
					union_count += p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
					intersect_count += min([p1dogs[d], p2dogs[d]])

				if union_count > 0:
					jaccard_index = intersect_count/union_count
					if jaccard_index >= threshold:
						outf_handle.write(p1 + '\t' + p2 + '\t' + str(jaccard_index*100.0) + '\n')
		outf_handle.close()
	except:
		msg = 'Issue performing pairwise assessment between proteins based on DOGs to determine protein-resolution ortholog groups.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')	

def splitNJT(input):
	og, cog_list_file, cog_dist_file, tre_file, spl_file, logObject, threads = input
	try:
		all_prots = set([])
		with open(cog_list_file) as oclf:
			for line in oclf:
				line = line.strip()
				all_prots.add(line)

		all_dogs = set([])
		for p in all_prots:
			for d in protein_dogs[p]:
				all_dogs.add(d)

		prot_dog_vectors = defaultdict(list)
		for p in all_prots:
			for d in sorted(all_dogs):
				prot_dog_vectors[p].append(protein_dogs[p][d])

		prot_dog_vectors_np = {}
		for p in prot_dog_vectors:
			prot_dog_vectors_np[p] = np.array(prot_dog_vectors[p])

		naming = {}
		outf = open(cog_dist_file, 'w')
		outf.write(str(len(prot_dog_vectors_np.keys())) + '\n')
		for p1i, p1 in enumerate(sorted(prot_dog_vectors_np)):
			p1_dist = []
			for p2 in sorted(prot_dog_vectors_np):
				dist = distance.cosine(prot_dog_vectors_np[p1], prot_dog_vectors_np[p2])
				p1_dist.append('{:.10f}'.format(dist))

			len_id = len(str(p1i))
			assert(10 > len_id)
			naming[str(p1i)] = p1
			name = str(p1i) + (' '*(10-len_id))
			outf.write(name + ' '.join(p1_dist) + '\n')
		outf.close()

		fastme_cmd = ['fastme', '-i', cog_dist_file, '-o', tre_file, '-T', str(threads)]
		runCmd(fastme_cmd, logObject, check_files=[tre_file])

		t = Tree(tre_file)

		samples_with_og = set([])
		leafs = set([])
		for n in t.traverse('postorder'):
			if n.is_leaf():
				n.name = naming[n.name]
				s = n.name.split('|')[0]
				samples_with_og.add(s)
				leafs.add(n.name)

		if len(leafs) < 3: return
		if len(samples_with_og) == 1: return

		R = t.get_midpoint_outgroup()
		t.set_outgroup(R)
		sp = recursive_splitting(t, samples_with_og, are_proteins=True)
		sp_refined = furtherSplitOutlierArtifactGroups(t, sp)
		sp_further_split = furtherSplitDisjointDOGPartitions(t, sp_refined)
		spl_outf = open(spl_file, 'w')
		for spi in sp_further_split:
			spl_outf.write(' '.join(sorted(spi)) + '\n')
		spl_outf.close()
		return
	except:
		msg = 'Issue with splitting protein ortholog group %s - based on domain cosine distances neighbor-joining tree.' % og
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(msg)
		sys.exit(1)

def determineProteinOrthogroup(dogs_file, protein_clustering_dir, ogs_file, logObject, dj=0.25, threads=1):
	# Add description to this this function
	"""
	Description:
	This function determines the protein ortholog groups based on the domain ortholog groups.
	********************************************************************************************************************	
	Parameters:
	- dogs_file: The input matrix file describing the membership of domain ortholog groups.
	- protein_clustering_dir: The directory where temporary files pertaining to protein clustering files will be stored.
	- ogs_file: The output file where the protein ortholog groups will be written.
	- logObject: The logger object for logging messages.
	- dj: The Jaccard index threshold for determining the similarity between two proteins.
	- threads: The number of threads to use for parallel processing.
	********************************************************************************************************************
	"""
	try:	
		input_dir = protein_clustering_dir + 'Comparison_Listings/'
		pairwise_dir = protein_clustering_dir + 'Protein_Pairs_Based_on_DOGs/' 
		p_clust_list_dir = protein_clustering_dir + 'Protein_Coarse_Cluster_Listings/'
		p_dist_dir = protein_clustering_dir + 'FASTME_Inputs/'
		p_tre_dir = protein_clustering_dir + 'Protein_DOG_Distance_Trees/'
		p_split_dir = protein_clustering_dir + 'Protein_Splitting/'

		setupReadyDirectory([input_dir, pairwise_dir, p_clust_list_dir, p_dist_dir, 
							 p_tre_dir, p_split_dir])
		pairwise_file = protein_clustering_dir + 'Protein_Pairs_Based_on_DOGs.txt'
		clusters_file = protein_clustering_dir + 'Protein_Clusters_Based_on_DOGs.txt'
		outf_handle = open(ogs_file, 'w')

		global single_copy_dogs
		global largely_idr_dogs
		global protein_dogs
		global dog_conservation

		num_batches = threads
		batch_handles = {}
		for batch in range(0, num_batches):
			batch_input_file = input_dir + str(batch) + '.txt'
			batch_handles[batch] = open(batch_input_file, 'w')

		all_proteins = set([])
		sample_count = None
		samples = []
		pair_count = 0

		with open(dogs_file) as odf:
			for i, line in enumerate(odf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0: 
					sample_count = len(ls[1:])
					samples = ls[1:]
					outf_handle.write(line + '\n')
				else:
					dog = ls[0]
					sc_flag = True
					sample_with = 0
					dog_lts = set([])
					tot = 0 
					idr = 0 
					for lts in ls[1:]:
						if ',' in lts:
							sc_flag = False
						for lt in lts.split(','):
							lt = lt.strip()
							if lt == '': continue
							sample_with += 1
							prot_id = '|'.join(lt.split('|')[:2])
							all_proteins.add(prot_id)
							protein_dogs[prot_id][dog] += 1
							dog_lts.add(prot_id)
							tot += 1
							if lt.split('|')[2] == 'inter-domain_region':
								idr += 1
					idr_prop = idr/float(tot)
					if idr_prop >= 0.8:
						largely_idr_dogs.add(dog)

					if sample_with == 1:
						sc_flag = False
					
					dog_conservation[dog] = sample_with/float(sample_count)
					if sc_flag:
						single_copy_dogs.add(dog)
					for j, p1 in enumerate(sorted(dog_lts)):
						for k, p2 in enumerate(sorted(dog_lts)):
							if j >= k: continue
							batch = pair_count % num_batches
							batch_input_handle = batch_handles[batch]
							batch_input_handle.write(p1 + '\t' + p2 + '\n')
							pair_count += 1

		pairwise_assessment_inputs = []
		for batch in range(0, threads):
			batch_handle = batch_handles[batch]
			batch_handle.close()

			input_listing_file = input_dir + str(batch) + '.txt'
			result_file = pairwise_dir + str(batch) + '.txt'
			pairwise_assessment_inputs.append([input_listing_file, result_file, dj])
		
		p = multiprocessing.Pool(num_batches)
		for _ in tqdm.tqdm(p.imap_unordered(runSetOfProteinComparisons, pairwise_assessment_inputs), total=num_batches):
			pass
		p.close()

		os.system('find %s -maxdepth 1 -type f | xargs cat >> %s' % (pairwise_dir, pairwise_file))

		#clust_cmd = ['slclust', '<', pairwise_file, '>', clusters_file]
		clust_cmd = ['mcl', pairwise_file, '--abc', '-I', '1.2', '-o', clusters_file, '-te', str(threads)]

		try:
			subprocess.call(' '.join(clust_cmd), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, executable='/bin/bash')
			assert (os.path.isfile(clusters_file))
			logObject.info('Successfully ran: %s' % ' '.join(clust_cmd))
		except Exception as e:
			logObject.error('Had an issue running concatenation: %s' % ' '.join(clust_cmd))
			sys.stderr.write('Had an issue running concatentation: %s\n' % ' '.join(clust_cmd))
			logObject.error(e)
			sys.exit(1)
		
		mcl_lists = []
		with open(clusters_file) as ocf:
			for line in ocf:
				line = line.strip()
				ls = line.split()
				mcl_lists.append(sorted(ls))

		protein_og_clusters = []
		large_protein_og_clusters = []
		paired_proteins = set([])
		split_nj_trees_input = []
		large_split_nj_trees_input = []
		for i, ls in enumerate(sorted(mcl_lists)):
			for p in ls: paired_proteins.add(p)
			sample_og_counts = defaultdict(int)
			for p in ls:
				s = p.split('|')[0]
				sample_og_counts[s] += 1
			og_counts = []
			for s in sample_og_counts:
				og_counts.append(sample_og_counts[s])
			max_og_count = max(og_counts)
			if max_og_count >= 2 and len(sample_og_counts) >= 2 and len(ls) >= 4:
				og_uniq_id = 'CoarseOG_' + str(i) 
				cog_list_file = p_clust_list_dir + og_uniq_id + '.txt'
				cog_dist_file = p_dist_dir + og_uniq_id + '.phylip'
				og_split_file = p_split_dir + og_uniq_id + '.txt'
				og_tre_file = p_tre_dir + og_uniq_id + '.tre'
				cl_handle = open(cog_list_file, 'w')
				for p in ls:
					cl_handle.write(p + '\n')
				cl_handle.close()
				if len(ls) > 400:
					large_split_nj_trees_input.append([og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, 
														og_split_file, logObject, threads])
				else:
					split_nj_trees_input.append([og_uniq_id, cog_list_file, cog_dist_file, og_tre_file, 
													og_split_file, logObject, 1])
				large_protein_og_clusters.append(ls)
			else:
				protein_og_clusters.append(ls)

		p = multiprocessing.Pool(1)
		for _ in tqdm.tqdm(p.imap_unordered(splitNJT, large_split_nj_trees_input), total=len(large_split_nj_trees_input)):
			pass
		p.close()

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(splitNJT, split_nj_trees_input), total=len(split_nj_trees_input)):
			pass
		p.close()

		accounted_for_in_splitting = set([])
		for f in os.listdir(p_split_dir):
			split_listing_file = p_split_dir + f 
			with open(split_listing_file) as oslf:
				for line in oslf:
					line = line.strip()
					ls = line.split()
					for p in ls:
						accounted_for_in_splitting.add(p)
					protein_og_clusters.append(ls)

		for pog_full in large_protein_og_clusters:
			missing = []
			for p in pog_full:
				if not p in accounted_for_in_splitting:
					missing.append(p)
			if len(missing) > 0:
				protein_og_clusters.append(missing)

		for prot in all_proteins:
			if not prot in paired_proteins:
				protein_og_clusters.append([prot])

		for i, c in enumerate(sorted(protein_og_clusters)):
			samp_lts = defaultdict(list)
			for p in c:
				s = p.split('|')[0]
				samp_lts[s].append(p.split('|')[1])
			og_id = generate_og_name(i)
			printlist = [og_id]
			for s in samples:
				printlist.append(', '.join(sorted(samp_lts[s])))
			outf_handle.write('\t'.join(printlist) + '\n')

		outf_handle.close()
	except:
		msg = 'Issues splitting domain ortholog groups from OrthoFinder using phylogenetic processing.'
		logObject.error(msg)
		sys.stderr.write(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def determine_tree_score(intree, all_og_samples, are_proteins=False):
	try:
		sample_lts = defaultdict(set)
		for n in intree.traverse('postorder'):
			if n.is_leaf():
				s = n.name.split('|')[0]
				if are_proteins:
					sample_lts[s].add(n.name)
				else:
					sample_lts[s].add('|'.join(n.name.split('|')[:-2]))

		curr_score = 0
		for s in all_og_samples:
			curr_score += abs(len(sample_lts[s]) - 1)

		return(curr_score)
	except:
		msg = 'Issues splitting ortholog group tree using the recursive function.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)


def get_children(intree):
	try:
		children = set([])
		for n in intree.traverse('postorder'):
			if n.is_leaf():
				if n.name.strip() != '':
					children.add(n.name)
		return(children)
	except:
		msg = 'Issue getting children leaves from input tree.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def determine_tree_score_bl_sim(ptree, min_bl=0.0005):
	try:
		phylo_breadth = 0.0
		for n in ptree.traverse('postorder'):
			if n.is_root(): continue
			phylo_breadth += min_bl # minimum branch length in normal FastTree 2 after v2.1.7
		return(phylo_breadth)
	except:
		msg = 'Issue determining branch-length score for tree partitioning.'
		sys.stderr.write(msg + '\n')

def recursive_splitting(intree, all_og_samples, are_proteins=False):
	try:
		all_leaves = get_children(intree)
		bl_sum = determine_tree_score_bl(intree)
		if are_proteins:
			if bl_sum == 0.0: 
				return([all_leaves])
		else:
			threshold = determine_tree_score_bl_sim(intree)
			if bl_sum <= threshold:
				return([all_leaves])

		full_score = determine_tree_score(intree, all_og_samples, are_proteins=are_proteins)

		subtree_info = []
		min_score = 1e100
		samples = set([])
		for n in intree.traverse('preorder'):
			if n.is_leaf(): 
				samples.add(n.name.split('|')[0])
			else:
				subtree_score = determine_tree_score(n, all_og_samples, are_proteins=are_proteins)
				subtree_info.append([get_children(n), subtree_score])				
				if subtree_score < min_score:
					min_score = subtree_score

		if min_score >= full_score or len(samples) == 1:
			return([all_leaves])

		accounted_for = set([])
		best_partitions = []
		for i, sti in enumerate(sorted(subtree_info, key=itemgetter(1))):
			if sti[1] == min_score:
				if len(sti[0].intersection(accounted_for)) == 0:
					best_partitions.append(sti[0])
					accounted_for = accounted_for.union(sti[0])

		complement_leaves = all_leaves.difference(accounted_for)
		complement_samples = set([x.split('|')[0] for x in complement_leaves])
		complement_splitting = [complement_leaves]
		if len(complement_samples) >= 2:
			complement_tree = intree.copy()
			complement_tree.prune(list(complement_leaves), preserve_branch_length=True)
			complement_splitting = recursive_splitting(complement_tree, all_og_samples, are_proteins=are_proteins)
			del complement_tree
		result = best_partitions + complement_splitting
		return(result)
	except:
		msg = 'Issues splitting ortholog group tree using the recursive function.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def furtherSplitOutlierArtifactGroups(rooted_t, sps):
	"""
	Further split ortholog groups that may contain small proteins or be the result
	of outlier sequences being claded together after midpoint rooting of FastME tree.
	"""	
	try:
		all_proteins = get_children(rooted_t)

		updated_sp = []
		for sp in sps:
			singletons = set([])
			for p1 in sorted(sp):
				p1dogs = protein_dogs[p1]
				max_jacc_internal = 0.0
				max_jacc_external = 0.0
				for p2 in sorted(sp):
					if p1 == p2: continue
					
					p2dogs = protein_dogs[p2]

					union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))

					union_count = 0
					intersect_count = 0
					for d in union_dogs:
						if d in largely_idr_dogs: continue
						union_count += p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
						intersect_count += min([p1dogs[d], p2dogs[d]])

					if union_count > 0:
						jaccard_index = intersect_count/union_count
						if jaccard_index > max_jacc_internal:
							max_jacc_internal = jaccard_index

				for p2 in sorted(all_proteins):
					if p2 in sp: continue
					p2dogs = protein_dogs[p2]

					union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))

					union_count = 0
					intersect_count = 0
					for d in union_dogs:
						if d in largely_idr_dogs: continue
						union_count += p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
						intersect_count += min([p1dogs[d], p2dogs[d]])

					if union_count > 0:
						jaccard_index = intersect_count/union_count
						if jaccard_index > max_jacc_external:
							max_jacc_external = jaccard_index
	
				if max_jacc_internal <= max_jacc_external:
					singletons.add(p1)
					updated_sp.append(set([p1]))
			
			remaining_sp = set([])
			for p in sp:
				if p in singletons: continue
				remaining_sp.add(p)
			
			if len(remaining_sp) == 0: continue
			updated_sp.append(remaining_sp)

		return(updated_sp)
	except:
		msg = 'Issues refining domain ortholog groups based on phylogenetics.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def furtherSplitOutlierArtifactGroupBurstApproach(rooted_t, sps, dj=0.25, bt=0.5):
	"""
	Burst clades with very disconnected proteins that might just artificially be produced
	by FastME into singletons.
	"""	
	try:
		updated_sp = []
		
		for sp in sps:
			total_comparisons = 0
			meet_threshold = 0
			for i, p1 in enumerate(sorted(sp)):
				for j, p2 in enumerate(sorted(sp)):
					if i >= j: continue
					total_comparisons += 1

					p1dogs = protein_dogs[p1]
					p2dogs = protein_dogs[p2]
					
					intersect_dogs = (set(p1dogs.keys())).intersection(set(p2dogs.keys()))
					union_dogs = (set(p1dogs.keys())).union(set(p2dogs.keys()))
					sc_dogs = single_copy_dogs.intersection(intersect_dogs)

					threshold = dj
					if len(sc_dogs) >= 1:
						for sd in sc_dogs:
							sd_conservation = dog_conservation[sd]
							updated_threshold = dj - (dj*sd_conservation)
							if updated_threshold < threshold:
								threshold = updated_threshold

					union_count = 0
					intersect_count = 0
					for d in union_dogs:
						if d in largely_idr_dogs: continue
						union_count += p1dogs[d] + p2dogs[d] - min([p1dogs[d], p2dogs[d]])
						intersect_count += min([p1dogs[d], p2dogs[d]])

					if union_count > 0:
						jaccard_index = intersect_count/union_count
						if jaccard_index >= threshold:
							meet_threshold += 1

			if total_comparisons > 1 and ((total_comparisons-meet_threshold)/float(total_comparisons)) >= bt:
				for p in sp:
					updated_sp.append(set([p]))
			else:
				updated_sp.append(sp)

		return(updated_sp)
	except:
		msg = 'Issues refining domain ortholog groups based on phylogenetics.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def furtherSplitDisjointDOGPartitions(rooted_t, sps):
	"""
	Function to further split disjoint domain ortholog groups based on the phylogenetic tree.
	"""
	try:
		node_id = 1
		innernode_children = []
		for n in rooted_t.traverse('postorder'):
			if n.is_leaf():
				innernode_children.append([node_id, set([n.name]), 1])
			else:
				children = get_children(n)
				innernode_children.append([node_id, children, len(children)])
			node_id += 1

		updated_sp = []
		for sp in sps:
			monophyletic_flag = False
			for n in innernode_children:
				if sp == n[1]:
					monophyletic_flag = True
			if monophyletic_flag:
				updated_sp.append(sp)
			else:
				accounted_leafs = set([])
				for n in sorted(innernode_children, key=itemgetter(2), reverse=True):
					nc = n[1]
					if nc.issubset(sp) and len(nc.difference(sp)) == 0 and len(nc.intersection(accounted_leafs)) == 0:
						accounted_leafs = accounted_leafs.union(nc)
						updated_sp.append(nc)
				assert(sp == accounted_leafs)
		return updated_sp

	except:
		msg = 'Issues further refining domain ortholog groups based on phylogenetics.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)
		
def assess_job_intensity(faa_file):
	try:
		heavy_job = False
		seq_lens = []
		with open(faa_file) as off:
			for rec in SeqIO.parse(off, 'fasta'):
				seq_lens.append(len(str(rec.seq)))
		med_seq_len = statistics.median(seq_lens)
		seq_count = len(seq_lens)
		if seq_count >= 100 or med_seq_len >= 1500:
			heavy_job = True
		return(heavy_job)
	except:
		msg = 'Issues with assessing intensity of constructing alignment/phylogenies for ortholog group.'
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def mergeBackDOGPartitions(t, pw_dists, split_partitions, fixation_index_cutoff=0.5):	
	all_sps = set([])
	leaf_to_clade = {}
	clade_leaves = defaultdict(set)
	for i, sp in enumerate(split_partitions):
		for l in sp:
			clade_leaves[i].add(l)
			leaf_to_clade[l] = i
			all_sps.add(i)

	if len(clade_leaves) == 1:
		return(split_partitions)

	node_clade_samples = defaultdict(lambda: defaultdict(list))
	node_tot_children = defaultdict(int)
	for ni, n in enumerate(t.traverse('postorder')):
		if n.is_leaf(): continue
		name = 'Node_' + str(ni+1)
		n.name = name
		for l in n.traverse('postorder'):
			if l.is_leaf():
				c = leaf_to_clade[l.name]
				node_clade_samples[name][c].append(l.name)
				node_tot_children[name] += 1			

	merge_sps = []
	paired_sps = set([])
	for n in t.traverse('postorder'):
		if n.is_leaf(): continue
		if len(node_clade_samples[n.name]) == 1: continue
		within_dists = []
		for c in sorted(node_clade_samples[n.name]):
			for i, s1 in enumerate(sorted(node_clade_samples[n.name][c])):
				for j, s2 in enumerate(sorted(node_clade_samples[n.name][c])):
					if i >= j: continue
					key = tuple(sorted([s1, s2]))
					dist = pw_dists[key]
					within_dists.append(dist)

		between_dists = []
		for i, c1 in enumerate(sorted(node_clade_samples[n.name])):
			for j, c2 in enumerate(sorted(node_clade_samples[n.name])):
				if i >= j: continue
				for s1 in node_clade_samples[n.name][c1]:
					for s2 in node_clade_samples[n.name][c2]:
						key = tuple(sorted([s1, s2]))
						dist = pw_dists[key]
						between_dists.append(dist)
		
		if len(within_dists) == 0: continue
		pi_within = statistics.mean(within_dists)
		pi_between = statistics.mean(between_dists)
		fix_index = (pi_between - pi_within)/pi_between

		if fix_index > fixation_index_cutoff: continue
		for i, c1 in enumerate(sorted(node_clade_samples[n.name])):
			for j, c2 in enumerate(sorted(node_clade_samples[n.name])):	
				if i >= j: continue
				merge_sps.append([c1, c2])
				paired_sps.add(c1)
				paired_sps.add(c2)

	if len(merge_sps) > 0:
		merged_sp_ids = single_linkage_cluster(merge_sps, all_sps, paired_sps)
		merged_sp_listing = []
		for msp in merged_sp_ids:
			msl = []
			for sp in msp:
				for l in clade_leaves[sp]:
					msl.append(l)
			merged_sp_listing.append(msl)
		return(merged_sp_listing)
	else:
		return(split_partitions)

def determine_tree_score_bl(ptree):
	try:
		phylo_breadth = 0.0
		for n in ptree.traverse('postorder'):
			if n.is_root(): continue
			phylo_breadth += n.dist
		return(phylo_breadth)
	except:
		msg = 'Issue determining branch-length score for tree partitioning.'
		sys.stderr.write(msg + '\n')

def pairwise_dist(inputs):
	d, l1, l2 = inputs
	dist = tree_obj.get_distance(l1, l2)
	key = tuple(sorted([l1, l2]))
	d[key] = dist

def splitDOGs(inputs):
	og, tre_file, spl_file, skip_merge_back_flag, rooting_seeds, fixation_index_cutoff, threads, logObject = inputs
	rooting_seeds = max([rooting_seeds, 1])
	try:		
		t = Tree(tre_file)
		samples_with_og = set([])
		leafs = set([])
		for n in t.traverse('postorder'):
			if n.is_leaf():
				s = n.name.split('|')[0]
				samples_with_og.add(s)
				leafs.add(n.name)

		pw_dists = None
		if threads > 1:
			with multiprocessing.Manager() as manager:
				d = manager.dict()
				pairs = []
				for i, l1 in enumerate(leafs):
					for j, l2 in enumerate(leafs):
						if i < j:
							pairs.append([d, l1, l2])

				with manager.Pool(threads) as pool:
					pool.map(pairwise_dist, pairs)

				pw_dists = dict(d)

			for l in leafs:
				pw_dists[tuple([l, l])] = 0.0

		else:
			pw_dists = {}
			for i, l1 in enumerate(leafs):
				for j, l2 in enumerate(leafs):
					if i >= j: continue
					dist = t.get_distance(l1, l2)
					key = tuple(sorted([l1, l2]))
					pw_dists[key] = dist

		all_rooting_partitions = []
		possible_nodes_for_rooting = set([])
		for node_id, n in enumerate(t.traverse('preorder')):
			if not n.is_leaf():
				n.name = 'node_' + str(node_id+1)	
			possible_nodes_for_rooting.add(n.name)

		if len(possible_nodes_for_rooting) > (rooting_seeds-1):
			possible_nodes_for_rooting = set(random.sample(list(possible_nodes_for_rooting), rooting_seeds-1))

		if len(possible_nodes_for_rooting) > 0:
			for n in t.traverse('preorder'):
				if not n.name in possible_nodes_for_rooting: continue
				rooted_t = copy.deepcopy(t)
				try:
					if n.name is rooted_t.name:
						pass
					else:
						rooted_t.set_outgroup(n.name)
				except:
					sys.stderr.write(traceback.format_exc() + '\n')
					continue

				sp = recursive_splitting(rooted_t, samples_with_og)
				sp_further_split = furtherSplitDisjointDOGPartitions(rooted_t, sp)
				sp_merge = sp_further_split
				if not skip_merge_back_flag:
					sp_merge = mergeBackDOGPartitions(rooted_t, pw_dists, sp_further_split, fixation_index_cutoff=fixation_index_cutoff)

				pscore_sum = 0 
				pscore_bl_sum = 0.0
				for p in sp_merge:
					ptree = copy.deepcopy(rooted_t)
					ptree.prune(p, preserve_branch_length=True)
					ptree_samples = set([x.split('|')[0] for x in p])
					pscore = determine_tree_score(ptree, ptree_samples)
					pscore_bl = determine_tree_score_bl(ptree)
					del ptree
					pscore_sum += pscore
					pscore_bl_sum += pscore_bl
				all_rooting_partitions.append([sp_merge, pscore_sum, len(sp_merge), pscore_bl_sum])			
				del rooted_t

		mp_t = copy.deepcopy(t)
		R = mp_t.get_midpoint_outgroup()
		mp_t.set_outgroup(R)
		sp = recursive_splitting(mp_t, samples_with_og)
		sp_further_split = furtherSplitDisjointDOGPartitions(mp_t, sp)
		sp_merge = sp_further_split
		if not skip_merge_back_flag:
			sp_merge = mergeBackDOGPartitions(mp_t, pw_dists, sp_further_split, fixation_index_cutoff=fixation_index_cutoff)
			
		pscore_sum = 0
		pscore_bl_sum = 0.0
		for p in sp_merge:
			ptree = copy.deepcopy(mp_t)
			ptree.prune(p, preserve_branch_length=True)
			ptree_samples = set([x.split('|')[0] for x in p])
			pscore = determine_tree_score(ptree, ptree_samples)
			pscore_bl = determine_tree_score_bl(ptree)
			del ptree
			pscore_sum += pscore
			pscore_bl_sum += pscore_bl
		all_rooting_partitions.append([sp_merge, pscore_sum, len(sp_merge), pscore_bl_sum])
		del mp_t

		if len(all_rooting_partitions) > 0:
			spl_outf = open(spl_file, 'w')
			for i, sp in enumerate(sorted(all_rooting_partitions, key=itemgetter(1,2,3))):
				if i == 0:
					for it, spi in enumerate(sp[0]):
						spog = og + '_' + str(it)
						for dom in spi:
							spl_outf.write(spog + '\t' + dom + '\n')
			spl_outf.close()

		return
	except:
		msg = 'Issue with splitting ortholog group %s - based on phylo from full MSA.' % og
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(msg)
		sys.exit(1)
		msg = 'Issue determining median sequence divergence from phylogeny for unsplit DOG %s' % dog
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		logObject.error(msg)

def resolveOrthogroupsUsingPhylogenetics(orthofinder_fasta_dir, orthofinder_tsv_file, orthofinder_tsv_singletons_file, resdog_dir, result_file, logObject, use_super5=True, exhaustive_rooting=False, skip_merge_back_flag=False, fixation_index_cutoff=0.25, rooting_seeds=100, trimal_options='-strict -keepseqs', threads=1):
	try:
		msa_dir = resdog_dir + 'Protein_MSAs/'
		trim_dir = resdog_dir + 'Protein_MSAs_Trimmed/'
		tre_dir = resdog_dir + 'Protein_Trees/'
		spl_full_dir = resdog_dir + 'Split_Protein_Listings/'
		setupReadyDirectory([msa_dir, trim_dir, tre_dir, spl_full_dir])

		muscle_cmds = []
		trimal_cmds = []
		for f in os.listdir(orthofinder_fasta_dir):
			in_faa = orthofinder_fasta_dir + f
			msa_file = msa_dir + '.'.join(f.split('.')[:-1]) + '.msa.faa'
			trim_file = trim_dir + '.'.join(f.split('.')[:-1]) + '.msa.trimmed.faa'
			
			heavy_job = assess_job_intensity(in_faa)

			if heavy_job:
				msa_cmd = ['muscle', '-super5', in_faa, '-output', msa_file, '-threads', str(threads), '-perturb', '12345']
				if not use_super5:
					msa_cmd = ['muscle', '-align', in_faa, '-output', msa_file, '-threads', str(threads), '-perturb', '12345']
				trim_cmd = ['trimal', '-in', msa_file, '-out', trim_file, trimal_options, logObject]
				trimal_cmds.append(trim_cmd)
				runCmd(msa_cmd, logObject)
			else:
				msa_cmd = ['muscle', '-super5', in_faa, '-output', msa_file, '-threads', '1', '-perturb', '12345', logObject]
				if not use_super5:
					msa_cmd = ['muscle', '-align', in_faa, '-output', msa_file, '-threads', '1', '-perturb', '12345', logObject]
				trim_cmd = ['trimal', '-in', msa_file, '-out', trim_file, trimal_options, logObject]
				trimal_cmds.append(trim_cmd)
				muscle_cmds.append(msa_cmd)

		msg = 'Running multiple-sequence alignments using MUSCLE5.'
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, muscle_cmds), total=len(muscle_cmds)):
			pass
		p.close()

		msg = 'Running alignment trimming using trimal.'
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, trimal_cmds), total=len(trimal_cmds)):
			pass
		p.close()

		fasttree_cmds = []
		trim_msa_dogs = set([])
		for f in os.listdir(orthofinder_fasta_dir):
			dog = '.'.join(f.split('.')[:-1])
			msa_file = msa_dir + dog + '.msa.faa'
			trim_file = trim_dir + dog + '.msa.trimmed.faa'
			tre_file = tre_dir + dog + '.tre'
		
			trim_num_sites = 0
			if os.path.isfile(trim_file):
				with open(trim_file) as otf:
					for i, rec in enumerate(SeqIO.parse(otf, 'fasta')):
						if i == 0:
							trim_num_sites = len(str(rec.seq))
			
			if trim_num_sites >= 10:
				tre_cmd = ['fasttree', '-out', tre_file, trim_file, logObject]
				fasttree_cmds.append(tre_cmd)
				trim_msa_dogs.add(dog)
			else:
				tre_cmd = ['fasttree', '-out', tre_file, msa_file, logObject]
				fasttree_cmds.append(tre_cmd)

		msg = 'Running phylogeny constructions using FastTree 2.'
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, fasttree_cmds), total=len(fasttree_cmds)):
			pass
		p.close()

		samples = []
		with open(orthofinder_tsv_file) as otf:
			for i, line in enumerate(otf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0: 
					samples = ['.ccds'.join(x.split('.ccds')[:-1]) for x in ls[1:]]
		
		split_inputs = []
		for f in os.listdir(tre_dir):
			dog = '.tre'.join(f.split('.tre')[:-1])
			tre_file = tre_dir + f
			spl_full_file = spl_full_dir + dog + '.txt'
			t = Tree(tre_file)
			if len(t.get_leaves()) < 500:
				split_inputs.append([dog, tre_file, spl_full_file, skip_merge_back_flag, rooting_seeds, fixation_index_cutoff, 1, logObject])
			else:
				global tree_obj
				tree_obj = Tree(tre_file)
				splitDOGs([dog, tre_file, spl_full_file, skip_merge_back_flag, rooting_seeds, fixation_index_cutoff, threads, logObject])
				tree_obj = None

		msg = 'Using phylogenetics to split coarse domain resolution ortholog groups.'
		sys.stderr.write(msg + '\n')
		logObject.info(msg)
		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(splitDOGs, split_inputs), total=len(split_inputs)):
			pass
		p.close()

		outf_handle = open(result_file, 'w')
		samples = []
		with open(orthofinder_tsv_file) as otf:
			for i, line in enumerate(otf):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0: 
					samples = ['.ccds'.join(x.split('.ccds')[:-1]) for x in ls[1:]]
					outf_handle.write('OG/Sample\t' + '\t'.join(samples) + '\n')
				else:
					dog = ls[0]
					dog_spl_file = spl_full_dir + dog + '.txt'
					if os.path.isfile(dog_spl_file):						
						cluster_prots = defaultdict(set)
						with open(dog_spl_file) as odsf:
							for line in odsf:
								line = line.strip()
								cluster_id, prot = line.split('\t')
								cluster_prots[cluster_id].add(prot)

						for cid in cluster_prots:
							printlist = [cid]
							for lts in ls[1:]:
								lts_retained = []
								for lt in lts.split(','):
									lt = lt.strip()
									if lt in cluster_prots[cid]:
										lts_retained.append(lt)
								printlist.append(', '.join(lts_retained))
							outf_handle.write('\t'.join(printlist) + '\n')
					else:
						outf_handle.write(line + '\n')
						
		with open(orthofinder_tsv_singletons_file) as ootsf:
			for i, line in enumerate(ootsf):
				line = line.strip('\n')
				if i == 0: continue
				outf_handle.write(line + '\n')
		outf_handle.close()	
		
	except:
		msg = 'Issues splitting domain ortholog groups from OrthoFinder using phylogenetic processing.'
		logObject.error(msg)
		sys.stderr.write(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def single_linkage_cluster(pairs, all_lts, paired_lts):
	try:
		"""	
		Solution for single-linkage clustering taken from mimomu's repsonse in the stackoverflow page:
		https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements?lq=1
		"""
		L = pairs
		LL = set(itertools.chain.from_iterable(L))
		for each in LL:
			components = [x for x in L if each in x]
			for i in components:
				L.remove(i)
			L += [list(set(itertools.chain.from_iterable(components)))]

		for lt in all_lts:
			if not lt in paired_lts:
				L.append([lt])

		return (L)
	except:
		msg = 'Issue running single linkage clustering!'
		sys.stderr.wrtie(msg + '\n')
		sys.stderr.write(traceback.format_exc() + '\n')
		sys.exit(1)

def combineOrthoFinderResults(orthofinder_tsv_file, orthofinder_tsv_singletons_file, result_file, logObject):
	try:
		outf_handle = open(result_file, 'w')
		with open(orthofinder_tsv_file) as ootf:
			for i, line in enumerate(ootf):
				line = line.strip('\n')
				if i == 0:
					ls = line.split('\t')
					samples = ['.ccds'.join(x.split('.ccds')[:-1]) for x in ls[1:]]
					outf_handle.write('OG/Sample\t' + '\t'.join(samples) + '\n')
				else:
					outf_handle.write(line + '\n')

		with open(orthofinder_tsv_singletons_file) as ootsf:
			for i, line in enumerate(ootsf):
				line = line.strip('\n')
				if i == 0: continue
				outf_handle.write(line + '\n')
		outf_handle.close()	
	except Exception as e:	
		msg = 'Issues combining ortholog groups and singleton ortholog groups into a single file.'
		logObject.error(msg)
		sys.stderr.write(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)


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
	#print(resource.getrlimit(resource.RLIMIT_AS))


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

def processGenomesUsingProdigal(sample_genomes, prodigal_outdir, logObject,
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

			prodigal_cmd = ['runProdigalAndMakeInputsForBofasa.py', '-i', sample_assembly, '-s', sample, '-gcm', gene_calling_method,
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

	except Exception as e:
		logObject.error(
			"Problem with creating commands for running prodigal via script runProdigalAndMakeProperGenbank.py. Exiting now ...")
		logObject.error(traceback.format_exc())
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
		
def processGenomesAsGenbanks(sample_genomes, gp_dir, logObject, threads=1, locus_tag_length=3,
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
	- threads: The number of threads to use.
	- locus_tag_length: The length of the locus tags to generate.
	- rename_locus_tags: Whether to rename locus tags.
	********************************************************************************************************************
	Returns:
	- sample_genomes_updated: Dictionary mapping sample names to paths of final/processed sample GenBanks.
	********************************************************************************************************************
	"""

	process_cmds = []
	try:
		alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
		possible_locustags = sorted(list(
			set([''.join(list(x)) for x in list(itertools.product(alphabet, repeat=locus_tag_length))])))

		for i, sample in enumerate(sorted(sample_genomes)):
			sample_locus_tag = possible_locustags[i]
			sample_genbank = sample_genomes[sample]
			process_cmd = ['processNCBIGenBankAndCreateInputs.py', '-i', sample_genbank, '-s', sample, '-o', gp_dir]
			if rename_locus_tags:
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
			faa_file = gp_dir + sample + '.faa'
			bed_file = gp_dir + sample + '.coords.bed'
			fna_file = gp_dir + sample + '.fna'
			try:
				assert (os.path.isfile(faa_file) and os.path.isfile(bed_file) and os.path.isfile(fna_file))
				assert (os.path.getsize(faa_file) > 0 and os.path.getsize(bed_file) > 0 and os.path.getsize(fna_file) > 0)
			except AssertionError:
				if os.path.isfile(faa_file):
					os.system('rm -f ' + faa_file)
				if os.path.isfile(bed_file):
					os.system('rm -f ' + bed_file)
				if os.path.isfile(fna_file):
					os.system('rm -f ' + fna_file)
				msg = "Unable to validate successful genbank reformatting/predicted-proteome creation for sample %s\n" % sample
				sys.stderr.write(msg + '\n')
				logObject.warning(msg)
				pass

		sys.stdout.write('Successfully processed %s genomes!\n' % successfully_processed)
	
	except Exception as e:
		logObject.error("Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ...")
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
		output_handle.write('\t'.join(['Protein', 'OG', 'Near scaffold edge?', 'Length', 'Upstream genes', 'Downstream genes', 'Upstream OGs', 'Downstream OGs']) + '\n')
		for scaffold in scaffold_features:
			scaffold_features_sorted = sorted(scaffold_features[scaffold], key=itemgetter(0))
			max_features = len(scaffold_features[scaffold])
			for cds_index, cds in enumerate(scaffold_features_sorted):
				cds_start = cds[0]
				cds_end = cds[1]
				cds_name = cds[2]
				if not cds_name in gene_to_og: continue
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

def determineOrthologGroupContexts(bofasa_prep_dir, og_context_info_file, surround_info_dir, og_tsv_file, logObject, surrounding_bp=10000, threads=1):
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
		isfinder_file = bofasa_prep_dir + 'Sample_IS_Element_Proteins.txt'
		plasmid_file = bofasa_prep_dir + 'Sample_Plasmid_Proteins.txt'
		phage_file = bofasa_prep_dir + 'Sample_Phage_Proteins.txt'

		ise_set = set([])
		plasmid_set = set([])
		phage_set = set([])

		genomad_flag = False
		if os.path.isfile(isfinder_file):
			with open(isfinder_file) as oif:
				for line in oif:
					line = line.strip()
					ls = line.split('\t')
					ise_set.add(ls[1])
		
		if os.path.isfile(plasmid_file):
			with open(plasmid_file) as opf:
				for line in opf:
					line = line.strip()
					ls = line.split('\t')
					plasmid_set.add(ls[1])
			genomad_flag = True
		
		if os.path.isfile(phage_file):
			with open(phage_file) as opf:
				for line in opf:
					line = line.strip()
					ls = line.split('\t')
					phage_set.add(ls[1])
			genomad_flag = True
		
		og_genes = defaultdict(set)
		og_samples = defaultdict(set)
		gene_to_og = {}
		if os.path.isfile(og_tsv_file):
			samples = []
			with open(og_tsv_file) as oot:
				for i, line in enumerate(oot):
					line = line.strip('\n')
					ls = line.split('\t')
					if i == 0:
						samples = ls[1:]
					else:
						og = ls[0]
						for j, gs in enumerate(ls[1:]):
							sample = samples[j]
							for g in gs.split(','):
								g = g.strip()
								if g != '':
									og_genes[og].add(g)
									og_samples[og].add(sample)
									gene_to_og[g] = og

		try:
			assert(len(og_genes) > 0)
		except Exception as e:
			msg = "Difficulties parsing input orthogroup results in the file: %s" % og_tsv_file
			logObject.error(msg)
			logObject.error(traceback.format_exc())
			sys.stderr.write(msg + '\n')
			sys.stderr.write(traceback.format_exc() + '\n')

		listing_file = bofasa_prep_dir + 'Info_on_Input_Genome_Files.txt'
		assert(os.path.isfile(listing_file))
		genome_params = []
		with open(listing_file) as olf:
			for i, line in enumerate(olf):
				if i == 0: continue
				line = line.strip()
				sample, ccds_proteome_file, proteome_file, coords_file, genome_file = line.split('\t')
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
		og_proteins = defaultdict(list)
		for genome_info in genome_params:
			sample, coords_file, output_file, gene_to_og, og_genes, surrounding_bp, logObject = genome_info
			with open(output_file) as oof:
				for i, line in enumerate(oof):
					if i == 0: continue
					line = line.strip('\n')
					protein, og, nse, length, upstream_genes, downstream_genes, upstream_ogs, downstream_ogs = line.split('\t')
					og_gene_lengths[og].append(float(length))
					og_proteins[og].append(protein)
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
		og_context_info_handle.write('\t'.join(['OG', 'Median OG length (bp)', 'Percentage contexts near scaffold edge', 'Number of genomes with OG', 
												'Number of protein in OG', 'Context conservation score', 'Context conservation score - complete contexts', 
												'Context entropy score', 'Context entropy score - complete contexts', 'Number of distinct neighbor OGs', 
												'Number of distinct OGs from complete contexts', 'Avg. number of distinct neighbor OGs', 
												'Avg. number of distinct neighbor OGs from complete contexts', 
												'Percentage instances on plasmid (based on geNomad annotation)', 
												'Percentage instances on phage (based on geNomad annotation)', 
												'Percentage homologous to IS-elements (based on ISfinder database)', 'Instances', 'Contexts']) + '\n')
		for og in sorted(og_contexts):
			median_length = og_gene_lengths[og][0]
			if len(og_gene_lengths) > 1:
				median_length = statistics.median(og_gene_lengths[og])
			num_samples = len(og_samples[og])
			num_contexts = len(og_contexts[og])
			nse_perc = round(100.0*(og_contexts_nses[og]/float(num_contexts)),2)
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
			
			plasmid_count = 0
			phage_count = 0
			is_count = 0
			for p in og_proteins[og]:
				if p in phage_set:
					phage_count += 1
				if p in plasmid_set:
					plasmid_count += 1
				if p in ise_set:
					is_count += 1
				
			plasmid_per = 100.0*(plasmid_count / float(num_contexts))
			phage_per = 100.0*(phage_count / float(num_contexts))
			ise_per = 100.0*(is_count / float(num_contexts))

			if not genomad_flag:
				plasmid_per = 'NA'
				phage_per = 'NA'
					
			og_context_info_handle.write('\t'.join([str(x) for x in [og, round(median_length,2), nse_perc, num_samples, num_contexts, context_var_score, 
										 context_var_score_complete, context_entropy, context_entropy_complete, total_nog, total_nog_complete,
										 avg_nog, avg_nog_complete, plasmid_per, phage_per, ise_per, '; '.join(og_proteins[og]), '; '.join(og_contexts[og])]]) + '\n')
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

def createProteinAlignments(bofasa_prep_dir, resulting_ogs_file, prot_dir, prot_algn_dir, logObject, use_super5=True, threads=1):
	"""
	Description:
	This function creates protein alignments from a directory of protein sequences.

	Function originally developed for zol.
	*******************************************************************************************************************
	Parameters:
	- bofasa_prep_dir: Input directory for bofasa generated by bofasa_prep
	- resulting_ogs_file: Resulting orthogroups file.
	- prot_dir: A directory containing protein sequences.
	- prot_algn_dir: A directory to write protein alignments.
	- logObject: A logging object.
	- use_super5: Whether to use the SUPER5 algorithm for MUSCLE alignment.
	- threads: The number of threads to use for alignment.
	"""
	try:
		prot_to_og = defaultdict(dict)
		samples = []
		with open(resulting_ogs_file) as orof:
			for i, line in enumerate(orof):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					samples = ls[1:]
				og = ls[0]
				for j, lts in enumerate(ls[1:]):
					s = samples[j]
					for lt in lts.split(', '):
						if lt != '':
							prot_to_og[s][lt] = og
		
		proteome_dir = bofasa_prep_dir + 'Genome_Processing/Proteomes/'
		for f in os.listdir(proteome_dir):
			if not f.endswith('.faa'): continue
			s = '.faa'.join(f.split('.faa')[:-1])
			with open(proteome_dir + f) as opf:
				for rec in SeqIO.parse(opf, 'fasta'):
					lt = rec.id
					if lt in prot_to_og[s]:
						og = prot_to_og[s][lt]
						outf = prot_dir + og + '.faa'
						outfh = open(outf, 'a+')
						outfh.write('>' + s + '|' + lt + '\n' + str(rec.seq) + '\n')
						outfh.close()

		msa_cmds = []
		for pf in os.listdir(prot_dir):
			prefix = '.faa'.join(pf.split('.faa')[:-1])
			prot_file = prot_dir + pf
			prot_algn_file = prot_algn_dir + prefix + '.msa.faa'
			heavy_job = assess_job_intensity(prot_file)

			if heavy_job:
				msa_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-threads', str(threads), '-perturb', '12345']
				if not use_super5:
					msa_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-threads', str(threads), '-perturb', '12345']
				runCmd(msa_cmd, logObject)
			else:
				msa_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-threads', '1', '-perturb', '12345', logObject]
				if not use_super5:
					msa_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-threads', '1', '-perturb', '12345', logObject]
				msa_cmds.append(msa_cmd)

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, msa_cmds), total=len(msa_cmds)):
			pass
		p.close()
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
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, hmmbuild_cmds), total=len(hmmbuild_cmds)):
			pass
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, hmmemit_cmds), total=len(hmmemit_cmds)):
			pass
		p.map(multiProcess, hmmemit_cmds)
		p.close()
	except Exception as e:
		sys.stderr.write('Issues with creating profile HMMs and consensus sequences.\n')
		logObject.error('Issues with creating profile HMMs and consensus sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def concatenateConsensusAlignment(og_cons_dir, concatenated_consensus_seqs_file, logObject):
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
		concatenated_consensus_seqs_handle.close()

	except Exception as e:
		sys.stderr.write('Issues with concatenating consensus sequences of multi-sequence ortholog groups and singleton protein sequences.\n')
		logObject.error('Issues with concatenating consensus sequences of multi-sequence ortholog groups and singleton protein sequences.')
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createNearSCCResolvedDomainProteinAlignments(bofasa_prep_dir, resulting_dogs_file, rdog_seqs_dir, 
												 rdog_algn_dir, rdog_trim_dir, merged_core_genome_file, 
												 logObject, use_super5=True, near_scc_prop=0.80, threads=1, 
												 trimal_options='-strict -keepseqs', allow_mge=False):
	try:
		isfinder_file = bofasa_prep_dir + 'Sample_IS_Element_Proteins.txt'
		plasmid_file = bofasa_prep_dir + 'Sample_Plasmid_Proteins.txt'
		phage_file = bofasa_prep_dir + 'Sample_Phage_Proteins.txt'

		mge_set = set([])
		
		if os.path.isfile(isfinder_file):
			with open(isfinder_file) as oif:
				for line in oif:
					line = line.strip()
					ls = line.split('\t')
					mge_set.add(ls[1])
		
		if os.path.isfile(plasmid_file):
			with open(plasmid_file) as opf:
				for line in opf:
					line = line.strip()
					ls = line.split('\t')
					mge_set.add(ls[1])

		if os.path.isfile(phage_file):
			with open(phage_file) as opf:
				for line in opf:
					line = line.strip()
					ls = line.split('\t')
					mge_set.add(ls[1])

		prot_to_og = {}
		samples = []
		with open(resulting_dogs_file) as orof:
			for i, line in enumerate(orof):
				line = line.strip('\n')
				ls = line.split('\t')
				if i == 0:
					samples = ls[1:]
					continue
				og = ls[0]
				# TODO: correct the following prot needs to be checked not og 
				og_has_mge_prot = False
				for lts in ls[1:]:
					for lt in lts.split(','):
						lt = lt.strip()
						if lt == '': continue
						if lt.split('|')[1] in mge_set:
							og_has_mge_prot = True

				if not allow_mge and og_has_mge_prot: continue

				samples_with_sc = set([])
				for j, lts in enumerate(ls[1:]):
					s = samples[j]		
					if ',' in lts: continue
					for lt in lts.split(', '):	
						lt = lt.strip()					
						if lt != '':
							samples_with_sc.add(lt)

				if float(len(samples_with_sc)/len(samples)) >= near_scc_prop:
					for lt in samples_with_sc:
						prot_to_og[lt] = og

		proteome_dir = bofasa_prep_dir + 'Domain_and_Interdomain_FASTAs/'
		for f in os.listdir(proteome_dir):
			if not f.endswith('.faa'): continue
			s = '.faa'.join(f.split('.faa')[:-1])
			with open(proteome_dir + f) as opf:
				for rec in SeqIO.parse(opf, 'fasta'):
					lt = rec.id
					if lt in prot_to_og:
						og = prot_to_og[lt]
						outf = rdog_seqs_dir + og + '.faa'
						outfh = open(outf, 'a+')
						outfh.write('>' + lt + '\n' + str(rec.seq) + '\n')
						outfh.close()

		msa_trim_cmds = []
		for pf in os.listdir(rdog_seqs_dir):
			prefix = '.faa'.join(pf.split('.faa')[:-1])
			prot_file = rdog_seqs_dir + pf
			prot_algn_file = rdog_algn_dir + prefix + '.msa.faa'
			prot_algn_trim_file = rdog_trim_dir + prefix + '.trimmed.msa.faa'
			heavy_job = assess_job_intensity(prot_file)

			if heavy_job:
				msa_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-threads', str(threads), '-perturb', '12345']
				if not use_super5:
					msa_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-threads', str(threads), '-perturb', '12345']
				runCmd(msa_cmd, logObject)
				trimal_cmd = ['trimal', '-in', prot_algn_file, '-out', prot_algn_trim_file, trimal_options, logObject]
				msa_trim_cmds.append(trimal_cmd)
			else:
				msa_cmd = ['muscle', '-super5', prot_file, '-output', prot_algn_file, '-threads', '1', '-perturb', '12345']
				if not use_super5:
					msa_cmd = ['muscle', '-align', prot_file, '-output', prot_algn_file, '-threads', '1', '-perturb', '12345']
				trimal_cmd = ['trimal', '-in', prot_algn_file, '-out', prot_algn_trim_file, trimal_options, logObject]
				msa_trim_cmd = msa_cmd + [';'] + trimal_cmd 
				msa_trim_cmds.append(msa_trim_cmd)

		p = multiprocessing.Pool(threads)
		for _ in tqdm.tqdm(p.imap_unordered(multiProcess, msa_trim_cmds), total=len(msa_trim_cmds)):
			pass
		p.close()

		sample_seqs = defaultdict(lambda: "")
		for f in os.listdir(rdog_trim_dir):
			samples_accounted = set([])
			seqlen = 0
			with open(rdog_trim_dir + f) as ortf:
				for rec in SeqIO.parse(ortf, 'fasta'):
					seq = str(rec.seq)
					seqlen = len(seq)
					if seqlen == 0: continue 
					s = rec.id.split('|')[0]
					samples_accounted.add(s)
					sample_seqs[s] += seq

			for s in samples:
				if not s in samples_accounted:
					sample_seqs[s] += ('-'*seqlen)

		aln_handle = open(merged_core_genome_file, 'w')
		for s in sample_seqs:
			aln_handle.write('>' + s + '\n' + sample_seqs[s] + '\n')
		aln_handle.close()
	except Exception as e:
		msg = 'Issues with create core genome alignment(s).'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createFinalReport(bofasa_prep_dir, og_context_info_file, final_result_file, logObject):
	"""
	"""
	try:
		plasmid_file = bofasa_prep_dir + 'Sample_Plasmid_Proteins.txt'
		phage_file = bofasa_prep_dir + 'Sample_Phage_Proteins.txt'

		genomad_flag = False
		if os.path.isfile(plasmid_file) or os.path.isfile(phage_file):
			genomad_flag = True

		# Generate Excel spreadsheet
		writer = pd.ExcelWriter(final_result_file, engine='xlsxwriter')
		workbook = writer.book
		dd_sheet = workbook.add_worksheet('Data Dictionary')
		dd_sheet.write(0, 0, 'Data Dictionary describing columns of "bofasa Results" spreadsheet can be found below and on bofasa\'s Wiki page at:')
		dd_sheet.write(1, 0, 'https://github.com/raufs/bofasa')


		wrap_format = workbook.add_format({'text_wrap': True, 'valign': 'vcenter', 'align': 'center', 'border': 1})
		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#FFFFFF', 'border': 1})

		numeric_columns = set(['Median OG length (bp)', 'Percentage contexts near scaffold edge', 'Number of genomes with OG', 
							   'Number of protein in OG', 'Context conservation score', 'Context conservation score - complete contexts', 
							   'Context entropy score', 'Context entropy score - complete contexts', 'Number of distinct neighbor OGs', 
							   'Number of distinct OGs from complete contexts', 'Avg. number of distinct neighbor OGs', 
							   'Avg. number of distinct neighbor OGs from complete contexts', 
							   'Percentage instances on plasmid (based on geNomad annotation)', 
							   'Percentage instances on phage (based on geNomad annotation)', 
							   'Percentage homologous to IS-elements (based on ISfinder database)'])

		warn_format = workbook.add_format({'bg_color': '#bf241f', 'bold': True, 'font_color': '#FFFFFF'})
		na_format = workbook.add_format({'font_color': '#a6a6a6', 'bg_color': '#FFFFFF', 'italic': True})
		header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#D7E4BC', 'border': 1})

		results_df = loadTableInPandaDataFrame(og_context_info_file, numeric_columns, cut_last_columns=2)
		num_rows = results_df.shape[0]+1
		results_df.to_excel(writer, sheet_name='bofasa Results', index=False, na_rep="NA")
		
		worksheet =  writer.sheets['bofasa Results']
		worksheet.conditional_format('A2:BA' + str(num_rows), {'type': 'cell', 'criteria': '==', 'value': '"NA"', 'format': na_format})
		worksheet.conditional_format('A1:BA1', {'type': 'cell', 'criteria': '!=', 'value': 'NA', 'format': header_format})

		max_values = defaultdict(lambda: 0.0)
		with open(og_context_info_file) as ocif:
			for i, line in enumerate(ocif):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				num_genomes = float(ls[3])
				num_proteins = float(ls[4])

				if num_genomes > max_values['num_genomes']:
					max_values['num_genomes'] = num_genomes
				if num_proteins > max_values['num_proteins']:
					max_values['num_proteins'] = num_proteins

				if ls[5] != 'NA':
					context_var_score = float(ls[5])
					if context_var_score > max_values['context_var_score']:
						max_values['context_var_score'] = context_var_score
				if ls[6] != 'NA':
					context_var_score_comp = float(ls[6])
					if context_var_score_comp > max_values['context_var_score_comp']:
						max_values['context_var_score_comp'] = context_var_score_comp
				if ls[7] != 'NA':
					context_ent_score = float(ls[7])
					if context_ent_score > max_values['context_ent_score']:
						max_values['context_ent_score_comp'] = context_ent_score
				if ls[8] != 'NA':
					context_ent_score_comp = float(ls[8])
					if context_ent_score_comp > max_values['context_ent_score_comp']:
						max_values['context_ent_score_comp'] = context_ent_score_comp

		# median OG length
		worksheet.conditional_format('B2:B' + str(num_rows), {'type': '2_color_scale', 'min_color': "#a9cafc", 'max_color': "#736991", "min_value": 0, "max_value": 2500, 'min_type': 'num', 'max_type': 'num'})

		# percentage instances near scaffold edge
		worksheet.conditional_format('C2:C' + str(num_rows), {'type': '2_color_scale', 'min_color': "#ffffff", 'max_color': "#ed9393", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

		# num genomes with OG
		worksheet.conditional_format('D2:D' + str(num_rows), {'type': '2_color_scale', 'min_color': "#f2c6f7", 'max_color': "#b07fb5", "min_value": 0.0, "max_value": max_values['num_genomes'], 'min_type': 'num', 'max_type': 'num'})

		# num proteins with OG
		worksheet.conditional_format('E2:E' + str(num_rows), {'type': '2_color_scale', 'min_color': "#e7cdf7", 'max_color': "#a186b3", "min_value": 0.0, "max_value": max_values['num_proteins'], 'min_type': 'num', 'max_type': 'num'})

		# context variability score
		worksheet.conditional_format('F2:F' + str(num_rows), {'type': '2_color_scale', 'min_color': "#e6f5ab", 'max_color': "#a4b36b", "min_value": 0.0, "max_value": max_values['context_var_score'], 'min_type': 'num', 'max_type': 'num'})
		
		# context variability score - complete
		worksheet.conditional_format('G2:G' + str(num_rows), {'type': '2_color_scale', 'min_color': "#e6f5ab", 'max_color': "#a4b36b", "min_value": 0.0, "max_value": max_values['context_var_score_comp'], 'min_type': 'num', 'max_type': 'num'})
	
		# context variability score
		worksheet.conditional_format('H2:H' + str(num_rows), {'type': '2_color_scale', 'min_color': "#b3e3d6", 'max_color': "#6aa192", "min_value": 0.0, "max_value": max_values['context_ent_score'], 'min_type': 'num', 'max_type': 'num'})
		
		# context variability score - complete
		worksheet.conditional_format('I2:I' + str(num_rows), {'type': '2_color_scale', 'min_color': "#b3e3d6", 'max_color': "#6aa192", "min_value": 0.0, "max_value": max_values['context_ent_score_comp'], 'min_type': 'num', 'max_type': 'num'})

		if genomad_flag:
			# percentage on plasmid
			worksheet.conditional_format('N2:N' + str(num_rows), {'type': '2_color_scale', 'min_color': "#ffffff", 'max_color': "#ed9393", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})
			
			# percentage on phage
			worksheet.conditional_format('O2:O' + str(num_rows), {'type': '2_color_scale', 'min_color': "#ffffff", 'max_color': "#ed9393", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

		# percentage homologous to IS-elements
		worksheet.conditional_format('P2:P' + str(num_rows), {'type': '2_color_scale', 'min_color': "#ffffff", 'max_color': "#ed9393", "min_value": 0.0, "max_value": 1.0, 'min_type': 'num', 'max_type': 'num'})

		worksheet.autofilter('A1:BA' + str(num_rows))
		workbook.close()
	except Exception as e:
		msg = 'Issues with create final XLSX report.'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def createFinalVisual(bofasa_prep_dir, tmp_result_file, og_context_info_file, final_result_plot, logObject):
	"""
	"""
	try:

		isfinder_file = bofasa_prep_dir + 'Sample_IS_Element_Proteins.txt'
		plasmid_file = bofasa_prep_dir + 'Sample_Plasmid_Proteins.txt'
		phage_file = bofasa_prep_dir + 'Sample_Phage_Proteins.txt'

		mge_set = set([])

		if os.path.isfile(isfinder_file):
			with open(isfinder_file) as oif:
				for line in oif:
					line = line.strip()
					ls = line.split('\t')
					mge_set.add(ls[1])
		
		if os.path.isfile(plasmid_file):
			with open(plasmid_file) as opf:
				for line in opf:
					line = line.strip()
					ls = line.split('\t')
					mge_set.add(ls[1])
		
		if os.path.isfile(phage_file):
			with open(phage_file) as opf:
				for line in opf:
					line = line.strip()
					ls = line.split('\t')
					mge_set.add(ls[1])
			
		tmp_file_header = ['OG', 'Number of protein in OG', 'Context entropy score', 'Over 50%% of protein instances homologous to IS-element or on plasmid or phage']
		outf_handle = open(tmp_result_file, 'w')
		outf_handle.write('\t'.join(tmp_file_header) + '\n')
		with open(og_context_info_file) as oocif:
			for i, line in enumerate(oocif):
				if i == 0: continue
				line = line.strip()
				ls = line.split('\t')
				mge_related = 'No'
				tot = 0
				mge = 0
				for p in ls[-2].split('; '):
					tot += 1
					if p in mge_set:
						mge += 1
				if mge/tot >= 0.5:
					mge_related = 'Yes'
				outf_handle.write('\t'.join([ls[0], ls[4], ls[7], mge_related]) + '\n')
		outf_handle.close()
				
		numeric_columns = set(['Number of protein in OG', 'Context entropy score'])
		simple_df = loadTableInPandaDataFrame(tmp_result_file, numeric_columns)
		fig = px.scatter(simple_df, x="Number of protein in OG", y="Context entropy score",
						color="Over 50%% of protein instances homologous to IS-element or on plasmid or phage", 
						marginal_x="histogram", marginal_y="histogram")
		fig.write_html(final_result_plot)
	except Exception as e:
		msg = 'Issues with create final HTML report.'
		sys.stderr.write(msg + '\n')
		logObject.error(msg)
		sys.stderr.write(str(e) + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)

def loadTableInPandaDataFrame(input_file, numeric_columns, cut_last_columns=None):
	"""
	Description:
	This function formats reads a TSV file and stores it as a pandas dataframe. Note, last two columns are skiped.
	********************************************************************************************************************
	Parameters:
	- input_file: The input TSV file, with first row corresponding to the header.
	- numeric_columns: Set of column names which should have numeric data.
	********************************************************************************************************************
	Returns:
	- panda_df: A pandas DataFrame object reprsentation of the input TSV file.
	********************************************************************************************************************
	"""
	panda_df = None
	try:
		data = []
		with open(input_file) as oif:
			for line in oif:
				line = line.strip('\n')
				ls = line.split('\t')
				if cut_last_columns != None:
					ls = ls[:(0-cut_last_columns)]
				data.append(ls)

		panda_dict = {}
		for ls in zip(*data):
			key = ls[0]
			cast_vals = ls[1:]
			if key in numeric_columns:
				cast_vals = []
				for val in ls[1:]:
					cast_vals.append(castToNumeric(val))
			panda_dict[key] = cast_vals
		panda_df = pd.DataFrame(panda_dict)

	except Exception as e:
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)
	return panda_df

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
					start = int(start)
					end = int(end)
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
				outf.write(sample + '\t' + p + '\t' + ', '.join(best_hits_by_bitscore[p][0]) + str(statistics.mean(best_hits_by_bitscore[p][1])) + '\n')
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
		- dom_coord_file: The domain coordinates file.
		- pfam_db_file: The Pfam HMM DB file.
		- pfam_z: The Pfam record count - for accurate E-value estimation.
		- minimal_length: The minimum length in amino acids for a domain matching or intra-domain region to be kept and
						  tagged as a chopped CDS feature [Default is 20].
		- logObject: A logging object.
		- threads: The number of threads to use [Default is 1].
	********************************************************************************************************************
	"""
	prot_file, ccds_prot_file, dom_coord_file, pfam_db_file, pfam_z, minimal_length, logObject, threads, skip_domain_splitting = inputs
	try:
		sample = '.'.join(prot_file.split('/')[-1].split('.')[:-1])
		if not skip_domain_splitting:
			# align Pfam domains and remove overlap similar to BiG-SCAPE
			alphabet = pyhmmer.easel.Alphabet.amino()
			sequences = []
			with pyhmmer.easel.SequenceFile(prot_file, digital=True, alphabet=alphabet) as seq_file:
				sequences = list(seq_file)

			target_dom_hits = defaultdict(list)
			with pyhmmer.plan7.HMMFile(pfam_db_file) as hmm_file:
				for hits in pyhmmer.hmmsearch(hmm_file, sequences, bit_cutoffs="trusted", Z=int(pfam_z), cpus=threads):
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
					if len(overlap_coords)/float(end-start+1) >= 0.1 or len(overlap_coords) >= minimal_length: continue
					accounted_coords = accounted_coords.union(set(range(start, end+1)))
					breakpoints[tg].append(start)
					breakpoints[tg].append(end+1)
					dom_start_names[tg + '|' + str(start)] = tg + '|' + dom_name + '|' + str(tg_dom_name_iter[dom_name]+1) 
					tg_dom_name_iter[dom_name] += 1

			cpf_handle = open(ccds_prot_file, 'w')
			dcf_handle = open(dom_coord_file, 'w')
			dcf_handle.write('Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n')
			with open(prot_file) as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					tg = rec.id
					tg_seq = str(rec.seq)
					prev_end_coord = 1
					tg_interdomain_index = 1
					if not tg in breakpoints and len(tg_seq) >= minimal_length:
						dn = sample + '|' + tg + '|full_protein|1'
						cpf_handle.write('>' + dn + '\n' + str(tg_seq) + '\n')
						dcf_handle.write('\t'.join([sample, tg, dn.split('|')[2], dn.split('|')[3], str(prev_end_coord), str(len(tg_seq))]) + '\n')
					else:
						for tg_seq_chunk in split_by_idx(tg_seq, ([0] + sorted(breakpoints[tg]))):
							if tg_seq_chunk.strip() == '': continue
							end_coord = prev_end_coord + len(tg_seq_chunk) - 1
							if len(tg_seq_chunk) >= minimal_length:
								dn = sample + '|' + dom_start_names[tg + '|' + str(prev_end_coord-1)]
								if dom_start_names[tg + '|' + str(prev_end_coord-1)] == 'NA':
									dn = sample + '|' + tg + '|inter-domain_region|' + str(tg_interdomain_index)
									tg_interdomain_index += 1
								cpf_handle.write('>' + dn + '\n' + str(tg_seq_chunk) + '\n')
								dcf_handle.write('\t'.join([sample, tg, dn.split('|')[2], dn.split('|')[3], str(prev_end_coord), str(end_coord)]) + '\n')
							prev_end_coord = end_coord + 1
			cpf_handle.close()
			dcf_handle.close()
		else:
			cpf_handle = open(ccds_prot_file, 'w')
			dcf_handle = open(dom_coord_file, 'w')
			dcf_handle.write('Sample\tProtein\tAnnotation\tAnnotation_Iterator\tStart\tEnd\n')
			with open(prot_file) as ocf:
				for rec in SeqIO.parse(ocf, 'fasta'):
					tg = rec.id
					tg_seq = str(rec.seq)
					tg_interdomain_index = 1
					if len(tg_seq) >= minimal_length:
						dn = sample + '|' + tg + '|full_protein|1'
						cpf_handle.write('>' + dn + '\n' + str(tg_seq) + '\n')
						dcf_handle.write('\t'.join([sample, tg, dn.split('|')[2], dn.split('|')[3], '1', str(len(tg_seq))]) + '\n')
			cpf_handle.close()
			dcf_handle.close()
	except:
		msg = 'An issue occurred with creating chopped up version of proteome file %s.' % prot_file
		logObject.error(msg)
		sys.stderr.write(traceback.format_exc())
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.exit(1)

def annotateAndSplitProteinsUsingPfam(sample_proteomes, split_proteins_dir, domain_coords_dir, domain_coord_info_file, logObject, minimal_length=20, threads=1, skip_domain_splitting=False):
	"""	
	Description:
	Create chopped CDS GenBank files from regular GenBank input with CDS features.
	********************************************************************************************************************
	Parameters:
	- sample_proteomes: Dictionary mapping sample names (keys) to proteome file paths (values).
	- split_proteins_dir: Directory where to write resulting cCDS proteomes.
	- domain_coords_dir: Directory where to write resulting domain coordinates.
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
			domain_coord_file = domain_coords_dir + sample + '.domain_coords.txt'
			prot_file = sample_proteomes[sample]
			prot_mod_inputs.append([prot_file, ccds_prot_file, domain_coord_file, pfam_hmm_path, pfam_z, minimal_length, logObject, threads, skip_domain_splitting])

		msg = "Creating domain-chopped up version of GenBank files for %d gene clusters" % len(prot_mod_inputs) 
		logObject.info(msg)
		sys.stdout.write(msg + '\n')

		p = multiprocessing.Pool(1)
		for _ in tqdm.tqdm(p.imap_unordered(createChoppedProteomes, prot_mod_inputs), total=len(prot_mod_inputs)):
			pass
		p.close()

		dci_handle = open(domain_coord_info_file, 'w')
		sample_ccds_proteomes = {}
		for f in os.listdir(domain_coords_dir):
			sample = f.split('.domain_coords.txt')[0]
			ccds_prot_file = split_proteins_dir + sample + '.ccds.faa'
			dom_coord_file = domain_coords_dir + f
			sample_ccds_proteomes[sample] = ccds_prot_file
			with open(dom_coord_file) as odf:
				for i, line in enumerate(odf):
					if i == 0: continue
					dci_handle.write(line)
		dci_handle.close()
		return(sample_ccds_proteomes)
	except:
		msg = 'An issue occurred with creating chopped up proteomes.' 
		logObject.error(msg)
		logObject.error(traceback.format_exc())
		sys.stderr.write(msg + '\n')
		sys.stderr.write(traceback.format_exc())
		sys.exit(1)