from collections import defaultdict, namedtuple, Counter
from itertools import combinations
import numpy as np
import scipy.stats
import sys
from os import listdir
import json
import argparse
import random

parser = argparse.ArgumentParser(description='Calculate inheritance pvalues for sibpairs.')
parser.add_argument('phase_dir', type=str, help='Directory with phase data.')
parser.add_argument('ped_file', type=str, help='Ped file for data.')
parser.add_argument('identicals_file', type=str, help='File containing identical twins.')
parser.add_argument('chrom', type=str, help='Chromosome.')

args = parser.parse_args()

def pull_phenotype_ped(ped_file):
	sample_to_sex = dict()
	sample_to_affected = dict()
	parents_to_children = defaultdict(list)
	with open(ped_file, 'r') as f:
		for line in f:
			pieces = line.strip().split('\t')
			sample_to_sex[pieces[1]] = pieces[4]
			sample_to_affected[pieces[1]] = pieces[5]
			if pieces[2] != '0' and pieces[3] != '0':
				parents_to_children[(pieces[3], pieces[2])].append(pieces[1])
	return sample_to_affected, sample_to_sex, parents_to_children

Sibpair = namedtuple('Sibpair', ['family', 'sibling1', 'sibling2', 'mom', 'dad', 'phase_dir', 'num_affected', 'num_males'])
def pull_sibpairs(phase_dir, sample_to_affected, sample_to_sex, parents_to_children, identicals_file):

	# pull identicals
	leave_out = set()
	with open(identicals_file, 'r') as f:
		for line in f:
			pieces = line.strip().split('\t')
			sibling1, sibling2 = pieces[:2]
			leave_out.add((sibling1, sibling2))
			leave_out.add((sibling2, sibling1))


	# pull sibpairs with phase data
	sibpair_to_famkey = dict()
	child_has_phase_data = set()
	family_to_inds = dict()
	for filename in listdir(phase_dir):
		if filename.endswith('.phased.txt'):
			family_key = filename[:-11]
			try:
				with open('%s/%s' % (phase_dir, filename), 'r')  as f:
					header = next(f).strip().split('\t')
					individuals = [header[i][:-4] for i in range(5, len(header)-3, 2)]
					family_to_inds[family_key] = individuals
					for sibling1, sibling2 in combinations(individuals[2:], 2):
						if (sibling1, sibling2) not in leave_out:
							sibpair_to_famkey[(sibling1, sibling2)] = family_key
							sibpair_to_famkey[(sibling2, sibling1)] = family_key
							child_has_phase_data.add(sibling1)
							child_has_phase_data.add(sibling2)
			except StopIteration:
				pass



	def form_sibpair(sibling1, sibling2, mom, dad):
		return Sibpair(sibpair_to_famkey[(sibling1, sibling2)], sibling1, sibling2, mom, dad,
			phase_dir,
			int(sample_to_affected[sibling1]=='2')+int(sample_to_affected[sibling2]=='2'),
			int(sample_to_sex[sibling1]=='1')+int(sample_to_sex[sibling2]=='1'))

	# pull sibpairs from families
	sibpairs = []
	for (mom, dad), children in parents_to_children.items():
		for sibling1, sibling2 in combinations([x for x in children if x in sample_to_affected], 2):
			if (sibling1, sibling2) in sibpair_to_famkey:
				sibpairs.append(form_sibpair(sibling1, sibling2, mom, dad))

		#aff_children = [x for x in children if x in sample_to_affected and sample_to_affected[x]=='2' and x in child_has_phase_data]
		#typ_children = [x for x in children if x in sample_to_affected and sample_to_affected[x]=='1' and x in child_has_phase_data]

		## may be no need
		## we need to keep sibpairs independent so we choose an affected proband and a typical proband
		## and we form sibpairs with the other siblings
		#
		#if len(aff_children)>0:
		#	aff_proband = random.choice(aff_children)
		#else:
		#	aff_proband = None
		#if len(typ_children)>0:
		#	typ_proband = random.choice(typ_children)
		#else:
		#	typ_proband = None
		#			
		#for sibling in [x for x in aff_children+typ_children if x != aff_proband and x != typ_proband]:
		#	if aff_proband is not None and (aff_proband, sibling) in sibpair_to_famkey:
		#		sibpairs.append(form_sibpair(aff_proband, sibling, mom, dad))
		#	if typ_proband is not None and (typ_proband, sibling) in sibpair_to_famkey:
		#		sibpairs.append(form_sibpair(typ_proband, sibling, mom, dad))
		#
		#
		## form sibpair between probands
		#if aff_proband is not None and typ_proband is not None and (aff_proband, typ_proband) in sibpair_to_famkey:
		#	sibpairs.append(form_sibpair(aff_proband, typ_proband, mom, dad))
	sibpairs = sorted(sibpairs)

	assert len(sibpairs) == len(set(sibpairs)) # should have no duplicates
	return sibpairs, family_to_inds

def pull_intervals(chrom, sibpairs, family_to_inds, mat_pat):
	positions = set([0])

	for sibpair_index, sibpair in enumerate(sibpairs):
		inds = family_to_inds[sibpair.family]
		sib1_ind_index, sib2_ind_index = inds.index(sibpair.sibling1), inds.index(sibpair.sibling2)
		sib1_mat_index, sib2_mat_index = 4+2*sib1_ind_index, 4+2*sib2_ind_index
		sib1_pat_index, sib2_pat_index = 5+2*sib1_ind_index, 5+2*sib2_ind_index

		if mat_pat == 'both':
			sib_phase_indices = [sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index]
		elif mat_pat == 'mat':
			sib_phase_indices = [sib1_mat_index, sib2_mat_index]
		elif mat_pat == 'pat':
			sib_phase_indices = [sib1_pat_index, sib2_pat_index]

		with open('%s/%s.phased.txt' % (sibpair.phase_dir, sibpair.family), 'r')  as f:
			next(f) # skip header

			prev_state = None
			for line in f:
				if line.startswith('chr%s\t'%chrom):
					pieces = line.strip().split('\t')
					
					start_pos, end_pos = [int(x) for x in pieces[-2:]]
					state = np.array([int(x) for x in pieces[1:-2]])

					if prev_state is None or not np.all(state[sib_phase_indices]==prev_state[sib_phase_indices]):
						prev_state = state
						positions.add(start_pos)
			positions.add(end_pos)

	positions = np.array(sorted(positions))

	# now remove positions that are too close to both neighbors
	#include = np.ones(positions.shape, dtype=bool)
	#too_close = positions[1:]-positions[:-1]<100000
	#include[np.where(too_close[:-1] & too_close[1:])[0]+1] = False
	#positions = positions[include]

	return positions

def pull_sibpair_matches(chrom, sibpairs, family_to_inds, interval_bins, mat_pat):
	sibpair_to_index = dict([((x.family, x.sibling1, x.sibling2), i) for i, x in enumerate(sibpairs)])
	interval_start_to_index = dict([(x, i) for i, x in enumerate(interval_bins)])

	# pull phase data
	# sibpair, interval, nomatch/match
	is_mat_match = -np.ones((len(sibpair_to_index), len(interval_bins)-1), dtype=int)
	is_pat_match = -np.ones((len(sibpair_to_index), len(interval_bins)-1), dtype=int)

	interval_lengths = interval_bins[1:]-interval_bins[:-1]

	for sibpair_index, sibpair in enumerate(sibpairs):
		flip = random.choice([0, 1])
		inds = family_to_inds[sibpair.family]
		sib1_ind_index, sib2_ind_index = inds.index(sibpair.sibling1), inds.index(sibpair.sibling2)
		sib1_mat_index, sib2_mat_index = 4+2*sib1_ind_index, 4+2*sib2_ind_index
		sib1_pat_index, sib2_pat_index = 5+2*sib1_ind_index, 5+2*sib2_ind_index
		
		if mat_pat == 'both':
			sib_phase_indices = [sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index]
		elif mat_pat == 'mat':
			sib_phase_indices = [sib1_mat_index, sib2_mat_index]
		elif mat_pat == 'pat':
			sib_phase_indices = [sib1_pat_index, sib2_pat_index]

		with open('%s/%s.phased.txt' % (sibpair.phase_dir, sibpair.family), 'r')  as f:
			next(f) # skip header

			current_state, current_start_pos, current_start_index = None, None, None
			for line in f:
				if line.startswith('chr%s\t'%chrom):
					pieces = line.strip().split('\t')
					
					start_pos, end_pos = [int(x) for x in pieces[-2:]]
					state = np.array([int(x) for x in pieces[1:-2]])

					if current_state is None:
						current_state, current_start_pos, current_start_index = state, start_pos, interval_start_to_index[start_pos]
					elif not np.all(state[sib_phase_indices]==current_state[sib_phase_indices]):
						# save data on current state
						end_index = interval_start_to_index[start_pos]
						assert np.sum(interval_lengths[current_start_index:end_index]) == start_pos - current_start_pos
						assert end_index > current_start_index
						
						if (mat_pat == 'both' or mat_pat == 'mat') and (current_state[sib1_mat_index] != -1) and (current_state[sib2_mat_index] != -1):
							is_mat_match[sibpair_index, current_start_index:end_index] = int(int(current_state[sib1_mat_index]==current_state[sib2_mat_index])==flip)
						if (mat_pat == 'both' or mat_pat == 'pat') and (current_state[sib1_pat_index] != -1) and (current_state[sib2_pat_index] != -1):
							is_pat_match[sibpair_index, current_start_index:end_index] = int(int(current_state[sib1_pat_index]==current_state[sib2_pat_index])==flip)

						# move on to next state
						current_state, current_start_pos, current_start_index = state, start_pos, end_index
			
			# last interval
			if current_state is not None:
				assert np.sum(interval_lengths[current_start_index:]) == end_pos - current_start_pos
							
				if (mat_pat == 'both' or mat_pat == 'mat') and (current_state[sib1_mat_index] != -1) and (current_state[sib2_mat_index] != -1):
					is_mat_match[sibpair_index, current_start_index:] = int(int(current_state[sib1_mat_index]==current_state[sib2_mat_index])==flip)
				if (mat_pat == 'both' or mat_pat == 'pat') and (current_state[sib1_pat_index] != -1) and (current_state[sib2_pat_index] != -1):
					is_pat_match[sibpair_index, current_start_index:] = int(int(current_state[sib1_pat_index]==current_state[sib2_pat_index])==flip)
	
	if mat_pat == 'mat':
		assert np.all(is_pat_match==-1)
	if mat_pat == 'pat':
		assert np.all(is_mat_match==-1)	

	return is_mat_match, is_pat_match

def calculate_pvalues(contingency):

	# pos
	pvalues = np.ones((contingency.shape[1],))
	for interval_index in range(contingency.shape[1]):
		try:
			pvalues[interval_index] = scipy.stats.binom_test(contingency[0, interval_index], contingency[1, interval_index],
					p=0.5, alternative='greater') 
		except:
			pass

	return pvalues


if __name__ == "__main__":

	# pull ped info
	sample_to_affected, sample_to_sex, parents_to_children = pull_phenotype_ped(args.ped_file,)

	# pull sibpairs
	sibpairs, family_to_inds = pull_sibpairs(args.phase_dir, sample_to_affected, sample_to_sex, parents_to_children, args.identicals_file)

	print('Overall')
	print('families', len(family_to_inds))
	print('sibpairs', len(sibpairs))
	print('num_affected', Counter([x.num_affected for x in sibpairs]))

	for num_affected in range(3):

		sibpairs_of_interest = [sp for sp in sibpairs if sp.num_affected == num_affected]
		print('Num affected', num_affected, 'sibpairs used', len(sibpairs_of_interest))

		for mat_pat in ['mat', 'pat', 'both']:

			interval_bins = pull_intervals(args.chrom, sibpairs_of_interest, family_to_inds, mat_pat)
			print('intervals', len(interval_bins)-1)

			is_mat_match, is_pat_match = pull_sibpair_matches(args.chrom, sibpairs_of_interest, family_to_inds, interval_bins, mat_pat)
				
			answer = 1 if (num_affected==0 or num_affected==2) else 0

			if mat_pat == 'mat':
				contingency = np.vstack((np.sum(is_mat_match==answer, axis=0)[np.newaxis, :], np.sum(is_mat_match!=-1, axis=0)[np.newaxis, :]))
			elif mat_pat == 'pat':
				contingency = np.vstack((np.sum(is_pat_match==answer, axis=0)[np.newaxis, :], np.sum(is_pat_match!=-1, axis=0)[np.newaxis, :]))
			elif mat_pat == 'both':
				contingency = np.vstack((np.sum(is_mat_match==answer, axis=0)[np.newaxis, :], np.sum(is_mat_match!=-1, axis=0)[np.newaxis, :]))
				contingency += np.vstack((np.sum(is_pat_match==answer, axis=0)[np.newaxis, :], np.sum(is_pat_match!=-1, axis=0)[np.newaxis, :]))
			else:
				raise Exception('mat/pat invalid')

			pvalues = calculate_pvalues(contingency)
			print('pvalues computed')

			print(pvalues.shape)
			np.save('%s/chr.%s.IST.baseline.pvalues.be.aff%d.%s' % (args.phase_dir, args.chrom, num_affected, mat_pat), pvalues)
			np.save('%s/chr.%s.IST.baseline.pvalues.be.aff%d.%s.contingency' % (args.phase_dir, args.chrom, num_affected, mat_pat), contingency)
			np.save('%s/chr.%s.IST.baseline.pvalues.be.aff%d.%s.regions' % (args.phase_dir, args.chrom, num_affected, mat_pat), interval_bins)
			print('results saved to %s/chr.%s.IST.baseline.pvalues.be.aff%d.%s' % (args.phase_dir, args.chrom, num_affected, mat_pat))


