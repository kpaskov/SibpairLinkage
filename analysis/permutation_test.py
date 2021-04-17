import numpy as np
import scipy.linalg
from collections import defaultdict, namedtuple, Counter
from itertools import combinations
import numpy as np
import scipy.stats
import sys
from os import listdir
import json
import random

# dataset_name = 'spark+ancestry'
# ped_files = ['../DATA/spark/spark.ped.quads.ped',
# '../DATA/ancestry/ancestry.ped.quads.ped'
# ]
# phase_dirs = ['../PhasingFamilies/phased_spark_quads',
# '../PhasingFamilies/phased_ancestry_quads_nohap_liftover38'
# ]
# identicals_files = ['../PhasingFamilies/sibpair_similarity/spark_quads_identicals.txt',
# '../PhasingFamilies/sibpair_similarity/ancestry_quads_identicals.txt'
# ]
# num_trials = 1000
# num_affected = 2
# interval_chrom, interval_start_pos, interval_end_pos = None, None, None

#dataset_name = 'spark'
#ped_files = ['../DATA/spark/spark.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_spark_quads']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/spark_quads_identicals.txt']
#num_trials = 1000
#num_affected = 0
#interval_chrom, interval_start_pos, interval_end_pos = "10", 123796859, 126917556


#dataset_name = 'ihart_nopass'
#ped_files = ['../DATA/ihart.ms2/ihart.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_ihart.ms2_quads_nopass']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/ihart.ms2_quads_identicals.txt']
#num_trials = 1000
#num_affected = 0
#interval_chrom, interval_start_pos, interval_end_pos = "10", 123796859, 126917556
#num_males = 2

dataset_name = 'mssng_common'
ped_files = ['../DATA/mssng/mssng.ped.quads.ped']
phase_dirs = ['../PhasingFamilies/phased_mssng_quads_common']
identicals_files = ['../PhasingFamilies/sibpair_similarity/mssng_quads_identicals.txt']
num_trials = 1000
num_affected = 0
interval_chrom, interval_start_pos, interval_end_pos = "10", 123796859, 126917556
num_males = None

#dataset_name = 'ssc_mf'
#ped_files = ['../DATA/ssc.hg38/ssc.ped']
#phase_dirs = ['../PhasingFamilies/phased_ssc.hg38']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/ssc_identicals.txt']
#num_trials = 1000
#num_affected = 1
#interval_chrom, interval_start_pos, interval_end_pos = "10", 123796859, 126917556
#num_males = 1

def pull_phenotype_ped(ped_file):
	sample_to_sex = dict()
	sample_to_affected = dict()
	parents_to_children = defaultdict(list)
	with open(ped_file, 'r') as f:
		for line in f:
			pieces = line.strip().split('\t')
			sample_to_sex[pieces[1]] = pieces[4]
			sample_to_sex[pieces[3]] = '2'
			sample_to_sex[pieces[2]] = '1'
			sample_to_affected[pieces[1]] = pieces[5]
			parents_to_children[(pieces[0], pieces[3], pieces[2])].append(pieces[1])
	return sample_to_affected, sample_to_sex, parents_to_children


Sibpair = namedtuple('Sibpair', ['family', 'sibling1', 'sibling2', 'mom', 'dad', 'phase_dir', 
                                 'num_affected', 'num_males', 
                                 'sibling1_aff', 'sibling2_aff', 'sibling1_male', 'sibling2_male'])
def pull_sibpairs(phase_dir, sample_to_affected, sample_to_sex, parents_to_children, identicals_file):

	# pull identicals
	leave_out = set()
	with open(identicals_file, 'r') as f:
		for line in f:
			pieces = line.strip().split('\t')
			for sibling1, sibling2 in combinations(pieces, 2):
				leave_out.add((sibling1, sibling2))
				leave_out.add((sibling2, sibling1))


	# pull sibpairs with phase data
	sibpair_has_phase_data = set()
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
							sibpair_has_phase_data.add((sibling1, sibling2))
							sibpair_has_phase_data.add((sibling2, sibling1))
			except StopIteration:
				pass



	def form_sibpair(family, sibling1, sibling2, mom, dad):
		return Sibpair(family, sibling1, sibling2, mom, dad,
			phase_dir,
			int(sample_to_affected[sibling1]=='2')+int(sample_to_affected[sibling2]=='2'),
			int(sample_to_sex[sibling1]=='1')+int(sample_to_sex[sibling2]=='1'),
			sample_to_affected[sibling1]=='2', sample_to_affected[sibling2]=='2',
			sample_to_sex[sibling1]=='1', sample_to_sex[sibling2]=='1')
	# pull sibpairs from families
	sibpairs = []
	for (family, mom, dad), children in parents_to_children.items():
		for sibling1, sibling2 in combinations([x for x in children if x in sample_to_affected], 2):
			if (sibling1, sibling2) in sibpair_has_phase_data:
				sibpairs.append(form_sibpair(family, sibling1, sibling2, mom, dad))
	sibpairs = sorted(sibpairs)

	for i in range(len(sibpairs)-1):
		if sibpairs[i] == sibpairs[i+1]:
			print(sibpairs[i])

	assert len(sibpairs) == len(set(sibpairs)) # should have no duplicates
	return sibpairs, family_to_inds


def pull_intervals(chrom, sibpairs, family_to_inds, interval_start_pos=None, interval_end_pos=None):

	positions = set()
	if interval_start_pos is not None:
		positions.add(interval_start_pos)
	else:
		positions.add(0)
	if interval_end_pos is not None:
		positions.add(interval_end_pos)

	for sibpair_index, sibpair in enumerate(sibpairs):
		inds = family_to_inds[sibpair.family]
		sib1_ind_index, sib2_ind_index = inds.index(sibpair.sibling1), inds.index(sibpair.sibling2)
		sib1_mat_index, sib2_mat_index = 4+2*sib1_ind_index, 4+2*sib2_ind_index
		sib1_pat_index, sib2_pat_index = 5+2*sib1_ind_index, 5+2*sib2_ind_index
		sib_phase_indices = [sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index]
		

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
						if ((interval_start_pos is None) or start_pos >= interval_start_pos) and ((interval_end_pos is None) or start_pos <= interval_end_pos):
							positions.add(start_pos)
						if ((interval_start_pos is None) or end_pos >= interval_start_pos) and ((interval_end_pos is None) or end_pos <= interval_end_pos):
							positions.add(end_pos)

	return np.array(sorted(positions))

def pull_sibpair_matches(sibpairs, family_to_inds, chroms, interval_bins, interval_chrom=None, interval_start_pos=None, interval_end_pos=None):
	sibpair_to_index = dict([((x.family, x.sibling1, x.sibling2), i) for i, x in enumerate(sibpairs)])
	interval_start_to_index = dict([((chrom, x), i) for i, (chrom, x) in enumerate(zip(chroms, interval_bins[0, :]))])

	# pull phase data
	# sibpair, interval
	#is_mat_match = np.zeros((len(sibpair_to_index), len(interval_bins)-1), dtype=int)
	#is_pat_match = np.zeros((len(sibpair_to_index), len(interval_bins)-1), dtype=int)
	is_mat_match, is_pat_match = [], []

	def update_matches(sibpair_index, chrom, start_pos, end_pos, state, sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index, interval_start_pos, interval_end_pos):
		if interval_start_pos is not None or interval_end_pos is not None:
			start_pos = np.clip(start_pos, interval_start_pos, interval_end_pos)
			end_pos = np.clip(end_pos, interval_start_pos, interval_end_pos)

		start_index, end_index = interval_start_to_index[(chrom, start_pos)], interval_start_to_index[(chrom, end_pos)]
						
		if state[sib1_mat_index]==-1 or state[sib2_mat_index]==-1:
			pass
		elif state[sib1_mat_index]==state[sib2_mat_index]:
			#is_mat_match[sibpair_index, start_index:end_index] = 1
			is_mat_match.append((sibpair_index, start_index, end_index, 1))
		else:
			#is_mat_match[sibpair_index, start_index:end_index] = -1
			is_mat_match.append((sibpair_index, start_index, end_index, -1))
		if state[sib1_pat_index]==-1 or state[sib2_pat_index]==-1:
			pass
		elif state[sib1_pat_index]==state[sib2_pat_index]:
			#is_pat_match[sibpair_index, start_index:end_index] = 1
			is_pat_match.append((sibpair_index, start_index, end_index, 1))
		else:
			#is_pat_match[sibpair_index, start_index:end_index] = -1
			is_pat_match.append((sibpair_index, start_index, end_index, -1))

	for sibpair_index, sibpair in enumerate(sibpairs):
		inds = family_to_inds[sibpair.family]
		sib1_ind_index, sib2_ind_index = inds.index(sibpair.sibling1), inds.index(sibpair.sibling2)
		sib1_mat_index, sib2_mat_index = 4+2*sib1_ind_index, 4+2*sib2_ind_index
		sib1_pat_index, sib2_pat_index = 5+2*sib1_ind_index, 5+2*sib2_ind_index
		sib_phase_indices = [sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index]
		
		with open('%s/%s.phased.txt' % (sibpair.phase_dir, sibpair.family), 'r')  as f:
			next(f) # skip header

			prev_chrom = None
			prev_state = None
			prev_start, prev_end = None, None
			for line in f:
				pieces = line.strip().split('\t')
				
				try:
					chrom = int(pieces[0][3:])
					start_pos, end_pos = [int(x) for x in pieces[-2:]]
					state = np.array([int(x) for x in pieces[1:-2]])


					if (prev_state is None) or (prev_chrom != chrom) or (not np.all(state[sib_phase_indices]==prev_state[sib_phase_indices])):
						# first save prev state
						if prev_state is not None and (interval_chrom is None or interval_chrom==str(prev_chrom)):
							update_matches(sibpair_index, prev_chrom, prev_start, prev_end, prev_state, 
								sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index,
								interval_start_pos, interval_end_pos)
						prev_state = state
						prev_chrom, prev_start, prev_end = chrom, start_pos, end_pos
					else:
						prev_end = end_pos
				except ValueError:
					pass
					

	return is_mat_match, is_pat_match


# pull ped info
sibpairs, family_to_inds = [], dict()
for ped_file, phase_dir, identicals_file in zip(ped_files, phase_dirs, identicals_files):
    sample_to_affected, sample_to_sex, parents_to_children= pull_phenotype_ped(ped_file)

    # pull sibpairs
    sps, f_to_i = pull_sibpairs(phase_dir, sample_to_affected, sample_to_sex, parents_to_children, identicals_file)
    sibpairs.extend(sps)
    family_to_inds.update(f_to_i)

if num_males is not None:
	sibpairs = [x for x in sibpairs if x.num_males==num_males]
    
print('Overall')
print('families', len(family_to_inds))
print('sibpairs', len(sibpairs))
print('num_affected', Counter([x.num_affected for x in sibpairs]))

chroms, interval_bins_starts, interval_bins_ends = [], [], []
if interval_chrom is not None:
	interval_bins = pull_intervals(interval_chrom, sibpairs, family_to_inds, interval_start_pos, interval_end_pos)
	print(interval_chrom, 'intervals', len(interval_bins))

	chroms.append(int(interval_chrom)*np.ones((interval_bins.shape[0],)))
	interval_bins_starts.append(interval_bins)
	interval_bins_ends.append(np.hstack((interval_bins[1:], [interval_bins[-1]])))
else:
	for chrom in [str(x) for x in range(1, 23)]:
		interval_bins = pull_intervals(chrom, sibpairs, family_to_inds, interval_start_pos, interval_end_pos)
		print(chrom, 'intervals', len(interval_bins))

		
		chroms.append(int(chrom)*np.ones((interval_bins.shape[0],)))
		interval_bins_starts.append(interval_bins)
		interval_bins_ends.append(np.hstack((interval_bins[1:], [interval_bins[-1]])))

interval_bins = np.vstack((np.hstack(interval_bins_starts), np.hstack(interval_bins_ends)))
chroms = np.hstack(chroms)

sibpairs_of_interest = [sp for sp in sibpairs if sp.num_affected == num_affected]
print('Num affected', num_affected, 'sibpairs used', len(sibpairs_of_interest))


is_mat_match, is_pat_match = pull_sibpair_matches(sibpairs_of_interest, family_to_inds, chroms, interval_bins, interval_chrom, interval_start_pos, interval_end_pos)
print(len(is_mat_match), len(is_pat_match))
#np.save('permutation_tests/%s.%d.is_mat_match.npy' % (dataset_name, num_affected), is_mat_match)
#np.save('permutation_tests/%s.%d.is_pat_match.npy' % (dataset_name, num_affected), is_pat_match)

# need to take into account sibling structure across quads
individuals = sorted(set([x.sibling1 for x in sibpairs_of_interest] + [x.sibling2 for x in sibpairs_of_interest]))
ind_to_index = dict([(x, i) for i, x in enumerate(individuals)])
sibling1_indices = np.array([ind_to_index[x.sibling1] for x in sibpairs_of_interest])
sibling2_indices = np.array([ind_to_index[x.sibling2] for x in sibpairs_of_interest])

A = np.random.randint(0, high=2, size=(num_trials, len(individuals), 2))
X1 = (A[:, sibling1_indices, 0] == A[:, sibling2_indices, 0]).astype(int)
X2 = (A[:, sibling1_indices, 1] == A[:, sibling2_indices, 1]).astype(int)
#X = np.hstack((X1, X2))
#Y = np.vstack((is_mat_match, is_pat_match))

X1[X1==0] = -1
X2[X2==0] = -1

#print(X.shape)

if interval_chrom is not None:
	dataset_name += '.chr%s' % interval_chrom
if interval_start_pos is not None or interval_end_pos is not None:
	dataset_name += '.%d-%d' % (interval_start_pos, interval_end_pos)

XdotY = np.zeros((num_trials, interval_bins.shape[1]), dtype=int)
actual = np.zeros((interval_bins.shape[1],), dtype=int)
actual_IBD = np.zeros((interval_bins.shape[1],), dtype=int)
actual_noIBD = np.zeros((interval_bins.shape[1],), dtype=int)

for i, (sibpair_index, start_index, end_index, value) in enumerate(is_mat_match):
	#print(XdotY.shape, XdotY[:, (offset+start_index):(offset+end_index)].shape, np.tile((X1[:, sibpair_index]*value)[:, np.newaxis], (1, end_index-start_index)).shape)
	XdotY[:, start_index:end_index] += np.tile((X1[:, sibpair_index]*value)[:, np.newaxis], (1, end_index-start_index))
	actual[start_index:end_index] += value
	if value==1:
		actual_IBD[start_index:end_index] += 1
	elif value==-1:
		actual_noIBD[start_index:end_index] += 1

	if i%10000==0:
		print(i, '/', len(is_mat_match))
print('mat done')

for i, (sibpair_index, start_index, end_index, value) in enumerate(is_pat_match):
	XdotY[:, start_index:end_index] += np.tile((X2[:, sibpair_index]*value)[:, np.newaxis], (1, end_index-start_index))
	actual[start_index:end_index] += value
	if value==1:
		actual_IBD[start_index:end_index] += 1
	elif value==-1:
		actual_noIBD[start_index:end_index] += 1

	if i%10000==0:
		print(i, '/', len(is_pat_match))
print('pat done')

if num_affected==1:
	z = np.min(XdotY, axis=1)
else:
	z = np.max(XdotY, axis=1)

np.save('permutation_tests/%s.%d.npy' % (dataset_name, num_affected), z)
np.save('permutation_tests/%s.%d.actual.npy' % (dataset_name, num_affected), actual)
np.save('permutation_tests/%s.%d.actual_IBD.npy' % (dataset_name, num_affected), actual_IBD)
np.save('permutation_tests/%s.%d.actual_noIBD.npy' % (dataset_name, num_affected), actual_noIBD)
np.save('permutation_tests/%s.%d.chroms.npy' % (dataset_name, num_affected), chroms)
np.save('permutation_tests/%s.%d.intervals.npy' % (dataset_name, num_affected), interval_bins)



