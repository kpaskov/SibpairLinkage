import numpy as np
import scipy.linalg
from collections import defaultdict, namedtuple, Counter
from itertools import combinations
import numpy as np
import scipy.stats as stats
import sys
from os import listdir
import json
import random
import csv

# dataset_name = 'spark+ancestry'
# ped_files = ['../DATA/spark/spark.ped.quads.ped',
# '../DATA/ancestry/ancestry.ped.quads.ped']
# phase_dirs = ['../PhasingFamilies/phased_spark20190423_quads',
# '../PhasingFamilies/phased_ancestry_quads_38']
# identicals_files = ['../PhasingFamilies/sibpair_similarity/spark_quads_identicals.txt',
# '../PhasingFamilies/sibpair_similarity/ancestry_quads_identicals.txt']
# num_trials = 1000
# num_affected = 2
# interval_chrom, interval_start_pos, interval_end_pos = None, None, None
# num_males = None

#dataset_name = 'spark'
#ped_files = ['../DATA/spark/sparkfam.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_spark_array_quads']
#identicals_files = ['../PhasingFamilies/phased_spark_array_quads/identicals.txt']
#num_trials = 1000
#interval_chrom, interval_start_pos, interval_end_pos = None, None, None
#num_males = None
#na = 0

dataset_name = 'spark'
ped_files = ['../DATA/spark/sparkfam.ped.quads.ped']
phase_dirs = ['../PhasingFamilies/phased_spark_wes1_array_quads_del',
              '../PhasingFamilies/phased_spark_wes2_array_quads_del',
              '../PhasingFamilies/phased_spark_wes3_array_quads_del',
              '../PhasingFamilies/phased_spark_wgs1_b01_array_quads_del',
              '../PhasingFamilies/phased_spark_wgs1_b02_array_quads_del',
              '../PhasingFamilies/phased_spark_wgs2_array_quads_del',
              '../PhasingFamilies/phased_spark_wgs3_array_quads_del',
              ]


interval_chrom, interval_start_pos, interval_end_pos = None, None, None
#dataset_name = 'ihart.chip'
#ped_files = ['../DATA/ihart.chip/ihart.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_ihart.chip_quads38_del']
#dataset_name = 'ancestry'
#ped_files = ['../DATA/ancestry/ancestry.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_ancestry_quads_del']
num_trials = 100
#interval_chrom, interval_start_pos, interval_end_pos = '22', 48443005, 48564598
#interval_chrom, interval_start_pos, interval_end_pos = '8', 72897465-1000000, 73361654+1000000
#interval_chrom, interval_start_pos, interval_end_pos = '17', 6426749-1000000, 6978790+1000000
#interval_chrom, interval_start_pos, interval_end_pos = '10', 125067164-1000000, 126635114+1000000
na = 3
flip = False


#individuals_of_interest = set()
#if family_type is not None:
#	with open('../PhasingFamilies/phenotypes/spark/individuals.csv', 'r') as f:
#		reader = csv.reader(f)
#		for pieces in reader:
#			if pieces[2] == family_type:
#				individuals_of_interest.add(pieces[0])
#
#dataset_name = 'spark.exome_missing_parent'
#ped_files = ['../DATA/spark/spark.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_spark.exome_quads_missing_parent']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/spark_quads_identicals.txt']
#num_trials = 1000
#num_affected = 0
#interval_chrom, interval_start_pos, interval_end_pos = "10", 123796859, 126917556
#num_males = None

#dataset_name = 'ihart'
#ped_files = ['../DATA/ihart.ms2/ihart.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_ihart.ms2_quads']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/ihart.ms2_quads_identicals.txt']
#num_trials = 100
#na = 2
#flip = True
#interval_chrom, interval_start_pos, interval_end_pos = None, None, None

#dataset_name = 'mssng_common'
#ped_files = ['../DATA/mssng/mssng.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_mssng_quads_common']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/mssng_quads_identicals.txt']
#num_trials = 1000
#num_affected = 0
#interval_chrom, interval_start_pos, interval_end_pos = "10", 123796859, 126917556
#num_males = None

#dataset_name = 'ssc'
#ped_files = ['../DATA/ssc.hg38/ssc.ped']
#phase_dirs = ['../PhasingFamilies/phased_ssc.hg38']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/ssc_identicals.txt']
#num_trials = 1000
#num_affected = 0
#interval_chrom, interval_start_pos, interval_end_pos = "10", 125067164, 126635114
#num_males = None

#dataset_name = 'ancestry'
#ped_files = ['../DATA/ancestry/ancestry.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_ancestry_quads']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/ancestry_quads_identicals.txt']
#num_trials = 100
#interval_chrom, interval_start_pos, interval_end_pos = "10", 125067164, 126635114
#interval_chrom, interval_start_pos, interval_end_pos = "14", 89209011, 89741920
#interval_chrom, interval_start_pos, interval_end_pos = "17", 6721825, 6908211
#flip = False
#na = 0

#dataset_name = 'mssng'
#ped_files = ['../DATA/mssng/mssng.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_mssng_quads']
#identicals_files = ['../PhasingFamilies/sibpair_similarity/mssng_quads_identicals.txt']
#num_trials = 1000
#num_affected = 1
#interval_chrom, interval_start_pos, interval_end_pos = "10", 125067164, 126635114
#num_males = None

# pull phenotype data
sample_to_affected, sample_to_sex, parents_to_children = dict(), dict(), defaultdict(set)
for ped_file in ped_files:
	with open(ped_file, 'r') as f:
		for line in f:
			pieces = line.strip().split('\t')
			sample_to_sex[pieces[1]] = pieces[4]
			sample_to_sex[pieces[3]] = '2'
			sample_to_sex[pieces[2]] = '1'
			sample_to_affected[pieces[1]] = pieces[5]
			parents_to_children[(pieces[0], pieces[3], pieces[2])].add(pieces[1])
print('families', len(parents_to_children))

# pull sibpairs
Sibpair = namedtuple('Sibpair', ['family', 'sibling1', 'sibling2', 'mom', 'dad', 'phase_dir', 
                                 'num_affected', 'num_males', 
                                 'sibling1_aff', 'sibling2_aff', 'sibling1_male', 'sibling2_male'])
def form_sibpair(family, sibling1, sibling2, mom, dad):
	return Sibpair(family, sibling1, sibling2, mom, dad,
			phase_dir,
			int(sample_to_affected[sibling1]=='2')+int(sample_to_affected[sibling2]=='2'),
			int(sample_to_sex[sibling1]=='1')+int(sample_to_sex[sibling2]=='1'),
			sample_to_affected[sibling1]=='2', sample_to_affected[sibling2]=='2',
			sample_to_sex[sibling1]=='1', sample_to_sex[sibling2]=='1')

def header_to_inds(header):
	header = header.strip().split('\t')
	return [header[i][:-4] for i in range(5, len(header)-3, 2)]

sibpairs, family_to_inds = dict(), dict()
for phase_dir in phase_dirs:
	# pull identicals
	leave_out = set()
	with open('%s/identicals.txt' % phase_dir, 'r') as f:
		for line in f:
			pieces = line.strip().split('\t')
			sibling1, sibling2 = pieces[:2]
			leave_out.add((sibling1, sibling2))
			leave_out.add((sibling2, sibling1))

	# pull sibpairs with phase data
	sibpair_has_phase_data = set()
	for filename in listdir(phase_dir):
		if filename.endswith('.phased.txt'):
			family_key = filename[:-11]
			try:
				with open('%s/%s' % (phase_dir, filename), 'r')  as f:
					header = next(f)
					individuals = header_to_inds(header)
					for sibling1, sibling2 in combinations(individuals[2:], 2):
						if (sibling1, sibling2) not in leave_out:
							sibpair_has_phase_data.add((sibling1, sibling2))
							sibpair_has_phase_data.add((sibling2, sibling1))
			except StopIteration:
				pass

	# pull sibpairs from families
	for (family, mom, dad), children in parents_to_children.items():
		for sibling1, sibling2 in combinations([x for x in children if x in sample_to_affected], 2):
			if (sibling1, sibling2) in sibpair_has_phase_data:
				sibpair = form_sibpair(family, sibling1, sibling2, mom, dad)
				sibpairs[(sibpair.family, sibpair.sibling1, sibpair.sibling2)] = sibpair

sibpairs = sorted(sibpairs.values())

print('Overall')
print('families', len(set([x.family.split('.')[0] for x in sibpairs])))
print('sibpairs', len(sibpairs))
print('num_affected', Counter([x.num_affected for x in sibpairs]))

if na ==3:
	pass
else:
	sibpairs = [x for x in sibpairs if x.num_affected==na]
num_sibpairs = len(sibpairs)

print('Overall')
print('families', len(set([x.family.split('.')[0] for x in sibpairs])))
print('sibpairs', len(sibpairs))
print('num_affected', Counter([x.num_affected for x in sibpairs]))

with open('permutation_tests/%s.%d.%ssibpairs.json' % (dataset_name, na, 'flip.' if flip else ''), 'w+') as f:
	json.dump(sibpairs, f)

def apply_interval_filter(chrom, start_pos, end_pos):
	if interval_start_pos is not None or interval_end_pos is not None:
		start_pos = np.clip(start_pos, interval_start_pos, interval_end_pos)
		end_pos = np.clip(end_pos, interval_start_pos, interval_end_pos)
	is_ok = (interval_chrom is None or interval_chrom == chrom) and (end_pos-start_pos>0)
	return is_ok, start_pos, end_pos

def process_phase_file(sibpair):
	with open('%s/%s.phased.txt' % (sibpair.phase_dir, sibpair.family), 'r')  as f:
		header = next(f) # skip header
		inds = header_to_inds(header)
		sib1_ind_index, sib2_ind_index = inds.index(sibpair.sibling1), inds.index(sibpair.sibling2)
		sib1_mat_index, sib2_mat_index = 4+(2*sib1_ind_index), 4+(2*sib2_ind_index)
		sib1_pat_index, sib2_pat_index = 5+(2*sib1_ind_index), 5+(2*sib2_ind_index)
		sib_phase_indices = [sib1_mat_index, sib2_mat_index, sib1_pat_index, sib2_pat_index]

		current_chrom, current_start_pos, current_end_pos, current_state = None, None, None, None
		for line in f:
			pieces = line.strip().split('\t')
			chrom = pieces[0][3:]
			start_pos, end_pos = [int(x) for x in pieces[-2:]]
			state = np.array([int(x) for x in pieces[1:-2]])[sib_phase_indices]

			if current_chrom is None:
				current_chrom, current_start_pos, current_end_pos, current_state = chrom, start_pos, end_pos, state
			elif current_chrom != chrom or np.any(current_state != state):
				is_ok, current_start_pos, current_end_pos = apply_interval_filter(current_chrom, current_start_pos, current_end_pos)
				if is_ok:
					yield current_chrom, current_start_pos, current_end_pos, current_state
				current_chrom, current_start_pos, current_end_pos, current_state = chrom, start_pos, end_pos, state
			else:
				current_end_pos = end_pos
		is_ok, current_start_pos, current_end_pos = apply_interval_filter(current_chrom, current_start_pos, current_end_pos)
		if is_ok:
			yield current_chrom, current_start_pos, current_end_pos, current_state

# pull intervals
positions = set()
for sibpair in sibpairs:
	#print(sibpair)
	for chrom, start_pos, end_pos, state in process_phase_file(sibpair):
		#print(chrom, start_pos, end_pos, state)
		positions.add((chrom, start_pos))
		positions.add((chrom, end_pos))

positions = sorted(positions, key=lambda x: (int(x[0]), x[1]) if x[0].isdigit() else x)
chroms, interval_starts, interval_ends = [], [], []
prev_chrom, prev_pos = None, None
for c, p in positions:
	if prev_chrom is not None and prev_chrom == c:
		chroms.append(c)
		interval_starts.append(prev_pos)
		interval_ends.append(p)
	prev_chrom, prev_pos = c, p


interval_starts = np.array(interval_starts)
interval_ends = np.array(interval_ends)
num_intervals = len(interval_starts)
print('intervals', num_intervals)
#print(interval_starts)
#print(interval_ends)
#print(interval_ends-interval_starts)


# pull sibpair IBD

# sibpair, interval
is_mat_match = np.zeros((num_sibpairs, num_intervals), dtype=int)
is_pat_match = np.zeros((num_sibpairs, num_intervals), dtype=int)

interval_start_to_index = dict([((chrom, x), i) for i, (chrom, x) in enumerate(zip(chroms, interval_starts))])
interval_end_to_index = dict([((chrom, x), i) for i, (chrom, x) in enumerate(zip(chroms, interval_ends))])

for sibpair_index, sibpair in enumerate(sibpairs):
	for chrom, start_pos, end_pos, state in process_phase_file(sibpair):
		start_index, end_index = interval_start_to_index[(chrom, start_pos)], interval_end_to_index[(chrom, end_pos)]+1

		if state[0]==-1 or state[1]==-1:
			pass
		elif state[0]==state[1]:
			is_mat_match[sibpair_index, start_index:end_index] = 1
		else:
			is_mat_match[sibpair_index, start_index:end_index] = -1
			
		if state[2]==-1 or state[3]==-1:
				pass
		elif state[2]==state[3]:
			is_pat_match[sibpair_index, start_index:end_index] = 1
		else:
			is_pat_match[sibpair_index, start_index:end_index] = -1


is_ok = interval_ends - interval_starts > 1
interval_starts = interval_starts[is_ok]
interval_ends = interval_ends[is_ok]
chroms = np.array([int(c) for c in chroms])[is_ok]
is_mat_match = is_mat_match[:, is_ok]
is_pat_match = is_pat_match[:, is_ok]
num_intervals = np.sum(is_ok)

np.save('permutation_tests/%s.%d.%sis_mat_match.npy' % (dataset_name, na, 'flip.' if flip else ''), is_mat_match)
np.save('permutation_tests/%s.%d.%sis_pat_match.npy' % (dataset_name, na, 'flip.' if flip else ''), is_pat_match)



print(is_mat_match.shape, is_pat_match.shape)

# take into account sibling structure across quads
individuals = sorted(set([x.sibling1 for x in sibpairs] + [x.sibling2 for x in sibpairs]))
ind_to_index = dict([(x, i) for i, x in enumerate(individuals)])
sibling1_indices = np.array([ind_to_index[x.sibling1] for x in sibpairs])
sibling2_indices = np.array([ind_to_index[x.sibling2] for x in sibpairs])

A = np.random.randint(0, high=2, size=(num_trials+1, len(individuals), 2))
X1 = (A[:, sibling1_indices, 0] == A[:, sibling2_indices, 0]).astype(int)
X2 = (A[:, sibling1_indices, 1] == A[:, sibling2_indices, 1]).astype(int)

# randomly flip IBD in sibpairs
#X1 = np.random.randint(0, high=2, size=(num_trials+1, len(sibpairs)))
#X2 = np.random.randint(0, high=2, size=(num_trials+1, len(sibpairs)))

X1[X1==0] = -1
X2[X2==0] = -1

# first entry is actual IBD relationships
X1[0, :] = 1
X2[0, :] = 1

print('ready')

if interval_chrom is not None:
	dataset_name += '.chr%s' % interval_chrom
if interval_start_pos is not None or interval_end_pos is not None:
	dataset_name += '.%d-%d' % (interval_start_pos, interval_end_pos)

# trial, interval, mat/pat
rand_pvalue = np.zeros((num_trials+1, num_intervals, 3), dtype=int)

print(na, 'mat')
rand_pvalue[:, :, 0] = X1.dot(is_mat_match)
print(na, 'pat')
rand_pvalue[:, :, 1] = X2.dot(is_pat_match)
print(na, 'both')
rand_pvalue[:, :, 2] = rand_pvalue[:, :, 0]+rand_pvalue[:, :, 1]


# we expect to see less IBD sharing between discordant sibpairs
if flip:
	rand_pvalue = -rand_pvalue

# -------------------- implementing Westfall-Young max T stepdown procedure

# indices are sorted along interval axis from interval with most IBD sharing
# to least IBD sharing
final_pvalues = np.zeros((num_intervals, 3))
for is_mat in range(3):

	orig_indices = np.flip(np.argsort(rand_pvalue[0, :, is_mat]))

	max_t_k = np.zeros((num_trials+1, num_intervals+1))
	max_t_k[:, -1] = np.min(rand_pvalue[:, :, is_mat], axis=1)
	for i, j in list(reversed(list(enumerate(orig_indices)))):
		max_t_k[:, i] = np.maximum(max_t_k[:, i+1], rand_pvalue[:, j, is_mat])
	max_t_k = max_t_k[:, :-1]

	#max_t_k = np.flip(np.sort(rand_pvalue[:, :, is_mat], axis=1), axis=1)
	
	assert np.all(max_t_k[0, :] == rand_pvalue[0, orig_indices, is_mat])

	# calculate pi(j)
	pvalues = np.sum(max_t_k[1:, :] >= np.tile(max_t_k[0, :], (num_trials, 1)), axis=0)/num_trials
	pvalues = np.array([np.max(pvalues[:(i+1)]) for i in np.arange(pvalues.shape[0])])
	final_pvalues[orig_indices, is_mat] = pvalues

np.save('permutation_tests/%s.%d.%snpy' % (dataset_name, na, 'flip.' if flip else ''), final_pvalues)
np.save('permutation_tests/%s.%d.%schroms.npy' % (dataset_name, na, 'flip.' if flip else ''), chroms)
np.save('permutation_tests/%s.%d.%sintervals.npy' % (dataset_name, na, 'flip.' if flip else ''), np.array([interval_starts, interval_ends]))



