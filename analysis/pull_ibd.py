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
import sys

dataset = sys.argv[1]
interval_chrom, interval_start_pos, interval_end_pos = None, None, None

output_file = dataset
if interval_chrom is not None:
	output_file+= '.chr%s' % interval_chrom
if interval_start_pos is not None or interval_end_pos is not None:
	output_file += '.%d-%d' % (interval_start_pos, interval_end_pos)

def header_to_inds(header):
	header = header.strip().split('\t')
	return [header[i][:-4] for i in range(5, len(header)-3, 2)]

with open('../PhasingFamilies/recomb_%s/sibpairs.json' % dataset, 'r') as f:
	sibpairs = json.load(f)

print('Overall')
print('families', len(set([x['family'].split('.')[0] for x in sibpairs])))
print('sibpairs', len(sibpairs))
num_sibpairs = len(sibpairs)

def apply_interval_filter(chrom, start_pos, end_pos):
	if interval_start_pos is not None or interval_end_pos is not None:
		start_pos = np.clip(start_pos, interval_start_pos, interval_end_pos)
		end_pos = np.clip(end_pos, interval_start_pos, interval_end_pos)
	is_ok = (interval_chrom is None or interval_chrom == chrom) and (end_pos-start_pos>0)
	return is_ok, start_pos, end_pos

def process_phase_file(sibpair):
	with open('../PhasingFamilies/%s/%s.phased.txt' % (sibpair['phase_dir'], sibpair['family']), 'r')  as f:
		header = next(f) # skip header
		inds = header_to_inds(header)
		sib1_ind_index, sib2_ind_index = inds.index(sibpair['sibling1']), inds.index(sibpair['sibling2'])
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
		if chrom != 'X':
			#print(chrom, start_pos, end_pos, state)
			positions.add((chrom, start_pos))
			positions.add((chrom, end_pos))

positions = sorted(positions, key=lambda x: (int(x[0]), x[1]))
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

print(is_mat_match.shape, is_pat_match.shape)

np.save('permutation_tests/phen.%s.chroms.npy' % output_file, chroms)
np.save('permutation_tests/phen.%s.intervals.npy' % output_file, np.array([interval_starts, interval_ends]))
np.save('permutation_tests/phen.%s.mat_ibd.npy' % output_file, is_mat_match)
np.save('permutation_tests/phen.%s.pat_ibd.npy' % output_file, is_pat_match)



