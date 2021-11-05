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
if len(sys.argv)>2 and sys.argv[2]=='--conservative':
	conservative = True
else:
	conservative = False
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
			if conservative and pieces[-3] == '1':
				state[:] = -1

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
with open('../PhasingFamilies/data/cytoband38.txt', 'r') as f:
	next(f) # skip header
	for line in f:
		pieces = line.strip().split('\t')
		chrom = pieces[0][3:]
		if chrom in [str(x) for x in range(1, 23)]:
			start_pos, end_pos = int(pieces[1]), int(pieces[2])
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

chrom_to_starts = defaultdict(list)
chrom_to_ends = defaultdict(list)
chrom_to_indices = defaultdict(list)
for i, (chrom, start, end) in enumerate(zip(chroms, interval_starts, interval_ends)):
	chrom_to_starts[chrom].append(start)
	chrom_to_ends[chrom].append(end)
	chrom_to_indices[chrom].append(i)

for sibpair_index, sibpair in enumerate(sibpairs):
	ibd_bp_mat = np.zeros((num_intervals,))
	noibd_bp_mat = np.zeros((num_intervals,))
	ibd_bp_pat = np.zeros((num_intervals,))
	noibd_bp_pat = np.zeros((num_intervals,))

	current_interval_index = 0
	for chrom, start_pos, end_pos, state in process_phase_file(sibpair):
		if chrom != 'X':
			if state[0]==-1 or state[1]==-1:
				pass
			elif state[0]==state[1]:
				starts, ends, indices = chrom_to_starts[chrom], chrom_to_ends[chrom], chrom_to_indices[chrom]
				start_index = np.searchsorted(starts, start_pos)
				end_index = np.searchsorted(starts, end_pos)
					
				ibd_bp_mat[indices[(start_index-1):(end_index+1)]] += np.clip(np.minimum(end_pos, ends[(start_index-1):(end_index+1)])-np.maximum(start_pos, starts[(start_index-1):(end_index+1)]), 0, None)
			else:
				starts, ends, indices = chrom_to_starts[chrom], chrom_to_ends[chrom], chrom_to_indices[chrom]
				start_index = np.searchsorted(starts, start_pos)
				end_index = np.searchsorted(starts, end_pos)

				noibd_bp_mat[indices[(start_index-1):(end_index+1)]] += np.clip(np.minimum(end_pos, ends[(start_index-1):(end_index+1)])-np.maximum(start_pos, starts[(start_index-1):(end_index+1)]), 0, None)

				
			if state[2]==-1 or state[3]==-1:
					pass
			elif state[2]==state[3]:
				starts, ends, indices = chrom_to_starts[chrom], chrom_to_ends[chrom], chrom_to_indices[chrom]
				start_index = np.searchsorted(starts, start_pos)
				end_index = np.searchsorted(starts, end_pos)
					
				ibd_bp_pat[indices[(start_index-1):(end_index+1)]] += np.clip(np.minimum(end_pos, ends[(start_index-1):(end_index+1)])-np.maximum(start_pos, starts[(start_index-1):(end_index+1)]), 0, None)

			else:
				starts, ends, indices = chrom_to_starts[chrom], chrom_to_ends[chrom], chrom_to_indices[chrom]
				start_index = np.searchsorted(starts, start_pos)
				end_index = np.searchsorted(starts, end_pos)

				noibd_bp_pat[indices[(start_index-1):(end_index+1)]] += np.clip(np.minimum(end_pos, ends[(start_index-1):(end_index+1)])-np.maximum(start_pos, starts[(start_index-1):(end_index+1)]), 0, None)

	is_mat_match[sibpair_index, ibd_bp_mat>noibd_bp_mat] = 1
	is_mat_match[sibpair_index, ibd_bp_mat<noibd_bp_mat] = -1
	is_pat_match[sibpair_index, ibd_bp_pat>noibd_bp_pat] = 1
	is_pat_match[sibpair_index, ibd_bp_pat<noibd_bp_pat] = -1


print(np.sum(is_mat_match), np.sum(is_pat_match))

np.save('permutation_tests/phen.%s.cytoband.chroms.npy' % output_file, chroms)
np.save('permutation_tests/phen.%s.cytoband.intervals.npy' % output_file, np.array([interval_starts, interval_ends]))
np.save('permutation_tests/phen.%s.cytoband.mat_ibd.npy' % output_file, is_mat_match)
np.save('permutation_tests/phen.%s.cytoband.pat_ibd.npy' % output_file, is_pat_match)



