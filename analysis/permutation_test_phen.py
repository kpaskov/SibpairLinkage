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

dataset = 'spark'
num_trials = 1000
interval_chrom, interval_start_pos, interval_end_pos = None, None, None
phen_index = int(sys.argv[1])

output_file = dataset
if interval_chrom is not None:
	output_file+= '.chr%s' % interval_chrom
if interval_start_pos is not None or interval_end_pos is not None:
	output_file += '.%d-%d' % (interval_start_pos, interval_end_pos)

#phase_dirs = ['../PhasingFamilies/phased_spark_wes1_array_quads_del',
#              '../PhasingFamilies/phased_spark_wes2_array_quads_del',
#              '../PhasingFamilies/phased_spark_wes3_array_quads_del',
#              '../PhasingFamilies/phased_spark_wgs1_b01_array_quads_del',
#              '../PhasingFamilies/phased_spark_wgs1_b02_array_quads_del',
#              '../PhasingFamilies/phased_spark_wgs2_array_quads_del',
#              '../PhasingFamilies/phased_spark_wgs3_array_quads_del',
#              ]

#dataset = 'ssc.hg38'
#ped_files = ['../DATA/ssc.hg38/ssc.ped']
#phase_dirs = ['../PhasingFamilies/phased_ssc.hg38_del',
#              ]

#dataset = 'mssng'
#ped_files = ['../DATA/mssng/mssng.ped.quads.ped']
#phase_dirs = ['../PhasingFamilies/phased_mssng_quads',
#              ]

#interval_chrom, interval_start_pos, interval_end_pos = '8', 72897465, 73361654
#interval_chrom, interval_start_pos, interval_end_pos = '17', 6426749, 6978790
#interval_chrom, interval_start_pos, interval_end_pos = '10', 125067164, 126635114

with open('../PhasingFamilies/recomb_%s/sibpairs.json' % dataset, 'r') as f:
	sibpairs = json.load(f)

is_mat_match = np.load('permutation_tests/phen.%s.mat_ibd.npy' % output_file)
is_pat_match = np.load('permutation_tests/phen.%s.pat_ibd.npy' % output_file)

# take into account sibling structure across quads
individuals = sorted(set([x['sibling1'] for x in sibpairs] + [x['sibling2'] for x in sibpairs]))
ind_to_index = dict([(x, i) for i, x in enumerate(individuals)])
sibling1_indices = np.array([ind_to_index[x['sibling1']] for x in sibpairs])
sibling2_indices = np.array([ind_to_index[x['sibling2']] for x in sibpairs])

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

print('SCQ', phen_index)
sample_to_affected = dict()
with open('../PhasingFamilies/phenotypes/spark_v5/spark_v5-scq-prep.csv', 'r') as f:
	reader = csv.reader(f)
	for pieces in reader:
		phen = pieces[13+phen_index]
		if phen=='1.0' or phen=='0.0':
			sample_to_affected[pieces[2]] = '1' if phen =='1.0' else '0'


aut_aut_na_response = [0]*2 + [2]*6 + [0] + [2]*9 + [0]*22

num_affected = np.array([-1 if (x['sibling1'] not in sample_to_affected or x['sibling2'] not in sample_to_affected) else int(sample_to_affected[x['sibling1']])+int(sample_to_affected[x['sibling2']]) for x in sibpairs])
#num_affected = np.array([-1 if (x.family + '.p1' not in sample_to_affected or x.family + '.s1' not in sample_to_affected) else int(sample_to_affected[x.family + '.p1'])+int(sample_to_affected[x.family + '.s1']) for x in sibpairs])
print(Counter(num_affected))
na = aut_aut_na_response[phen_index]
print(na, np.sum(num_affected==na), end=' ')


is_match_reduced, reduced_inverse = np.unique(np.vstack((is_mat_match[num_affected==na, :],
                                                                 is_pat_match[num_affected==na, :])), axis=1, return_inverse=True)
print(is_mat_match.shape[1], is_match_reduced.shape)

			
num_intervals = is_match_reduced.shape[1]

all_pvalues_reduced = np.zeros((num_intervals, ))

# trial, interval, mat/pat
rand_pvalue = np.hstack((X1[:, num_affected==na], X2[:, num_affected==na])).dot(is_match_reduced)

# -------------------- implementing Westfall-Young max T stepdown procedure

# indices are sorted along interval axis from interval with most IBD sharing
# to least IBD sharing

orig_indices = np.flip(np.argsort(rand_pvalue[0, :]), axis=0)

max_t_k = np.zeros((num_trials+1, num_intervals+1))
max_t_k[:, -1] = np.min(rand_pvalue, axis=1)
for i, j in list(reversed(list(enumerate(orig_indices)))):
	max_t_k[:, i] = np.maximum(max_t_k[:, i+1], rand_pvalue[:, j])
max_t_k = max_t_k[:, :-1]
				
assert np.all(max_t_k[0, :] == rand_pvalue[0, orig_indices])

# calculate pi(j)
pvalues = np.sum(max_t_k[1:, :] >= np.tile(max_t_k[0, :], (num_trials, 1)), axis=0)/num_trials
pvalues = np.array([np.max(pvalues[:(i+1)]) for i in np.arange(num_intervals)])
all_pvalues_reduced[orig_indices] = pvalues

all_pvalues = all_pvalues_reduced[reduced_inverse]
print(all_pvalues.shape)

print(np.min(all_pvalues))

np.save('permutation_tests/scq%d.%s.npy' % (phen_index+1, output_file), all_pvalues)



