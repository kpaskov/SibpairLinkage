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
num_trials = 1000
interval_chrom, interval_start_pos, interval_end_pos = None, None, None

output_file = dataset
if interval_chrom is not None:
	output_file+= '.chr%s' % interval_chrom
if interval_start_pos is not None or interval_end_pos is not None:
	output_file += '.%d-%d' % (interval_start_pos, interval_end_pos)


with open('../PhasingFamilies/recomb_%s/sibpairs.json' % dataset, 'r') as f:
	sibpairs = json.load(f)

is_mat_match = np.load('permutation_tests/phen.%s.mat_ibd.npy' % output_file)
is_pat_match = np.load('permutation_tests/phen.%s.pat_ibd.npy' % output_file)

# simple binomial test, ignoring sibling structure
ibd = np.zeros((2, is_mat_match.shape[1]))
ibd[0, :] = np.sum(is_mat_match==1, axis=0) + np.sum(is_pat_match==1, axis=0)
ibd[1, :] = np.sum(is_mat_match==-1, axis=0) + np.sum(is_pat_match==-1, axis=0)

ibd_reduced, ibd_inverse = np.unique(ibd, axis=1, return_inverse=True)
print(ibd.shape[1], ibd_reduced.shape)

			
num_intervals = is_match_reduced.shape[1]

pvalues = np.ones((ibd_reduced.shape[1],))
for i in range(ibd_reduced.shape[1]):
	pvalues[i] = stats.binom_test(ibd_reduced[:, i])


np.save('permutation_tests/binombaseline.%s.npy' % output_file, pvalues[reduced_inverse])



