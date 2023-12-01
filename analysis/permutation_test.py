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
import argparse

import sys
sys.path.append('../PhasingFamilies')
sys.path.append('../PhasingFamilies/phase')
from phase.input_output import PhaseData


parser = argparse.ArgumentParser(description='Pull crossovers from phasing output.')
parser.add_argument('dataset_name', type=str, help='Name of test.')
parser.add_argument('data_dir', type=str, help='Directory of genotype data for the cohort in .npy format.')
parser.add_argument('ped_file', type=str, help='.ped file')
parser.add_argument('sibpair_type', type=int, help='Sibpair type to use in test. 0 indicates control-control sibpairs, 1 indicates control-affected sibpairs, 2 indicates affected-affected sibpairs.')
parser.add_argument('--num_trials', type=int, default=1000, help='num permutation trials to run')
parser.add_argument('--interval', type=str, default=None, help='interval to consider')
parser.add_argument('--num_males', type=int, default=None, help='restrict sibpairs by sex')
args = parser.parse_args()

dataset_name = args.dataset_name

if args.num_males is not None:
    dataset_name += '.nm%d' % args.num_males
if args.interval is not None:
    dataset_name += '.%s' % args.interval.replace(':', '.')

# pull phenotype data
sample_to_affected, sample_to_sex = dict(), dict()
with open(args.ped_file, 'r') as f:
    for line in f:
        pieces = line.strip().split('\t')
        sample_to_sex[pieces[1]] = pieces[4]
        sample_to_sex[pieces[3]] = '2'
        sample_to_sex[pieces[2]] = '1'
        sample_to_affected[pieces[1]] = pieces[5]

# pull sibpairs
phase_data = PhaseData(args.data_dir)

sibpairs = phase_data.get_sibpairs()
print('sibpairs', len(sibpairs))

is_fully_phased = np.array([x['is_fully_phased'] for x in sibpairs])
is_identical = np.array([x['is_identical'] for x in sibpairs])
is_ibd_outlier = np.array([x['is_ibd_outlier'] for x in sibpairs], dtype=bool)
is_co_outlier = np.array([x['is_crossover_outlier'] for x in sibpairs], dtype=bool)

print('not fully phased', np.sum(~is_fully_phased))
print('identicals', np.sum(is_identical))
print('ibd outliers', np.sum(is_ibd_outlier))
print('crossover outliers', np.sum(is_co_outlier))

sibpairs = [x for x in phase_data.get_sibpairs() if x['is_fully_phased'] and \
                                                    not x['is_identical'] and \
                                                    not x['is_ibd_outlier'] and \
                                                    not x['is_crossover_outlier']]

if args.num_males is not None:
    sibpairs = [x for x in sibpairs if int(sample_to_sex[x['sibling1']]=='1')+int(sample_to_sex[x['sibling2']]=='1')==args.num_males]


family_chrom_to_cos = defaultdict(list)
for co in phase_data.get_crossovers():
    family_chrom_to_cos[(co['family'], co['chrom'])].append(co)

for gc in phase_data.get_crossovers():
    family_chrom_to_cos[(gc['family'], gc['chrom'])].append(gc)

for sibpair in sibpairs:
    sibpair['num_affected'] = int(sample_to_affected[sibpair['sibling1']]=='2') + int(sample_to_affected[sibpair['sibling2']]=='2')

if args.sibpair_type == 4:
    pass
elif args.sibpair_type == 3:
    sibpairs = [x for x in sibpairs if x['num_affected']>0]
else:
    sibpairs = [x for x in sibpairs if x['num_affected']==args.sibpair_type]
num_sibpairs = len(sibpairs)


print('Overall')
print('families', len(set([x['family'].split('.')[0] for x in sibpairs])))
print('sibpairs', len(sibpairs))
#print('num_affected', Counter([x['num_affected'] for x in sibpairs]))

def parse_interval(x):
    return x.split(':')[0], int(x.split(':')[1].split('-')[0]), int(x.split(':')[1].split('-')[1])

def process_phase_file(sibpair):
    inds = phase_data.get_phase_info(sibpair['family'])['individuals']
    chroms, starts, ends, mat_phases, pat_phases, is_htss = phase_data.parse_phase_file_into_arrays(sibpair['family'])
    sib1_ind_index, sib2_ind_index = inds.index(sibpair['sibling1']), inds.index(sibpair['sibling2'])
       
    mat_match = (mat_phases[sib1_ind_index, :]==mat_phases[sib2_ind_index, :]).astype(int)
    mat_match[mat_match==0] = -1
    mat_match[mat_phases[sib1_ind_index, :]==-1] = 0
    mat_match[mat_phases[sib2_ind_index, :]==-1] = 0
    #mat_match[is_htss] = 0
        
    pat_match = (pat_phases[sib1_ind_index, :]==pat_phases[sib2_ind_index, :]).astype(int)
    pat_match[pat_match==0] = -1
    pat_match[pat_phases[sib1_ind_index, :]==-1] = 0
    pat_match[pat_phases[sib2_ind_index, :]==-1] = 0
    #pat_match[is_htss] = 0
    
    # prune interval
    if (args.interval is not None):
        is_selected = np.zeros((len(chroms),), dtype=bool)

        chrom, start_pos, end_pos = parse_interval(args.interval)
        is_in_chrom = np.array([c==chrom for c in chroms], dtype=bool)
        is_selected_chrom = is_in_chrom & (np.minimum(ends, end_pos)-np.maximum(starts, start_pos)>0)
        is_selected[is_selected_chrom] = True
        starts[is_selected_chrom] = np.maximum(starts[is_selected_chrom], start_pos)
        ends[is_selected_chrom] = np.minimum(ends[is_selected_chrom], end_pos)

        chroms = [chroms[i] for i in np.where(is_selected)[0]]
        starts = starts[is_selected]
        ends = ends[is_selected]
        mat_match = mat_match[is_selected]
        pat_match = pat_match[is_selected]


    # we don't know what's going on inside crossovers
    for chrom in set(chroms):
        is_in_chrom = np.array([c==chrom for c in chroms])
        for co in family_chrom_to_cos[(sibpair['family'], chrom)]:
            if co['is_mat']:
                mat_match[is_in_chrom & (np.minimum(ends, co['end_pos'])-np.maximum(starts, co['start_pos']+1)>0)] = 0
            else:
                pat_match[is_in_chrom & (np.minimum(ends, co['end_pos'])-np.maximum(starts, co['start_pos']+1)>0)] = 0

        
    chrom_breaks = np.array([c1!=c2 for c1, c2 in zip(chroms[:-1], chroms[1:])])
    
    #print(len(chroms))
    #for p in zip(chroms, starts, ends, mat_match, pat_match):
    #    print(p)
                
    # collapse
    num_intervals = len(chroms)
    if num_intervals>1:
        change_index = np.where(chrom_breaks | (mat_match[:-1] != mat_match[1:]) | (pat_match[:-1] != pat_match[1:]))[0]
        chroms = [chroms[i] for i in np.hstack((change_index, num_intervals-1))]
        starts = starts[np.hstack(([0], change_index+1))]
        ends = ends[np.hstack((change_index, num_intervals-1))]
        mat_match = mat_match[np.hstack((change_index, num_intervals-1))]
        pat_match = pat_match[np.hstack((change_index, num_intervals-1))]
    
    #print(len(chroms))
    #for p in zip(chroms, starts, ends, mat_match, pat_match):
    #    print(p)
    #print(chroms, starts, ends)
    
    return chroms, starts, ends, mat_match, pat_match
       


# pull intervals
positions = set()
for sibpair in sibpairs:
    chroms, starts, ends, mat_match, pat_match = process_phase_file(sibpair)
    for chrom, start_pos, end_pos in zip(chroms, starts, ends):
        positions.add((chrom, start_pos))
        positions.add((chrom, end_pos))

positions = sorted(positions, key=lambda x: (int(x[0]), x[1]) if x[0].isdigit() else (23, x[1]))
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

with open('permutation_tests/%s.%d.intervals.json' % (dataset_name, args.sibpair_type), 'w+') as f:
    json.dump([{'index': int(i),
                'chrom': chroms[i],
                'start_pos': int(interval_starts[i]),
                'end_pos': int(interval_ends[i]),
                } for i in range(num_intervals)], 
        f, indent=4)

# pull sibpair IBD

# sibpair, interval
is_mat_match = np.zeros((num_sibpairs, num_intervals), dtype=np.int8)
is_pat_match = np.zeros((num_sibpairs, num_intervals), dtype=np.int8)

interval_start_to_index = dict([((chrom, x), i) for i, (chrom, x) in enumerate(zip(chroms, interval_starts))])
interval_end_to_index = dict([((chrom, x), i) for i, (chrom, x) in enumerate(zip(chroms, interval_ends))])

for sibpair_index, sibpair in enumerate(sibpairs):
    for chrom, start_pos, end_pos, mat_match, pat_match in zip(*process_phase_file(sibpair)):
        start_index, end_index = interval_start_to_index[(chrom, start_pos)], interval_end_to_index[(chrom, end_pos)]+1
        is_mat_match[sibpair_index, start_index:end_index] = mat_match
        is_pat_match[sibpair_index, start_index:end_index] = pat_match

print('mat match', np.sum(is_mat_match==1, axis=0))
print('mat dont match', np.sum(is_mat_match==-1, axis=0))
print('pat match', np.sum(is_pat_match==1, axis=0))
print('pat dont match', np.sum(is_pat_match==-1, axis=0))
print(is_mat_match.shape, is_pat_match.shape)

u = np.arange(num_intervals)

if (args.interval is not None):
    chrom, start_pos, end_pos = parse_interval(args.interval)
    is_in_chrom = np.array([c==chrom for c in chroms], dtype=bool)
    is_selected = is_in_chrom & (np.minimum(interval_ends, end_pos)-np.maximum(interval_starts, start_pos)>0)

    u = u[is_selected]

print(interval_starts[u])

num_intervals = len(u)
print(num_intervals)

np.save('permutation_tests/%s.%d.mat_match.npy' % (dataset_name, args.sibpair_type), np.sum(is_mat_match==1, axis=0))
np.save('permutation_tests/%s.%d.pat_match.npy' % (dataset_name, args.sibpair_type), np.sum(is_pat_match==1, axis=0))

np.save('permutation_tests/%s.%d.is_mat_match.npy' % (dataset_name, args.sibpair_type), is_mat_match)
np.save('permutation_tests/%s.%d.is_pat_match.npy' % (dataset_name, args.sibpair_type), is_pat_match)

# take into account sibling structure across quads
individuals = sorted(set([x['sibling1'] for x in sibpairs] + [x['sibling2'] for x in sibpairs]))
ind_to_index = dict([(x, i) for i, x in enumerate(individuals)])
sibling1_indices = np.array([ind_to_index[x['sibling1']] for x in sibpairs])
sibling2_indices = np.array([ind_to_index[x['sibling2']] for x in sibpairs])

A = np.random.randint(0, high=2, size=(args.num_trials+1, len(individuals), 2))
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

# now make X1 and X2 the IBD-phen relationships
is_discordant = np.array([sibpair['num_affected']==1 for sibpair in sibpairs])
X1[:, is_discordant] = -X1[:, is_discordant]
X2[:, is_discordant] = -X2[:, is_discordant]

print('ready')

# trial, interval
print('mat')
mat_rand_pvalue = X1.dot(is_mat_match)
print('pat')
pat_rand_pvalue = X2.dot(is_pat_match)


# -------------------- implementing Westfall-Young max T stepdown procedure

# indices are sorted along interval axis from interval with most IBD sharing
# to least IBD sharing
def stepdown(rand_pvalue):
    final_pvalues = np.zeros((num_intervals, ))
    
    orig_indices = np.flip(np.argsort(rand_pvalue[0, :]), axis=0)

    max_t_k = np.zeros((args.num_trials+1, num_intervals+1))
    max_t_k[:, -1] = np.min(rand_pvalue, axis=1)
    for i, j in list(reversed(list(enumerate(orig_indices)))):
        max_t_k[:, i] = np.maximum(max_t_k[:, i+1], rand_pvalue[:, j])
    max_t_k = max_t_k[:, :-1]

    #max_t_k = np.flip(np.sort(rand_pvalue[:, :, is_mat], axis=1), axis=1)
    
    assert np.all(max_t_k[0, :] == rand_pvalue[0, orig_indices])

    # calculate pi(j)
    pvalues = np.sum(max_t_k[1:, :] >= np.tile(max_t_k[0, :], (args.num_trials, 1)), axis=0)/args.num_trials
    pvalues = np.array([np.max(pvalues[:(i+1)]) for i in np.arange(pvalues.shape[0])])
    final_pvalues[orig_indices] = pvalues
    return final_pvalues

mat_final_pvalues = stepdown(mat_rand_pvalue)
pat_final_pvalues = stepdown(pat_rand_pvalue)

print(np.min(mat_final_pvalues), np.min(pat_final_pvalues))

np.save('permutation_tests/%s.%d.mat.npy' % (dataset_name, args.sibpair_type), mat_final_pvalues)
np.save('permutation_tests/%s.%d.pat.npy' % (dataset_name, args.sibpair_type), pat_final_pvalues)



