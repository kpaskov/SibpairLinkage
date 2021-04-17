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


ped_file = '../DATA/mssng/mssng.ped.quads.ped'
sample_file = '../DATA/mssng/genotypes/samples.json'
ad_file = '../SBSE/depth/mssng_chr10.npy'

with open(sample_file, 'r') as f:
	samples = json.load(f)
sample_to_index = dict([(x, i) for i, x in enumerate(samples)])

parents_to_children = defaultdict(list)
with open(ped_file, 'r') as f:
	for line in f:
		pieces = line.strip().split('\t')
		if pieces[1] in sample_to_index and pieces[2] in sample_to_index and pieces[3] in sample_to_index:
			parents_to_children[(pieces[0], pieces[3], pieces[2])].append(pieces[1])
parents_to_children = dict([(k, v) for k, v in parents_to_children.items() if len(v)==4])

print('families', len(parents_to_children))

ad = np.load(ad_file)

for (famkey, mom, dad), inds in parents_to_children.items():
	indices = [sample_to_index[x] for x in inds]

	np.save('../SBSE/depth/%s_chr10' % famkey, ad[indices, :, :])
