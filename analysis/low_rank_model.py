import numpy as np
import matplotlib.pyplot as plt
import csv
import json
import cvxpy as cp
import sys

dataset = sys.argv[1]
k = int(sys.argv[2]) # 5
lambda_X = float(sys.argv[3]) # 0.05
lambda_Y_G = float(sys.argv[4]) # 0.01
lambda_Y_P = float(sys.argv[5]) # 0.0

print('k', k, 'lambda_X', lambda_X, 'lambda_Y_G', lambda_Y_G, 'lambda_Y_P', lambda_Y_P)

with open('../PhasingFamilies/recomb_%s/sibpairs.json' % dataset.split('.')[0], 'r') as f:
    sibpairs = json.load(f)
print('sibpairs', len(sibpairs))

mat_ibd = np.load('permutation_tests/phen.%s.mat_ibd.npy' % dataset)
pat_ibd = np.load('permutation_tests/phen.%s.pat_ibd.npy' % dataset)
G_orig = mat_ibd + pat_ibd
G_orig[G_orig==0] = 1
G_orig[G_orig==-2] = 0
G_orig[(mat_ibd==0) | (pat_ibd==0)] = -1

sample_to_phen = dict()
with open('../PhasingFamilies/phenotypes/spark_v5/spark_v5-scq-prep.csv', 'r') as f:
    reader = csv.reader(f)
    for pieces in reader:
        phen = pieces[13:53]
        sample_to_phen[pieces[2]] = np.array([1 if x=='1.0' else -1 if x=='0.0' else 0 for x in phen])
            
print('phen', len(sample_to_phen))

P_orig = np.zeros((len(sibpairs), 40))
for i, sibpair in enumerate(sibpairs):
    if sibpair['sibling1'] in sample_to_phen and sibpair['sibling2'] in sample_to_phen:
        phen1 = sample_to_phen[sibpair['sibling1']]
        phen2 = sample_to_phen[sibpair['sibling2']]
        is_not_missing = (phen1 != 0) & (phen2 != 0)
        P_orig[i, is_not_missing] = np.sign(np.multiply(phen1[is_not_missing], phen2[is_not_missing]))


no_pheno = np.all(P_orig==0, axis=1)


G = G_orig[~no_pheno, :]#[:500, :]
P = P_orig[~no_pheno, :]#[:500, :]
print('G', G.shape, 'P', P.shape)

num_sibpairs = G.shape[0] # num sibpairs
num_features_G = G.shape[1] # num features
num_features_P = P.shape[1] # num features
print(num_sibpairs, num_features_G, num_features_P)

sibpair_is_masked_G = np.random.random((num_sibpairs,))<0.05
is_masked_G = np.zeros((num_sibpairs, num_features_G), dtype=bool)
is_masked_G[sibpair_is_masked_G, :] = True
is_masked_G = is_masked_G & (G!=-1)

sibpair_is_masked_P = (np.random.random((num_sibpairs,))<0.05) & ~sibpair_is_masked_G
is_masked_P = np.zeros((num_sibpairs, num_features_P), dtype=bool)
is_masked_P[sibpair_is_masked_P, :] = True
is_masked_P = is_masked_P & (P!=0)

is_missing_G = (G==-1) | is_masked_G
is_missing_P = (P==0) | is_masked_P
print('G', np.sum(is_masked_G)/(num_sibpairs*num_features_G), 
      np.sum(G==-1)/(num_sibpairs*num_features_G))
print('P', np.sum(is_masked_P)/(num_sibpairs*num_features_P), 
      np.sum(P==0)/(num_sibpairs*num_features_P))
np.save('low_rank_models/is_masked_G.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), is_masked_G)
np.save('low_rank_models/is_masked_P.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), is_masked_P)
np.save('low_rank_models/is_missing_G.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), is_missing_G)
np.save('low_rank_models/is_missing_P.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), is_missing_P)


P[P==-1] = 0


def fitX(Y_G, Y_P, k):
    X = cp.Variable((num_sibpairs, k))

    G_est = cp.hstack((np.ones((num_sibpairs, 1)), X)) @ Y_G
    P_est = cp.hstack((np.ones((num_sibpairs, 1)), X)) @ Y_P
    log_likelihood = cp.sum(
        (cp.multiply(G, G_est) - 2*cp.logistic(G_est))[~is_missing_G]
    ) + \
    cp.sum(
        (cp.multiply(P, P_est) - cp.logistic(P_est))[~is_missing_P]
    )
    
    problem = cp.Problem(cp.Maximize(log_likelihood - lambda_X * cp.pnorm(X, 1)),
                        [])
    
    result = problem.solve(solver='MOSEK', verbose=True)
    return np.hstack((np.ones((num_sibpairs, 1)), X.value))

def fitY(X, k):
    Y_G = cp.Variable((k+1, num_features_G))
    Y_P = cp.Variable((k+1, num_features_P))

    G_est = X @ Y_G
    P_est = X @ Y_P
    log_likelihood = cp.sum(
        (cp.multiply(G, G_est) - 2*cp.logistic(G_est))[~is_missing_G]
    ) + \
    cp.sum(
        (cp.multiply(P, P_est) - cp.logistic(P_est))[~is_missing_P]
    )
    
    problem = cp.Problem(cp.Maximize(log_likelihood - lambda_Y_G * cp.pnorm(Y_G[1:, :], 1) - lambda_Y_P * cp.pnorm(Y_P[1:, :], 1)),
                        [Y_G[1:, :]>=0, Y_P[1:, :]>=0])
    
    result = problem.solve(solver='MOSEK')#, verbose=True)
    return Y_G.value, Y_P.value

X = np.hstack((np.ones((num_sibpairs, 1)), np.random.randn(num_sibpairs, k)))

for i in range(3):
    print('fitY')
    Y_G, Y_P = fitY(X, 5)
    print('fitX')
    X = fitX(Y_G, Y_P, 5)

    np.save('low_rank_models/X.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), X)
    np.save('low_rank_models/Y_G.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), Y_G)
    np.save('low_rank_models/Y_P.%s.%d.%0.2f.%0.2f.%0.2f' % (dataset, k, lambda_X, lambda_Y_G, lambda_Y_P), Y_P)



    