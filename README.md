# SibpairLinkage

## Purpose
This project contains code for conducting the Genome-wide sibling-pair linkage test, as published in

Paskov K, Chrisman B, Stockham N, McNealis M, Dunlap K, Jung JY, Wall DP. Genome-wide sibling-pair linkage test reveals parent-of-origin-specific risk regions for autism spectrum disorder. Under review.

## Input and output
This code starts with sibling IBD data calculated by the https://github.com/kpaskov/PhasingFamilies project. It produces parent of origin pvalues for every genomic region in the dataset, defined by crossover endpoints in all sibling pairs available in

```
permutation_tests/[dataset_name].[sibpair_type].intervals.json
permutation_tests/[dataset_name].[sibpair_type].mat.npy
permutation_tests/[dataset_name].[sibpair_type].pat.npy
```

The `[dataset_name].[sibpair_type].intervals.json` contains the genomic regions included in the test in .json format. Each region is defined by a chromosome and start and end positions.

The `[dataset_name].[sibpair_type].mat.npy` file contains adjusted pvalues for maternal IBD sharing.

The `[dataset_name].[sibpair_type].mat.npy` file contains adjusted pvalues for paternal IBD sharing.

## Instructions for running code

### 1. Start by calculating sibling pair IBD. 
using the PhasingFamilies model available at https://github.com/kpaskov/PhasingFamilies

### 2. Run the Genome-wide sibling-pair linkage test.
The algorithm calculates parent-of-origin specific pvalues throughout the genome. 

```
python analysis/permutation_test.py [dataset_name]  [data_dir] [ped_file] [sibpair_type]
```

The required input arguments are
- `dataset_name` a unique name for the dataset, will be used in naming the output files.
- `data_dir` the data directory containing IBD information, produced during Step 1.
- `ped_file` the .ped file describing family relationships
- `sibpair_type` a numerical value indicating the affected/unaffected status of sibling pairs to be analyzed. A value of 2 indicates that the analysis should be run on affected-affected sibling pairs (the most common case). A value of 1 indicates that the analysis should be run on affected-unaffected sibling pairs. A value of 0 indicates that the analysis should be run on unaffected-unaffected sibling pairs.

The script has options
- `--num_trials` Number of permutations to run in order to calculate pvalue. Default value is 1000.
- `--interval` Restrict the analysis to a genomic interval in format chr1:1000000-2000000
- `--num_males [num_males]` Restrict the analysis by sibpair sex. A value of 0 indicates that the analysis should be run only on female-female sibling pairs. A value of 1 indicates that the analysis should be run only on male-female sibling pairs. A value of 2 indicates that the analysis should be run only on male-male sibling pairs.

The example below runs the analysis on all affected-affected sibling pairs in a dataset named ancestry.

```
python3 analysis/permutation_test.py ancestry ../DATA/ancestry ../DATA/ancestry/ancestry.ped.quads.ped 2
 ```
 

