#!/bin/bash                        
#                                                             
#                                                                   
#SBATCH --job-name=phen                                              
#SBATCH --output=logs/phen%a.out                                      
#SBATCH --error=logs/phen%a.err                                         
#SBATCH -p dpwall                                                    
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/SibpairLinkage     
#SBATCH -t 5:00:00                                                      
#SBATCH --mem=128G
#SBATCH --array=0-0

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36


srun python3 analysis/permutation_test_phen_comb.py spark 1000 $SLURM_ARRAY_TASK_ID

#srun python3 analysis/permutation_test_phen_comb.py spark 10000 7 36
