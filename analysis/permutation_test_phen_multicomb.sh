#!/bin/bash                        
#                                                             
#                                                                   
#SBATCH --job-name=comb                                              
#SBATCH --output=logs/comb.out                                      
#SBATCH --error=logs/comb.err                                         
#SBATCH -p dpwall                                                    
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/SibpairLinkage     
#SBATCH -t 5:00:00                                                      
#SBATCH --mem=128G

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36


#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 4 9 23 29 31

#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 7 36

#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 1 30 18

#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 4 31 29

#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 4 23 29

#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 4 23 31

#srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 4 23 9

srun python3 analysis/permutation_test_phen_multicomb.py spark 10000 9 23
