#!/bin/bash                        
#                                                             
#                                                                   
#SBATCH --job-name=ibd                                              
#SBATCH --output=logs/ibd.out                                      
#SBATCH --error=logs/ibd.err                                         
#SBATCH -p dpwall                                                    
#SBATCH -D /oak/stanford/groups/dpwall/users/kpaskov/SibpairLinkage     
#SBATCH -t 5:00:00                                                      
#SBATCH --mem=128G                                                     

module load py-numpy/1.14.3_py36
module load py-scipy/1.1.0_py36


#srun python3 analysis/pull_ibd.py spark
srun python3 analysis/pull_ibd.py ssc.hg38
#srun python3 analysis/pull_ibd.py ihart.ms2
