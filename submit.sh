#!/bin/bash -l

#SBATCH -J E1-O3-Re
#SBATCH -o nohup.out
#SBATCH -p epyc
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=64gb
#SBATCH --time=168:00:00


./modRS ../SIDM60-E1-O3/output snapshot 000 
