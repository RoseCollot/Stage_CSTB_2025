#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry

/home/collot/.conda/envs/stage/bin/python3 /home/collot/stage/code/orthologs/write_fasta.py