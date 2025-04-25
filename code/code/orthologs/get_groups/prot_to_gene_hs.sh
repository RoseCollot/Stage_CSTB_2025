#!/bin/bash
#SBATCH --job-name=gene_name_orth
#SBATCH --output=gene_name_orth_output.out
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=vega

/home/collot/.conda/envs/stage/bin/python3 /home/collot/stage/code/orthologs/prot_to_gene_hs.py