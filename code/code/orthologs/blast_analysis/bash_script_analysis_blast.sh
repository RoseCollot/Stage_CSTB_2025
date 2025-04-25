#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry

/home/collot/.conda/envs/stage/bin/python3 /home/collot/stage_git/code/orthologs/script_blast_analysis.py stage/collot/collot/out_ortho/blast/out/${1}.out /home/collot/stage/collot/collot/out_ortho/fasta_homo_sapiens/${1}.fasta stage/collot/collot/out_ortho/all_hits_tsv/${1}.tsv