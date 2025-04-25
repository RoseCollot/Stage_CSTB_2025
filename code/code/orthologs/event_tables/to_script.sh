#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry

/home/collot/.conda/envs/stage/bin/python3 /home/collot/stage_git/code/orthologs/events_tables.py /home/collot/stage/collot/collot/out_ortho/species_tables/${1}.tsv /home/collot/stage/collot/collot/out_ortho/species_df_gene/df${1}