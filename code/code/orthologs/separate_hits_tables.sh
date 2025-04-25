#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry

awk -F $'\t' -v specie_id="${1}" '$6 == specie_id' /home/collot/stage/collot/collot/out_ortho/hits_table.tsv > /home/collot/stage/collot/collot/out_ortho/species_tables/${1}.tsv