#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry


for file in /home/collot/stage/collot/collot/out_ortho/one_to_one_events/species_tables/* ; do 
    id=$(basename "$file" .tsv) 
    sbatch --job-name=one_to_one${id} --output=job_output_${id} /home/collot/stage/code/orthologs/get_one_to_one.sh ${id}
done






