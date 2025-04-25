#!/bin/bash
#--job-name=dependancy
#--output=job_output_dependancy
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry


for file in /home/collot/stage/collot/collot/out_ortho/species_tables/* ; do 
    id=$(basename "$file" .tsv) 
    sbatch --job-name=table${id} --output=job_output_${id} /home/collot/stage/code/orthologs/event_tables/to_script.sh $id
done 