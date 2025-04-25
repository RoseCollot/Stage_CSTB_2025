#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=vega


for file in /home/collot/stage/collot/out_ortho/x* ; do 
    id=$(basename "$file") 
    sbatch --job-name=oto_hs_${id} --output=job_output_${id} /home/collot/stage/code/orthologs/orthologs_oto_proteins.sh $id
    done 
