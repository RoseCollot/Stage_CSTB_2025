#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry



for file in /home/collot/stage/collot/PrimateData/ncbi_dataset/data/GCF* ; do 
    id=$(echo $file |rev |cut -d '/' -f1|rev) 
    sbatch --job-name=df_${id} --output=job_output_${id} /home/collot/stage/code/statistics/batch1.sh $id
done 






