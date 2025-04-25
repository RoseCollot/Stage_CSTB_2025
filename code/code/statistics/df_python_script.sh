#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry

/home/collot/.conda/envs/stage/bin/python3  /home/collot/stage_git/code/script_df.py /home/collot/stage/collot/PrimateData/ncbi_dataset/data/${1}/genomic.gff /home/collot/stage/collot/collot/out_stats/databases/db_${1} /home/collot/stage/collot/collot/out_stats/output_dataframes/${1}/