#!/bin/bash

for file in /home/collot/stage/collot/collot/out_ortho/fasta_homo_sapiens/* ; do 
    id=$(basename "$file" .fasta) 
    echo "$id" 
    sbatch --job-name=${id} --output=stage/collot/collot/out_ortho/out_job_blast_analysis/${id}.out /home/collot/stage/code/orthologs/blast_analysis/bash_script_analysis_blast.sh ${id}
done 

#while IFS= read -r file ; do 
#    id=$(echo "$file") 
#    echo "$id" 
#    sbatch --job-name=blast_analysis_${id} --output=stage/collot/collot/out_ortho/out_job_blast_analysis/${id}.out /home/collot/stage/code/orthologs/bash_script_analysis_blast.sh ${id}
#done < errors.txt