#!/bin/bash

while IFS= read -r file ; do 
    id=$(echo "$file")
    sbatch --job-name=table${id} --output=job_output_${id} /home/collot/stage/code/orthologs/separate_hits_tables.sh ${id}
done < stage/collot/collot/out_stats/tax_id.txt 
