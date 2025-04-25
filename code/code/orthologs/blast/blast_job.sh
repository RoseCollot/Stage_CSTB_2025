#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --output=blast.out
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=larry

for file in /home/collot/stage/collot/collot/out_ortho/fasta_homo_sapiens/* ; do 
    id=$(basename "$file" .fasta) 
    echo "$id" 
    makeblastdb -in /home/collot/stage/collot/collot/out_ortho/fasta_homo_sapiens/${id}.fasta -dbtype nucl -out /home/collot/stage/collot/collot/out_ortho/blast/db/blastdb_${id}
    blastn -db /home/collot/stage/collot/collot/out_ortho/blast/db/blastdb_${id} -query /home/collot/stage/collot/collot/out_ortho/fasta_homo_sapiens/${id}.fasta -out /home/collot/stage/collot/collot/out_ortho/blast/out/${id}.out  -outfmt 6 
done 