#!/bin/bash
#SBATCH --job-name=HS_orth
#SBATCH --output=HS_orth_output.out
#SBATCH --ntasks=1
#SBATCH --time=2-00:00:00
#SBATCH --partition=lab
#SBATCH --nodelist=vega

touch /home/collot/stage/collot/out_ortho/prot_oto_homo_sapiens${1}.tsv
IFS=$'\n' read -d '' -r -a proteins < /home/collot/stage/collot/out_ortho/${1}
for prot in "${proteins[@]}"; do
    oto_list=""
    while read -r col1 col2 col3 ; do 
            if [[ "$col2" == *"$prot"* ]]; then
                ortho=$col3
                oto_list="${oto_list}${ortho},"
            elif [[ "$col3" == *"$prot"* ]]; then
                ortho=$col2
                oto_list="${oto_list}${ortho},"
            fi
    done < <(zcat /data/kress/RefseqPrimates/res.tsv.gz | grep "OTO.*$prot" )
    row="${prot}\t${oto_list}"
    row="${row%,}" 
    echo -e "$row" >> /home/collot/stage/collot/out_ortho/prot_oto_homo_sapiens${1}.tsv
done
