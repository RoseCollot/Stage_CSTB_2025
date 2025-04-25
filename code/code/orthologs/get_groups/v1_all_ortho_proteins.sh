#!/bin/bash

my_gene="XP_039336287.1" 
group_number=-1
declare -A group_map
while read -r line ; do 
    if [[ "$line" == SUMMARY* ]]; then
        group_number=-1  
        group_map=()
        gene_groups=()
        echo "summary"
        continue
    fi
    if [[ "$line" == GROUP* ]]; then
        ((group_number++))
        group_map["$group_number"]=$(echo "$line" | cut -f2- | tr '\t' ',')
        if [[ "$line" == *"$my_gene"* ]]; then
            gene_groups+=("$group_number")
            echo "$group_map["$group_number"]"
            echo "group $group_number"
        fi
    fi
    if [[ "$line" == OTO* && "$line" == *"$my_gene"* ]]; then
        echo "my gene is here"
        col2=$(echo "$line" | awk '{print $2}')
        col3=$(echo "$line" | awk '{print $3}')
        if [[ "$col2" == *"$my_gene"* ]]; then
            ortho=$col3
            echo "$ortho"
            oto_list="${oto_list}${ortho},"
        elif [[ "$col3" == *"$my_gene"* ]]; then
            ortho=$col2
            echo "$ortho"
            oto_list="${oto_list}${ortho},"
        fi
        echo "oto $oto_list"
    fi
    if [[ "$line" == OTM* ]]; then
        for group in "${gene_groups[@]}"; do
            if [[ "$line" == *$'\t'"$group"* ]]; then
                echo "my group is in OTM"
                otm_list="$(echo "$line" | cut -f2 | tr '\t' ','),"
                echo "otm group $otm_list"
            fi
        done
        if [[ "$line" == *"$my_gene"* ]]; then 
            referenced_group=$(echo "$line" | awk '{print $3}')
            otm_list="${otm_list}${group_map["$referenced_group"]},"
            echo "otm gene $otm_list"
        fi 
    fi
    if [[ "$line" == MTM* ]]; then 
        for group in "${gene_groups[@]}"; do
            if [[ "$line" == *$'\t'"$group"* ]]; then
                echo "my group is in mtm"
                col2=$(echo "$line" | awk '{print $2}')
                col3=$(echo "$line" | awk '{print $3}')                
                if [[ "$col2" == *"$group"* ]]; then
                    ortho_group=$col3
                    mtm_list="${mtm_list}${group_map["$ortho_group"]},"
                    echo "col 3 $ortho_group"
                elif [[ "$col3" == *"$group"* ]]; then
                    ortho_group=$col2
                    mtm_list="${mtm_list}${group_map["$ortho_group"]},"
                    echo "col 2 $ortho_group"
                fi
                mtm_list="${mtm_list/$my_gene}"
                echo "mtm $mtm_list"
            fi
        done
    fi
done < <(zcat /data/kress/RefseqPrimates/res.tsv.gz)
row="${my_gene}\t${oto_list}${otm_list}${mtm_list}"
row="${row%,}" 
echo "$row" >> /home/collot/stage/collot/data/genes_orthologs.tsv
