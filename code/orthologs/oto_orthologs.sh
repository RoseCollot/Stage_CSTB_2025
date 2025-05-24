#creation of the output file 
touch "output.tsv"

#get the names of the proteins
IFS=$'\n' read -d '' -r -a proteins < "list of proteins"
for prot in "${proteins[@]}"; do    #loop through each protein
    oto_list=""
    while read -r col1 col2 col3 ; do 
            if [[ "$col2" == *"$prot"* ]]; then
                ortho=$col3
                oto_list="${oto_list}${ortho},"
            elif [[ "$col3" == *"$prot"* ]]; then
                ortho=$col2
                oto_list="${oto_list}${ortho},"
            fi
    done < <(zcat "OrthoInspector output" | grep "OTO.*$prot" )
    row="${prot}\t${oto_list}"
    row="${row%,}" 
    echo -e "$row" >> "output.tsv"
done