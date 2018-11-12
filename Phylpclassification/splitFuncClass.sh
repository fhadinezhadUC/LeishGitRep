#!/bin/bash

array=(A R N D C Q E G H I L K M X F P Z S T W Y V)
echo "Array size: ${#array[*]}"

for item in ${array[*]}
do
    printf "   %s\n" $item
    cat EditedCovea.fasta | fasgrep "_${item}$" > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput/${item}.fasta"
    cat "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput/${item}.fasta" | fasconvert -o clustalw > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput/clustalW/TryTryp_${item}.aln"

done

