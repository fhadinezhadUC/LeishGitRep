#!/bin/bash

array=(A R N D C Q E G H I L K M X F P Z S T W Y V)
echo "Array size: ${#array[*]}"

for item in ${array[*]}
do
    printf "   %s\n" $item
    cat HomoTryTryp_EditedCovea.fasta | fasgrep "_${item}$" > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfmHomoCinput/TryTryp/${item}.fasta"
    # for HOMOC genes
    cat HomoTryTryp_EditedCovea.fasta | fasgrep "HOMO${item}" > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfmHomoCinput/HomoC/Homo${item}.fasta"

    cat "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfmHomoCinput/TryTryp/${item}.fasta" | fasconvert -o clustalw > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfmHomoCinput/TryTryp/clustalW/TryTryp_${item}.aln"
    cat "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfmHomoCinput/HomoC/Homo${item}.fasta" | fasconvert -o clustalw > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfmHomoCinput/HomoC/clustalW/Homo_${item}.aln"

done

