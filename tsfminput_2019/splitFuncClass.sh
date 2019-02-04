#!/bin/bash

#functions we have from homo
#A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  X  Y 
# we need to go through all the files in this path and split the functional class

filenames=$(ls /home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/*.fasta)

array=(A R N D C Q E G H I L K M X F P S T W Y V)

for ref in $filenames; do
    ref2=${ref#"/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/"}
    ref3=${ref2%".fasta"}
    cd  /home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/
    mkdir "$ref3"
    printf " %s\n" $ref3
mkdir "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/$ref3/clustalW"
for item in ${array[*]}
do
    printf "   %s\n" $item
    cat "$ref" | fasgrep "_${item}$" > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/$ref3/${item}.fasta"
    cat "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/$ref3/${item}.fasta" | fasconvert -o clustalw > "/home/fatemeh/Leishmania_Aug2018/phyloclassification/tsfminput_2019/tRNADBClustered/$ref3/clustalW/${ref3}_${item}.aln"
  

done
done
