#!/bin/bash
#Script for dowloading the fasta and gff file and extracting the "*_Genome.fasta" files
if [[ $# -eq 0 ]]
	then 
	path=(.)
	else
	path="$1"
fi

wget -A fasta,gff -m -p -e robots=off -E -k -K -np -P "$path" http://tritrypdb.org/common/downloads/Current_Release/ 
# the output will be in directory given as the first argument of this script
# there are pages that need agreement so wont be downloaded with the above link alone. 
# checking each directory, if there was a fasta directory that was empty, take the parent folder and make the new url to download the fasta file 
for i in $(find "$path" -name "fasta" -type d); do
 if [ -z "$(ls -A $i)" ]; then 
 pardir=$(dirname "$i")
 parentfolder=$(echo "$pardir" | awk 'BEGIN { FS="/" } {print $NF}')
 mainpath="http://tritrypdb.org/common/downloads/Current_Release"
 mypath="${mainpath}/${parentfolder}/fasta/data/"
 wget -A fasta,gff -m -p -e robots=off -E -k -K -np -P "$path" "$mypath"
 echo "$mypath"
 fi
done
# after downloading the files, extract all the file with the pattern "*_Genome.fasta" and redirect them to the folder "$path"
find "$path" -name "*_Genome.fasta" -exec cp {} "$path" \;
#find /home/fatemeh/Leishmania_Aug2018/tritrypdb.org/ -name "*_Genome.fasta" -exec cp {} /home/fatemeh/Leishmania_Aug2018/tritrypdb.org/ \;

