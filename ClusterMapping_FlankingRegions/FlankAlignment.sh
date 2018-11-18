#!/bin/bash
# remove the first line of the file flankingpositions.txt
indexqryfile=$1
Leftoutputfile=$2 
Rightoutputfile=$3
ref=$4
Flankingref=$5
outputfolder=$6
#Alignment () {
while read line 
do
start=0
end=0
seqname=$(echo "$line" | awk -F  "," '{print $4}')
start=$(($(echo "$line" | awk -F  "," '{print $5}')+0))
end=$(($(echo "$line" | awk -F  "," '{print $6}')+0))
name=$(echo "$seqname" | tr -d '"')

samtools faidx "$ref" $name:$start-$end > "./${outputfolder}/flankingseq.fasta"
bowtie2 -x "./ReferenceIndex/${indexqryfile}" -f "./${outputfolder}/flankingseq.fasta" -S "./${outputfolder}/aligned.sam"
samtools view -S -b "./${outputfolder}/aligned.sam" > "./${outputfolder}/aligned.bam"
samtools view "./${outputfolder}/aligned.bam" > "./${outputfolder}/aligned.txt"
leftcoordinate=$(awk '{print $3,$4}' < "./${outputfolder}/aligned.txt")
echo "left"
echo "$leftcoordinate"
leftcoordinatearr[$i]="$leftcoordinate"

# do the same thing for the right flanking region 
start=0
end=0
seqname=$(echo "$line" | awk -F  "," '{print $4}')
start=$(($(echo "$line" | awk -F  "," '{print $7}')+0))
end=$(($(echo "$line" | awk -F  "," '{print $8}')+0))
name=$(echo "$seqname" | tr -d '"')

samtools faidx "$ref" $name:$start-$end > "./${outputfolder}/flankingseq.fasta"
bowtie2 -x "./ReferenceIndex/${indexqryfile}" -f "./${outputfolder}/flankingseq.fasta" -S "./${outputfolder}/aligned.sam"
samtools view -S -b "./${outputfolder}/aligned.sam" > "./${outputfolder}/aligned.bam"
samtools view "./${outputfolder}/aligned.bam" > "./${outputfolder}/aligned.txt"
rightcoordinate=$(awk '{print $3,$4}' < "./${outputfolder}/aligned.txt")
echo 	"right"
echo "$rightcoordinate"
rightcoordinatearr[$i]="$rightcoordinate"        
i=$((i+1))
done < "$Flankingref"

printf "%s\n" "${rightcoordinatearr[@]}" >  "./${outputfolder}/$Rightoutputfile"
printf "%s\n" "${leftcoordinatearr[@]}" > "./${outputfolder}/$Leftoutputfile"

#}


