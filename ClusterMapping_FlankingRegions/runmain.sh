#!/bin/bash
reffilenames=$(ls ./ReferenceIndex/*.fasta)
yourfilenames=$(ls ./GenomesFlankingRegions/*.txt)

for ref in $reffilenames; do

ref2=${ref#"./ReferenceIndex/TriTrypDB-38_"}
ref3=${ref2%"_Genome.fasta"}
Flankingref="./GenomesFlankingRegions/FlankingR_${ref3}.txt"
outputfolder="Alignment_${ref3}"
mkdir "${outputfolder}"
for FILE in $yourfilenames; do
        FILE2=${FILE#"./GenomesFlankingRegions/FlankingR_"}
	FILE3=${FILE2%".txt"}
        indexqryfile="TriTrypDB-38_${FILE3}_Genome"
        Leftoutputfile="${ref3}_LeftF_AlTo_${FILE3}.txt"
        Rightoutputfile="${ref3}_RightF_AlTo_${FILE3}.txt"
        FlankAlignment.sh "$indexqryfile"  "$Leftoutputfile" "$Rightoutputfile" "$ref" "$Flankingref" "$outputfolder" 
	 
done

done
