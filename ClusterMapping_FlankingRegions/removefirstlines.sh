#!/bin/bash

yourfilenames=$(ls ./GenomesFlankingRegions/*.txt)
for FILE in $yourfilenames; do
        #FILE1=${FILE#"./GenomesFlankingRegions/"}
        sed -i '1d' "${FILE}"
done

