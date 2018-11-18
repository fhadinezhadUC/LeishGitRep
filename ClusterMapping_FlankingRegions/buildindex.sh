#!/bin/bash

yourfilenames=$(ls *.fasta)
for FILE in $yourfilenames; do
	FILE2=${FILE%".fasta"}
bowtie2-build "${FILE}" "${FILE2}"
done

