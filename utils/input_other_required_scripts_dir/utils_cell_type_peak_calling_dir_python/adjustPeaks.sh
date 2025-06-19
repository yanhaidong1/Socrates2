#!/bin/bash

# vars
name=$1
peaks=$( ls *summits.bed )

# iterate
for i in ${peaks[@]}; do

	echo "$i"

	id=$( echo "$i" | cut -d'.' -f1-2 )	

	python normalize_score.py $i \
	| python expand_and_annotate_peaks.py > "$id.temp"

done

# merge peaks
cat *.temp | sort -k1,1 -k2,2n - > $1.temp2

# actually merge
bedtools merge -i $1.temp2 -c 5 -o collapse | python selectNonOverlapping.py - > $1.unique500bpPeaks.bed
rm *.temp *.temp2
