#!/bin/bash

##this version we will set the arguments and delele the submission part

# vars
tissue=$1
threads=$2
clusters=$3
bam=$4
ref=$5
g_size=$6
qval=$7
targetclust=$8

#input=''

#perl -ne 'chomp;if($_=~/^cell/){next;}else{print"$_\n";}' $clusters > $clusters.new

grep -v '^cell' "$clusters" > "$clusters.new"

# make directory
if [ ! -d bigwigs ];then
	mkdir bigwigs
else
	rm -rf bigwigs
	mkdir bigwigs
fi

if [ ! -d bedgraphs ]; then
	mkdir bedgraphs
else
	rm -rf bedgraphs
	mkdir bedgraphs
fi

#################################################
# function to merge cells given cluster and BAM #
#################################################
mergeCells(){

	# vars
	clust=$1
	bamf=$2
	clustf=$3
	tissue=$4
	targetclust=$5

	# prep cluster data file
	#perl -ne 'chomp;if($_=~/^cellID/){next;}else{print"$_\n";}' $clustf > $clustf.$clust

    grep -v '^cellID' "$clustf" > "${clustf}.${clust}"


	##the cluster is in the 13th column since we add the singlet information
	##now we replace the 13 to the targetclust information
	cut -f1,$targetclust $clustf.$clust > $clustf.$clust.temp
	mv $clustf.$clust.temp $clustf.$clust

	# select cells from correct cluster
	awk -F'\t' -v cluster=$clust '$2==cluster' $clustf.$clust | cut -f1 - > cluster.$clust.bc_IDs.txt
	grep -Ff cluster.$clust.bc_IDs.txt $bamf | cut -f1-4 - > cluster_$clust.$tissue.pool.bed

	# make pseudoreps
	numcells=$( wc -l < cluster.$clust.bc_IDs.txt )
	repnum=$( echo "$numcells / 2" | bc -l )
	repnum=$(printf %.0f $repnum)
	echo " - splitting cluster into 2 reps each with $repnum cells ..."
	shuf cluster.$clust.bc_IDs.txt > cluster.$clust.bc_IDs.shuffled.txt
	head -n $repnum cluster.$clust.bc_IDs.shuffled.txt > cluster.$clust.bc_IDs.rep1.txt
	tail -n $repnum cluster.$clust.bc_IDs.shuffled.txt > cluster.$clust.bc_IDs.rep2.txt
	grep -Ff cluster.$clust.bc_IDs.rep1.txt cluster_$clust.$tissue.pool.bed > cluster_$clust.$tissue.rep1.bed
	grep -Ff cluster.$clust.bc_IDs.rep2.txt cluster_$clust.$tissue.pool.bed > cluster_$clust.$tissue.rep2.bed

	# clean temp files
	#rm cluster.$clust.bc_IDs.txt
	rm $clustf.$clust
	rm cluster.$clust.bc_IDs.shuffled.txt
	rm cluster.$clust.bc_IDs.rep1.txt
	rm cluster.$clust.bc_IDs.rep2.txt

}
export -f mergeCells


#############################
# iterate over all clusters #
#############################

# function
iterateClusters(){

	# load parameters
	i=$1
	bam=$2
	clusters=$3
	tissue=$4
	ref=$5
	#input=$6
	g_size=$6
	qval=$7
	targetclust=$8

	# make directory
	if [ ! -d $PWD/$tissue.$i ]; then
		mkdir $PWD/$tissue.$i
	fi

	# verbose
	echo "merging cells from cluster $i ..."
	mergeCells $i $bam $clusters $tissue $targetclust

	# call peaks
	echo "calling peaks with macs2 cluster $i ..."
	readdepth=$( wc -l < cluster_$i.$tissue.pool.bed)
	echo "total reads = $readdepth"
	macs2 callpeak -t cluster_$i.$tissue.pool.bed \
		-f BED \
        	-g $g_size \
                --nomodel \
                --keep-dup all \
                --extsize 150 \
                --shift -75 \
		--qvalue $qval \
                --outdir $PWD/$tissue.$i \
		--bdg \
                -n cluster.$i.macs2

	# un-corrected
	sort -k1,1 -k2,2n $PWD/$tissue.$i/cluster.$i.macs2_treat_pileup.bdg | python cleanBED.py $ref $readdepth - > cluster.$i.macs2_treat_pileup.clean.bdg
	bedGraphToBigWig cluster.$i.macs2_treat_pileup.clean.bdg $ref cluster_$i.$tissue.raw.bw
	mv cluster.$i.macs2_treat_pileup.clean.bdg cluster_$i.$tissue.raw.bdg

	# clean
#	rm $PWD/$tissue.$i/cluster.$i.macs2_treat_pileup.bdg
	rm $PWD/$tissue.$i/cluster.$i.macs2_control_lambda.bdg


	# call peaks rep1
        echo "calling peaks with macs2 cluster $i ..."
        readdepth=$( wc -l < cluster_$i.$tissue.rep1.bed)
        echo "total reads rep1 = $readdepth"
        macs2 callpeak -t cluster_$i.$tissue.rep1.bed \
                -f BED \
                -g $g_size \
                --nomodel \
                --keep-dup all \
                --extsize 150 \
                --shift -75 \
                --qvalue $qval \
                --outdir $PWD/$tissue.$i \
                -n cluster.$i.macs2.rep1

	# call peaks rep2
        echo "calling peaks with macs2 cluster $i ..."
        readdepth=$( wc -l < cluster_$i.$tissue.rep2.bed)
        echo "total reads rep2 = $readdepth"
        macs2 callpeak -t cluster_$i.$tissue.rep2.bed \
                -f BED \
                -g $g_size \
                --nomodel \
                --keep-dup all \
                --extsize 150 \
                --shift -75 \
                --qvalue $qval \
                --outdir $PWD/$tissue.$i \
                -n cluster.$i.macs2.rep2

	# check overlaps
	bedtools intersect -a $PWD/$tissue.$i/cluster.$i.macs2_peaks.narrowPeak -b $PWD/$tissue.$i/cluster.$i.macs2.rep1_peaks.narrowPeak -u \
		| bedtools intersect -a - -b $PWD/$tissue.$i/cluster.$i.macs2.rep2_peaks.narrowPeak -u > $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak

	# get summits for reproducible peaks
	bedtools intersect -a $PWD/$tissue.$i/cluster.$i.macs2_summits.bed -b $PWD/$tissue.$i/cluster.$i.reproducible.narrowPeak -u > $PWD/cluster.$i.reproducible_summits.bed

	# move bed files to cluster directory
	mv cluster_$i.$tissue.pool.bed cluster_$i.$tissue.rep1.bed cluster_$i.$tissue.rep2.bed $PWD/$tissue.$i
	mv cluster_$i.$tissue.raw.bdg $PWD/bedgraphs
	mv cluster_$i.$tissue.raw.bw $PWD/bigwigs
}
export -f iterateClusters

# iterate over clusters
parallel -j $threads iterateClusters {1} $bam $clusters.new $tissue $ref $g_size $qval $targetclust ::: $( cut -f$targetclust $clusters.new | sort -k1,1n | uniq )
#iterateClusters 1 $bam $clusters.new $tissue $ref $input ::: $( cut -f11 $clusters.new | sort -k1,1n | uniq )


# remove temp1
rm $clusters.new

# create merged set of peaks
bash ./adjustPeaks.sh $tissue
