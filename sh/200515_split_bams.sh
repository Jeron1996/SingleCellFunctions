#!/bin/bash -e

id_table=$1
samtool_bam=$2
output_dir=$3
sample=$4

#Following three lines isolate the cluster number from the current text file
suffix="${id_table##*[0-9]}"
number="${id_table%"$suffix"}"
number="${number##*[!-0-9]}"
echo $id_table
echo $number
output=${output_dir}cluster${number}/${sample}clust${number}.bam
#Extracts header from original bam file and uses that as the start for the new bam file
samtools view -H ${samtool_bam} > $output
#Extracts all cell barcodes belonging to cluster X and pastes them below the header of the original bam file
samtools view ${samtool_bam} | grep -w -f $id_table >> $output
