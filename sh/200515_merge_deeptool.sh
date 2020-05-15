#!/bin/bash -e

cluster=$1
outputDir=$2

out=${outputDir}${cluster}/merged_${cluster}.bam
bams="$(ls -d ${outputDir}${cluster}/*.bam )"
echo ${out}
echo ${bams}

#merge bam files and index
samtools merge ${out} ${bams}
samtools index ${out}

#Create result folder for bigwig files and convert bam to bigwig
RESULTS_FOLDER=${outputDir}${cluster}/bigwig/
mkdir -p ${RESULTS_FOLDER}
echo ${RESULTS_FOLDER}
CMD="bamCoverage -b ${out} -o ${RESULTS_FOLDER}merge_${cluster}.bw \
--extendReads 300 \
--normalizeUsing CPM \
--ignoreDuplicates \
--ignoreForNormalization chrY chrM\
--numberOfProcessors max "
eval $CMD

#Create result folder for BedGraph file and convert bam to BedGraph
results_bedgraph="${outputDir}${cluster}/BedGraph/"
mkdir -p ${results_bedgraph}
# generate bedgraph file for each sample
# Normalize using Counts per million
CMD_bed="bamCoverage -b ${out} -o ${results_bedgraph}merge_${cluster}.bdg \
--binSize 1
--outFileFormat "bedgraph" \
--normalizeUsing CPM \
--ignoreForNormalization chrY chrM"
eval $CMD_bed
