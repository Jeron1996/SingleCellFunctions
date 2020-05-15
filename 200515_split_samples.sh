sample=$1
sampleDir=$2
output_dir=$3
echo ${sample}
echo $sampleDir
echo $output_dir

bams=$sampleDir$sample/
echo $bams
bam=$(ls -d $bams*.bam)
echo $bam

#parallel run the script that splits all cells in a bam file accoriding to their cluster identity
ls ${bams}ids* | parallel --gnu -j0 -k bash /share/ScratchGeneral/jerven/ATAC_seq/Final_test/200515_split_bams.sh {} $bam ${output_dir} $sample
