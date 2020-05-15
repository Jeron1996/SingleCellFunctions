#!/bin/bash -e

source /etc/profile.d/modules.sh
export PATH=/home/yolsan/python-2.7.9-yolsan/bin:$PATH
export MODULEPATH=/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib:$MODULEPATH
module load gi/samtools/1.2
module load yolsan/python/2.7.9
source /home/yolsan/python-2.7.9-yolsan/bin/activate

#Each sample should have its own directory containing the sequencing bam file (.bam), a txt file with the cluster information (.txt) and txt files for every cluster ID with the corresponding cell barcodes (.txt; filename should start with "ids_" and end with "clust*.txt")
sampleDir="/share/ScratchGeneral/jerven/ATAC_seq/Final_test/bam_files/"
samples=$(ls ${sampleDir})

#Script needed to start the parallel to split bam files for each sample
script1="/share/ScratchGeneral/jerven/ATAC_seq/Final_test/200515_split_samples.sh"

#Script needed to start the second parallel, to merge all bam files off one cluster and make bigwig/BedGraph files
script2="/share/ScratchGeneral/jerven/ATAC_seq/Final_test/200515_merge_deeptool.sh"
#Make a directory for all the split bam files. Change "n" and "N" according to the lowest and highest cluster number. Also change "j" and "J" at the end of the script
#This was written for clusters IDs which are distringuished by numbers. I don't know if it also works if cluster IDs are words
output_dir=${sampleDir}../split_by_cluster/
mkdir ${output_dir}
n=0
N=11
while [ "$n" -le "$N" ]; do
  mkdir ${output_dir}cluster$n
  n=$(($n + 1))
done

#parallel run two script. test1.sh will start a second parallel script that distributes cell barcodes to their cluster ID. This should all run in parallel, so it may be a lot quicker than doing it in loops
ls ${sampleDir} | parallel --gnu -j0 -k bash ${script1} {} $sampleDir $output_dir

#parallel run for every cluster --> merge indiviudal bam files, index them + generate bigwig file and bedgraph file of every merged bam file
ls ${output_dir} | parallel --gnu -j0 -k bash ${script2} {} $output_dir
