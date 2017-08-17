#!/bin/sh


#sampleIDlist => Need to be defined and must correspond to the names of the samples 
#otuIDlist => Need to be defined and must correspond to the names of isolates 
#16s_trim_v5v6_KG_DepExp_Derep.fa => is a fasta file, needs to be defined and must be trimmed

#Uncompress MiSeq files
gzip -d  001_forward_reads.fastq.gz
gzip -d  001_reverse_reads.fastq.gz
gzip -d  001_barcodes.fastq.gz

#Add barcodes Seq to forward and reverse reads 
$path_to_script/add_barcode.py 001_forward_reads.fastq 001_reverse_reads.fastq 001_barcodes.fastq &

#Assemble Paired ends 
pandaseq -N -o 80 -q 30 -F -d rbfkms -l 344 -L 389 -f 001_forward_reads_tagged.fastq -r 001_reverse_reads_tagged.fastq 1> assembled.fastq 2> stat.txt &

#Name assembled read by IdSample
$path_to_script/NameBySample.py assembled.fastq 002_mapping.txt  fasta -

#cat  ../Rawseq/001_mapping.txt | sed -e 's/\t.*//g' > idsamples.txt
#cat  list | sed -e 's/.fa//g' > idsamples
#grep ">" 16s_exp.fa | sed -e 's/ .*//g;s/>//g;'> outIDs

#Split the seqs fasta file by Idsample names for mapping
bsub -q multicore40 -R "rusage[mem=20000]" $path_to_scripts/sortByHeaderStr.py ../seqs_pandaseq.fasta ../idsamples.txt

#Aligment of reads to Reference 16S
java -Xmx40g -jar /home/hassani/RDPTools/AlignmentTools.jar pairwise-knn -k 1 -m global -t 3 -p 0 -o outputalig/$output \
$QUERY $path/16s_trim_v5v6_KG_DepExp_Derep.fa

#OTUtable 
$path_to_script/getSum.py sampleIDlist otuIDlist cutoff pairwiseKNNfiles > tab.txt


