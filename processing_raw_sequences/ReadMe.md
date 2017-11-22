#customized scripts used for generating CTU count table in Hassani et al., 
#M. Amine Hassani - ahassani@bot.uni-kiel.de
#uncompress fastq files, the files can be downloaded from >>> need to specify permanent link at mpipz <<< 
gzip -d  001_forward_reads.fastq.gz
gzip -d  001_reverse_reads.fastq.gz
gzip -d  001_barcodes.fastq.gz

#add barcode sequences to forward and reverse reads (custom script - Benli Chai)
$path/add_barcode.py 001_forward_reads.fastq 001_reverse_reads.fastq 001_barcodes.fastq

#reads assembly - pre-requist pandaseq (2.8v) 
pandaseq -N -o 80 -q 30 -F -d rbfkms -l 344 -L 389 -f 001_forward_reads_tagged.fastq -r 001_reverse_reads_tagged.fastq 1> assembled.fastq 2> stat.txt

#change read Ids by sample Ids (custom script - Benli Chai)
$path/NameBySample.py assembled.fastq 001_mapping.txt  fasta -

#split reads by sample Ids (custom script - Benli Chai)
$path/sortByHeaderStr.py assembled.fastqSampleName.fasta id_samples.txt

#mapping reads to 16s rRNA reference sequences for full community - prerequisites RDPTools
java -Xmx40g -jar /home/hassani/RDPTools/AlignmentTools.jar pairwise-knn -k 1 -m global -t 3 -p 0 -o [$path/output_$QUERY.al] [$QUERY] 16s_trim_v5v6_KG_DepExp_Derep.fa

#mapping reads to 16s rRNA reference sequences for -13HS community - prerequisites RDPTools
java -Xmx40g -jar /home/hassani/RDPTools/AlignmentTools.jar pairwise-knn -k 1 -m global -t 3 -p 0 -o [$path/output_$QUERY.al] [$QUERY] 16s_trim_v5v6_KG_DepExp_Derep_CC.fa

#mapping reads to 16s rRNA reference sequences for -13HS community - prerequisites
java -Xmx40g -jar /home/hassani/RDPTools/AlignmentTools.jar pairwise-knn -k 1 -m global -t 3 -p 0 -o [$path/output_$QUERY.al] [$QUERY] 16s_trim_v5v6_KG_DepExp_Derep_CS.fa

#generating OTU table (custom script - Benli Chai) - used cutoff 1
$path/getSum.py sampleidlist ctuidlist cutoff output_$QUERY.al > count_table.txt
