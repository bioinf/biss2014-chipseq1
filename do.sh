#!/usr/bin/bash

INFILE="data/K562_IFNg30_IRF1_rep1.fastq"
BASE=`basename "$INFILE" .fastq`
mkdir work

################## Trimming #################

# quality control
fastqc $INFILE

# trimming by quality
fastq_quality_trimmer -v -t 28 -l 36 -i $INFILE -o work/$BASE.qtrimmed.fastq

# quality control
fastqc work/$BASE.qtrimmed.fastq

# trimming by length
fastx_trimmer -v -f 1 -l 75 -i $INFILE -o work/$BASE.ltrimmed.fastq

# quality control
fastqc work/$BASE.ltrimmed.fastq


################## Mapping ##################

# building reference
(
    cd data; 
    mkdir hg19;
    cd hg19; 
    tar -xzf ../chromFa.tar.gz
    cat `ls *.fa | grep -P "chr(\d{1,2}|X|Y).fa"` > hg19.fasta
    bowtie2-build hg19.fasta hg19
)

# mapping

bowtie2 -p 4 -x data/hg19/hg19 -U $INFILE -S work/$BASE.sam
#bowtie2 -p 4 -x data/hg19/hg19 -U work/$BASE.qtrimmed.fastq -S wor/$BASE.ltrimmed.sam

samtools view -bS work/$BASE.sam -o work/$BASE.bam
samtools sort work/$BASE.bam work/$BASE.sorted
samtools index work/$BASE.sorted.bam

# :ToDo: what -F 4 means?
samtools view -b -F 4 -q 1 work/$BASE.sorted.bam -o work/$BASE.filtered.bam
samtools index work/$BASE.filtered.bam

macs2 -t work/$BASE.filtered.bam -g 1e9 --name=work/IRF1 --format=BAM --tsize=36 --mfold=10,30



