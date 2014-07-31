#!/usr/bin/bash

INFILE="data/K562_IFNg30_IRF1_rep1.fastq"
BASE=`basename "$INFILE" .fastq`
mkdir work

################## Trimming #################

# quality control
fastqc $INFILE

# trimming by quality
fastq_quality_trimmer -Q 64 -v -t 28 -l 30 -i $INFILE -o work/$BASE.qtrimmed.fastq |& ts | tee work/$BASE.qtrimmed.log

# quality control
fastqc work/$BASE.qtrimmed.fastq

# trimming by length
#fastx_trimmer -v -f 1 -l 75 -i $INFILE -o work/$BASE.ltrimmed.fastq

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

bowtie2 -p 4 -x data/hg19/hg19 -U $INFILE -S work/$BASE.sam |& ts | tee work/$BASE.sam.log
bowtie2 -p 4 -x data/hg19/hg19 -U work/$BASE.qtrimmed.fastq -S work/$BASE.qtrimmed.sam |& ts | tee work/$BASE.qtrimmed.sam.log
bowtie2 -p 4 -t -k 1 -x data/hg19/hg19 -U work/$BASE.qtrimmed.fastq -S work/$BASE.qtrimmed.nomm.sam |& ts | tee work/$BASE.qtrimmed.nomm.sam.log
bowtie2 -p 4 -t -k 1 -x data/hg19/hg19 -U work/$BASE.qtrimmed.fastq -S work/$BASE.qtrimmed.nomm.sam |& ts | tee work/$BASE.qtrimmed.nomm.sam.log

function sam_sort {
    SAM=$1
    BAM=${1/sam/bam}
    SORTED=${1/sam/sorted}
    samtools view -bS $1 -o $BAM
    samtools sort $BAM $SORTED
    samtools index $SORTED.bam
}

sam_sort work/$BASE.sam
sam_sort work/$BASE.qtrimmed.sam

function sam_filter {
    BAM=$1
    OPTIONS="$2"
    FILTERED=${BAM/bam/filtered.bam}

    samtools view -b $OPTIONS $BAM -o $FILTERED
    samtools index $FILTERED
}

sam_filter work/$BASE.sorted.bam "-F 4 -q 1"
sam_filter work/$BASE.qtrimmed.sorted.bam "-q 30"

#macs2 -t work/$BASE.filtered.bam -g 1e9 --name=work/IRF1 --format=BAM --tsize=36 --mfold=10,30
macs2 -t work/$BASE.filtered.bam -g hs --name=work/IRF1.2 --format=BAM --tsize=36 --mfold=10,30 |& ts | tee work/$BASE.IRF1.2.log
macs2 -t work/$BASE.qtrimmed.sorted.filtered.bam -g hs --name=work/IRF1.3 --format=BAM --tsize=36 --mfold=10,30 |& ts | tee work/$BASE.IRF1.3.log
mv work/IRF1.3_treat_pvalue.{bdg,bedgraph}
igvtools toTDF ./work/IRF1.3_treat_pvalue.{bedgraph,tdf} hg19
mv work/IRF1.3_treat_qvalue.{bdg,bedgraph}
igvtools toTDF ./work/IRF1.3_treat_qvalue.{bedgraph,tdf} hg19



