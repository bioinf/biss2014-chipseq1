#!/usr/bin/bash

# quality control
fastqc test3.fastq

# trimming
fastq_quality_trimmer -Q 33 -v -t 28 -l 36 -i test3.fastq -o test3.trimmed.fastq

# quality control
fastqc test3.trimmed.fastq
