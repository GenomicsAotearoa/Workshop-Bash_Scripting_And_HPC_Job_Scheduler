#!/bin/bash 

module purge
module load BWA/0.7.18-GCC-12.3.0
module load SAMtools/1.22-GCC-12.3.0
module load BCFtools/1.22-GCC-12.3.0
module load HISAT2/2.2.1-gompi-2023a
module load Subread/2.0.7-GCC-12.3.0


echo "Loaded modules BWA, SAMtools, BCFtools, HISAT2, Subread"
