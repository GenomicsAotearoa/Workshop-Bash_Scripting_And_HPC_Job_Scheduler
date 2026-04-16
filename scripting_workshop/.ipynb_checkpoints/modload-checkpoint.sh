#!/bin/bash 

module purge
module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.13-GCC-9.2.0
module load BCFtools/1.13-GCC-9.2.0
module load HISAT2/2.2.0-gimkl-2020a
module load Subread/2.0.0-GCC-9.2.0


echo "Loaded modules BWA, SAMtools, BCFtools,HISAT2,Subread"
