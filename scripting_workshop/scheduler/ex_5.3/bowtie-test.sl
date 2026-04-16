#!/bin/bash -e

#SBATCH --account		nesi02659
#SBATCH --job-name 		Demo_Bowtie2
#SBATCH --cpus-per-task 	2
#SBATCH --time 			00:03:00
#SBATCH --mem 			1G
#SBATCH --output 		slurmout/bowtie2Test-%j.out
#SBATCH --error 		slurmout/bowtie2-%j.err
#SBATCH --mail-type		ALL
#SBATCH --mail-user		yourname@email.com

##====Loading_Modules===============##
module purge
module load Bowtie2/2.5.4-GCC-12.3.0
module load SAMtools/1.22-GCC-12.3.0
module load BCFtools/1.22-GCC-12.3.0

##====Path_variables===============##
#We will be using $PWD represent current working directory

echo "$PWD"

#indexing a reference genome
bowtie2-build $PWD/input_data/reference/lambda_virus.fa lambda_virus

#Aligning examples reads
bowtie2 -x lambda_virus -U $PWD/input_data/reads/reads_1.fq -S eg1.SAMtools

#Paired-end example
bowtie2 -x lambda_virus -1 $PWD/input_data/reads/reads_1.fq -2 $PWD/input_data/reads/reads_2.fq -S eg2.sam

#Local alignment example
bowtie2 --local -x lambda_virus -U $PWD/input_data/reads/longreads.fq -S eg3.sam

#Convert the SAM to BAM:
samtools view -bS eg2.sam > eg2.bam

#Convert the BAM file to a sorted BAM file
samtools sort eg2.bam > eg2.sorted.bam

#We now have a sorted BAM file called `eg2.sorted.bam`. Sorted BAM is a useful
#format because the alignments are (a) compressed, which is convenient for
#long-term storage, and (b) sorted, which is conveneint for variant discovery.
#To generate variant calls in VCF format, run
bcftools mpileup -f $PWD/input_data/reference/lambda_virus.fa eg2.sorted.bam | bcftools call -mv -Oz -o eg2.raw.bcf
