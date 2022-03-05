# Variant Calling Workflow
This material is extracted from the Genomics Data Carpentry Lesson
## Aim
- To understand the steps to perform variant calling then overall put all these steps into a script.
## Objectives
- Understand and perform the steps involved in variant calling.
- Describe the types of data formats encountered during variant calling.
- Use command line tools to perform variant calling.

### Assumptions
- You have already performed trimming and filtering of your reads and saved in a directory called trimmed_reads.
- You have a reference genome saved in a directory called ref_genome.

In this workshop, we have already trimmed the reads and downloaded the reference genome for you.
First, it is always good to verify where we are:

```bash

$ cd ~

$ pwd
/home/[your_username]
# good I am ready to work

```

Checking to make sure we have the directory and files for the workshop.

```bash

$ ls
scripting_workshop ...

```
>hint : If you do not have the workshop directory, you can copy it using the command: `cp -r  /nesi/project/nesi02659/scripting_workshop/ ~`  

```bash
$ cd scripting_workshop

$ ls
ref_genome  trimmed_reads 
```

## Alignment to a reference genome
First we need to create directories for the results that will be generated as part of this workflow. We can do this in a single line of code, because mkdir can accept multiple new directory names as input.

```bash
$ mkdir -p results/sam results/bam results/bcf results/vcf
```
### Index the reference genome
Our first step is to index the reference genome for use by BWA. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. Indexing the reference only has to be run once. The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different tool for alignment.

Since we are working on the NeSI HPC, we need to search and load the package before we start using it.
- More on packages will be discussed in the HPC and Slurm section

Search
```bash
$ module spider bwa
```

and then load 
```bash
$ module purge
$ module load BWA/0.7.17-GCC-9.2.0
```

indexing the genome
```bash
$ bwa index ref_genome/ecoli_rel606.fasta
[bwa_index] Pack FASTA... 0.03 sec
[bwa_index] Construct BWT for the packed sequence...
[bwa_index] 1.04 seconds elapse.
[bwa_index] Update BWT... 0.03 sec
[bwa_index] Pack forward-only FASTA... 0.02 sec
[bwa_index] Construct SA from BWT and Occ... 0.57 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index ref_genome/ecoli_rel606.fasta
[main] Real time: 2.462 sec; CPU: 1.702 sec
```

