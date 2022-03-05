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
```

```bash
$ ls
ref_genome  trimmed_reads 
```

## Alignment to a reference genome
First we need to create directories for the results that will be generated as part of this workflow. We can do this in a single line of code, because mkdir can accept multiple new directory names as input.

```bash
$ mkdir -p results/sam results/bam results/bcf results/vcf
```

Since we are working on the NeSI HPC, we need to search and load the package before we start using it.

Search

```bash
$ module spider fastqc
```

and then load 

```bash
$ module purge
$ module load FastQC/0.11.9
```
