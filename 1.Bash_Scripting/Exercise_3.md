# RNA-seq Mapping And Count Data Workflow
This material is extracted from the [RNA-seq workshop](https://github.com/GenomicsAotearoa/RNA-seq-workshop) Lesson
## Aim
- To develop a pipeline that does mapping and count the number of reads that mapped then overall put all these steps into a script.
## Objectives
- Understand and perform the steps involved in RNA-seq mapping and read count.
- Use command line tools to run the pipeline.

### Assumptions
- You have already performed trimming and filtering of your reads and saved in a directory called trimmed_reads.
- You have a reference genome saved in a directory called ref_genome.

In this workshop, we have already trimmed the reads and downloaded the reference genome for you.
First, it is always good to verify where we are:

```bash

$ cd ~

$ pwd
/home/[your_username]/scripting_workshop
# good I am ready to work

```

Checking to make sure we have the directory and files for the workshop.

```bash

$ ls
rna_seq  variant_calling

```
>hint : If you do not have the workshop directory, you can copy it using the command: `cp -r  /nesi/project/nesi02659/scripting_workshop/ ~`  

```bash
$ cd rna_seq

$ ls
ref_genome  trimmed_reads 
```


