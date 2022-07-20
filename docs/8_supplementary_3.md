# S3 : Solutions

<p style="text-align:left;">
    <b><a class="btn" href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/workshop_material/7_supplementary_2.html" style="background: var(--bs-green);font-weight:bold">&laquo;7. Supplementary 2</a></b>
</p>

??? success "Exercise 5.4 😬"	

    ```bash

    #!/bin/bash -e

    #SBATCH --account		nesi02659
    #SBATCH --job-name 	    variant_calling_workflow
    #SBATCH --cpus-per-task 2
    #SBATCH --time 			00:15:00
    #SBATCH --mem 			4G
    #SBATCH --output 		slurmout/variant_calling-%j.out
    #SBATCH --error 		slurmout/variant_calling-%j.err
    #SBATCH --mail-type		END
    #SBATCH --mail-user		myemail@email.co.nz

    # Load all the required modules
    module purge
    module load BWA/0.7.17-GCC-9.2.0
    module load SAMtools/1.13-GCC-9.2.0
    module load BCFtools/1.13-GCC-9.2.0

    echo "$PWD"

    # create the results directories
    mkdir -p results/sam results/bam results/bcf results/vcf

    # indexing the genome
    genome=~/scripting_workshop/variant_calling/ref_genome/ecoli_rel606.fasta
    trimmed=~/scripting_workshop/variant_calling/trimmed_reads
    bwa index $genome

    # create a loop that map reads to the genome, sort the bam files and call variants
    for fq1 in ${trimmed}/*_1.trim.sub.fastq
     do
        echo "working with file $fq1"

        base=$(basename $fq1 _1.trim.sub.fastq)
        echo "base name is $base"

        # setting the variables
        fq1=${trimmed}/${base}_1.trim.sub.fastq
        fq2=${trimmed}/${base}_2.trim.sub.fastq
        sam=results/sam/${base}.aligned.sam
        bam=results/bam/${base}.aligned.bam
        sorted_bam=results/bam/${base}.aligned.sorted.bam
        raw_bcf=results/bcf/${base}_raw.bcf
        variants=results/vcf/${base}_variants.vcf
        final_variants=results/vcf/${base}_final_variants.vcf

        # running the analysis steps
        bwa mem $genome $fq1 $fq2 > $sam
        samtools view -S -b $sam > $bam
        samtools sort -o $sorted_bam $bam
        samtools index $sorted_bam
        bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
        bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
        vcfutils.pl varFilter $variants > $final_variants

    done

    echo "DONE"
    ```



??? success "Exercise 5.5 😬"	


    ```bash
    #!/bin/bash -e

    #SBATCH --account		nesi02659
    #SBATCH --job-name 	rna-seq_workflow
    #SBATCH --cpus-per-task 	2
    #SBATCH --time 			00:15:00
    #SBATCH --mem 			4G
    #SBATCH --output 		rna-seq_workflow-%j.out
    #SBATCH --error 		rna-seq_workflow-%j.err
    #SBATCH --mail-type		END
    #SBATCH --mail-user		myemail@email.org.nz


    echo "$PWD"

    mkdir -p ~/scripting_workshop/scheduler/ex_5.5/{Mapping,Counts} && cd ~/scripting_workshop/scheduler/ex_5.5/
    cp -r /nesi/project/nesi02659/scripting_workshop/rna_seq/* ./

    module purge
    module load HISAT2/2.2.0-gimkl-2020a
    module load SAMtools/1.10-GCC-9.2.0
    module load Subread/2.0.0-GCC-9.2.0

    echo $PWD

    #index file
    hisat2-build -p 4 -f $PWD/ref_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa $PWD/ref_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel


    #Mapping Samples to the reference genome

    for filename in $PWD/trimmed_reads/*
      do
        base=$(basename ${filename} .fastq)
        hisat2 -p 4 -x $PWD/ref_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel -U $filename -S $PWD/Mapping/${base}.sam --summary-file $PWD/Mapping/${base}_summary.txt
      done

    #Convert SAMfiles to BAM

    for filename in $PWD/Mapping/*.sam
      do
       base=$(basename ${filename} .sam)
       samtools view -S -b ${filename} -o $PWD/Mapping/${base}.bam
      done

    #Sort BAM files

    for filename in $PWD/Mapping/*.bam
      do
        base=$(basename ${filename} .bam)
        samtools sort -o $PWD/Mapping/${base}_sorted.bam ${filename}
      done

    #count how many reads aligned to each genome feature (exon).

    featureCounts -a $PWD/ref_genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf -o $PWD/Counts/yeast_counts.txt -T 2 -t exon -g gene_id $PWD/Mapping/*sorted.bam
```

{% endcapture %}

{% include exercise.html title="es3dot2" content=es3dot2%}


<p align="center"><b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/">Back to homepage</a></b></p>