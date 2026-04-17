# S3 : Solutions

??? circle-check "Exercise 3.1 - RNA-Seq Mapping and Count Data 😺"

    Copy the text below and save it as a bash script `rnaseq.sh`

    ```bash linenums="1"
    #!/bin/bash -e

    # Jane Doe
    # 16 April 2026
    #This script uses relative paths. It has to be executed from rna_seq parent directory OR edit paths accordingly. 

    # Load required modules and setup the environment
    module purge
    module load HISAT2/2.2.1-gompi-2023a
    module load SAMtools/1.22-GCC-12.3.0
    module load Subread/2.0.7-GCC-12.3.0

    #Print the current working directory
    echo "${PWD}"

    #Create results directories. We are moving away from "Mapping" and "Counts"
    mkdir -p results/{sam,bam,counts}

    # Index genome
    hisat2-build -p 2 -f ref_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa \
                         ref_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel

    # Align to indexed genome
    for filename in trimmed_reads/*.fastq
    do
          # Extract base name
          base=$(basename ${filename} .fastq)

          # Align to the reference genome
          hisat2 -p 2 -x ref_genome/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel -U $filename \
          -S results/sam/${base}.sam --summary-file results/sam/${base}.summary.txt

          # Convert SAM to BAM
          samtools view -S -b results/sam/${base}.sam | samtools sort -o results/bam/${base}_sorted.bam

          # Extract stats for mapping
          samtools flagstat results/bam/${base}_sorted.bam > results/bam/${base}_mapstat.txt
    done

    # count how many reads aligned to each genome feature (exon).
    featureCounts -a ref_genome/Saccharomyces_cerevisiae.R64-1-1.99.gtf \
                  -o results/counts/yeast_counts.txt -T 4 -t exon -g gene_id \
                     results/bam/${base}_sorted.bam

    ```

??? circle-check "Exercise 5.4 - Variant calling workflow slurm script 📜"	
    
    Copy the text below and save it as a slurm script `variant-calling.sl`

    ```bash linenums="1"
    #!/bin/bash -e

    #SBATCH --account		nesi02659
    #SBATCH --job-name 	    variant_calling_workflow
    #SBATCH --cpus-per-task 2
    #SBATCH --time 			00:15:00
    #SBATCH --mem 			4G
    #SBATCH --output 		variant_calling-%j.out
    #SBATCH --error 		variant_calling-%j.err
    #SBATCH --mail-type		END
    #SBATCH --mail-user		myemail@email.co.nz # remember to update with your own email!

    # Load all the required modules
    module purge
    module load BWA/0.7.18-GCC-12.3.0
    module load SAMtools/1.22-GCC-12.3.0
    module load BCFtools/1.22-GCC-12.3.0


    echo "$PWD"

    # create the results directories
    mkdir -p results/sam results/bam results/bcf results/vcf

    # indexing the genome
    genome=~/scripting_workshop/variant_calling/ref_genome/ecoli_rel606.fasta # you may need to explicity set the path here!
    trimmed=~/scripting_workshop/variant_calling/trimmed_reads # you may need to explicity set the path here!
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


??? circle-check "Exercise 5.5 - RNA-seq workflow slurm script 📜"

    Copy the text below and save it as a slurm script `rnaseq.sl`

    ```bash linenums="1"
    #!/bin/bash -e

    #SBATCH --account       nesi02659
    #SBATCH --job-name      rna-seq_workflow
    #SBATCH --cpus-per-task 4
    #SBATCH --time          00:15:00
    #SBATCH --mem           4G
    #SBATCH --output        rna-seq_workflow-%j.out
    #SBATCH --error         rna-seq_workflow-%j.err
    #SBATCH --mail-type     END
    #SBATCH --mail-user     myemail@email.org.nz

    # load modules
    module purge
    module load HISAT2/2.2.1-gompi-2023a
    module load SAMtools/1.22-GCC-12.3.0
    module load Subread/2.0.7-GCC-12.3.0

    # create results directories
    mkdir -p ex_5.5/{Mapping,Counts} && cd ex_5.5

    # set variables
    genomedir=~/scripting_workshop/rna_seq/ref_genome  # update path as needed!
    trimmeddir=~/scripting_workshop/rna_seq/trimmed_reads  # update path as needed!

    # index genome
    hisat2-build -p 4 -f ${genomedir}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa ${genomedir}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel

    # map reads, convert to BAM and sort
    for filename in ${trimmeddir}/*
      do
        base=$(basename ${filename} .fastq)
        hisat2 -p 4 -x ${genomedir}/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel \
        -U $filename -S Mapping/${base}.sam \
        --summary-file Mapping/${base}_summary.txt
        samtools view -S -b Mapping/${base}.sam -o Mapping/${base}.bam
        samtools sort -o Mapping/${base}_sorted.bam Mapping/${base}.bam
        samtools flagstat Mapping/${base}_sorted.bam > Mapping/${base}_mapstat.txt
      done

    # count reads per feature
    featureCounts -a ${genomedir}/Saccharomyces_cerevisiae.R64-1-1.99.gtf \
                  -o Counts/yeast_counts.txt \
                  -T 4 -t exon -g gene_id \
                  Mapping/*sorted.bam

    echo "All done!"
    ```



<p align="center"><b><a href="https://genomicsaotearoa.github.io/Workshop-Bash_Scripting_And_HPC_Job_Scheduler/">Back to homepage</a></b></p>