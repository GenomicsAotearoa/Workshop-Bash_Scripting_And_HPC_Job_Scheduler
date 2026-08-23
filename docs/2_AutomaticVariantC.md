# Automating a Variant Calling Workflow

!!! clipboard-list "Aim"

    - Put all the steps from the previous lesson into a script.

### Variant calling workflow


!!! rectangle-list "Remember our variant calling workflow has the following steps:"

    - Index the reference genome for use by bwa and samtools.
    - Align reads to reference genome.
    - Convert the format of the alignment to sorted BAM, with some intermediate steps.
    - Calculate the read coverage of positions in the genome.
    - Detect the single nucleotide variants (SNVs).
    - Filter and report the SNVs in VCF (variant calling format).


!!! terminal-2 "Let's start with creating a new directory as our script working space and copy all the required resources."

    ```bash
    pwd
    ```
    **Output** - `/home/$USER/scripting_workshop`
    ```bash
    mkdir script_workspace
    ```
    ```bash
    cd script_workspace
    ```
    ```bash
    #Don't forget the  . at the end of the line
    cp -r ../variant_calling/ref_genome ../variant_calling/trimmed_reads .
    ```
    ```bash
    ls
    ```

    **Output** - `ref_genome trimmed_reads `
    


!!! terminal-2 "Now we are ready to start building the script."

    ```bash
    nano variant_calling.sh
    ```
    <br>
    ??? jupyter "Prefer Jupyter *file* over `nano` ?"
        You are welcome to choose Jupyter interactive *file* option over `nano`. If this is your preferred option, below are the instructions on how to **create** a file

        1. **Make sure** to navigate the left file explorer to correct *path*
        2. <KBD>Right click</KBD> on the explorer to open the menu and choose `New file`
        3. Rename the file as `variant_calling.sh`

        ![image](./nesi_images/open_text_file.png)


??? tip "Variant calling command line steps"
    
    In the order they were run in the previous step-by-step lesson:

    ```bash

    mkdir -p results/sam results/bam results/bcf results/vcf

    module purge
    module load BWA/0.7.18-GCC-12.3.0

    bwa index ref_genome/ecoli_rel606.fasta

    bwa mem ref_genome/ecoli_rel606.fasta trimmed_reads/SRR2584866_1.trim.sub.fastq trimmed_reads/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam

    module load SAMtools/1.22-GCC-12.3.0

    samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam

    samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam

    samtools flagstat results/bam/SRR2584866.aligned.sorted.bam 

    module load BCFtools/1.22-GCC-12.3.0

    bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf -f ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam

    bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf

    vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf > results/vcf/SRR2584866_final_variants.vcf
    ```

!!! terminal-2 "In the text editor, type the commands"

    ```bash linenums="1"
    #!/bin/bash 
    
    # This script runs the variant calling pipeline from mapping to vcf.

    # Jane Doe
    echo $(date)

    # Setting -e flag tells bash to exit immediately if a command exits with a non-zero (error) status
    set -e
    
    # Load all the required modules
    module purge
    module load BWA/0.7.18-GCC-12.3.0
    module load SAMtools/1.22-GCC-12.3.0
    module load BCFtools/1.22-GCC-12.3.0
    
    # create the results directories
    mkdir -p results/sam results/bam results/bcf results/vcf
    
    # indexing the genome
    genome=ref_genome/ecoli_rel606.fasta
    bwa index $genome
    
    # create a loop that map reads to the genome, sort the bam files and call variants
    for file in trimmed_reads/*_1.trim.sub.fastq
        do
        echo "working with file $file"
    
        base=$(basename $file _1.trim.sub.fastq)
        echo "base name is $base"
    
       # setting the variables
       fq1=trimmed_reads/${base}_1.trim.sub.fastq
       fq2=trimmed_reads/${base}_2.trim.sub.fastq
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
      bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
      bcftools call --ploidy 1 -m -v -o $variants $raw_bcf
      vcfutils.pl varFilter $variants > $final_variants
    
    done
    ```

    <br>

    ??? tip "Optional: add flagstat output too"

        ```bash
        # setting the variables - add AFTER sorted_bam variable
        flagstat=results/bam/${base}.aligned.sorted.bam.summary

        # running the analysis steps  - add AFTER samtools sort command
        samtools flagstat $sorted_bam > $flagstat
        ```

    ??? hand-holding-dollar "Shell variables"
        A variable is a character string to which we assign a value. The value assigned could be a number, text, filename, device, or any other type of data. A variable is nothing more than a pointer to the actual data. The shell enables you to create, assign, and delete variables.

    !!! bell "What is `#!/bin/bash` ?"
        * This line is called a [**shebang**.](https://en.wikipedia.org/wiki/Shebang_%28Unix%29) It's not a comment, even though it starts with `#`.
        * It tells the operating system which interpreter to use to run the rest of the file. `#!/bin/bash` means "run this script with bash."
        * It **must be the very first line** of the file; no blank lines or comments before it, or it won't be recognised.
        * If you run your script as `bash script.sh`, the shebang line is ignored and bash is used as the interpreter, but if you run a directly executed script *i.e.,* `./script.sh` the shebang line is needed to call the interpreter. To be safe, always have the shebang line in your scripts.  
        * `#!/bin/bash` assumes bash lives at `/bin/bash`. This is generally the case, but if bash lives somewhere else you should specify the correct path instead. 



!!! terminal-2 "Running the script (Running the script)"

    ```bash
    bash ./variant_calling.sh
    ```


!!! tip "Adding executable permissions"

    The way the script is written means we have to indicate which program to use whenever we are running it. 
    To be able to run without calling bash, we need to change the script permissions.

    !!! terminal "script"
    
        ```bash 
        ls -l variant_calling.sh 
        ```

        **Output** `-rw-rw----+ 1 userID nesi02659 1411 Apr 16 13:53 variant_calling.sh`
        
        ```bash
        chmod u+x variant_calling.sh
        ```
        ```bash
        ls -l variant_calling.sh 
        ```

        **Output** - `-rwxrw----+ 1 userID nesi02659 1411 Apr 16 13:53 variant_calling.sh`

        - note colour change on the script filename
        
    Now we can execute the script without calling bash

    !!! terminal "script"
        ```bash
        ./variant_calling.sh
        ```

In the [Next Lesson](3_RNAseq.md), we will do a similar pipeline, but with RNAseq data. 



