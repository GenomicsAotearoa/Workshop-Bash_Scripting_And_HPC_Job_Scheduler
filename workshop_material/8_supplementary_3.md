## Solutions

### Exercise 5.4 😬	
{% capture es3dot1 %}

```bash

#!/bin/bash -e

#SBATCH --account		nesi02659
#SBATCH --job-name 	    variant_calling_workflow
#SBATCH --cpus-per-task 	2
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

{% endcapture %}

{% include exercise.html title="es3dot1" content=es3dot1%}