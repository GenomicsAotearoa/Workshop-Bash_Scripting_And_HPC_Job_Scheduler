# Automating a Variant Calling Workflow

## Aim
- Put all the steps from the previous lesson into a script.

Remember our variant calling workflow has the following steps:
- Index the reference genome for use by bwa and samtools.
- Align reads to reference genome.
- Convert the format of the alignment to sorted BAM, with some intermediate steps.
- Calculate the read coverage of positions in the genome.
- Detect the single nucleotide variants (SNVs).
- Filter and report the SNVs in VCF (variant calling format).

