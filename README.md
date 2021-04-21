# analytical_methods
Excersice for the subject: Analytical methods in cancer and population genomics and transcriptomics

At first I would like download the reference human genome:

```console
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz
```

At first I would like to align tumor data to the reference genome:
- first of all -> I creating BWA index for reference genome data:
```console
bwa index -p chrX chrX.fa 
```

At the second I tumour data and wild type data to the reference genome:
```console
bwa mem -M chrX tu.r1.fq tu.r2.fq 2> logs/bwa.err > results/tu.sam
bwa mem -M chrX wt.r1.fq wt.r2.fq 2> logs/bwa.err > results/wt.sam
```

Covert a SAM file to a BAM file
```console
samtools view -O BAM -o tu.bam tu.sam
samtools view -O BAM -o wt.bam wt.sam
```

Sort and index the BAM file
```console
samtools sort -T temp -O bam -o tu.sorted.bam tu.bam
samtools sort -T temp -O bam -o wt.sorted.bam wt.bam
samtools index tu.sorted.bam
samtools index wt.sorted.bam
```

Mark PCR Duplicates / to remove PCR suplicates that may have been introduces during the library construction stage
```console
picard-tools MarkDuplicates I=tu.sorted.bam O=tu.markdup.bam M=tu.metrics.txt
picard-tools MarkDuplicates I=wt.sorted.bam O=wt.markdup.bam M=wt.metrics.txt

samtools index tu.markdup.bam
samtools index wt.markdup.bam
```

