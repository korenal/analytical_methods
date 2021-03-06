# Cancer Genomics Data Analysis Exercise

At first I would like download the reference human genome:

```console
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz
```

### 1. Aligning reads with BWA-MEM:
Creating BWA-MEM index on reference genome data:
```console
bwa index -p chrX chrX.fa 
```

Alignment tumor and reference data & control and reference data:
```console
bwa mem -M chrX tu.r1.fq tu.r2.fq 2> logs/bwa.err > results/tu.sam
bwa mem -M chrX wt.r1.fq wt.r2.fq 2> logs/bwa2.err > results/wt.sam
```

Covert a SAM file to a BAM file:
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

Mark PCR Duplicates - to remove PCR duplicates that may have been introduces during the library construction stage:
```console
picard-tools MarkDuplicates I=tu.sorted.bam O=tu.markdup.bam M=tu.metrics.txt
picard-tools MarkDuplicates I=wt.sorted.bam O=wt.markdup.bam M=wt.metrics.txt

samtools index tu.markdup.bam
samtools index wt.markdup.bam
```

Extract data for chromosome X from 20Mbp to 40Mbp:
```console
samtools view -b tu.markdup.bam chrX:20000000-40000000 > tu.output.bam
samtools view -b wt.markdup.bam chrX:20000000-40000000 > wt.output.bam
samtools index tu.output.bam
samtools index wt.output.bam
```

### 2. Generating a read-depth plot:
Computing the read depth at each position:
```console
samtools depth tu.output.bam > tu.coverage
samtools depth wt.output.bam > wt.coverage
```
Plot the data in R studio:
```r
library(reshape)
tu.chrX <- read.table("C:/Users/lucie/OneDrive/Dokumenty/tu.coverage", header=FALSE, sep='\t', na.strings="NA", dec=".", strip.white=TRUE)
wt.chrX <- read.table("C:/Users/lucie/OneDrive/Dokumenty/wt.coverage", header=FALSE, sep='\t', na.strings="NA", dec=".", strip.white=TRUE)

tu.chrX <-rename(tu.chrX,c(V1="Chr", V2="locus", V3="depth")) 
wt.chrX <-rename(wt.chrX,c(V1="Chr", V2="locus", V3="depth"))

library(lattice, pos=10) 
xyplot(depth ~ locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), tu.chrX, main="depth by tumor data")

library(lattice, pos=10) 
xyplot(depth ~ locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), scales=list(x=list(relation='same'), y=list(relation='same')), wt.chrX, main="Read-depth plot of control data", col = "red")
```

![tumor - read-depth graph](https://user-images.githubusercontent.com/53037749/115723040-dd113a80-a37f-11eb-83e2-2e6b7b9c1442.png)
![cotrol - OK](https://user-images.githubusercontent.com/53037749/115750409-6e8ca680-a398-11eb-904e-50317c7bdadd.png)
