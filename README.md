# analytical_methods
Excersice for the subject: Analytical methods in cancer and population genomics and transcriptomics

At first I would like download the reference human genome:

```console
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
```

At first I would like to align tumor data to the reference genome:
- first of all -> I creating BWA index for reference genome data:
```console
bwa index -p GRCh37 GRCh37_latest_genomic.fna
```
