# Kea Sex Chromosomes

Kea do not have annotated sex chromosomes in reference genome assembly. Therefore, to identify scaffolds of sex chromosomal origin in the fragmented Kea reference assembly, we need to align the scaffolds to the kakapo (Strigops habroptilus) reference genome, which has both sex chromosomes (Z & W) assembled (GenBank accession: GCF_004027225.2).

### Kea reference genome index file

Created a .fai file (fasta index file) from the kea reference genome. We get the length of scaffolds from the second column of the kea_ref_genome.fa.fai
```
#!/bin/sh
module load SAMtools
cd ../ref_genome/
samtools faidx kea_ref_genome.fasta
less kea_ref_genome.fasta.fai 
```

## Complete pipeline:

Followed pipeline in github file [Identifying_Sex_Chromosome.txt](https://github.com/akstubbs/keaGBS/blob/main/Identifying_Sex_Chromosomes.txt)

Made a final list of kea scaffold names that align with Z & W sex chromosomes.

Now will exclude all SNPs found for these scaffolds from vcf file. 

## Removing kea sex chromosomes.

vcftools --not-chr ??
