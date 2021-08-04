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

Created a bed file for the scaffolds that are on the sex chromosomes.

```
#!/bin/sh
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ../../kea_ref_genome.fasta.fai > kea_genome.bed
sort kea_genome.bed >> kea_genome.bed
join kea_ZWscaffold_names.txt kea_genome.bed > sex_chr.bed
```

Excluded the sex chromosome scaffolds in the bed file using VCFtools

```
#!/bin/sh
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
vcftools --vcf filtered_dp3_34_md80.vcf --exclude-bed sex_chr.bed --recode
```
Output below. Left with 22310 possible sites.
```
VCFtools - 0.1.15
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
        --vcf filtered_dp3_34_md80.vcf
        --out sex_chr_rm
        --recode
        --exclude-bed sex_chr.bed

After filtering, kept 187 out of 187 Individuals
Outputting VCF file...
        Read 2452 BED file entries.
After filtering, kept 22310 out of a possible 22944 Sites
Run Time = 8.00 seconds
```
