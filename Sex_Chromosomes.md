# Kea Sex Chromosomes

Kea do not have annotated sex chromosomes in reference genome assembly. Therefore we want to blast the reference genome of Kea to the reference genome of Kakapo, which do have annotated sex chromosomes.

Kea - align to kākāpō

Z: https://www.ncbi.nlm.nih.gov/nuccore/CM013763.2

W: https://www.ncbi.nlm.nih.gov/nuccore/CM013773.2


Ran blastn as a job.
```
module load BLAST
cd ../source_files_kea/

mkdir kakapo
cd kakapo

makeblastdb -in kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -dbtype nucl

blastn -query ../ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -outfmt 7 -num_threads 8 > blast_results.txt
```

Kakapo sex chromosome alignment:

Z?: NC_044301.1

W?: NC_044302.1


Create a list of kea genome scaffolds that are potentially sex chromosomes in the kea genome. 

```
cd uoo03341/source_files_kea/kakapo/
less blast_results.txt

grep NC_044301.1 blast_results.txt | cut -f 1-2 | sort | uniq > potential_sex_chr.scaffolds.txt

grep NC_044302.1 blast_results.txt | cut -f 1-2 | sort | uniq > potential_sex_chr.scaffolds.txt
```
Create a .fai file (fasta index file) from the kea reference genome. We get the length of scaffolds from the second column of the kea_ref_genome.fa.fai

```
module load SAMtools
cd ../ref_genome/
samtools faidx kea_ref_genome.fasta
less kea_ref_genome.fasta.fai 
```

Combine file one and two, to find the proportion of your genome that is attributed to sex-chromosomes by blast.


... 
potential_sex_chr.scaffolds.txt -> will find sex chromosomes in kea genome
lengths of scaffolds from the second column of kea_ref_genome.fa.fai file

could do this by extracting total scaffold lengths in second column of the kea_ref_genome.fa.fai file, and extracting total of potential sex chromosome scaffold lengths ?? 

1. List of sex chr from potential_sex_chr.scaffolds.txt file
2. For each of thoese sex chromosomes, find corresponding length in kea_ref_genome.fa.fai file

Sum those lengths.

vcftools --not-chr ??
