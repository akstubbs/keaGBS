# Kea Sex Chromosomes

Kea do not have annotated sex chromosomes in reference genome assembly. Therefore we want to blast the reference genome of Kea to the reference genome of Kakapo, which do have annotated sex chromosomes.

Kea - align to kākāpō

W: https://www.ncbi.nlm.nih.gov/nuccore/CM013773.2

Z: https://www.ncbi.nlm.nih.gov/nuccore/CM013763.2


Ran blastn as a job.
```
#!/bin/sh
module load BLAST
cd ../source_files_kea/

mkdir kakapo_sex_chr
cd kakapo_sex_chr

makeblastdb -in kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -dbtype nucl

blastn -query ../kea_ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -outfmt 7 -num_threads 8 -max_target_seqs 3 > blast_results.txt
```
## Kakapo sex chromosome alignment:

W: NC_044301.1

Z: NC_044302.1

Create a list of kea genome scaffolds that are potentially sex chromosomes in the kea genome. 

```
#!/bin/sh
cd uoo03341/source_files_kea/kakapo/
less blast_results.txt

grep NC_044301.1 blast_results.txt | cut -f 1-2 | sort | uniq > potential_sex_chrW.scaffolds.txt

grep NC_044302.1 blast_results.txt | cut -f 1-2 | sort | uniq > potential_sex_chrZ.scaffolds.txt
```
Create a .fai file (fasta index file) from the kea reference genome. We get the length of scaffolds from the second column of the kea_ref_genome.fa.fai

```
#!/bin/sh
module load SAMtools
cd ../ref_genome/
samtools faidx kea_ref_genome.fasta
less kea_ref_genome.fasta.fai 
```

Combine file one (W) and two (Z) by concatenating. Then sort and remove duplicate parts of kea genome scaffolds.
```
cat potential_sex_chrW.scaffolds.txt potential_sex_chrZ.scaffolds.txt > potential_sex_chr_concatenated.scaffolds.txt
sort -u -k1,1 potential_sex_chr_concatenated.scaffolds.txt > potential_sex_chrsort.scaffolds.txt
```
```-u``` for unique, and ```-k1,1``` for the key field 1 (the filed with the kea genome scaffold).


Now need to find the proportion of kea genome that is attributed to sex-chromosomes by blast.

Lengths of scaffolds are from the second column of the kea_ref_genome.fa.fai file.

Used ```awk``` to get sum of the total scaffold lengths:

```
module load bioawk
awk '{sum+=$2;} END {print "Total of 2nd Column:" sum}2' ../kea_ref_genome/kea_ref_genome.fasta.fai 
```
```
Total of 2nd Column:1053559886
```

Sort the kea reference genome, and join it together with the potential kea sex chromosome scaffolds - to then determine the lengths of the potential sex chromosomes.
```
sort ../kea_ref_genome/kea_ref_genome.fasta.fai > kea_ref_genome_sorted.fasta.fai
join potential_sex_chrsort.scaffolds.txt kea_ref_genome_sorted.fasta.fai > potential_sex_chr_withlength.txt

less potential_sex_chr_withlength.txt 
awk '{sum+=$3;} END {print "Total of 3nd Column:" sum}3' potential_sex_chr_withlength.txt 
```
``` 
Total of 3nd Column:657976221 
```

Potential sex chromosome in genome =  657976221/1053559886 = 62.5%

---
---
**Map kakapo sex chromosomes to kea genome**

The proportion of the genome that was potentially sex chromosomes was really high (~65%, but should be ~10%).

Therefore, instead of mapping whole kea to whole kakapo - Map the kakapo Z and W chromosome to the whole kea genome. The kea becomes the database.

```
#!/bin/sh
module load BLAST

cd kakapo_sex_chr
makeblastdb -in ../kea_ref_genome/kea_ref_genome.fasta -dbtype nucl
```
```
#!/bin/sh
mkdir kakapo_w_chr kakapo_z_chr
```
Copy in the kakapo sex chromosomes from NCBI.

**W chromosome:**
```
#!/bin/sh
curl -L -o kakapo_w_chr/kakapo_w_chr.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Strigops_habroptila/all_assembly_versions/GCF_004027225.2_bStrHab1.2.pri/GCF_004027225.2_bStrHab1.2.pri_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chrW.fna.gz

gunzip kakapo_w_chr/kakapo_w_chr.fasta.gz
```
**Z chromosome:**
```
#!/bin/sh
curl -L -o kakapo_z_chr/kakapo_z_chr.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Strigops_habroptila/all_assembly_versions/GCF_004027225.2_bStrHab1.2.pri/GCF_004027225.2_bStrHab1.2.pri_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chrZ.fna.gz

gunzip kakapo_z_chr/kakapo_z_chr.fasta.gz
```

Ran w and z blastn as jobs:
```
#!/bin/sh
blastn -query kakapo_w_chr/kakapo_w_chr.fasta -db ../kea_ref_genome/kea_ref_genome.fasta -outfmt 6 -num_threads 10 -max_target_seqs 100000 > blast_results_kakapo_w_to_kea.txt

blastn -query kakapo_z_chr/kakapo_z_chr.fasta -db ../kea_ref_genome/kea_ref_genome.fasta -outfmt 6 -num_threads 10 -max_target_seqs 100000 > blast_results_kakapo_z_to_kea.txt
```


... 
~~potential_sex_chr.scaffolds.txt -> will find sex chromosomes in kea genome
lengths of scaffolds from the second column of kea_ref_genome.fa.fai file~~

could do this by ~~extracting total scaffold lengths in second column of the kea_ref_genome.fa.fai file~~, and extracting total of potential sex chromosome scaffold lengths ??

1. List of sex chr from potential_sex_chr.scaffolds.txt file
2. For each of thoese sex chromosomes, find corresponding length in kea_ref_genome.fa.fai file
3. use potential txt file to extract lines of scaffold in whole kea file.
4. repeat awk to get sum of scaffold lengths that are potential sex chromosomes.

Sum those lengths.

vcftools --not-chr ??
