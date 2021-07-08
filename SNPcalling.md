# PRACTICE SNPcalling.md

The raw data is available on the High Capacity Storage of Otago University. 

Contact akstubbs.nz@gmail.com for access

## Quality control

The data is single end __ across two lanes. For this practice, only lane 1 was used. 
The adapter content and the barcodes were subset to understand the structure before being used on whole raw data file:

```
#!/bin/sh
# in source_files_kea

ls
less Kea_conversion_table.txt

zcat SQ1609_CD82MANXX_s_6_fastq.txt.gz | head -n 1000000 > lane2_sample.fq
zcat SQ1630_CD9F4ANXX_s_7_fastq.txt.gz | head -n 1000000 > lane3_sample.fq
```

Then I check the quality of the sequeniging run using fastqc. 

```
#!/bin/sh
module load FastQC

fastqc *fq
```
Looking at the html output, sequence quality is okay but quite some error in first and last few bases.

I will first remove the last bases in each lane as there is a lot of adapter contamination. 
As we are using a reference genome for kea, we do not need to have a common read length. 
Shorter reads are removed by setting minimum number of bp to 40. 

The following uses the whole raw data file. Unzipped first using zcat.

```
#!/bin/sh
module load cutadapt

cd source_files_kea
zcat SQ1609_CD82MANXX_s_6_fastq.txt.gz > SQ1609_1.fastq
zcat SQ1630_CD9F4ANXX_s_7_fastq.txt.gz > SQ1630_2.fastq

cutadapt -j 2 -a ACCGAGATCGGAAGAGC -m 40 -o trimmed_SQ1609_1.fastq SQ1609_1.fastq 
cutadapt -j 2 -a ACCGAGATCGGAAGAGC -m 40 -o trimmed_SQ1630_2.fastq SQ1630_2.fastq 
```
```
=== Summary SQ1609 ===

Total reads processed:             261,124,595
Reads with adapters:               144,422,194 (55.3%)
Reads that were too short:          26,004,723 (10.0%)
Reads written (passing filters):   235,119,872 (90.0%)

Total basepairs processed: 26,135,975,036 bp
Total written (filtered):  18,588,542,969 bp (71.1%)

=== Summary SQ1630 ===

Total reads processed:             259,007,851
Reads with adapters:               127,761,396 (49.3%)
Reads that were too short:          19,101,143 (7.4%)
Reads written (passing filters):   239,906,708 (92.6%)

Total basepairs processed: 25,938,225,049 bp
Total written (filtered):  19,534,525,392 bp (75.3%)
```

## Reference Genome
Created a folder for kea reference genome:
```
#!/bin/sh
cd 
cd source_files_kea/
mkdir -p ref_genome
```
Copied kea reference genome into ref_genome folder and unzipped file

```
#!/bin/sh
curl -L -o ref_genome/kea_ref_genome.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/696/875/GCA_000696875.1_ASM69687v1/GCA_000696875.1_ASM69687v1_genomic.fna.gz

gunzip ref_genome/kea_ref_genome.fasta.gz
```
Load BWA and index reference

```
#!/bin/sh
module load BWA

bwa index kea_ref_genome.fasta
```

## SNP calling

### Demultiplexing

First I extract barcodes from the .key file (SQ1609.txt) of the sequencing platform.

```
#!/bin/sh
##Key provided in folder
cat SQ1609.txt | grep -E "K38" | cut -f 3-4 > barcodes_1.txt
cat SQ1630.txt | grep -E "K38" | cut -f 3-4 > barcodes_2.txt
```
Then I create different folders to deal with samples sequenced across the two different lanes before concatenating them together.

```
#!/bin/sh
mkdir raw1 samples1 raw2 samples2 samples_concatenated

cd ../raw1
ln -s ../source_files_kea/trimmed_SQ1609_1.fastq

cd ../raw2
ln -s ../source_files_kea/trimmed_SQ1630_2.fastq
```

I am now ready to run process_radtags to demultiplex.

Then we will have clean sample files and ready for alignment.

```
#!/bin/sh
module load Stacks

process_radtags -p raw1/ -o samples1/ -b source_files_kea/barcodes_1.txt  --renz_1 pstI --renz_2 mspI -r -c -q --inline_null

process_radtags -p raw2/ -o samples2/ -b source_files_kea/barcodes_2.txt  --renz_1 pstI --renz_2 mspI -r -c -q --inline_null
```

Looks good, keeping ~98% of reads **

```
Outputing details to log: 'samples1/process_radtags.raw1.log'

235119872 total sequences
  4456537 barcode not found drops (1.9%)
   394714 low quality read drops (0.2%)
   682381 RAD cutsite not found drops (0.3%)
229586240 retained reads (97.6%)

Outputing details to log: 'samples2/process_radtags.raw2.log'

239906708 total sequences
  3833078 barcode not found drops (1.6%)
   275458 low quality read drops (0.1%)
   627642 RAD cutsite not found drops (0.3%)
235170530 retained reads (98.0%)

```
Concatenate samples from the two lanes using python:

```
#python
import os
#
for sample in os.listdir ("samples2"):
  print(sample)
  if sample.endswith("gz") and not sample in os.listdir("samples1"):
    raise Exception

for sample in os.listdir ("samples2"):
  if sample.endswith("gz"):
    os.system("cat samples1/"+sample+" samples2/"+sample+" > samples_concatenated/"+sample)
```

Two samples, K38552 and K38547, were merged as they are the same bird.
```
zcat K38552.fq.gz K38547.fq.gz | gzip -c > K38552_47.fq.gz
```

### Alignment and variant calling

First alignment for every sample using BWA. Here is one example command:

```
#!/bin/sh
bwa mem -t 8 $bwa_db $src/${sample}.fq.gz | samtools view -b | samtools sort --threads 4 > ${sample}.bam
```

The complete list of commands is in [realign.sh](realign.sh):

Ran script in realign.sh

```
#!/bin/sh
bash realign.sh
```

I created a population map file for the three haplotype regions of wild kea (North, Central, South).

```
#!/bin/sh
nano popmap_NCS?.txt
<sample file prefix><tab><population ID>

<K38*><tab><north/central/south>

```
I created a new folder for the output reference map:

```
#!/bin/sh
mkdir output_refmap_NCS?
```

Then I run refmap from Stacks quick run to identify low quality individuals:
**run ref_map.pl as a job**

```
#!/bin/sh
ref_map.pl --samples samples_concatenated/ --popmap popmap_NCS?.txt -T 8 -o output_refmap_NCS?/
```

Now check samples and output files for low quality individuals with low sample numbers. 
Remove these, the blank (if it hasn't already been removed?) and any misidentified individuals?

### Re-run ref_map populations again cleanly:

```
module load Stacks

mkdir output_pop_NCS?

populations -P output_refmap_NCS?/ -O output_pop_NCS?/ -M popmap_NCS?.txt -p 1 --write-single-snp --plink --vcf --verbose

Removed 

```
*Link to output files here* 

## Remove individuals with high missing data

VCFtools to look at missing data per individual. Make a new directory for copying VCF file after running populations
```
mkdir filtering
cp output_pop_NCS?/populations.snps.vcf filtering/
```

Run vcftools to output a file that contains missingness on an individual level.

```
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
vcftools --vcf populations.snps.vcf --missing-indv
```
Downloaded ".imiss" file and opened in Excel, sorted individuals by the most missing data (F_MISS). Created a scatterplot to visualise data (F_MISS x N_GENOTYPES-FILTERED)

Four individuals >0.70 F_MISS threshold:

K38612  K38699  K38682  K38689

Vcftools run again to output a new VCF file, which removes individuals with >0.70 missing data. New VCF file renamed appropriately.
```
vcftools --vcf populations.snps.vcf --remove-indv K38612 --remove-indv K38699 --remove-indv K38682 --remove-indv K38689 --recode
mv out.recode.vcf removed_missing.vcf
```

## Filtering


## Structure stuff ??

Created a whitelist.txt file to generate a list of 1000 random loci. 
```
cat output_refmap_NCS/populations.sumstats.tsv | grep -v "#" | cut -f 1 | sort | uniq | shuf | head -n 1000 > whitelist.txt
```

Tried to load structure but missing a mainparams file...??
