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
mkdir raw1P samples1P

module load cutadapt

cd source_files_kea
zcat SQ1609_CD82MANXX_s_6_fastq.txt.gz > SQ1609_1.fastq

cutadapt -j 2 -a AGATCGGAAGAGC -m 40 -o trimmed_SQ1609_1.fastq SQ1609_1.fastq 
##storage file error appears - re run cutadapt and scripts from demultiplexing process_radtags step onwards (so all sequences completed for all samples)

cd ../raw1P
ln -s ../source_files_kea/trimmed_SQ1609_1.fastq

fastqc trimmed_SQ1609_1.fastq #don't need to check this
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
```
Then I create different folders to deal with samples sequenced across the two different lanes before concatenating them together.


I am now ready to run process_radtags to demultiplex.

Then we will have clean sample files and ready for alignment.

```
#!/bin/sh
process_radtags -p raw1P/ -o samples/ -b source_files_kea/barcodes_1.txt  --renz_1 pstI --renz_2 mspI -r -c -q --inline_null
```

Looks good, keeping nearly 98% of reads

```
Outputing details to log: 'samples1P/process_radtags.raw1P.log'

172273235 total sequences
  3030216 barcode not found drops (1.8%)
   268364 low quality read drops (0.2%)
   443076 RAD cutsite not found drops (0.3%)
168531579 retained reads (97.8%)
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

I created a population map file for the three haplotype regions of wild kea (North, Central, South)

```
#!/bin/sh
nano popmap_2_NCS.txt
<sample file prefix><tab><population ID>

<K38*><tab><north/central/south>

```
I created a new folder for the output reference map:

```
#!/bin/sh
cd practice/
mkdir output_refmap_NCS
```

Then I run refmap from Stacks quick run to identify low quality individuals:
**run ref_map.pl as a job

```
#!/bin/sh
module load Stacks
ref_map.pl --samples ./practice/alignment --popmap popmap_2_NCS.txt -T 8  -o ./practice/output_refmap_NCS
```

Now check samples and output files for low quality individuals with low sample numbers. 
Remove these, the blank (if it hasn't already been removed?) and any misidentified individuals?

Then run ref_map populations again cleanly: - for different population filters

```

```
*Link to output files here* 

## Re-filtering populations individually

*For population-level analyses*

