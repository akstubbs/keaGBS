# PRACTICE SNPcalling.md

The raw data is available on the High Capacity Storage of Otago University. 

Contact akstubbs.nz@gmail.com for access

## Quality control

The data is single end __ across two lanes. 

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
Looking at the html output, sequence quality is okay but quite some error in lane 2 first few bases and last few.

I will first remove the last bases in each lane as there is a lot of adapter contamination. 
As we are using a reference genome for kea, we do not need to have a common read length. 
Shorter reads are removed by setting minimum number of bp to 40. 

```
#!/bin/sh
module load cutadapt

cutadapt -a AGATCGGAAGAGC -m 40 -o trim_practice1.fastq SQ1609_CD82MANXX_s_6_fastq.txt.gz 

fastqc trim_practice1.fastq
```
#### *Just specifically for practice*

*I then created a practice folder specifically for this code.*

```
#!/bin/sh
mkdir practice
```
*The created fastqc files for practicing code were moved into the practice folder:*

```
#!/bin/sh
mv lane* practice/
mv trim* practice/
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
cat SQ1609.txt | grep -E "K38" | cut -f 3-4 > barcodes_practice1.txt
```
Then I create different folders to deal with samples sequenced across the two different lanes before concatenating them together.

//

*below code for when have both data plates - will need to change above for actual analysis*
```
#!/bin/sh
mkdir raw1 samples1 raw2 samples2 samples_concatenated
cd raw1
ln -s ../source_files/trim_plate_1.fastq
cd ..

cd raw2
ln -s ../source_files/trim_plate_2.fastq
cd ..
```
//
I am now ready to run process_radtags to demultiplex.

```
#!/bin/sh
process_radtags -p raw1/ -o ./samples1/ -b source_files_kea/barcodes_lane1.txt --renz_1 pstI --renz_2 mspI -r -c -q --inline_null

process_radtags -p raw2/ -o ./samples2/ -b source_files_kea/barcodes_lane1.txt --renz_1 pstI --renz_2 mspI -r -c -q --inline_null (inline??**check the --inline code)
```

^^for the above need to include correct process_radtags for both lanes/plates, THEN concatenate samples. 

```
#python
import os
#
for sample in os.listdir ("samples2"):
  print sample
  if sample.endswith("gz") and not sample in os.listdir("samples1"):
    raise Exception

for sample in os.listdir ("samples2"):
  print sample
  if sample.endswith("gz") and not sample in os.listdir("samples1"):
    os.system("cat samples1/"+sample+" samples2/"+sample+" > samples_concatenated/"+sample)
```

Then we will have clean sample files and ready for alignment.

#### *Demultiplexing Just specifically for practice* 

*Instead of working inside created folders (i.e. raw2, samples2) this was done in the practice folder*

```
#!/bin/sh
cd uoo03341/practice/
process_radtags -p practice/ -o practice/ -b source_files_kea/barcodes_practice1.txt  --renz_1 pstI --renz_2 mspI -r -c -q --inline_null
```

*Looks good, keeping nearly 97% of reads*

```
Outputing details to log: 'practice/process_radtags.practice.log'

297351 total sequences
  9019 barcode not found drops (3.0%)
   441 low quality read drops (0.1%)
   736 RAD cutsite not found drops (0.2%)
287155 retained reads (96.6%)
```

### Alignment and variant calling

First alignment for every sample using BWA. Here is one example command:

```
#!/bin/sh
bwa mem -t 8 $bwa_db $src/${sample}.fq.gz | samtools view -b | samtools sort --threads 4 > ${sample}.bam
```

The complete list of commands is in [realign.sh](/keaGBS/blob/main/realign.sh):

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

```
#!/bin/sh
module load Stacks
ref_map.pl --samples ./practice/alignment --popmap popmap_2_NCS.txt -T 8  -o ./practice/output_refmap_NCS
```

Now check samples and output files for low quality individuals with low sample numbers. 
Remove these, the blank (if it hasn't already been removed?) and any misidentified individuals?

Then run ref_map again cleanly:

```

```
*Link to output files here* 

## Re-filtering populations individually

*For population-level analyses*

