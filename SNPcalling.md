# SNPcalling.md

The raw data is available on the High Capacity Storage of Otago University. 

Contact akstubbs.nz@gmail.com for access.

## Quality control

The data is single end __ across two lanes. 

- Plate 1 / lane1 = SQ1609...6_fastq.txt.gz
- Plate 2 / lane2 = SQ1630...7_fastq.txt.gz

The adapter content and the barcodes were subset to understand the structure before being used on whole raw data file:

```
#!/bin/sh
# in source_files_kea

ls
less Kea_conversion_table_plate1.txt

zcat SQ1609_CD82MANXX_s_6_fastq.txt.gz | head -n 1000000 > lane1_sample.fq
zcat SQ1630_CD9F4ANXX_s_7_fastq.txt.gz | head -n 1000000 > lane2_sample.fq
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
mkdir -p kea_ref_genome
```
Copied kea reference genome into kea_ref_genome folder and unzipped file

```
#!/bin/sh
curl -L -o kea_ref_genome/kea_ref_genome.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/696/875/GCA_000696875.1_ASM69687v1/GCA_000696875.1_ASM69687v1_genomic.fna.gz

gunzip kea_ref_genome/kea_ref_genome.fasta.gz
```
Load BWA and index reference

```
#!/bin/sh
module load BWA

bwa index kea_ref_genome.fasta
```

## SNP calling

### Demultiplexing

First I extract barcodes from the .key file (SQ1609.txt and SQ1630.txt) of the sequencing platform.

```
#!/bin/sh
##Key provided in folder
cat SQ1609.txt | grep -E "K38" | cut -f 3-4 > barcodes_plate_1.txt
cat SQ1630.txt | grep -E "K38" | cut -f 3-4 > barcodes_plate_2.txt
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

process_radtags -p raw1/ -o samples1/ -b source_files_kea/barcodes_plate_1.txt  --renz_1 pstI --renz_2 mspI -r -c -q --inline_null

process_radtags -p raw2/ -o samples2/ -b source_files_kea/barcodes_plate_2.txt  --renz_1 pstI --renz_2 mspI -r -c -q --inline_null
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

Make directory for bam files to be kept in:

```
#!/bin/sh
mkdir bam_files_ALLsamples
```

First alignment for every sample using BWA. Here is one example command:

```
#!/bin/sh
bwa mem -t 8 $bwa_db $src/${sample}.fq.gz | samtools view -b | samtools sort --threads 4 > ${sample}.bam
```

The complete list of commands is in [realign.sh](realign.sh):

Ran script in realign.sh

```
#!/bin/sh
cd bam_files_ALLsamples
bash ../realign.sh
```

I created a population map file for the three haplotype regions of wild kea (North, Central, South) and the captive population (unknown haplotype region).

```
#!/bin/sh
cd nano popmap_haplotype.txt
<sample file prefix><tab><population ID>

<K38*><tab><north/central/south/?>
```
I created a new folder for the output reference map:

```
#!/bin/sh
mkdir output_refmap_haplotype
```

Then I run refmap from Stacks quick run to identify low quality individuals:
**run ref_map.pl as a job**

```
#!/bin/sh
ref_map.pl --samples bam_files_ALLsamples/ --popmap popmap_haplotype.txt -T 8 -o output_refmap_haplotype/
```

Now check samples and output files for low quality individuals with low sample numbers. 
Remove these, the blank (if it hasn't already been removed?) and any misidentified individuals?

### Re-run ref_map populations again cleanly:

```
#!/bin/sh
module load Stacks

mkdir output_pop_haplotype

populations -P output_refmap_haplotype/ -O output_pop_haplotype/ -M popmap_haplotype.txt -p 1 --write-single-snp --plink --vcf --verbose

Removed ...
```
*Link to output files here* 

## Remove individuals with high missing data

VCFtools to look at missing data per individual. Make a new directory for copying VCF file after running populations
```
#!/bin/sh
mkdir filtering_haplotype
cp output_pop_haplotype/populations.snps.vcf filtering_haplotype/
```

Run vcftools to output a file that contains missingness on an individual level.

```
#!/bin/sh
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
vcftools --vcf populations.snps.vcf --missing-indv
```
Downloaded ".imiss" file and opened in Excel, sorted individuals by the most missing data (F_MISS). Created a scatterplot to visualise data (F_MISS x N_GENOTYPES-FILTERED)

Four individuals >0.70 F_MISS threshold:

K38612  K38699  K38682  K38689

  **We can remove these using the following code. However - this was not done, so we can check the filtering for all individuals. Individuals can be removed later on, after checking other filtering.**

  *Vcftools run again to output a new VCF file, which removes individuals with >0.70 missing data. New VCF file renamed appropriately. And file moved into appropriate new folder for future analysis to take place and file storage.*

```
#!/bin/sh
vcftools --vcf populations.snps.vcf --remove-indv K38612 --remove-indv K38699 --remove-indv K38682 --remove-indv K38689 --recode
mv out.recode.vcf removed_4ind_70-FMISS.vcf

mkdir filtering_70-FMISS_rm4ind
mv removed_4ind_70-FMISS.vcf filtering_70-FMISS_rm4ind/
```

## Filtering

VCFtools to output the mean depth per individual:

```
#!/bin/sh
vcftools --vcf populations.snps.vcf --depth
```
Generated 'out.idepth' file. Individual with highest depth, and doubled it to obtain maximum depth parameter. ```34```

Remove SNPs with depth <2 and >34, and 80% missing data

```
#!/bin/sh
vcftools --vcf ../populations.snps.vcf --minDP 2 --maxDP 34 --max-missing 0.80 --recode
After filtering, kept 27652 out of a possible 92797 Sites
```
Rename new vcf to suitable name:

```
#!/bin/sh
mv out.recode.vcf filtered_dp2_34_md80.vcf
```
Tried for different levels of filtering - depth <2 or <3, and 70% or 80% missing data :

```
#!/bin/sh
vcftools --vcf ../populations.snps.vcf --minDP 3 --maxDP 34 --max-missing 0.80 --recode
After filtering, kept 22944 out of a possible 92797 Sites
mv out.recode.vcf filtered_dp3_34_md80.vcf

vcftools --vcf ../populations.snps.vcf --minDP 3 --maxDP 34 --max-missing 0.70 --recode
After filtering, kept 27857 out of a possible 92797 Sites
mv out.recode.vcf filtered_dp3_34_md70.vcf

vcftools --vcf ../populations.snps.vcf --minDP 2 --maxDP 34 --max-missing 0.70 --recode
After filtering, kept 32616 out of a possible 92797 Sites
mv out.recode.vcf filtered_dp2_34_md70.vcf
```
Number of SNPs kept compared. 

Stringent filtering (3-34, 0.8 = 22944), less stringent (3-34, 0.7 = 27857). Difference ~6000.

Stringent filtering (2-34, 0.8 = 27652), less stringent (2-34, 0.7 = 32616). Difference ~5000.

**Stringent filtering 3min 34max 0.8 missing ** -> still keeps 23,000 SNPs

## PCA

### Plink

Load Plink. Give Plink the vcf (MinDP 3, missingness 0.80) and covert it into PLINK binary files (.bim,.bed,.fam). 
```
#!/bin/sh
module load PLINK/1.09b6.16

plink --vcf filtered_dp3_34_md80.vcf --allow-extra-chr --double-id --make-bed --out pca_filtered_dp3_34_md80
```
Create pca output files:
```
#!/bin/sh
plink --bfile pca_filtered_dp3_34_md80 --allow-extra-chr --pca --out pca_filtered_dp3_34_md80
```
Downloaded two output files:
- pca_filtered_dp3_34_md80.eigenvec
- pca_filtered_dp3_34_md80.eigenval

Plotting completed in R.

### Plotting PCA - detect population structure

Install appropriate R packages:

```
install.packages("tidyverse")
install.packages("readr")

library(tidyverse)
library(readr)
```

#### Base R PCA with sample names

**PC1xPC2**

```
eigenvec_table <- read.table('pca_filtered_dp3_34_md80.eigenvec')
plot(eigenvec_table[,3], eigenvec_table[,4], xlab = "PC1", ylab = "PC2", col = "seagreen", cex=0.50)
```

PCA with sample names in text for PC1xPC2
```
plot(eigenvec_table[,3], eigenvec_table[,4], xlab = "PC1", ylab = "PC2", col = "seagreen", cex=0)
text(eigenvec_table[,3], eigenvec_table[,4], eigenvec_table[,2], cex=0.50, xlim=c(-0.25,0.4))
```

<img src="https://user-images.githubusercontent.com/85653223/125231404-c818c780-e32e-11eb-82a0-6025e73fa5ff.png" width=75% height=75%>


**PC2xPC3**

```
plot(eigenvec_table[,4], eigenvec_table[,5], xlab = "PC2", ylab = "PC3", col = "seagreen", cex=0)
text(eigenvec_table[,4], eigenvec_table[,5], eigenvec_table[,2], cex=0.50, xlim=c(-0.25,0.4))
```
PCA with sample names in text for PC2xPC3

```
plot(eigenvec_table[,5], eigenvec_table[,6], xlab = "PC3", ylab = "PC4", col = "seagreen", cex=0)
text(eigenvec_table[,5], eigenvec_table[,6], eigenvec_table[,2], cex=0.35, xlim=c(-0.25,0.4))
```

### PCA with ggplot, & variances per PC

Read in data

```
pca <- read_table2("./pca_filtered_dp3_34_md80.eigenvec", col_names = FALSE)
eigenval <- scan("./pca_filtered_dp3_34_md80.eigenval")
```
Cleaning data up: remove nuisance column (had 2 col of ind names FID/ID)
```
pca <- pca[,-1]

#Set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

#Convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
```
Plot of PC variances
```
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_classic()
```

<img src="https://user-images.githubusercontent.com/85653223/125231985-d3b8be00-e32f-11eb-93ee-2f9e3d6fe06b.png" width=50% height=50%>


**plot PC1 x PC2 with variance**

I wanted to be able to identify the samples by their haplotype region (plus unknown captive) as mentioned above for the population map (ref_map and populations). 

Therefore, I extracted the pop_map file as a .csv file, and converted the V2 column (haplotypes) into factor variables that can be used to change the colour of the plot by haplotype group.

```
data = read.csv(file="pop test.csv", sep=',', header=F)
Haplotype <- as.factor(data$V2)
```

Final code for plot:

```
b <- ggplot(pca, aes(PC1, PC2, color = Haplotype)) + geom_point(size = 1)
b <- b + coord_equal() + theme_classic()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
```

<img src="https://user-images.githubusercontent.com/85653223/125230662-5db35780-e32d-11eb-9c13-a6379009bf2f.png" width=70% height=70%>


**plot PC2 x PC3**
```
c <- ggplot(pca, aes(PC2, PC3, color = Haplotype)) + geom_point(size = 1)
c <- c + coord_equal() + theme_classic()
c + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
```
**plot PC3 x PC4**
```
d <- ggplot(pca, aes(PC3, PC4, color = Haplotype)) + geom_point(size = 1)
d <- d + coord_equal() + theme_classic()
d + xlab(paste0("PC3 (", signif(pve$pve[3], 3), "%)")) + ylab(paste0("PC4 (", signif(pve$pve[4], 3), "%)"))
```

### PCA plots excluding captive kea samples

Created file containing all captive kea sample names: kea_sample_IDs_captive.txt

Re-run Plink with flag --remove-fam to exclude samples from captive kea file.

```
plink --vcf filtered_dp3_34_md80.vcf --remove-fam ../kea_sample_IDs_captive.txt --allow-extra-chr --double-id --make-bed --out pca_filtered_dp3_34_md80_NoCaptiveInd
```

Some of the output from the plink log - showing that all the captive individuals (N=48) were successfully removed. 139 individuals remaining.
``` 
187 people (0 males, 0 females, 187 ambiguous) loaded from .fam.
Ambiguous sex IDs written to pca_filtered_dp3_34_md80_NoCaptiveInd.nosex .
--remove-fam: 139 people remaining.
```

Create pca output files:
```
#!/bin/sh
plink --bfile pca_filtered_dp3_34_md80_NoCaptiveInd --allow-extra-chr --pca --out pca_filtered_dp3_34_md80_NoCaptiveInd
```
Downloaded two output files:
- pca_filtered_dp3_34_md80_NoCaptiveInd.eigenvec
- pca_filtered_dp3_34_md80_NoCaptiveInd.eigenval

Plotting completed in R.


## Structure stuff ??

Created a whitelist.txt file to generate a list of 1000 random loci. 
```
#!/bin/sh
cat output_refmap_haplotype/populations.sumstats.tsv | grep -v "#" | cut -f 1 | sort | uniq | shuf | head -n 1000 > whitelist.txt
```

Tried to load structure but missing a mainparams file...??
