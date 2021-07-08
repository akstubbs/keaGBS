# Kea Sex Chromosomes

cd ../source_files_kea/
ls
makeblastdb -in kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -dbtype nucl
module load blast
module spide blast
module load BLAST
makeblastdb -in kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -dbtype nucl
mkdir kakapo
ls
mv kakapo_gen* kakapo/
ls
cd kakapo/
ls
ls ../ref_genome/
blastx -query ../ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -outfmt 7 > blast_results.txt
head kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna
ls
cp kakapo_genome_GCF_004027225-Copy1.1_bStrHab1_v1.p_genomic.fna kakapo_genome_GCF_004027225.fa
ls
head kakapo_genome_GCF_004027225.fa 
blastx -query ../ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.fa -outfmt 7 > blast_results.txt
ls
blastx -query ../ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.fa -outfmt 7 > blast_results.txt
blastn -query ../ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.fa -outfmt 7 > blast_results.txt
makeblastdb -in kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -dbtype nucl
blastn -query ../ref_genome/kea_ref_genome.fasta -db kakapo_genome_GCF_004027225.1_bStrHab1_v1.p_genomic.fna -outfmt 7 > blast_results.txt

//

cd 
cd uoo03341/source_files_kea/kakapo/
ls
less blast_results.txt 
grep NC044302.1 blast_results.txt
less blast_results.txt 
grep NC044301.1 blast_results.txt
less blast_results.txt 
grep NC_044301.1 blast_results.txt
grep NC_044301.1 blast_results.txt | cut -f 1-2
grep NC_044301.1 blast_results.txt | cut -f 1-2 | sort | uniq
module load SAMtools
samtools faidx ../ref_genome/kea_ref_genome.fasta
ls
cd ../ref_genome/
ls
less kea_ref_genome.fasta.fai 
grep JJRH01000001.1 kea_ref_genome.fasta -A 10
x="TACCAGCACAATTAGCAGGCAGGtttcaagggtattcaaaggccatccaAAACCTTTAGAACTCTCAAATGGTACTGTAA
ATAGCATGGTGGGGAAGCGGGGGGTATATGGGAATATCTCCTTCATAGGTTGACCcccaaagggaaggaaaaagaaagat
gtttaatTGCTAAAGAATTCCATTATATGTGTCCCAAAGTATAGAGGCAATGACCACTCTGAGAACACGTGCCAGCTCAG
TTTCATaatcaatgagtctattgtatta"
wc -c $x
grep [a-zA-Z] -c TACCAGCACAATTAGCAGGCAGGtttcaagggtattcaaaggccatccaAAACCTTTAGAACTCTCAAATGGTACTGTAA
ATAGCATGGTGGGGAAGCGGGGGGTATATGGGAATATCTCCTTCATAGGTTGACCcccaaagggaaggaaaaagaaagat
gtttaatTGCTAAAGAATTCCATTATATGTGTCCCAAAGTATAGAGGCAATGACCACTCTGAGAACACGTGCCAGCTCAG
TTTCATaatcaatgagtctatt
grep [a-zA-Z] -c $x
echo "TACCAGCACAATTAGCAGGCAGGtttcaagggtattcaaaggccatccaAAACCTTTAGAACTCTCAAATGGTACTGTAA
ATAGCATGGTGGGGAAGCGGGGGGTATATGGGAATATCTCCTTCATAGGTTGACCcccaaagggaaggaaaaagaaagat
gtttaatTGCTAAAGAATTCCATTATATGTGTCCCAAAGTATAGAGGCAATGACCACTCTGAGAACACGTGCCAGCTCAG
TTTCATaatcaatgagtctattgtatta" > test
wc -c test
less kea_ref_genome.fasta.fai 
grep NC_044301.1 blast_results.txt | cut -f 1-2 | sort | uniq > potential_sex_chr.scaffolds.txt
cd ../kakapo/
grep NC_044301.1 blast_results.txt | cut -f 1-2 | sort | uniq > potential_sex_chr.scaffolds.txt
