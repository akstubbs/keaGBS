#!/bin/bash -e
#SBATCH -t 12:00:00
#SBATCH --account=uoo03341
#SBATCH --job-name=LastZ
#SBATCH --partition=large
#SBATCH --mem=3000MB
#SBATCH --hint=nomultithread

### This script will run lastz to align to genomes with each other. 
### The analysis is run separately for each chromosome
### Stefanie Grosser, 23/03/2020, University of Otago, Dunedin
### To run script: sbatch alignGenomes_LastZ.sh ${query} ${target_path} ${target} ${out_path} ${query_species} ${target_species}


echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
date
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo


## Load required modules
module load LASTZ/1.02.00-gimkl-2018b


## Import arguments from command line
query_chromosome_file=$1 	# A file containing the names of the individual chromosome fasta files
target_path=$2				# Path to the genome file against which to align the chromosomes
target=$3						# Genome file
out_path=$4					# Output path for analysis
query_species=$5			# Species of the query sequence
target_species=$6			# Species of the target sequence

chr_name=$(basename "$query_chromosome_file" ".fa" | awk -F "_" '{print $NF}')

## The alala reference genome contains lower case bases which are not marking low complexity regions but stem from the assembly process with falcon
## To avoid that LastZ recognises these lower case letters include the [unmask] action following the file name.

## Run lastZ

echo "lastz ${query_chromosome_file}[unmask] ${target_path}${target}[unmask]... --output=${out_path}LastZalignment_${query_species}_${target_species}_${chr_name}.fasta.txt"
 lastz ${query_chromosome_file}[unmask] ${target_path}${target}[unmask] \
        --masking=254 \
        --hspthresh=4500 \
        --gappedthresh=300 \
        --ydrop=15000 \
        --gapped \
        --chain \
        --seed=12of19 \
        --notransition \
        --matchcount=5000 \
        --format=general:name1,start1,end1,name2,start2,end2 \
        --ambiguous=n \
        --ambiguous=iupac \
        --output=${out_path}LastZalignment_${query_species}_${target_species}_${chr_name}.txt

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
date
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo