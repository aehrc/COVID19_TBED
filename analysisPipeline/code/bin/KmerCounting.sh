#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --output=kmer_analysis-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=laurence.wilson@csiro.au

## Load HPC modules
module load amazon-corretto/8.212.04.2
export _JAVA_OPTIONS="-Xmx32g -XX:+UseSerialGC"
module load python/3.7.2

##### Not sure whether the above lines are still required.
##### Assuming that the master script runs as a batch script, then probably not?
##### Testing is required to see whether the modules they are needed.
##### I've overestimated the batch job resource requirement
##### but it ensures that the script runs to completion without errors.

## **********
## *** Kmer counting script directory
## *** THIS SHOULD NOT BE CHANGED
## **********
KMER_DIR=/home/wil9cq/scratch1/COVID19_TBED/analysisPipeline/code/bin

## Input file and output directory
## Read them as arguments from the master script
fastaFile=$1
dupIdsFile=$2
outputDir=$3
kmerFreqDir=freqs

## **********
## *** Script Setup
## **********
## *** Fix up the format of some FASTA headers
## Remove trailing whitespaces
sed 's/^> \+/>/g' ${fastaFile} > tmp.fa
sed -i 's/^>\(.*\) \+$/>\1/g' tmp.fa
sed -i 's/hCoV-19\///g' tmp.fa

## *** Rename some FASTA headers
sed -i 's/BetaCoV\///g' tmp.fa
sed -i 's/ SARS-CoV-2.*//g' tmp.fa
sed -i 's/\/[0-9]\+|/|/g' tmp.fa
sed -i 's/ /_/g' tmp.fa
sed -i 's/\/2020$//g' tmp.fa
sed -i 's/\/2020_/_/g' tmp.fa

## Run Kmer counting
python ${KMER_DIR}/calculate_kmer_frequencies.py \
    split \
    -f tmp.fa \
    -k 10 \
    -n \
    -o ${outputDir}/${kmerFreqDir}

## Run Kmer analysis
python ${KMER_DIR}/analyse.py \
    ${kmerFreqDir} \
    ${dupIdsFile} \
    ${outputDir}

## Clean-up
rm tmp.fa


