#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --mem=60GB
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20

module load amazon-corretto/8.212.04.2
export _JAVA_OPTIONS="-Xmx32g -XX:+UseSerialGC"
module load R/3.6.1
module load python/3.7.2
module load mafft

INPUT=$1

DATE=$2 # This is used to name the files
SEQS=data/${DATE}.fasta #

# merge GISAID sequences with AAHL and EVAg sequences
cat $1 data/AAHL_sequences/* data/EVAg_download/ > $2

# Sanitize seq file
# Suzanne's webcrawler script seems to add <unkown description> to the end of each header which messes downstream stuff up
sed -i 's|<unknown description>||g' $SEQS

OUTDIR=results/${DATE}
mkdir $OUTDIR # this is where all output files will be stored

mkdir $OUTDIR/sequences
mkdir $OUTDIR/alignments
mkdir $OUTDIR/figures
mkdir $OUTDIR/other

FILTERED=${OUTDIR}/sequences/${DATE}.filtered.fastq   # sequences filtered for N content
ALIGNMENT=${OUTDIR}/alignments/${DATE}.mafftAlignment
GAPFREQ=${OUTDIR}/other/${DATE}.gapFrequncies.txt # file contain the gap frequncy at each position based on $ALIGNMENT
GAPFULL=${OUTDIR}/figures/${DATE}.gapCoverage.wholeGenome.pdf # pdf of gap coverage across full genome
GAP5=${OUTDIR}/figures/${DATE}.gapCoverage.5prime.pdf # pdf of gap coverage at 5'
GAP3=${OUTDIR}/figures/${DATE}.gapCoverage.3prime.pdf # pdf of gap coverage at 3'
TRIMALIGN=${OUTDIR}/alignments/${DATE}.trimmed.mafftAlignment
TRIMSEQS=${OUTDIR}/sequences/${DATE}.trimmed.fastq # trimmed sequences
NODUPALIGN=${OUTDIR}/alignments/${DATE}.trimmed.noDup.mafftAlignment
NODUPALIGNREC=${OUTDIR}/other/${DATE}.trimmed.noDup.mafftAlignment.identicalIDs.txt
NODUPSEQS=${OUTDIR}/sequences/${DATE}.trimmed.noDup.fasta # no duplicate sequences based on $TRIMSEQS
NODUPSEQSREC=${OUTDIR}/other/${DATE}.trimmed.noDup.fasta.identicalIDs.txt # record of identical seq ID's
CONSENSUSSEQ=${OUTDIR}/sequences/${DATE}.consensusSequence.fasta
CONFILSEQ=${OUTDIR}/sequences/${DATE}.filteredByConsensus.fasta
KMERSIG=${OUTDIR}/${DATE}_kmerSigs # kmer signatures based on $NODUPSEQS
MUTFREQUNC=${OUTDIR}/other/${DATE}.uncondensedMutationFreq.txt
MUTFREQCON=${OUTDIR}/other/${DATE}.condensedMutationFreq.txt

# Filter sequences by N content
perl code/bin/filterSeqsByNfreq.pl --input $SEQS --freq 1 --nlength 50 --output ${FILTERED}

# Run muscle alignment
# muscle -in $FILTERED -out $ALIGNMENT
mafft --thread 20 $FILTERED > $ALIGNMENT

# Count gap frequencies and plot
perl code/bin/countGapFrequencies.pl --input ${ALIGNMENT} --output $GAPFREQ
Rscript code/bin/plotGapCoverage.R $GAPFREQ $GAPFULL $GAP5 $GAP3

# Trim and sanatize alignment and sequences
perl code/bin/trimAndSanitize.pl --input $ALIGNMENT --freqs $GAPFREQ --alignout $TRIMALIGN --sequences $TRIMSEQS

# Collapse duplicates
perl code/bin/collapseIdenticalSequences.pl --input $TRIMALIGN --output $NODUPALIGN --record $NODUPALIGNREC
perl code/bin/collapseIdenticalSequences.pl --input $TRIMSEQS --output $NODUPSEQS --record $NODUPSEQSREC

# As a sanity check, $NODUPALIGNREC and $NODUPSEQSREC should have the same content

# Construct mutation frequency plots
python code/bin/extractMutationFrequency.py $NODUPALIGN ${OUTDIR}/other/${DATE}
Rscript code/bin/plotMutationFrequency.R $MUTFREQUNC $MUTFREQCON ${OUTDIR}/figures/${DATE}

# Construct consensus sequence and filter alignments
python code/bin/constructReference.py $NODUPALIGN $CONSENSUSSEQ
python code/bin/filterByHammingDistance.py --consensus $CONSENSUSSEQ --alignment $NODUPALIGN --threshold 30 --output $CONFILSEQ

# Construct kmer freqs
sbatch code/bin/KmerCounting.sh $CONFILSEQ $NODUPSEQSREC $KMERSIG

# submit script to create phylogeny
sbatch code/bin/estimateTree.q $NODUPALIGN ${DATE}
