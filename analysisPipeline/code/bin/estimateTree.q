#!/bin/bash
#SBATCH --ntasks 30 
#SBATCH --mem 60GB
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=4
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=laurence.wilson@csiro.au

module load raxml-ng/0.9.0
alignmentFile=$1 #the alignment file
originalAlignmentFilementFile=$alignmentFile

#remove spaces
sed -i 's/ /_/g' $alignmentFile

raxml-ng --parse --msa $alignmentFile --model GTR+G+I
#this may pruduce a compressed version of the alignment file
if [ -f "$alignmentFile.raxml.reduced.phy" ]; then
    alignmentFile=$alignmentFile.raxml.reduced.phy
    raxml-ng --parse --msa $alignmentFile --model GTR+G+I
fi

BalignmentFile=$alignmentFile.raxml.rba #binary multiple sequence alignment file

OUTDIR=results/${2}/${2}_trees
mkdir $OUTDIR

echo check 1

for i in `seq 1 10`; 
do
  srun -N 1 -n 1 --exclusive /apps/raxml-ng/0.9.0/bin/raxml-ng --search --msa $BalignmentFile --tree pars{100} --prefix ${OUTDIR}/CT$i --seed $RANDOM --threads 4 &
done

for i in `seq 11 20`; 
do
  srun -N 1 -n 1 --exclusive /apps/raxml-ng/0.9.0/bin/raxml-ng --search --msa $BalignmentFile --tree rand{100} --prefix ${OUTDIR}/CT$i --seed $RANDOM --threads 4 &
done

for i in `seq 1 10`
do
    srun -N 1 -n 1 --exclusive /apps/raxml-ng/0.9.0/bin/raxml-ng --bootstrap --msa $BalignmentFile --bs-trees 100 --prefix ${OUTDIR}/CB$i --seed $RANDOM --threads 4 &
done
wait

echo check 2

echo final fixup

cd ${OUTDIR}
#sorts through logs to find the best tree
bestTree=$(grep "Final LogLikelihood" CT*.raxml.log | sort -k 3 | awk -F"." 'NR==1{print $1;}')
bestTree=$bestTree.raxml.bestTree
#appends all bootstrapped trees into one file
cat CB*.raxml.bootstraps > allbootstraps
#finds bootstrap support of the best tree
raxml-ng --support --tree $bestTree --bs-trees allbootstraps --prefix CS --threads 1 --redo
#copies out tree with support values
cd ../../
cp ${OUTDIR}/CS.raxml.support $originalAlignmentFilementFile.nex

#because trees only bifurcate there is a minimum branch length between identical strains
#set this length to 0
sed -i 's/0.000001/0/g' $originalAlignmentFilementFile.nex


