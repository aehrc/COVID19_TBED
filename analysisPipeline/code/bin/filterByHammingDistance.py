#!/usr/bin/python

import sys, getopt
import argparse

def hamdist(str1, str2):
	diffs = 0
	for ch1, ch2 in zip(str1, str2):
		if ch1 != ch2:
			diffs += 1

	return diffs

def main(argv):
	parser = argparse.ArgumentParser(description='Filter sequence based on hemming distance from supplied reference')
	parser.add_argument('--consensus', help='Consensus sequence fasta file', required=True)
	parser.add_argument('--alignment', help='MAFFT alignment file for processing', required=True)
	parser.add_argument('--threshold', help='Threshold % for filtering out sequences', required=True)
	parser.add_argument('--output', help='Output for retained sequences', required=True)
	args = parser.parse_args()

	refname = ''
	refseq = ''
	with open(args.consensus, 'r') as f:
		for line in f:
			line = line.rstrip()
			if line[0] == ">":
				refname = line
			else:
				refseq = line

	total = int(len(refseq))

	with open(args.alignment, 'r') as f, open(args.output, 'w') as outfile:
		header = ''
		seq = ''

		for line in f:
			line = line.rstrip()
			if line[0] == ">":
				header = line
			else:
				seq = line

			if header and seq:
				diff = hamdist(refseq, seq)
				perc = (diff/total)*100

				if perc < int(args.threshold):
					outfile.write('{}\n{}\n'.format(header,seq.replace('-','')))

				header = ''
				seq = ''

if __name__ == "__main__":
	main(sys.argv[1:])

