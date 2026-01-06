#!/usr/bin/env python3

import argparse
from collections import defaultdict

def main(args):

	counts = defaultdict(int)
	genes = set()
	samples = set()

	for input_file in args.input_files:

		with open(input_file) as infile:
			header_key = {}
			header = infile.readline().rstrip().split('\t')
			for i, column in enumerate(header[1:]):
				header_key[i+1] = column.rstrip('.bam.C').rstrip('.bam.T')
				samples.add(column.rstrip('.bam.C').rstrip('.bam.T'))

			for line in infile:
				cur = line.rstrip().split('\t')
				genes.add(cur[0].strip('"'))
				for i, col in enumerate(cur[1:]):
					counts[(header_key[i+1], cur[0].strip('"'))] = int(col)

	output = ['gene\t%s' % '\t'.join(sample for sample in sorted(samples))]
	for gene in sorted(genes):
		output.append('%s\t%s' % (gene, '\t'.join(['%s' % counts[(sample, gene)] for sample in sorted(samples)])))

	with open(args.output, 'w') as out:
		out.write('\n'.join(output))

def parseArguments():
	parser = argparse.ArgumentParser(prog="prepare_comparisons", description='', usage='%(prog)s [options]')	
	required = parser.add_argument_group('Input')
	required.add_argument('-i', '--input-files', required=True, nargs='+', help='List of TE_Transcript count files.', metavar='', dest='input_files')
	required.add_argument('-o', '--output', required=True, help='Name of output file.', metavar='', dest='output')
	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)