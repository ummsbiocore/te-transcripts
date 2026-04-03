#!/usr/bin/env python3

import argparse

def main(args):

	output = []

	with open(args.input) as infile:
		header = infile.readline().rstrip().split('\t')
		output.append('gene\t%s' % (header[1].replace('_sorted.bam', '')))

		for line in infile:
			cur = line.rstrip().split('\t')
			output.append('%s\t%s' % (cur[0].strip('"'), cur[1]))

	with open(args.output, 'w') as out:
		out.write('\n'.join(output))

def parseArguments():
	parser = argparse.ArgumentParser(prog="clean_TE_count", description='', usage='%(prog)s [options]')	
	required = parser.add_argument_group('Input')
	required.add_argument('-i', '--input', required=True,  help='TE_Count count files.', metavar='', dest='input')
	required.add_argument('-o', '--output', required=True, help='Name of output file.', metavar='', dest='output')
	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)