import argparse
import gzip
import re
from tqdm import tqdm

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	'''
	alist.sort(key=natural_keys) sorts in human order
	http://nedbatchelder.com/blog/200712/human_sorting.html
	(See Toothy's implementation in the comments)
	'''
	return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def sort_input_files(file_list):
	return sorted(file_list, key = natural_keys)

def format_output(chromosome, pos1, pos2):
	return '\t'.join([chromosome, str(pos1), str(pos2)])

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, nargs = "+", help="input files")
	# ap.add_argument("-o", "--output", required=True, type=str, help="output file")

	# args = ap.parse_args()

	callable_files = sort_input_files(snakemake.input.callable_sites)
	# callable_files = sort_input_files(args.input)

	out_file = open(snakemake.output.full_bed, 'w')
	# out_file = open(args.output, 'w')

	for file in tqdm(callable_files, total = len(callable_files)):
		with gzip.open(file, 'rt') as chr_file:
			header = chr_file.readline()

			try:
				first_site = chr_file.readline().split('\t')
				prev_pos = starting_pos = pos = int(first_site[1])
				chrom = first_site[0]
			except IndexError:
				continue

			for line in chr_file:
				pos = int(line.split('\t')[1])
				if pos != prev_pos + 1:
					ending_pos = prev_pos
					print(format_output(chrom, starting_pos, ending_pos), file = out_file)
					starting_pos = pos

				prev_pos = pos

			print(format_output(chrom, starting_pos, pos), file = out_file)

	out_file.close()
