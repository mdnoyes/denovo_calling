import argparse
import gzip
import re
from tqdm import tqdm
from pybedtools import BedTool

def get_callable_bed(callable_file, chromosome_file, outfile):
	callable_bed = BedTool(callable_file)
	chromosome_bed = BedTool(chromosome_file)
	intersection = callable_bed.intersect(chromosome_bed)
	intersection.saveas(outfile)

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

def main(input_files, sex_bed, temp_file, output_file):
#	callable_files = sort_input_files(input_files)
#	temp_bed = open(temp_file, 'w')

#	for file in tqdm(callable_files, total = len(callable_files)):
#		with gzip.open(file, 'rt') as chr_file:
#			header = chr_file.readline()
#
#			try:
#				first_site = chr_file.readline().split('\t')
#				prev_pos = starting_pos = pos = int(first_site[1])
#				chrom = first_site[0]
#			except IndexError:
#				continue
#
#			for line in chr_file:
#				pos = int(line.split('\t')[1])
#				if pos != prev_pos + 1:
#					ending_pos = prev_pos
#					print(format_output(chrom, starting_pos, ending_pos), file = temp_bed)
#					starting_pos = pos
#
#				prev_pos = pos
#
#			print(format_output(chrom, starting_pos, pos), file = temp_bed)
#
#	temp_bed.close()

	get_callable_bed(temp_file, sex_bed, output_file)


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, nargs = "+", help="input files")
	ap.add_argument("-b", "--bed", required=True, type=str, help="output file")
	ap.add_argument("-t", "--temp", required=True, type=str, help="output file")
	ap.add_argument("-o", "--output", required=True, type=str, help="output file")
	
	args = ap.parse_args()
	
	main(args.input, args.bed, args.temp, args.output)
	#main(snakemake.input.callable_sites, snakemake.params.sex_bed, snakemake.output.temp_bed, snakemake.output.full_bed)


