import argparse
import numpy as np
from pysam import FastaFile
from collections import Counter

def get_var_data(file_line, header_list):
	var_info = file_line.rstrip().split('\t')

	return dict(zip(header_list, var_info))

def get_sequence(chr, pos, ref_fasta):
	start = pos - 6
	end = pos + 5

	return ref_fasta.fetch(chr, start=start, end=end).upper()

def check_homopolymer(sequence):
	start = 0
	window_size = 5

	base_counts = Counter(sequence[:window_size])

	while start < len(sequence) - window_size:
		if window_size in base_counts.values():
			return True, next(k for k,v in base_counts.items() if v > 0)

		base_counts[sequence[start]] -= 1
		base_counts[sequence[start + window_size]] += 1

		start += 1

	return False, None

def validate_homopolymer(var_data, subunit):
	ab = [float(var_data[f'{platform}_ab']) for platform in ['hifi', 'ont', 'illumina']]
	ab = np.nan_to_num(ab)

	if var_data['ref'] == subunit or var_data['alt'] == subunit:
		if ab[2] == 0:
			return False

		if int(var_data['others_with_allele_bq20']) > 0 or int(var_data['others_with_allele_bq10']) > 0:
			return False

	return True

def main(variant_file, fasta_file, output_file):
	reference_fasta = FastaFile(fasta_file)

	with open(variant_file, 'r') as variants, open(output_file, 'w') as outfile:
		header = variants.readline().rstrip().split('\t')
		header.append('in_homopolymer')
		print('\t'.join(header), file = outfile)
		for line in variants:
			variant = get_var_data(line, header[:-1])
			flanking_sequence = get_sequence(variant['chr'], int(variant['pos']), reference_fasta)
			in_homopolymer, subunit = check_homopolymer(flanking_sequence)

			if in_homopolymer and not validate_homopolymer(variant, subunit):
				continue

			out_data = line.rstrip() + f"\t{in_homopolymer}"
			print(out_data, file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="file of de novo variants")
	# ap.add_argument("-f", "--fasta", required=True, help="fasta file")
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()

	# main(args.variants, args.fasta, args.output)
	main(snakemake.input.variants, snakemake.params.fasta, snakemake.output.filtered_snvs)



