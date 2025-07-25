import argparse

def get_indel_size(variant_info):
	size = len(var_info[4]) - len(variant_info[3])
	return str(size)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	# ap.add_argument("-o", "--output", required=True, type=str, help="family id")
	# args = ap.parse_args()

	out_header = ['chr', 'pos', 'id', 'ref', 'alt', 'caller', 'indel_size']

	# with open(args.output, 'w') as outfile, open(args.input, 'r') as infile:
	with open(snakemake.output.indel_calls, 'w') as outfile, open(snakemake.input.vcf, 'r') as infile:
		print('\t'.join(out_header), file = outfile)
		infile.readline()
		for line in infile.readlines():
			var_data = line.rstrip().split('\t')
			var_info = var_data[:5] + [var_data[7].split(';')[-1].split('=')[-1]]
			indel_size = get_indel_size(var_info)
			var_info.append(indel_size)
			print('\t'.join(var_info), file = outfile)
