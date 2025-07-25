import argparse
import re

def main(vcf_files, sample, output_file):
	variant_dict = {}

	for filename in vcf_files:
		caller = re.split('\.|_', filename)[-2]
		with open(filename, 'r') as variants:
			for line in variants:
				var_info = line.rstrip().split('\t')
				var_id = (var_info[0], var_info[1], var_info[3], var_info[4])
				if var_id not in variant_dict:
					variant_dict[var_id] = [[caller], var_info[5:]]
				else:
					variant_dict[var_id][0].append(caller)

	all_variants = []

	for variant, info in variant_dict.items():
		info[1][2] = info[1][2] + ';CALLER=' + '_'.join(info[0])
		var_data = [variant[0],
			variant[1],
			f'{sample}_' + '_'.join(variant),
			variant[2],
			variant[3]] + \
			info[1]
		all_variants.append(var_data)

	autosomal_variants = [x for x in all_variants if \
		'X' not in x[0] and \
		'Y' not in x[0] and \
		'M' not in x[0]]

	sorted_variants = sorted(autosomal_variants, key=lambda x: [int(y) for y in [x[0][3:], x[1]]])

	with open(output_file, 'w') as out:
		for variant in sorted_variants:
			print('\t'.join(variant), file=out)

if __name__ == "__main__":
	# to run from the command line, use commented out lines
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-s", "--sample", required=True, help="de novo sample id")
	# ap.add_argument("-f", "--files", required=True, help="first input file path and name", nargs='+')
	# ap.add_argument("-o", "--output", required=True, help="path to output file")
	# args = ap.parse_args()
	#
	# main(args.files, args.sample, args.output)
	main(snakemake.input.vcfs, snakemake.params.sample, snakemake.output.merged_calls)
