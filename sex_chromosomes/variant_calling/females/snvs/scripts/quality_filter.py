import argparse

def get_gq(sample_data, idx):
	"""
	Fetches the genotype quality from a colon separated string of VCF data.
	:param sample_data: A string of one sample's data for one variant from a VCF.
	:param idx: Index of the GQ value.
	:returns: Integer value for genotype quality.
	"""
	return int(sample_data.split(':')[idx])

def check_child_gq(variant_data, idx, threshold):
	"""
	Determines if the child's genotype quality passes the threshold.
	:param variant_data: List of column values for a variant in VCF format.
	:param idx: Index of the GQ value in the genotype column.
	:param threshold: Minumum threshold required for genotype quality.
	:returns: True if child GQ above threshold, False if not.
	"""
	if get_gq(variant_data[11], idx) > threshold:
		return True

	return False

def main(variant_file, snv_output, indel_output):
	genotype_quality = 20

	snvs = open(snv_output, 'w')
	indels = open(indel_output, 'w')

	with open(variant_file, 'r') as variants:
		for line in variants:
			var_info = line.rstrip().split('\t')
			gq_idx = var_info[8].split(':').index('GQ')

			if check_child_gq(var_info, gq_idx, genotype_quality):
				if len(var_info[3]) != 1 or len(var_info[4]) != 1:
					if len(var_info[3]) == len(var_info[4]):
						raise ValueError("you have MNMs in your data. go back and normalize your VCF with bcftools please!")
					print(line.rstrip(), file=indels)
				else:
					print(line.rstrip(), file=snvs)

	snvs.close()
	indels.close()

if __name__ == "__main__":
	# to run from the command line, use commented out lines
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="input file")
	# ap.add_argument("-s", "--snv_out", required=True, help="snv output file")
	# ap.add_argument("-d", "--indel_out", required=True, help="indel output file")
	# args = ap.parse_args()

	# main(args.input, args.genotype_quality, args.snv_out, args.indel_out)
	main(snakemake.input.merged_calls, snakemake.output.snv_calls, snakemake.output.indel_calls)
