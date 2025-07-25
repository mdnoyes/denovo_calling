import argparse
import pandas as pd
from pybedtools import BedTool

def read_variant_file(variant_file):
	"""
	Import variant file as list of items in header and dictionary of variant data mapped to a bed interval.
	:param variant_file: Path to variant file.
	:return: Header list, dictionary of variant data mapped to bed intervals.
	"""
	variant_dict = {}
	with open(variant_file, 'r') as var_file:
		header = var_file.readline().rstrip().split('\t')
		for line in var_file:
			var_data = line.rstrip().split()
			pos = int(var_data[1])
			bed_info = (var_data[0], pos - 1, pos)
			variant_dict[bed_info] = var_data

	return header, variant_dict

def get_variants_in_trs(variant_dict, trf_file):
	"""
	Intersects variant bed intervals with TRF bed file.
	:param variant_dict: Dictionary of variant data mapped to bed intervals.
	:param trf_file: Path to bed file of TRF annotations.
	:return: List of bed intervals that overlap with TRF annotation.
	"""
	trf_bed = BedTool(trf_file)
	variant_bed = BedTool(list(variant_dict.keys()))
	tr_variants = variant_bed.intersect(trf_bed)
	tr_df = tr_variants.to_dataframe()
	return list(tr_df.itertuples(index=False, name=None))

def validate_tr_variant(variant_info, header):
	"""
	Validate variant in tandem repeat space based on the number of individuals that share the alternate allele.
	:param var_info: List of variant data from input file.
	:param header: List of header information from input file.
	:return: True if the variant is likely de novo, False if it is like an error.
	"""
	validation = variant_info[header.index('final_validation')]

	if variant_info[header.index('singleton')] == 'singleton':
		return True

	if variant_info[header.index('singleton')] == 'potential_recurrent_error' and validation == 'true_de_novo':
		share_high_qual_allele = int(variant_info[header.index('others_with_allele_bq20')])
		share_low_qual_allele = int(variant_info[header.index('others_with_allele_bq10')])
		if share_high_qual_allele < 2 and share_low_qual_allele < 3:
			return True

	return False

def validate_nontr_variant(variant_info, header):
	"""
	Validate variant not in tandem repeat space based on the number of individuals that share the alternate allele.
	:param var_info: List of variant data from input file.
	:param header: List of header information from input file.
	:return: True if the variant is likely de novo, False if it is like an error.
	"""
	share_high_qual_allele = int(variant_info[header.index('others_with_allele_bq20')])
	share_low_qual_allele = int(variant_info[header.index('others_with_allele_bq10')])

	if share_high_qual_allele > 0 or share_low_qual_allele > 1:
		return False

	return True

def main(variant_file, trf_bed, output_file):
	out_header, variants = read_variant_file(variant_file)
	out_header.append("in_tr")

	tr_variants = get_variants_in_trs(variants, trf_bed)

	with open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)
		for interval, var_info in variants.items():
			if interval in tr_variants:
				real_variant = validate_tr_variant(var_info, out_header)
				var_info.append('True')
			else:
				real_variant = validate_nontr_variant(var_info, out_header)
				var_info.append('False')

			if real_variant:
				print('\t'.join(var_info), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="file of de novo variants")
	# ap.add_argument("-t", "--trf", required=True, help="trf bed file")
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()
	#
	# main(args.variants, args.trf, args.output)
	main(snakemake.input.singletons, snakemake.params.trf_bed, snakemake.output.filtered_snvs)
