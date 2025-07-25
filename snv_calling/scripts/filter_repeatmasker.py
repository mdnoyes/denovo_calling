import argparse
import pandas as pd
import numpy as np
from pybedtools import BedTool
from typing import Tuple, List, Dict, Iterable, Optional

def read_variant_file(variant_file: str):
	"""
	Import variant file as list of items in header and dictionary of variant data mapped to a bed interval.
	:param variant_file: Path to variant file.
	:return: Header list, dictionary of variant data mapped to bed intervals.
	"""
	variant_dict = {}
	with open(variant_file, 'r') as var_file:
		header = var_file.readline().rstrip().split('\t')
		for line in var_file:
			var_data = line.rstrip().split('\t')
			pos = int(var_data[1])
			bed_info = (var_data[0], pos - 1, pos)
			variant_dict[bed_info] = var_data[:-1]

	return header[:-1], variant_dict

def get_variants_in_repeats(variant_dict: Dict[tuple[str, int, int], list[str]], repeat_file: str):
	"""
	Intersects variant bed intervals with TRF bed file.
	:param variant_dict: Dictionary of variant data mapped to bed intervals.
	:param repeat_file: Path to bed file of repeatmasker annotations.
	:return: List of bed intervals that overlap with repeatmasker annotation.
	"""
	repeat_bed = BedTool(repeat_file)
	variant_bed = BedTool(list(variant_dict.keys()))
	repeat_variants = variant_bed.intersect(repeat_bed)
	repeat_df = repeat_variants.to_dataframe()
	return list(repeat_df.itertuples(index=False, name=None))

def validate_repeat_variant(variant_info: list[str], header: list[str]):
	"""
	Validate variant in tandem repeat space based on the number of individuals that share the alternate allele.
	:param var_info: List of variant data from input file.
	:param header: List of header information from input file.
	:return: True if the variant is likely de novo, False if it is like an error.
	"""

	singleton = variant_info[header.index('singleton')] == "singleton"
	in_tr = variant_info[header.index('in_tr')] == "True"

	platforms = ['hifi', 'ont', 'illumina']
	ab = [float(variant_info[header.index(f'{platform}_ab')]) for platform in platforms]
	ab = np.nan_to_num(ab)
	mean_ab = np.mean(ab)

	if singleton:
		if in_tr:
			if mean_ab > 0.1:
				return True
			else:
				return False
		if mean_ab > 0.08:
			return True
		if 0 in ab and mean_ab > 0.04:
			return True


	return False


def main(variant_file, repeat_bed, output_file):
	out_header, variants = read_variant_file(variant_file)
	out_header.append("in_repeatmasker")

	repeat_variants = get_variants_in_repeats(variants, repeat_bed)

	with open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)
		for interval, var_info in variants.items():
			if interval in repeat_variants:
				real_variant = validate_repeat_variant(var_info, out_header)
				var_info.append('True')
			else:
				real_variant = True
				var_info.append('False')

			if real_variant:
				print('\t'.join(var_info), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="file of de novo variants")
	# ap.add_argument("-r", "--repeats", required=True, help="repeat bed file")
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()
	#
	# main(args.variants, args.repeats, args.output)
	main(snakemake.input.variants, snakemake.params.repeat_bed, snakemake.output.filtered_snvs)
