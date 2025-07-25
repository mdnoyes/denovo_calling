import argparse
import numpy as np
from typing import Tuple, List, Dict, Iterable, Optional
from collections import Counter
from collections import defaultdict

def import_sources(source_file: str):
	"""
	Import file of hifi data sources to dictionary of read names mapped to 'blood' or 'cells'.
	:param source_file: Path to tab-separated file with read prefix and source info. Only child data is necessary.
	:return: Dictionary of read prefix string mapped to 'blood' or 'cells'.
	"""
	source_dict = {}
	with open(source_file, 'r') as sources:
		header = sources.readline().rstrip().split('\t')
		sample = header.index('sample')
		source = header.index('source')
		prefix = header.index('readname_prefix')

		for line in sources:
			bam_info = line.rstrip().split('\t')
			source_dict[(bam_info[sample], bam_info[prefix])] = bam_info[source]

	return source_dict


def get_var_info(var: List[str]):
	"""
	Subset relevant variant info for the output file.
	:param var: List of variant data from the input variant file.
	:return: List of relevant info, including variant identifiers and caller.
	"""
	var_info = var[0:5]
	var_info.append(var[7].split('=')[-1])
	return var_info


def import_reads(read_file: str):
	"""
	Import sample's read_data file to a dictionary of sample read data mapped to variants.
	:param read_file: Path to read_data file.
	:return: Dict of dicts, for each variant key, there should be a dict of read data mapped to each family member.
	"""
	read_dict = defaultdict(lambda: defaultdict(list))

	with open(read_file, 'r') as reads:
		header = reads.readline().rstrip().split('\t')
		basequal_idx = header.index('base_quality')
		mapq_idx = header.index('mapq')

		for line in reads:
			read = line.rstrip().split('\t')
			var_id, sample = read[0:2]

			read_data = (read[2], read[3], read[basequal_idx], read[mapq_idx])
			read_dict[var_id][sample].append(read_data)

	return read_dict

def mapq_filter(mapq: str, threshold: int):
	"""
	Filter mapq based on threshold.
	:param mapq: Mapq value for a given read.
	:param threshold: Mapq threshold as defined for sequencing platform.
	:return: True if above threshold, False if below threshold.
	"""
	try:
		 mapq = int(mapq)
	except ValueError:
		 return False

	if mapq < threshold:
		return False

	return True


def base_qual_filter(base_qual: str):
	"""
	Filter base_qual based on threshold.
	:param mapq: Base quality value for a given read.
	:return: True if above threshold, False if below threshold.
	"""
	try:
		 base_qual = int(base_qual)
	except ValueError:
		return 0

	if base_qual >= 20:
		return 1

	if base_qual >= 10:
		return 0

	return -1


def filter_parent_alleles(read_list: List[str]):
	"""
	Filter parental reads based on mapq, and then assign alleles based on base quality.
	:param read_list: List of read data, including read name, allele, base quality, and mapq.
	:return: List of alleles with base quality under 20, and alleles with base quality 20 or above.
	"""
	mapq_threshold = mapq_thresholds[platform]
	bq_over_20_alleles = []
	bq_under_20_alleles  = []

	for read in read_list:
		if not mapq_filter(read[3], mapq_threshold):
			continue

		if base_qual_filter(read[2]) == 1:
			bq_over_20_alleles.append(read[1])
		elif base_qual_filter(read[2]) == 0:
			bq_under_20_alleles.append(read[1])

	return bq_over_20_alleles, bq_under_20_alleles


def determine_read_source(read_name: str, hifi_sources: Dict[str, str], sample_name: str):
	"""
	Determine if the sequencing read was derived from blood or cell line data.
	:param read_name: Name of sequencing read.
	:param hifi_sources: Dictionary of read prefix string mapped to 'blood' or 'cells'.
	:return: Either 'blood' or 'cells'.
	"""
	if platform == 'ont':
		return 'cells'

	if platform == 'illumina':
		return 'blood'

	read_prefix = read_name.split('/')[0]

	return hifi_sources[(sample_name, read_prefix)]


def filter_child_alleles(read_list: List[str], hifi_sources: Dict[str, str], sample_name: str):
	"""
	Filter parental reads based on mapq, and then assign alleles based on base quality.
	:param read_list: List of read data, including read name, allele, base quality, and mapq.
	:param hifi_sources: Dictionary of read prefix string mapped to 'blood' or 'cells'.
	:return: List of alleles from blood, and alleles from cell line.
	"""
	mapq_threshold = mapq_thresholds[platform]
	allele_dict = {'blood': [], 'cells': []}

	for read in read_list:
		if not mapq_filter(read[3], mapq_threshold):
			continue

		if not base_qual_filter(read[2]):
			continue

		source = determine_read_source(read[0], hifi_sources, sample_name)
		allele_dict[source].append(read[1])

	return allele_dict['blood'], allele_dict['cells']


def check_allele_count_threshold(allele_list: list, quality: str, alt_allele: str):
	"""
	Determine if alternate allele may have been inherited from list of parental alleles.
	:param allele_list: Either high or low quality parental allele list.
	:param quality: 'high' if base quality 20 or greater, 'low' if not.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: True if alternate allele count is below threshold, False if it is equal or above threshold.
	"""
	threshold = inheritance_thresholds[platform][quality]
	alt_allele_count = Counter(allele_list)[alt_allele]

	if alt_allele_count > threshold:
		return False

	return True


def determine_possible_inheritance(high_qual_alleles: list, low_qual_alleles: list, alt_allele: str):
	"""
	Count the high and low quality alternate alleles, and compare them to inheritance threshold.
	:param high_qual_alleles: List of alleles with base quality 20 or greater.
	:param low_qual_alleles: List of alleles with base quality under 20.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: Tuple of True/False if total allele count is greater than 0, True/False if high qual alts below threshold, True/False if low qual alts below threshold.
	"""

	depth = len(high_qual_alleles) + len(low_qual_alleles)
	if depth == 0:
		return (False, False, False)

	high_qual_pass = check_allele_count_threshold(high_qual_alleles, 'high', alt_allele)
	low_qual_pass = check_allele_count_threshold(high_qual_alleles + low_qual_alleles, 'low', alt_allele)

	return (True, high_qual_pass, low_qual_pass)


def calculate_depth_and_ab(allele_list: list, alt_allele: str):
	"""
	Calculate read depth and allele balance for blood or cell data.
	:param allele_list: List of alleles.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: Tuple of alt allele count, total read depth, and alt allele balance.
	"""
	depth = len(allele_list)
	if depth == 0:
		return (0, 0, np.nan)

	alt_count = Counter(allele_list)[alt_allele]
	ab = round(alt_count / depth, 5)

	return (alt_count, depth, ab)


def determine_false_positive(allele_stats: Tuple):
	"""
	Determine if variant is true, inherited, false, or missing data.
	:param allele_stats: Tuple with child alt allele count, read depth, and allele balance.
	:return: Tuple of True/False if depth is greater than 0, True/False if alt allele is represented.
	"""
	if allele_stats[1] == 0:
		return (False, False)

	alt_threshold = 0
	if platform == 'hifi':
		alt_threshold = 1

	if allele_stats[0] > alt_threshold:
		return (True, True)

	return (True, False)


def validate_variant(blood_summary: Tuple, cell_summary: Tuple, dad_summary: Tuple, mom_summary: Tuple):
	"""
	Determine if variant is true, inherited, false, or missing data.
	:param blood_summary: Tuple of True/False if child has blood data, True/False if child has alt allele in blood.
	:param cell_summary: Tuple of True/False if child has cell data, True/False if child has alt allele in cells.
	:param dad_summary: Tuple of True/False if dad has data, True/False if alt in high qual reads, True/False if alt in low qual reads.
	:param mom_summary: Tuple of True/False if mom has data, True/False if alt in high qual reads, True/False if alt in low qual reads.
	:return: Final validation of variant.
	"""
	if dad_summary[0] == False or mom_summary[0] == False:
		return 'missing_parent_data'

	if False in dad_summary[1:] or False in mom_summary[1:]:
		return 'inherited'

	if blood_summary[0] == False and cell_summary[0] == False:
		return 'missing_child_data'

	if platform == 'hifi':
		if blood_summary[1] == True:
			return 'true_de_novo'
		if cell_summary[1] == True:
			return 'cell_line_only'

	if platform == 'ont':
		if cell_summary[1] == True:
			return 'true_de_novo'

	if platform == 'illumina':
		if blood_summary[1] == True:
			return 'true_de_novo'

	return 'false_positive'


def get_parent_data(variant_reads, alt_allele):
	"""
	Count the number of alternate alleles on filtered reads and determine possible inheritance.
	:param variant_read_dict: List of reads containing variant from parent.
	:param alt_allele: Alternate allele from de novo mutation.
	:returns:  Tuple of True/False if parent has data, True/False if alt in high qual reads, True/False if alt in low qual reads.
	"""
	high_qual_alleles, low_qual_alleles = filter_parent_alleles(variant_reads)
	return determine_possible_inheritance(high_qual_alleles, low_qual_alleles, alt_allele)


def main(variant_file, read_file, hifi_sources, father, mother, child, output_file):
	if platform == 'hifi':
		out_header = ['chr', 'pos', 'id', 'ref', 'alt', 'callers', 'child_cell_filtered_count', 'child_cell_depth', 'child_cell_allele_balance', 'child_blood_filtered_count', 'child_blood_depth', 'child_blood_allele_balance', 'filtered_status']
	else:
		out_header = ['chr', 'pos', 'id', 'ref', 'alt', 'callers', 'child_filtered_count', 'child_depth', 'child_allele_balance', 'filtered_status']

	hifi_source_dict = import_sources(hifi_sources)
	variant_reads = import_reads(read_file)

	with open(variant_file, 'r') as variants, open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)
		for line in variants:
			var_data = line.rstrip().split('\t')
			id, ref, alt = var_data[2:5]

			variant_info = get_var_info(var_data)

			dad_inh_summary = get_parent_data(variant_reads[id][father], alt)
			mom_inh_summary = get_parent_data(variant_reads[id][mother], alt)

			child_blood_alleles, child_cell_alleles = filter_child_alleles(variant_reads[id][child], hifi_source_dict, child)
			blood_stats = calculate_depth_and_ab(child_blood_alleles, alt)
			blood_dnm = determine_false_positive(blood_stats)
			cell_stats = calculate_depth_and_ab(child_cell_alleles, alt)
			cell_dnm = determine_false_positive(cell_stats)

			validation = validate_variant(blood_dnm, cell_dnm, dad_inh_summary, mom_inh_summary)

			if platform == 'hifi':
				out_data = variant_info + [str(x) for x in cell_stats] + [str(x) for x in blood_stats] + [validation]
			elif platform == 'ont':
				out_data = variant_info + [str(x) for x in cell_stats] + [validation]
			elif platform == 'illumina':
				out_data = variant_info + [str(x) for x in blood_stats] + [validation]

			# if id == '13230_p1_chr1_207259371_C_G':
			# 	print([x for x in variant_reads[id][father] if x[1] == 'G'])
			# 	print(dad_inh_summary)
			# 	high_qual_alleles, low_qual_alleles = filter_parent_alleles(variant_reads[id][father])
			# 	print(high_qual_alleles)
			# 	print(low_qual_alleles)

			print('\t'.join(out_data), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-r", "--reads", required=True, help="File of read data")
	# ap.add_argument("-f", "--father", required=True, help="father ID", type=str)
	# ap.add_argument("-m", "--mother", required=True, help="mother ID", type=str)
	# ap.add_argument("-c", "--child", required=True, help="sample ID", type=str)
	# ap.add_argument("-s", "--source", required=False, help="read data source (blood/cells)", type=str)
	# ap.add_argument("-p", "--platform", required=True, help="Platform")
	# ap.add_argument("-o", "--output", required=True, help="Output file")
	# args = ap.parse_args()
	#
	# platform = args.platform
	platform = snakemake.params.platform

	inheritance_thresholds = {
		'hifi': {'high': 1, 'low': 2},
		'ont': {'high': 2, 'low': 3},
		'illumina': {'high': 1, 'low': 2}
	}

	mapq_thresholds = {
		'hifi': 59,
		'ont': 59,
		'illumina': 0
	}

	# main(args.variants, args.reads, args.source, args.father, args.mother, args.child, args.output)
	main(snakemake.input.snv_calls,
		snakemake.input.read_data,
		snakemake.params.hifi_sources,
		snakemake.params.father,
		snakemake.params.mother,
		snakemake.params.child,
		snakemake.output.platform_validation)
