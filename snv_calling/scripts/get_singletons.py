import argparse
import numpy as np
from typing import Tuple, List, Dict, Iterable, Optional
from collections import Counter
from collections import defaultdict

def make_sample_dict(parent_dict: Dict[str, str], child: str):
	"""
	Map each sample's ID to its relationship with the de novo child.
	:param parent_dict: Dictionary of every sample's parents.
	:param child: Child ID.
	:return: Dictionary of sample relationships mapped to their ID.
	"""
	fm, sm = child.split('_')
	dad = parent_dict[child][0]
	mom = parent_dict[child][1]

	sample_dict = {}
	sample_dict[child] = 'child'
	sample_dict[dad] = 'father'
	sample_dict[mom] = 'mother'

	for sample, parents in parent_dict.items():
		if sample in sample_dict:
			continue
		if parents == (dad, mom):
			sample_dict[sample] = 'sibling'
		else:
			sample_dict[sample] = 'unrelated'

	return sample_dict

def import_manifest(manifest: str, child: str):
	"""
	Read in sample manifest to dictionary of samples and their parents.
	:param manifest: Path to manifest file.
	:param child: De novo child's ID.
	:return: Dictionayr of parents mapped to their child's ID.
	"""
	parent_dict = {}
	with open(manifest, 'r') as  manifest_file:
		manifest_file.readline()
		for line in manifest_file:
			sample_data = line.split('\t')
			sample_id = sample_data[0] + "_" + sample_data[1]
			parent_dict[sample_id] = (sample_data[2], sample_data[3])

	return make_sample_dict(parent_dict, child)

def import_reads(readfile: str):
	"""
	Map read data from input file for each sample and variant.
	:param readfile: Path to file of read data.
	:return: Dictionary of read data mapped to a tuple of variant ID and sample name.
	"""
	read_dict = defaultdict(lambda: defaultdict(list))

	with open(readfile, 'r') as reads:
		header = reads.readline().rstrip().split('\t')
		basequal_idx = header.index('base_quality')
		mapq_idx = header.index('mapq')

		for line in reads:
			read = line.rstrip().split('\t')
			var_id, sample = read[0:2]

			read_data = (read[3], read[basequal_idx], read[mapq_idx])
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

def filter_sample_alleles(read_list: List[str]):
	"""
	Filter reads based on mapq, and then assign alleles based on base quality.
	:param read_list: List of read data, including read name, allele, base quality, and mapq.
	:return: List of alleles with base quality under 20, and alleles with base quality 20 or above.
	"""
	mapq_threshold = 59
	bq_over_20_alleles = []
	bq_under_20_alleles  = []

	for read in read_list:
		if not mapq_filter(read[2], mapq_threshold):
			continue

		if base_qual_filter(read[1]) == 1:
			bq_over_20_alleles.append(read[0])
		elif base_qual_filter(read[1]) == 0:
			bq_under_20_alleles.append(read[0])

	return bq_over_20_alleles, bq_under_20_alleles

def check_allele_count_threshold(allele_list: List[str], quality: str, alt_allele: str):
	"""
	Determine if alternate allele may have been inherited from list of parental alleles.
	:param allele_list: Either high or low quality parental allele list.
	:param quality: 'high' if base quality 20 or greater, 'low' if not.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: True if alternate allele count is below threshold, False if it is equal or above threshold.
	"""
	threshold = inheritance_thresholds[quality]
	alt_allele_count = Counter(allele_list)[alt_allele]

	if alt_allele_count > threshold:
		return False

	return True

def check_for_alternate_allele(high_qual_alleles: List[str], low_qual_alleles: List[str], alt_allele: str):
	"""
	Count the high and low quality alternate alleles, and compare them to inheritance threshold.
	:param high_qual_alleles: List of alleles with base quality 20 or greater.
	:param low_qual_alleles: List of alleles with base quality under 20.
	:param alt_allele: Str alt allele from de novo mutation.
	:return: Tuple of True/False if high qual alts below threshold, True/False if low qual alts below threshold.
	"""

	not_in_high_qual = check_allele_count_threshold(high_qual_alleles, 'high', alt_allele)
	not_in_low_qual = check_allele_count_threshold(high_qual_alleles + low_qual_alleles, 'low', alt_allele)

	return (not_in_high_qual, not_in_low_qual)

def get_sample_data(variant_reads: List[List[str]], alt_allele: str):
	"""
	Count the number of alternate alleles on filtered reads.
	:param variant_read_dict: List of reads containing variant from sample.
	:param alt_allele: Alternate allele from de novo mutation.
	:returns:  Tuple of True/False if alt in high qual reads, True/False if alt in low qual reads.
	"""
	high_qual_alleles, low_qual_alleles = filter_sample_alleles(variant_reads)
	return check_for_alternate_allele(high_qual_alleles, low_qual_alleles, alt_allele)

def get_samples_with_allele(absence_list: List[Tuple[bool, bool]]):
	"""
	Count the number of samples with the alternate allele.
	:param absence_list: List of tuples, True/False if alt in high qual reads, True/False if alt in low qual reads.
	:return: Number of samples with high quality allele, number of samples with low quality allele.
	"""
	high_quality_count = 0
	low_quality_count = 0

	for sample_absence in absence_list:
		if sample_absence[0] == False:
			high_quality_count += 1
		if sample_absence[1] == False:
			low_quality_count += 1

	return high_quality_count, low_quality_count


def get_status(allele_data: Dict[str, List[Tuple[bool, bool]]]):
	"""
	Determine how many samples share the alternate allele, and determine potential error status.
	:param allele_data: Dictionary of relationships mapped to allele absence.
	:return: Allele sharing status, number of samples with high qual allele, number with low qual allele.
	"""
	unrelated_high_qual, unrelated_low_qual = get_samples_with_allele(allele_data['unrelated'])

	if 'sibling' in allele_data.keys():
		sibling_high_qual, sibling_low_qual = get_samples_with_allele(allele_data['sibling'])
		if sibling_high_qual > 0 or sibling_low_qual > 0:
			return 'potential_parental_pzm', str(unrelated_high_qual + sibling_high_qual), str(unrelated_low_qual + sibling_low_qual)

	if unrelated_high_qual > 0 or unrelated_low_qual > 1:
		return 'potential_recurrent_error', str(unrelated_high_qual), str(unrelated_low_qual)

	return 'singleton', str(unrelated_high_qual), str(unrelated_low_qual)

def main(read_file, variant_file, manifest_file, child, output_file):
	reads = import_reads(read_file)
	samples = import_manifest(manifest_file, child)


	all_variants = []

	with open(variant_file, 'r') as variants, open(output_file, 'w') as outfile:
		var_header = variants.readline().rstrip().split('\t')
		out_header = var_header + ['singleton', 'others_with_allele_bq20', 'others_with_allele_bq10']
		print('\t'.join(out_header), file = outfile)

		for line in variants:
			var_data = line.rstrip().split('\t')
			var_id = var_data[2]
			alt = var_data[4]

			sample_allele_data = defaultdict(list)

			for sample_id, relationship in samples.items() :
				if relationship in ['father', 'mother', 'child']:
					continue

				try:
					sample_reads = reads[var_id][sample_id]
				except KeyError:
					sample_reads = []

				alt_allele_absence = get_sample_data(sample_reads, alt)
				sample_allele_data[relationship].append(alt_allele_absence)

			singleton_status, others_with_high_qual, others_with_low_qual = get_status(sample_allele_data)

			variant_info = var_data + [singleton_status, others_with_high_qual, others_with_low_qual]
			print('\t'.join(variant_info), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-m", "--manifest", required=True, help="manifest file")
	# ap.add_argument("-r", "--reads", required=True, help="File of read data")
	# ap.add_argument("-c", "--child", required=True, help="sample ID", type=str)
	# ap.add_argument("-o", "--output", required=True, help="Output file")
	# args = ap.parse_args()

	inheritance_thresholds = {'high': 1, 'low': 2}

	# main(args.reads, args.variants, args.manifest, args.child, args.output)
	main(snakemake.input.depth_stats, snakemake.input.dnms, snakemake.params.bam_manifest, snakemake.params.child, snakemake.output.singletons)
