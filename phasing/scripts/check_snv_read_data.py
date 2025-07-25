import pysam
import re
import pandas as pd
import sys
import argparse
import numpy as np
from tqdm import tqdm
from collections import Counter

def import_sources(source_file: str):
	"""
	Import file of hifi data sources to dictionary of read names mapped to 'blood' or 'cells'.
	:param source_file: Path to tab-separated file with read prefix and source info. Only child data is necessary.
	:return: Dictionary of read prefix string and sample name mapped to 'blood' or 'cells'.
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

def import_cell_line_only(sample_file: str):
	"""
	Import file of cell line only samples, save as list.
	:param source_file: File with each sample name on a new line.
	:return: Tuple of sample names.
	"""
	sample_list = []

	with open(sample_file, 'r') as samples:
		for sample in samples:
			sample_list.append(sample.rstrip())

	return tuple(sample_list)

def get_bam(sample, bams, sequence_platform):
	"""
	Retrieves the path to the sample's bam file.
	:param sample: Sample to find bam file for.
	:param bams: Manifest dataframe of bam files.
	:param sequence platform: illumina, ont, or hifi
	:returns: PySam AlignmentFile object of the sample's read data.
	"""
	fm, sm = sample.split("_")
	bam_path = bams.loc[(bams['sample'] == sm) & (bams['family'] == fm)][sequence_platform].values[0]
	bam = pysam.AlignmentFile(bam_path, "rc")
	return bam

def check_blood(read, source_dict, sample_name):
	read_prefix = read.query_name.split('/')[0]
	if source_dict[(sample_name, read_prefix)] == 'blood':
		return True

	return False

def get_idx(read, position):
	try:
		return read.get_reference_positions(full_length=True).index(position - 1)
	except ValueError:
		return None

def check_read_quality(read, var_idx):
	mapq = read.mapping_quality
	if mapq < 59:
		return False 
	
	# have to be annoying about this because some reads are length 0, which makes read.quer_qualities None
	try:
		dnm_quality = read.query_qualities[var_idx]
		if dnm_quality < 20:
			return False
	except TypeError:
		return False

	return True

def get_allele(read, var_idx, ref, alt):
	dnm_allele = read.query_sequence[var_idx]
	if dnm_allele == ref:
		return 'ref'
	elif dnm_allele == alt:
		return 'alt'
	else:
		return None

def import_tagging_snps(snp_file):
	tagging_snp_dict = {}
	with open(snp_file, 'r') as snps:
		header = snps.readline().split('\t')
		for snp in snps:
			snp = snp.rstrip().split('\t')
			dnm_id = snp[0]
			if dnm_id not in tagging_snp_dict.keys():
				tagging_snp_dict[dnm_id] = []
			snp_info = snp[1:]
			tagging_snp_dict[dnm_id].append(snp_info)

	return tagging_snp_dict


def get_tagging_snps(variant_id, variant_pos, snp_list):
	tagging_snp_dict = {}

	for snp in snp_list:
		snp_dict = {}
		snp_pos = int(snp[1])

		ref, alt, dad_gt, mom_gt = snp[2:6]
		if '1' in dad_gt and '1' not in mom_gt:
			snp_dict[ref] = 'mom'
			snp_dict[alt] = 'dad'
		elif '1' in mom_gt and '1' not in dad_gt:
			snp_dict[ref] = 'dad'
			snp_dict[alt] = 'mom'
		elif '0' in dad_gt and '0' not in mom_gt:
			snp_dict[ref] = 'dad'
			snp_dict[alt] = 'mom'
		elif '0' in mom_gt and '0' not in dad_gt:
			snp_dict[ref] = 'mom'
			snp_dict[alt] = 'dad'
		else:
			raise ValueError('there is a genotype you dont know')
		tagging_snp_dict[(snp[0], int(snp[1]))] = snp_dict

	return tagging_snp_dict

def check_snp(snp_position, read):
	snp_idx = get_idx(read, snp_position)

	if not snp_idx:
		return None

	base = read.query_sequence[snp_idx]
	quality = read.query_qualities[snp_idx]

	if quality < 20:
		return None

	return base

def check_tagging_snps(read, dnm_pos, tagging_snp_dict, window_size):
	haplotypes = {
		'dad': [],
		'mom': []
	}

	for snp, allele_dict in tagging_snp_dict.items():
		snp_chrom, snp_pos = snp
		snp_base = check_snp(snp_pos, read)

		try:
			haplotype = allele_dict[snp_base]
		except KeyError:
			continue

		distance_to_dnm = abs(snp_pos - dnm_pos)
		haplotypes[haplotype].append(distance_to_dnm)

	read_inheritance = assign_inheritance(haplotypes, window_size)
	return read_inheritance

def check_reads(dnm_info, bam, tagging_snp_dict, read_sources, window_size, platform, read_classes, sample_name, cell_line_samples):
	chrom = dnm_info['chr']
	pos = int(dnm_info['pos'])

	read_types = []

	for aligned_segment in bam.fetch(chrom, pos - 1, pos):
		if sample_name not in cell_line_samples:
			if platform == 'hifi' and not check_blood(aligned_segment, read_sources, sample_name):
				continue

		dnm_idx = get_idx(aligned_segment, pos)
		if not dnm_idx:
			continue
		if not check_read_quality(aligned_segment, dnm_idx):
			continue
		allele = get_allele(aligned_segment, dnm_idx, dnm_info['ref'], dnm_info['alt'])
		if not allele:
			continue

		inheritance = check_tagging_snps(aligned_segment, pos, tagging_snp_dict, window_size)
		if not inheritance:
			continue

		read_types.append(f'{inheritance}_{allele}')

	observed_reads = Counter(read_types)
	count_summary = summarize_reads(observed_reads, read_classes)

	return count_summary


def assign_inheritance(haplotype_dict, window_size):
	dist = window_size / 2

	# paternal alleles get scaled by -1, maternal alleles by +1
	inheritance_values = [-1 * x for x in haplotype_dict['dad']] + haplotype_dict['mom']

	if len(inheritance_values) == 0:
		return None

	inheritance_score = sum([dist / x for x in inheritance_values]) / sum([dist / abs(x) for x in inheritance_values])

	if inheritance_score > 0:
		return 'maternal'
	elif inheritance_score < 0:
		return 'paternal'

def summarize_reads(read_counts, read_classes):
	count_dict = {}
	for read_class in read_classes:
		try:
			count_dict[read_class] = read_counts[read_class]
		except KeyError:
			count_dict[read_class] = 0

	return count_dict

def phase(read_counts):
	if read_counts['paternal_alt'] == 0 and read_counts['maternal_alt'] > 0:
		return 'maternal'
	elif read_counts['maternal_alt'] == 0 and read_counts['paternal_alt'] > 0:
		return 'paternal'
	else:
		return 'cannot_determine'

def determine_pzm(inheritance, read_counts):
	if inheritance == 'paternal':
		if read_counts['paternal_ref'] > 1:
			var_type = 'pzm'
		else:
			var_type = 'dnm'
	elif inheritance == 'maternal':
		if read_counts['maternal_ref'] > 1:
			var_type = 'pzm'
		else:
			var_type = 'dnm'
	else:
		var_type = 'cannot_determine'

	return var_type

def main(denovo_file, tagging_snps, bam_manifest, child_id, platform, window_size, hifi_source_file, cell_line_file, output_file):
	hifi_sources = import_sources(hifi_source_file)
	all_tagging_snps = import_tagging_snps(tagging_snps)
	cell_line_samples = import_cell_line_only(cell_line_file)

	variant_df = pd.read_csv(denovo_file, sep='\t')

	bam_df = pd.read_csv(bam_manifest, sep='\t', dtype=str)

	child_bam = get_bam(child_id, bam_df, platform)

	out_columns =  variant_df.columns.tolist()[:11]

	read_classes = ['paternal_ref', 'paternal_alt', 'maternal_ref', 'maternal_alt']
	out_header = out_columns + read_classes + ['new_phase', 'new_type']

	outfile = open(output_file, 'w')
	print('\t'.join(out_header), file = outfile)

	for idx, row in tqdm(variant_df.iterrows(), total=len(variant_df)):
		out_info = [row[x] for x in out_columns]
		id = row['id']
		pos = int(row['pos'])

		try:
			tagging_snps = get_tagging_snps(id, pos, all_tagging_snps[id])
		except KeyError:
			out_info = out_info + ['NA', 'NA', 'NA', 'NA', 'cannot_determine', 'cannot_determine']
			print('\t'.join([str(x) for x in out_info]), file = outfile)
			continue

		read_summary = check_reads(row, child_bam, tagging_snps, hifi_sources, window_size, platform, read_classes, child_id, cell_line_samples)

		out_info = out_info + [read_summary[x] for x in read_classes]

		dnm_inheritance = phase(read_summary)
		out_info.append(dnm_inheritance)

		variant_type = determine_pzm(dnm_inheritance, read_summary)
		out_info.append(variant_type)

		print('\t'.join([str(x) for x in out_info]), file = outfile)

	outfile.close()

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	# ap.add_argument("-s", "--tagging_snps", required=True, help="file of tagging snps")
	# ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	# ap.add_argument("-c", "--child", required=True, type=str, help="sample id")
	# ap.add_argument("-hs", "--hifi_sources", required=True, type=str, help="hifi sources")
	# ap.add_argument("-w", "--window_size", required=True, type=int, help="window size")
	# ap.add_argument("-p", "--platform", required=True, type=str, help="sequencing platform")
	# ap.add_argument("-cl", "--cell_line", required=True, type=str, help="cell line only samples")
	# ap.add_argument("-o", "--output", required=True, type=str, help="output file")
	#
	# args = ap.parse_args()
	#
	# main(args.input, args.tagging_snps, args.bams, args.child, args.platform, args.window_size, args.hifi_sources, args.output)
	main(snakemake.input.denovo_file, snakemake.input.tagging_snps, snakemake.params.bam_manifest, snakemake.params.child, snakemake.params.platform, snakemake.params.window_size, snakemake.params.hifi_sources, snakemake.params.cell_line, snakemake.output.platform_haplotypes)
