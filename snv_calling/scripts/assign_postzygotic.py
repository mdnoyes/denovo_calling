import argparse
from scipy.stats import chi2_contingency
from scipy.stats import binomtest
from typing import Tuple, List, Dict

def get_allele_balance(var_data: List[str], index_dict: Dict[str, List[int]]):
	"""
	Pull out list of allele balance values from variant data.
	:param var_data: List of variant data from input file.
	:param index_dict: Dictionary of indexes for sequencing platform data.
	:return: List of allele balance values for each platform.
	"""
	ab_list = []

	for platform, indeces in index_dict.items():
		ab_list.append(var_data[indeces[2]])

	return ab_list


def remove_missing_data(num_platforms: int, contingency_table: List[List[int]]):
	"""
	Remove platforms with no read data (0 ref alleles, 0 alt alleles) from contingency_table.
	:param num_platforms: Int number of sequencing platforms in dataset.
	:param contingency_table: Contingency table, two lists of ref and alternate allele counts.
	:return: Contingency table with missing sequencing data removed.
	"""
	missing_data = []
	for x in range(0, num_platforms):
		if contingency_table[0][x] == 0 and contingency_table[1][x] == 0:
			missing_data.append(x)

	for x in missing_data:
		contingency_table[0].pop(x)
		contingency_table[1].pop(x)

	return contingency_table


def make_contingency_table(var_data: List[str], index_dict: Dict[str, List[int]]):
	"""
	Make contingency_table for allele count chi squared test.
	:param var_data: List of variant data from input file.
	:param index_dict: Dictionary of indexes for sequencing platform data.
	:return: List of alternate allele count list, reference allele count list.
	"""
	contingency_table = [[], []]

	for platform, indeces in index_dict.items():
		alt_allele_count = int(var_data[indeces[0]])
		ref_allele_count = int(var_data[indeces[1]]) - alt_allele_count
		contingency_table[0].append(alt_allele_count)
		contingency_table[1].append(ref_allele_count)

	contingency_table = remove_missing_data(len(index_dict), contingency_table)

	return contingency_table


def calculate_concordance(contingency_table: List[List[int]], num_platforms: int):
	"""
	Perfrom chi squared test to determine if allele balance is concordant across platforms.
	:param contingency_table: Contingency table, two lists of ref and alternate allele counts.
	:param num_platforms: Number of platforms represented in the table.
	:return: True if allele balance is concordant across platforms, False if not.
	"""
	if [0, 0, 0] in contingency_table or [0, 0] in contingency_table:
		return True

	if len(contingency_table[0]) != num_platforms:
		return False

	pval = chi2_contingency(contingency_table)[1]

	if pval >= 0.05:
		concordant_ab = True
	else:
		concordant_ab = False

	return concordant_ab


def assign_pzm(contingency_table: List[List[int]], concordant_ab: bool):
	"""
	Determine if variant is postzygotic using a binomial test to determin if allele balance is different from 0.5 across platforms.
	:param contingency_table: Contingency table, two lists of ref and alternate allele counts.
	:param concordant_ab: True if allele balance is concordant across platforms, False if not.
	:return: 'dnm' is AB is not different from 0.5, 'pzm' if it is significanlty less than 0.5.
	"""
	if concordant_ab == True:
		total_alt_alleles = sum(contingency_table[0])
		total_alleles = sum(contingency_table[1]) + total_alt_alleles
		net_ab = total_alt_alleles / total_alleles

		pval = binomtest(total_alt_alleles, total_alleles, 0.5).pvalue

		if pval < 0.05 and net_ab < 0.5:
			return 'pzm'

	return 'dnm'


def make_index_dict(header: List[str]):
	"""
	Generate dictionary of allele count, depth, and AB indeces for each sequencing platform.
	:param header: Header from input file.
	:return: Dictionary of data indeces mapped to sequencing platforms.
	"""
	platforms = [x.split('_')[0] for x in header if 'allele_count' in x]

	index_dict = {}
	for platform in platforms:
		index_dict[platform] = (header.index(f'{platform}_alt_allele_count'), \
			header.index(f'{platform}_dp'), \
			header.index(f'{platform}_ab'))

	return index_dict


def main(variant_file, output_file):
	with open(variant_file, 'r') as variants, open(output_file, 'w') as outfile:
		input_header = variants.readline().rstrip().split('\t')
		val_idx = input_header.index('final_validation')

		pf_index_dict = make_index_dict(input_header)

		out_header = input_header[:6] + \
			[f'{x}_ab' for x in pf_index_dict.keys()] + \
			['concordant_ab', 'variant_type'] +  \
			input_header[val_idx + 1: ] + \
			['validation']

		print('\t'.join(out_header), file=outfile)

		for line in variants:
			variant_data = line.rstrip().split('\t')
			data_table = make_contingency_table(variant_data, pf_index_dict)
			pass_ab_concordance = calculate_concordance(data_table, len(pf_index_dict))
			variant_type = assign_pzm(data_table, pass_ab_concordance)

			allele_balances = get_allele_balance(variant_data, pf_index_dict)

			var_output = variant_data[:6] + \
				allele_balances + \
				[str(pass_ab_concordance), variant_type] + \
				variant_data[val_idx + 1: ] + \
				[variant_data[val_idx]]


			print('\t'.join(var_output), file=outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-o", "--output", required=True, help="Output file")
	# args = ap.parse_args()
	#
	# main(args.variants, args.output)
	main(snakemake.input.filtered_snvs, snakemake.output.pzm_annotated)
