import argparse
import numpy as np
from typing import Tuple, List, Dict
from collections import Counter
from collections import defaultdict

class Variant:
	chrom: str
	pos: int
	id: str
	old_validation: str
	tandem_repeat: bool
	concordant_ab: bool
	relatives_with_allele: int
	allele_balance: float
	allele_counts: Dict[str, Dict[str, List[Tuple[int]]]]

	def __init__(
		self,
		variant_data: List[str],
		file_header: List[str],
		read_dict: Dict[str, Dict[str, Tuple]]
	):
		self.chrom, self.pos, self.id = variant_data[:3]
		self.pos = int(self.pos)
		self.old_validation = variant_data[file_header.index('final_validation')]
		self.tandem_repeat = eval(variant_data[file_header.index('in_tr')])
		self.is_singleton = variant_data[file_header.index('singleton')] == 'singleton'
		self.allele_balance = np.mean([float(variant_data[file_header.index(f'{pf}_ab')]) for pf in ['hifi', 'ont', 'illumina']])
		self.allele_counts = dict([(pf, read_dict[pf][self.id]) for pf in ['hifi', 'ont', 'illumina']])


	def __str__(self):
		return self.id

	def find_missing_support(self, platform: str) -> str:
		if platform == 'hifi':
			threshold = 0
		else:
			threshold = 1

		for platform in ['hifi', 'ont', 'illumina']:
			if self.allele_counts[platform]['mom'][1] > threshold:
				return 'inherited'

		if self.noise_in_parent('mom'):
			return 'too_noisy_to_call'

		if self.allele_counts[platform]['child'][1] > 0:
			return 'true_de_novo'

		if self.allele_balance >= 0.9:
				return 'true_de_novo'

		return self.old_validation

	def noise_in_parent(self, parent: str) -> bool:
		if self.allele_counts['illumina'][parent][1] > 3 or self.allele_counts['ont'][parent][1] > 3:
			return True

		if self.allele_counts['illumina'][parent][2] > 3 or self.allele_counts['ont'][parent][2] > 3:
			return True

		return False

	def relitigate(self) -> str:
		# tandem repeats require support across all three platforms, no noise in parents, and AB of at least 10%
		if self.tandem_repeat:
			if self.allele_counts['illumina']['child'][1] < 1 or self.allele_counts['ont']['child'][1] < 1:
				return 'tandem_repeat_error'

			if self.noise_in_parent('mom'):
				return 'tandem_repeat_error'

			if self.allele_balance < 0.9:
				return 'tandem_repeat_error'

		# variants with discordant AB must have AB of at least 5%
		if not self.allele_balance >= 0.9:
			return 'sequencing_error'

		if 'no_hifi' in self.old_validation:
			return self.find_missing_support('hifi')

		if 'no_ont' in self.old_validation:
			return self.find_missing_support('ont')

		if 'ilm_false_positive' in self.old_validation:
			return self.find_missing_support('illumina')

		return self.old_validation


def import_all_reads(read_files: List[str], sample_dict: Dict[str, str]):
	"""
	Import read data from all sequencing platforms into dictionary.
	:param read_files: List of paths so sequencing read data files.
	:param sample_dict: Sample IDs mapped to their relationship.
	:return: Dict of dict of dicts - read data dicts mapped to sequencing platform name.
	"""
	reads_dict = {}
	for read_file in read_files:
		platform = read_file.split('_')[-2].split('.')[0]
		platform_reads = import_reads(read_file, sample_dict)
		reads_dict[platform] = platform_reads

	return reads_dict

def import_reads(read_file: str, sample_dict: Dict[str, str]):
	"""
	Import sample's read data file to a dictionary of sample read data mapped to variants.
	:param read_file: Path to read_data file.
	:param sample_dict: Sample IDs mapped to their relationship.
	:return: Dict of dicts, for each variant key, there should be a dict of read data mapped to each family member.
	"""
	read_dict = defaultdict(lambda: defaultdict(list))

	with open(read_file, 'r') as reads:
		header = reads.readline().rstrip().split('\t')
		ref_idx = header.index('ref_count')
		alt_idx = header.index('alt_count')
		other_idx = header.index('other_count')

		for line in reads:
			read = line.rstrip().split('\t')
			var_id, sample = read[0:2]
			try:
				sample_name = sample_dict[sample]
			except KeyError:
				continue

			try:
				read_data = (int(read[ref_idx]), int(read[alt_idx]), int(read[other_idx]))
			except ValueError:
				read_data = (np.nan, np.nan, np.nan)

			read_dict[var_id][sample_name] = read_data

	return read_dict

def main(variant_file, read_file, father, mother, child, output_file):
	samples = {mother: 'mom', child: 'child'}
	all_reads = import_all_reads(read_file, samples)

	with open(variant_file, 'r') as variants, open(output_file, 'w') as outfile:
		header = variants.readline().rstrip().split('\t')
		print('\t'.join(header), file = outfile)
		for line in variants:
			variant_line = line.rstrip().split('\t')
			variant = Variant(variant_line, header, all_reads)
			new_validation = variant.relitigate()

			# if new_validation == 'true_de_novo':
			variant_line[-1] = new_validation
			print('\t'.join(variant_line), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="file of de novo variants")
	# ap.add_argument("-r", "--read_data", required=True, help="file of read_data", nargs = '+')
	# ap.add_argument("-o", "--output", required=True, help="output")
	# ap.add_argument("-f", "--father", required=True, help="father ID", type=str)
	# ap.add_argument("-m", "--mother", required=True, help="mother ID", type=str)
	# ap.add_argument("-c", "--child", required=True, help="sample ID", type=str)
	#
	# args = ap.parse_args()

	# main(args.variants, args.read_data, args.father, args.mother, args.child, args.output)
	main(snakemake.input.pzm_annotated, snakemake.input.read_data, snakemake.params.father, snakemake.params.mother, snakemake.params.child, snakemake.output.true_dnms)
