import argparse
import csv
from dataclasses import dataclass
from typing import Tuple

Supporting_reads = int
Total_reads = int

def make_sample_dict(parent_dict, child):
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
		elif sample in parent_dict[mom] or sample in parent_dict[dad]:
			sample_dict[sample] = 'grandparent'
		else:
			sample_dict[sample] = 'unrelated'

	return sample_dict

def import_manifest(manifest, child):
	parent_dict = {}
	males = []
	with open(manifest, 'r') as  manifest_file:
		manifest_file.readline()
		for line in manifest_file:
			sample_data = line.split('\t')
			sample_id = sample_data[0] + "_" + sample_data[1]
			if sample_data[2] == 'M':
				males.append(sample_id)
			dad = sample_data[3]
			mom = sample_data[4]
			parent_dict[sample_id] = (dad, mom)

	sample_dict = make_sample_dict(parent_dict, child)

	return {k: sample_dict[k] for k in males}

class SampleAlleles:
	variant_id: str
	sample_id: str
	role: str
	ont: Tuple[Supporting_reads, Total_reads]
	hifi: Tuple[Supporting_reads, Supporting_reads, Total_reads]
	ilm: Tuple[Supporting_reads, Total_reads]

	def __init__(self, variant_id: str, sample_id: str):
		self.variant_id = variant_id
		self.sample_id = sample_id

		self.ont = None
		self.hifi = None
		self.ilm = None

	def set_role(self, sample_role: str):
		self.role = sample_role

	def add_platform(self, platform: str, supporting_reads: str, cell_supporting_reads: str, total_reads: str):
		if platform == 'hifi':
			self.hifi = (int(supporting_reads), int(cell_supporting_reads), int(total_reads))
		elif platform == 'ont':
			self.ont = (int(supporting_reads), int(total_reads))
		elif platform == 'illumina':
			self.ilm = (int(supporting_reads), int(total_reads))

class VariantAlleles:
	variant_id: str
	dad: SampleAlleles
	mom: SampleAlleles
	child: SampleAlleles

	def __init__(self, variant_id:str):
		self.variant_id = variant_id
		self.dad = None
		# self.mom = None
		self.child = None

	def add_sample(self, sample_alleles: SampleAlleles):
		if sample_alleles.role == 'dad':
			self.dad = sample_alleles
		# elif sample_alleles.role == 'mom':
			# self.mom = sample_alleles
		elif sample_alleles.role == 'child':
			self.child = sample_alleles

def get_platform(file_name):
	file = file_name.split('/')[-1]
	platform = file.split('_')[2]
	return platform

def add_platform(platform, var_file):
	for line in var_file:
		yield platform, line

def check_sorted(variant_id, sample_id, variant_dict):
	if variant_dict['id'] != variant_id or variant_dict['sample'] != sample_id:
		raise ValueError("your files are not in sorted order!")

def get_sample_alleles(variant_id, sample_id, platform_data, sample_role):
	alleles = SampleAlleles(variant_id, sample_id)
	alleles.set_role(sample_role)
	for platform, var_dict in platform_data:
		check_sorted(variant_id, sample_id, var_dict)
		alleles.add_platform(platform, var_dict['supporting_reads'], var_dict['cell_line_supporting_reads'], var_dict['total_reads'])

	return alleles

def check_finish(file_list):
	for file in file_list:
		try:
			next(file)
		except StopIteration:
			pass
		else:
			raise ValueError("one of your files is longer than the other!")

def check_parent(read_tuple, platform):
	allele_depth = read_tuple[0]

	if platform in ('hifi', 'illumina'):
		threshold = 1
	else:
		threshold = 2

	if allele_depth < threshold:
		return True
	else:
		return False

def check_hifi_de_novo_allele(dad_alleles, child_alleles):
	blood_alleles = child_alleles[0]
	cell_alleles = child_alleles[1]

	if blood_alleles + cell_alleles < 2:
		return 'false_positive'

	dad = check_parent(dad_alleles, 'hifi')

	if dad:
		if blood_alleles > 0:
			return 'true_de_novo'
		else:
			return 'cell_line_only'
	else:
		return 'inherited'


# def check_de_novo_allele(dad_alleles, mom_alleles, child_alleles):
def check_de_novo_allele(dad_alleles, child_alleles, platform):
	dad = check_parent(dad_alleles, platform)

	if not dad:
		return 'inherited'

	if child_alleles[0] < 2:
		return 'false_positive'

	return 'true_de_novo'


def get_platform_validations(variant_alleles):
	hifi = check_hifi_de_novo_allele(variant_alleles.dad.hifi,
		# variant_alleles.mom.hifi,
		variant_alleles.child.hifi)
	ilm = check_de_novo_allele(variant_alleles.dad.ilm,
		# variant_alleles.mom.ilm,
		variant_alleles.child.ilm, 'illumina')
	ont = check_de_novo_allele(variant_alleles.dad.ont,
		# variant_alleles.mom.ont,
		variant_alleles.child.ont, 'ont')

	return hifi, ilm, ont


def validate_dnm(variant_alleles):
	hifi_val, illumina_val, ont_val = get_platform_validations(variant_alleles)

	if 'inherited' in [hifi_val, illumina_val, ont_val]:
		return 'inherited'

	if 'false_positive' in [hifi_val, illumina_val]:
		return 'false_positive'

	if illumina_val == 'true_de_novo' and hifi_val == 'true_de_novo':
		return 'true_de_novo'

	if hifi_val == 'cell_line_only':
		if illumina_val == 'true_de_novo':
			return 'true_de_novo'
		else:
			return 'cell_line_only'

	raise Exception(f"you missed this case: hifi {hifi_val} illumina {illumina_val}")

def add_to_sib_dict(variant_id, variant_data, sib_dict):
	if variant_id not in sib_dict.keys():
		sib_dict[variant_id] = {}

	for platform, var_dict in variant_data:
		allele_count = int(var_dict['supporting_reads'])

		if allele_count < 1:
			sib_dict[variant_id][platform] = 'true_de_novo'
		else:
			sib_dict[variant_id][platform] = 'inherited'

	return sib_dict

def check_sibling(sibling_validations):
	if 'inherited' in sibling_validations.values():
		return 'inherited'

	return 'true_de_novo'

def validate_indels(read_data_files, family_dict, sample_dict, out_columns):
	sibling_allele_counts = {}

	validation_dict = {}
	var_data = None
	for variants in zip(*read_data_files):
		variant = variants[0][1]['id']
		if var_data is None:
			var_data = VariantAlleles(variant)

		sample = variants[0][1]['sample']

		try:
			sample_role = family_dict[sample]
		except KeyError:
			if sample_dict[sample] == 'sibling':
				sibling_allele_counts = add_to_sib_dict(variant, variants, sibling_allele_counts)
			continue

		sample_alleles = get_sample_alleles(variant, sample, variants, sample_role)

		# Checking that the new sample has a matching variant
		var_data.add_sample(sample_alleles)

		# if var_data.dad and var_data.mom and var_data.child:
		if var_data.dad and var_data.child:
			validation = validate_dnm(var_data)
			# print(variant, str(variants[0][1]['indel_size']), validation)
			out_data = [variants[0][1][x] for x in out_columns[:7]] + [validation]
			validation_dict[variant] = out_data
			var_data = None
		else:
			continue

	return validation_dict, sibling_allele_counts


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="validation files", nargs='+')
	# ap.add_argument("-f", "--father", required=True, help="dad's id")
	# # ap.add_argument("-m", "--mother", required=True, help="mom's id")
	# ap.add_argument("-s", "--sample", required=True, help="sample's id")
	# ap.add_argument("-o", "--output", required=True, help="output")
	# ap.add_argument("-mf", "--manifest", required=True, help="output")
	#
	# args = ap.parse_args()

	# family_ids = {
	# 	args.father: 'dad',
	# 	# args.mother: 'mom',
	# 	args.sample: 'child'
	# }
	# variant_files = {get_platform(path): path for path in args.input}
	#
	# samples = import_manifest(args.manifest, args.sample)

	family_ids = {
		snakemake.params.father: 'dad',
		# snakemake.params.mother: 'mom',
		snakemake.params.sample: 'child'
	}
	samples = import_manifest(snakemake.input.manifest, snakemake.params.sample)
	variant_files = {get_platform(path): path for path in snakemake.input.read_summary}


	sibling = False
	if 'sibling' in samples.values():
		sibling = True


	files_with_platform = [add_platform(pf, csv.DictReader(open(path, 'r'), delimiter = '\t')) for pf, path in variant_files.items()]

	out_columns = ['chr', 'pos', 'id', 'ref', 'alt', 'caller', 'indel_size', 'validation']

	variant_dict, sibling_validation_dict = validate_indels(files_with_platform, family_ids, samples, out_columns)

	check_finish(files_with_platform)

	# with open(args.output, 'w') as outfile:
	with open(snakemake.output.validated_indels, 'w') as outfile:
		print('\t'.join(out_columns), file = outfile)
		for variant, out_data in variant_dict.items():
			if sibling and out_data[-1] == 'true_de_novo':
				if variant in sibling_validation_dict:
					out_data[-1] = check_sibling(sibling_validation_dict[variant])
			print('\t'.join(out_data), file = outfile)
