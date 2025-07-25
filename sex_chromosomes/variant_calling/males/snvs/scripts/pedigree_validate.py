import pandas as pd
import sys
import argparse
import numpy as np
from tqdm import tqdm
from typing import Tuple, List, Dict, Iterable
from dataclasses import dataclass

class FamilyTree:
	its_just_people: Dict[str, 'Person'] = {}

	def keys(self) -> Iterable[str]:
		return self.its_just_people.keys()

	def __contains__(self, sample_id: str):
		return sample_id in self.its_just_people

	def __getitem__(self, sample_id: str) -> 'Person':
		return self.its_just_people[sample_id]

	@staticmethod
	def grow_tree(relationships: Iterable[Tuple[str, str, str]]) -> 'FamilyTree':
		tree = FamilyTree()

		for relationship in relationships:
			sample_id, father_id, mother_id = relationship

			try:
				person = tree.its_just_people[sample_id]
			except KeyError:
				person = tree.its_just_people[sample_id] = Person(sample_id)

			if father_id != 'NA':
				try:
					father = tree.its_just_people[father_id]
				except KeyError:
					father = tree.its_just_people[father_id] = Person(father_id)
				person.father = father
				father.children.append(person)

			if mother_id != 'NA':
				try:
					mother = tree.its_just_people[mother_id]
				except KeyError:
					mother = tree.its_just_people[mother_id] = Person(mother_id)
				person.mother = mother
				mother.children.append(person)

		return tree


class Person:
	sample_id: str
	father: 'Person'
	mother: 'Person'
	children: List['Person']

	def __init__(self, sample_id: str, father: 'Person' = None, mother: 'Person' = None):
		if sample_id == 'NA':
			raise ValueError(f'invalid sample ID {sample_id}')

		self.sample_id = sample_id
		self.father = father
		self.mother = mother
		self.children = []

	def grandparents(self, lineage: str) -> List['Person']:
		parent = {
			'paternal': self.father,
			'maternal': self.mother
		}[lineage]

		return [parent.father, parent.mother] if parent else None

	def avunculi(self, lineage: str) -> List['Person']:
		parent = {
			'paternal': self.father,
			'maternal': self.mother
		}[lineage]

		return parent.siblings if parent else None

	@property
	def siblings(self) -> List['Person']:
		siblings = set()
		if self.father:
			siblings.update(self.father.children)
		if self.mother:
			siblings.intersection_update(self.mother.children)
		siblings.discard(self)
		return list(siblings) if siblings else None

	def __repr__(self) -> str:
		return self.sample_id


def read_manifest(manifest_file):
	all_families = {}
	with open(manifest_file, 'r') as manifest:
		header = manifest.readline().rstrip().split('\t')
		indeces = [
			header.index(col_name) for col_name in ('sample', 'father', 'mother')
		]
		for line in manifest:
			data = line.rstrip().split('\t')
			family = data[0]
			if (family not in all_families):
				all_families[family] = []
			sm, fa, mo = [f'{family}_{data[idx]}' for idx in indeces]
			if fa == f'{family}_':
				fa = 'NA'
			if mo == f'{family}_':
				mo = 'NA'
			all_families[family].append((sm, fa, mo))

	return all_families


def get_grouped_allele_counts(allele_string):
	if allele_string == 'NA':
		return None

	alleles = allele_string.split(';')[1:]
	if 'NA' in alleles:
		alleles = [x != 'NA' for x in alleles]

	if len(alleles) == 0:
		return None
	else:
		return [int(x) for x in alleles]

def check_group(allele_list, relationship):
	if relationship == 'grandparent':
		threshold = 1
	else:
		threshold = 2

	alt_allele_present = [x > threshold for x in allele_list]
	if sum(alt_allele_present) == 0:
		return True
	else:
		return False

def check_extended_family(parent_name, grandparent_counts, avuncular_counts, sibling_counts):
	explanation = None

	if sibling_counts:
		sibling_pass = check_group(sibling_counts, 'sibling')
	else:
		sibling_pass = True

	if grandparent_counts:
		grandparent_pass = check_group(grandparent_counts, 'grandparent')
	else:
		grandparent_pass = True

	if avuncular_counts:
		avuncular_pass = check_group(avuncular_counts, 'avunculi')
	else:
		avuncular_pass = True

	if grandparent_pass and sibling_pass:
		if avuncular_pass:
			# parent does not have the de novo allele and neither does anyone else in their family
			explanation = 'true_de_novo'
	elif grandparent_pass and not sibling_pass:
		if avuncular_pass:
			# parent, parent's siblings, and grandparents do not have the de novo allele but at least one of their siblings do
			if avuncular_counts:
				explanation = f"{parent_name}_parental_pzm"
			else:
				explanation = f"{parent_name}_parental_pzm_no_avunculi"
		else:
			explanation = f"{parent_name}_dropped_haplotype_potential_grandparent_pzm"
	elif not grandparent_pass:
		# parent does not have the de novo allele but their parents do, parental siblings do not matter in this case
		explanation = f'{parent_name}_dropped_haplotype'

	if explanation is None:
		raise ValueError("you have encountered a case not specified in your logic you loser")

	if explanation == 'true_de_novo':
		return True, explanation
	else:
		return False, explanation


def resolve_explanation(dad, mom):
	if '_parental_pzm_no_avunculi' in dad:
		return mom
	elif '_parental_pzm_no_avunculi' in mom:
		return dad
	if 'dropped_haplotype' in dad and 'dropped_haplotype' in mom:
		return 'segregating_site'

	return dad + "_" + mom

def validate_variant(var_data):
	# paternal_gps = get_grouped_allele_counts(var_data['paternal_gp'])
	maternal_gps = get_grouped_allele_counts(var_data['maternal_gps'])
	siblings = get_grouped_allele_counts(var_data['siblings'])
	# paternal_avunculi = get_grouped_allele_counts(var_data['paternal_avunculi'])
	maternal_avunculi = None

	mom_pass, mom_situation = check_extended_family('mom', maternal_gps, maternal_avunculi, siblings)
	# dad_pass, dad_situation = check_extended_family('dad', paternal_gps, paternal_avunculi, siblings)
	dad_pass = True
	dad_situation = 'irrelevant'

	if mom_pass and dad_pass:
		return 'true_de_novo', 'true_de_novo'
	elif mom_pass and not dad_pass:
		return 'inherited', dad_situation
	elif dad_pass and not mom_pass:
		return 'inherited', mom_situation
	else:
		return 'inherited', resolve_explanation(dad_situation, mom_situation)

def count_relatives(sample, site_id, read_data, relative_list):
	all_children = [x.sample_id for x in sample.children]
	rels = [x for x in relative_list if x != sample.sample_id]

	children_with_allele = []
	relatives_with_allele = []

	for sample in rels:
		try:
			count = int(read_data.loc[site_id][sample])
		except ValueError:
			count = 0
		if count > 0:
			if sample in all_children:
				children_with_allele.append(sample)
			else:
				relatives_with_allele.append(sample)

	if len(children_with_allele) > 0:
		kids = ';'.join(children_with_allele)
	else:
		if len(all_children) > 0:
			kids = 'none'
		else:
			kids = 'NA'

	if len(relatives_with_allele) > 0:
		fam = ';'.join(relatives_with_allele)
	else:
		fam = 'none'

	return kids, fam

def main(variant_file, sample_summary, family_reads, manifest_file, output_file):
	variants = pd.read_csv(variant_file, sep='\t')
	variants = variants.astype(str)
	all_read_data = pd.read_csv(sample_summary, sep='\t', keep_default_na=False, index_col = 0)
	family_data = pd.read_csv(family_reads, sep='\t', keep_default_na=False)

	families_dict = read_manifest(manifest_file)
	print(families_dict)

	pedigree = FamilyTree.grow_tree(families_dict['K1463'])
	samples = sorted(pedigree.keys())
	print(samples)

	output = open(output_file, 'w')
	header = variants.columns.to_list() +  ['pedigree_validation', 'origin', 'children_with_allele', 'relatives_with_allele']
	print('\t'.join(header), file = output)

	for idx, row in family_data.iterrows():
		validation, explanation = validate_variant(row)
		denovo_sample = '_'.join(row['site'].split('_')[:2])

		person = pedigree[denovo_sample]
		kids_with_variant, relatives_with_variant = count_relatives(person, row['site'], all_read_data, samples)

		var_info = variants.loc[variants['id'] == row['site']].values.flatten().tolist() + \
			[validation, explanation, kids_with_variant, relatives_with_variant]
		print('\t'.join(var_info), file = output)

	output.close()


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-m", "--manifest", required=True, help="Sample manifest with parental information")
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-a", "--all_samples", required=True, help="File of read data")
	# ap.add_argument("-f", "--family_data", required=True, help="File of family data")
	# ap.add_argument("-o", "--output", required=True, help="Output file")
	# args = ap.parse_args()
	#
	# main(args.variants, args.all_samples, args.family_data, args.manifest, args.output)
	main(snakemake.input.dnms, snakemake.input.summary, snakemake.input.family_reads, snakemake.params.bam_manifest, snakemake.output.singletons)
