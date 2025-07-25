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

def get_allele_counts(read_file):
	read_data = {}

	with open(read_file, 'r') as reads:
		header = reads.readline().rstrip().split('\t')
		id_idx, sm_idx, allele_idx, qual_idx, mapq_idx = [
			header.index(col_name) for col_name in ('id', 'sample', 'allele', 'base_quality', 'mapq')
		]

		current_site = current_sample = None
		for line in reads.readlines():
			read = line.rstrip().split('\t')
			if current_site != read[id_idx] or current_sample != read[sm_idx]:
				if current_site is not None:
					read_data[(current_site, current_sample)] = alt_allele_count
				current_site = read[id_idx]
				current_sample = read[sm_idx]
				alt_allele_count = 0

			if filter_read(read, qual_idx, mapq_idx) and check_alt_allele(read, id_idx, allele_idx):
				alt_allele_count += 1

	return read_data

def filter_read(read_data, base_qual, mapq):
	try:
		if int(read_data[base_qual]) >= 10 and int(read_data[mapq]) >= 57:
			return True
	except ValueError:
		pass

	return False

def check_alt_allele(read_data, site_id, allele):
	alt_allele = read_data[site_id].split('_')[-1]
	if read_data[allele] == alt_allele:
		return True
	else:
		return False

def main(variant_file, manifest_file, read_file, denovo_sample):
	variants = pd.read_csv(variant_file, sep='\t')
	families_dict = read_manifest(manifest_file)

	all_data = open(f'read_data/chrX_no_PAR/{denovo_sample}_all_sample_summary.tsv', 'w')
	family_specific = open(f'read_data/chrX_no_PAR/{denovo_sample}_family_specific_reads.tsv', 'w')

	for family in families_dict:
		pedigree = FamilyTree.grow_tree(families_dict[family])
		samples = sorted(pedigree.keys())

		alt_allele_counts = get_allele_counts(read_file)

		# Print to first file
		print('\t'.join(('site', *samples)), file = all_data)
		# Print to second file
		print('\t'.join(('site', 'siblings', 'children', 'paternal_gp', 'maternal_gps', 'paternal_avunculi',)), file = family_specific)
		for idx, row in variants.iterrows():
			site_id = row['id']

			# Print to first file
			print('\t'.join((site_id, *map(str, (alt_allele_counts.get((site_id, sample_id), 'NA') for sample_id in samples)))), file = all_data)

			person = pedigree[denovo_sample]
			relationships = [
				person.siblings,
				person.children,
				person.grandparents('paternal'),
				person.grandparents('maternal'),
				person.avunculi('paternal')
			]

			familial_allele_counts = []
			for rel in relationships:
				if rel is None:
					count = 'NA'
				elif isinstance(rel, list):
					counts = tuple(alt_allele_counts.get((site_id, r.sample_id), 'NA') for r in rel if r is not None)
					try:
						all_allele_counts = [count for count in counts if count != 'NA']
						avg = sum(all_allele_counts) / len(all_allele_counts)
						count = ';'.join(map(str, (avg, *counts)))
					except ZeroDivisionError:
						count = 'NA'
				else:
					count = str(alt_allele_counts.get((site_id, rel.sample_id), 'NA'))
				familial_allele_counts.append(count)

			# Print to second file
			print('\t'.join((site_id, *familial_allele_counts)), file = family_specific)

	all_data.close()
	family_specific.close()

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-m", "--manifest", required=True, help="Sample manifest with parental information")
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-r", "--reads", required=True, help="File of read data")
	# ap.add_argument("-s", "--denovo_sample", required=True, help="denovo_sample")
	# args = ap.parse_args()
	#
	# main(args.variants, args.manifest, args.reads, args.denovo_sample)
	main(snakemake.input.dnms, snakemake.params.bam_manifest, snakemake.input.read_data, snakemake.params.child)
