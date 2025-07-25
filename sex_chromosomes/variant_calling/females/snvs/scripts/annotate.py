from pybedtools import BedTool
from liftover import ChainFile
import pysam
from typing import Tuple, List, Dict, Iterable, Optional

def get_complements(complement_file: str):
	"""
	Read in triplet completment lookup table and save as dictionary.
	:param complement_file: Path to complement lookup table file.
	:return" Dictionary of triplets mapped to their complement.
	"""
	complements = {}
	with open(complement_file) as f:
		for line in f:
			key, val = line.rstrip().split('\t')
			complements[key] = val
			complements[val] = val

	return complements

def read_variant_file(variant_file: str):
	"""
	Read in de novo variants, save them to a dictionary with bed intervals.
	:param variant_file: Path to file of de novo mutations.
	:return: Header of variant file as list of strings.
	:return: Dictionary of variant info mapped to bed interval.
	"""
	variants = {}
	with open(variant_file, 'r') as var_file:
		header = var_file.readline().rstrip().split('\t')
		for line in var_file:
			var_data = line.rstrip().split()
			pos = int(var_data[1])
			bed_info = (var_data[0], pos - 1, pos)
			variants[bed_info] = {k:v for k,v in zip(header, var_data)}

	return header, variants

def get_variants_in_bed(test_bed: BedTool, bed_file: str):
	"""
	Find the overlap between the input variants and the bed regions.
	:param test_bed: BedTool object of intervals for your de novo variants.
	:param bed_file: Path to bed file of regions of interest.
	:return: List of variant intervals that overlap with the bed regions.
	"""
	ref_bed = BedTool(bed_file)
	intersected_variants = test_bed.intersect(ref_bed)
	variant_df = intersected_variants.to_dataframe()

	return list(variant_df.itertuples(index=False, name=None))

def add_regional_annotation(variant_interval: Tuple[str, int, int], variants_in_region: List[Tuple[str, int, int]]):
	"""
	Check whether a variant falls in an annotated region.
	:param variant_interval: Tuple of bed interval for the variant.
	:variants_in_region: List of bed intervals that fall in annotated region.
	:return: 'True' if variant in region, 'False' if not.
	"""
	if variant_interval in variants_in_region:
		return 'True'
	else:
		return 'False'

def lift(variant_info: List[str], converter: ChainFile):
	"""
	Get the coordinates of a variant in a different coordinate space.
	:param variant_info: List of variant data from input file.
	:param converter: ChainFile object from old to new coordinate space.
	:return: String of the new coordinate, or 'NA' if missing in the new reference.
	"""
	lifted_pos = converter[variant_info['chr']][int(variant_info['pos'])]
	if lifted_pos:
		return str(lifted_pos[0][1])
	else:
		return 'NA'

def make_annotation_dict(params, variants: BedTool):
	"""
	Create dictionary of regional annotations.
	:def params: Snakemake params list.
	:def variants: BedTool object of all the variants.
	:return: Dictionary of region names mapped to list of variant intervals that fall in them.
	"""
	annotation_dict = {}

	annotation_dict['seg_dup'] = get_variants_in_bed(variants, params.seg_dup_bed)
	annotation_dict['centromere'] = get_variants_in_bed(variants, params.centromere_bed)
	annotation_dict['lcr'] = get_variants_in_bed(variants, params.low_complexity_bed)
	annotation_dict['recent_repeat'] = get_variants_in_bed(variants, params.recent_repeat_bed)
	annotation_dict['refseq_exon'] = get_variants_in_bed(variants, params.refseq_bed)

	return annotation_dict

def transition_transversion(ref: str, alt: str):
	"""
	Assign variant as transition or transversion mutation.
	:param ref: Reference allele.
	:param alt: Alternate allele.
	:return: Transition or transversion.
	"""
	bases = set([ref, alt])
	if bases == set(['A', 'G']) or bases == set(['C', 'T']):
		return 'transition'
	else:
		return 'transversion'

def check_triplet(variant_info: List[str], ref_fasta: FastaFile, complement_lookup: Dict[str, str]):
	"""
	Create triplet name for variant, with 3' and 5' base converted into correct complement, and describes mutation type.
	:param variant_info: List of variant data from input file.
	:param ref_fasta: FastaFile object of reference genome sequence.
	:complement_lookup: Dictionary of complement sequence mapped to correct complement.
	:return: Dictionary of triplet sequence, whether the mutation is CpG>CpG, and substitution_type.
	"""
	ref = variant_info['ref']
	alt = variant_info['alt']
	pos1 = int(variant_info['pos']) - 2
	pos2 = int(variant_info['pos']) + 1
	seq = ref_fasta.fetch(variant_info['chr'], start=pos1, end=pos2).upper()
	triplet = "%s-%s"%(seq, alt)
	triplet = complement_lookup[triplet]

	if triplet[-4:-2] == 'CG':
		CpG = 'True'
	else:
		CpG = 'False'

	trans = transition_transversion(ref, alt)

	return {'tiplet': triplet, 'CpG': CpG, 'substitution_type': trans}

def main(variant_file, chain_file, fasta_file, complement_lookup, snakemake_params, output_file):
	out_header, variant_dict = read_variant_file(variant_file)
	variant_bed = BedTool(list(variant_dict.keys()))
	t2t_to_grch38 = ChainFile(chain_file, 'chm13v2', 'grch38')
	t2t_fasta = pysam.FastaFile(fasta_file)
	complement_dict = get_complements(complement_lookup)

	annotations  = make_annotation_dict(snakemake_params, variant_bed)

	out_header = out_header + list(annotations.keys()) + ['hg38_pos', 'triplet', 'CpG', 'substitution_type']

	with open(output_file, 'w') as outfile:
		print('\t'.join(out_header), file = outfile)

		for interval, var_info in variant_dict.items():
			for annotation, affected_variants in annotations.items():
				var_info[annotation] = add_regional_annotation(interval, affected_variants)
			var_info['hg38_pos'] = lift(var_info, t2t_to_grch38)
			var_info.update(check_triplet(var_info, t2t_fasta, complement_dict))

			print('\t'.join(var_info.values()), file = outfile)


if __name__ == "__main__":
	main(snakemake.input.variants, snakemake.params.chain_file, snakemake.params.t2t_fasta, snakemake.params.complement_lookup, snakemake.params, snakemake.output.annotated)
