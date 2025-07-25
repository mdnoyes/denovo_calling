import pandas as pd
import argparse
from operator import itemgetter
import gzip

def import_dnms(dnm_file):
	dnm_dict = {}
	with open(dnm_file, 'r') as infile:
		header = infile.readline().rstrip().split('\t')
		for line in infile:
			variant = line.rstrip().split('\t')
			dnm_dict[(variant[1], int(variant[2]))] = variant[5]

	return dnm_dict

def check_informative_site(dad_gt, mom_gt, child_gt):
	ref_genotype = ['0' , '0']
	alt_genotypes = [
		['1', '1'],
		['0', '1'],
		['1', '0']
	]

	if len(set(child_gt)) == 2:
		if dad_gt == ref_genotype and mom_gt in alt_genotypes:
			return True

		if mom_gt == ref_genotype and dad_gt in alt_genotypes:
			return True

	return False

def get_variant_info(variant):
	return variant[0], int(variant[1]), variant[3], variant[4]

def get_tagged_snps(chrom, pos, dnm_dict, window_size):
	tagged_dnms = [x for x in dnm_dict.keys() if x[0] == chrom and abs(x[1]-pos) < window_size]

	return tagged_dnms

def read_vcf(vcf_file, dnm_dict, window_size):
	tagging_list = []

	with open(vcf_file, 'r') as infile:
		for line in infile:
			variant = line.rstrip().split('\t')
			if line.startswith('#'):
				continue

			chrom, pos, ref, alt = get_variant_info(variant)
			if len(ref) != 1 or len(alt) != 1:
				 continue

			dad_gt, mom_gt, child_gt = [x[:3].replace('|', '/').split('/') for x in variant[9:]]
			if not check_informative_site(dad_gt, mom_gt, child_gt):
				continue

			tagged_dnms = get_tagged_snps(chrom, pos, dnm_dict, window_size)

			for dnm in tagged_dnms:
				id = dnm_dict[dnm]
				distance_to_dnm = str(pos - dnm[1])
				tagging_list.append([
					id, chrom, str(pos), ref, alt,
					'/'.join(dad_gt), '/'.join(mom_gt), '/'.join(child_gt),
					distance_to_dnm
					])

	return tagging_list

def main(dnm_file, vcf_file, window_size, output_file):
		window = window_size / 2

		dnms = import_dnms(dnm_file)
		tagging_snps = read_vcf(vcf_file, dnms, window)

		out_header = ['denovo_id', 'chrom', 'pos', 'ref', 'alt', 'dad_gt', 'mom_gt', 'child_gt', 'distance_to_dnm']

		with open(output_file, 'w') as outfile:
			print('\t'.join(out_header), file = outfile)
			for snp in tagging_snps:
				print('\t'.join(snp), file = outfile)


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-d", "--denovo", required=True, help="File of de novo variants")
	# ap.add_argument("-v", "--vcf", required=True, help="VCF of variants near DNMs")
	# ap.add_argument("-o", "--output", required=True, help="Output file name")
	# ap.add_argument("-w", "--window_size", required=False, help="Size of window to examine around variants of interest", default="40000")
	# args = ap.parse_args()

	# main(args.denovo, args.vcf, args.window_size, args.output)
	main(snakemake.input.dnms, snakemake.input.vcf, snakemake.params.window_size, snakemake.output.tagging_snps)
