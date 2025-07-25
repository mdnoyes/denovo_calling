from pybedtools import BedTool
import argparse
from tqdm import tqdm

def get_callable_bed(callable_bed, chromosome_file):
	chromosome_bed = BedTool(chromosome_file)
	return callable_bed.intersect(chromosome_bed)

def get_bed_size(bed):
	space = 0
	for interval in bed:
		start, end = [int(x) for x in interval[1:3]]
		size = 1 + end - start
		space += size
	return str(space)

def get_callable_space(callable_bed, region_file, sex_chromosome):
	region_bed = BedTool(region_file)
	callable_regions = callable_bed.intersect(region_bed)
	return get_bed_size(callable_regions), get_bed_size(region_bed.intersect(BedTool(sex_chromosome)))

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="callable genome bed")
	# ap.add_argument("-s", "--sample", required=True, type=str, help="sample name")
	# ap.add_argument("-o", "--output", required=True, type=str, help="output")
	# args = ap.parse_args()

	sample = snakemake.params.sample
	sex = snakemake.params.sex
	# sample = args.sample

	X_bed = 'bed/chrX_no_PARS.bed'
	Y_bed = 'bed/chrY_no_PARS.bed'

	file_dict = {
		'exons' : 'region_files/exons_1idx.bed',
		'acrocentric' : 'region_files/acros.bed',
		'all_centromeres' : 'region_files/centromeres.bed',
		'cen_hor' : 'region_files/cen_hor_1idx.bed',
		'cen_hsat' : 'region_files/cen_hsat_1idx.bed',
		'cen_monomeric' : 'region_files/cen_mon_1idx.bed',
		'centromeric_transition' : 'region_files/cen_transition_1idx.bed',
		'all_segdup' : 'region_files/seg_dup_1idx.bed',
		'segdup_90_95' : 'region_files/segdups_90_95_identity_1idx.bed',
		'segdup_95_98' : 'region_files/segdups_95_98_identity_1idx.bed',
		'segdup_98_99' : 'region_files/segdups_98_99_identity_1idx.bed',
		'segdup_over_99' : 'region_files/segdups_over_99_identity_1idx.bed'
	}

	out_header = ['sample', 'chromosome', 'region', 'callable_size', 'total_size']

	print('importing callable genome....')
	full_callable_bed = BedTool(snakemake.input.callable_bed)

	callable_chrX_bed = get_callable_bed(full_callable_bed, X_bed)
	callable_X, X_size = get_callable_space(full_callable_bed, X_bed, X_bed)

	if sex == 'male':
		callable_chrY_bed = get_callable_bed(full_callable_bed, Y_bed)
		callable_Y, Y_size = get_callable_space(full_callable_bed, Y_bed, Y_bed)


	with open(snakemake.output.callable_sizes, 'w') as outfile:
	# with open(args.output, 'w') as outfile:
		print('\t'.join([sample, 'chrX', 'chrX', callable_X, X_size]), file=outfile)
		if sex == 'male':
			print('\t'.join([sample, 'chrY', 'chrY', callable_Y, Y_size]), file=outfile)
		for region_name, file_path in tqdm(file_dict.items(), total = len(file_dict)):
			callable, total = get_callable_space(callable_chrX_bed, file_path, X_bed)
			print('\t'.join([sample, 'chrX', region_name, callable, total]), file = outfile)
			if sex == 'male':
				callable, total = get_callable_space(callable_chrY_bed, file_path, Y_bed)
				print('\t'.join([sample, 'chrY', region_name, callable, total]), file=outfile)
