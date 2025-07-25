import argparse
import gzip

def check_finish(file_list):
	for file in file_list:
		try:
			next(file)
		except StopIteration:
			pass
		else:
			raise ValueError("one of your files is longer than the other!")


def process_sample(sample_info):
	site_info = sample_info.rstrip().split('\t')
	gt_data = site_info[9].split(':')
	gt = gt_data[0].replace('|', '/')
	try:
		gq_idx = site_info[8].split(':').index('GQ')
		gq = gt_data[gq_idx]
	except ValueError:
		gq = 0

	return (gq, gt)

def filter_female_X(dad, mom, child):
	if dad[1] not in ['0/0', './.']:
		return False

	if mom[1] != '0/0':
		return False

	if int(child[0]) <= 20:
		return False

	return True

def filter_male_X(mom, child):
	if int(child[0]) <= 20:
		return False

	if mom[1] in ['0/0', './.', '1/1']:
		return True

	return False

def filter_male_Y(dad, child):
	if int(child[0]) <= 20:
		return False

	if dad[1] in ['0/0', '1/1']:
		return True

	return False

def male_X(mom_file, son_file, output_file):
	files = [gzip.open(path, 'rt') for path in (mom_file, son_file)]

	out_header = ['chr', 'pos', 'ref', 'mom_gt', 'child_gt']

	with gzip.open(output_file, 'wt') as outfile:
		print('\t'.join(out_header), file = outfile)

		for mom, child in zip(*files):
			if mom.startswith('#'):
				if not child.startswith('#'):
					raise ValueError('headers are of unequal length! figure yourself out')
				continue

			mom_info, child_info = map(process_sample, (mom, child))

			if filter_male_X(mom_info, child_info):
				site_info = mom.rstrip().split('\t')
				site_pos = site_info[:2]
				ref = site_info[3][0]
				out_info = site_pos + [ref, mom_info[1], child_info[1]]

				print('\t'.join(out_info), file = outfile)

	check_finish(files)

def male_Y(dad_file, son_file, output_file):
	files = [gzip.open(path, 'rt') for path in (dad_file, son_file)]

	out_header = ['chr', 'pos', 'ref', 'dad_gt', 'child_gt']

	with gzip.open(output_file, 'wt') as outfile:
	# with gzip.open(args.output, 'wt') as outfile:
		print('\t'.join(out_header), file = outfile)

		for dad, child in zip(*files):
			if dad.startswith('#'):
				if not child.startswith('#'):
					raise ValueError('headers are of unequal length! figure yourself out')
				continue

			dad_info, child_info = map(process_sample, (dad, child))

			if filter_male_Y(dad_info, child_info):
				site_info = dad.rstrip().split('\t')
				site_pos = site_info[:2]
				ref = site_info[3][0]
				out_info = site_pos + [ref, dad_info[1], child_info[1]]

				print('\t'.join(out_info), file = outfile)

	check_finish(files)

def female_X(dad_file, mom_file, daughter_file, output_file):
	files = [gzip.open(path, 'rt') for path in (dad_file, mom_file, daughter_file)]

	out_header = ['chr', 'pos', 'ref', 'dad_gt', 'mom_gt', 'child_gt']

	with gzip.open(output_file, 'wt') as outfile:
		print('\t'.join(out_header), file = outfile)

		for dad, mom, child in zip(*files):
			if dad.startswith('#'):
				if not mom.startswith('#') or not child.startswith('#'):
					raise ValueError('headers are of unequal length! figure yourself out')
				continue

			dad_info, mom_info, child_info = map(process_sample, (dad, mom, child))

			if filter_female_X(dad_info, mom_info, child_info):
				site_info = dad.rstrip().split('\t')
				site_pos = site_info[:2]
				ref = site_info[3][0]
				out_info = site_pos + [ref, dad_info[1], mom_info[1], child_info[1]]

				print('\t'.join(out_info), file = outfile)

	check_finish(files)

def main(sex, chromosome, dad_file, mom_file, child_file, output_file):
	if sex == 'female':
		print('starting female X')
		female_X(dad_file, mom_file, child_file, output_file)

	elif sex == 'male':
		if chromosome == 'chrX':
			print('starting male X')
			male_X(mom_file, child_file, output_file)
		elif chromosome == 'chrY':
			print('starting male Y')
			male_Y(dad_file, child_file, output_file)

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-c", "--chromosome", required=True, help="chromosome")
	# ap.add_argument("-x", "--sex", required=True, help="sex")
	# ap.add_argument("-f", "--father", required=True, help="dad vcf")
	# ap.add_argument("-m", "--mother", required=True, help="mom vcf")
	# ap.add_argument("-s", "--sample", required=True, help="sample vcf")
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()
	#
	#
	# main(args.sex, args.chromosome, args.father, args.mother, args.sample, args.output)
	main(snakemake.params.sex, snakemake.params.chrom, snakemake.input.father, snakemake.input.mother, snakemake.input.sample, snakemake.output.combined_data)
