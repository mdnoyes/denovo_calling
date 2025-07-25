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

def filter_site(dad, mom, child):
	if dad[1] != '0/0':
		return False

	if mom[1] != '0/0':
		return False

	if int(child[0]) <= 20:
		return False

	return True

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-f", "--father", required=True, help="dad vcf")
	# ap.add_argument("-m", "--mother", required=True, help="mom vcf")
	# ap.add_argument("-s", "--sample", required=True, help="sample vcf")
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()

	files = [gzip.open(path, 'rt') for path in (snakemake.input.father, snakemake.input.mother, snakemake.input.sample)]
	# files = [gzip.open(path, 'rt') for path in (args.father, args.mother, args.sample)]

	out_header = ['chr', 'pos', 'ref', 'dad_gt', 'mom_gt', 'child_gt']

	with gzip.open(snakemake.output.combined_data, 'wt') as outfile:
	# with gzip.open(args.output, 'wt') as outfile:
		print('\t'.join(out_header), file = outfile)

		for dad, mom, child in zip(*files):
			if dad.startswith('#'):
				if not mom.startswith('#') or not child.startswith('#'):
					raise ValueError('headers are of unequal length! figure yourself out')
				continue

			dad_info, mom_info, child_info = map(process_sample, (dad, mom, child))

			if filter_site(dad_info, mom_info, child_info):
				site_info = dad.rstrip().split('\t')
				site_pos = site_info[:2]
				ref = site_info[3][0]
				out_info = site_pos + [ref, dad_info[1], mom_info[1], child_info[1]]

				print('\t'.join(out_info), file = outfile)

	check_finish(files)
