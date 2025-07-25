import argparse
import gzip

def get_file_list(*files):
	assert len(set([len(file) for file in files])) == 1

	f_list = []
	for trio in zip(*files):
		interval = (int(x) for x in trio[0].split('.')[-3].split('_')[1:])
		f_list.append((*interval, *trio))

	return sorted(f_list)

def in_interval(position, int_tuple):
	if position >= int_tuple[0] and position <= int_tuple[1]:
		return True

	return False

def advance_interval(position, interval_index, f_list):
	next_interval = interval_index
	while not in_interval(position, f_list[next_interval]):
		next_interval += 1

	return next_interval

def open_files(f_list):
	files = []

	for file in f_list:
		f = gzip.open(file, 'rt')
		next(f)
		files.append(f)

	return files

def ref_base_to_idx(base):
	match base:
		case 'A':
			return 1
		case 'C':
			return 2
		case 'G':
			return 3
		case 'T':
			return 4

def idx_to_base(idx):
	match idx:
		case 1:
			return 'A'
		case 2:
			return 'C'
		case 3:
			return 'G'
		case 4:
			return 'T'

def depth_filter(dad_bases, mom_bases, child_bases):
	for bases in [dad_bases, mom_bases, child_bases]:
		dp = sum([int(x) for x in bases[1:]])
		if dp == 0:
			return False

	return True

def get_child_alt_alleles(ref_base, child_bases):
	bases = [1,2,3,4]
	bases.remove(ref_base)

	observed_bases = []
	for base in bases:
		if int(child_bases[base]) >= 2:
			observed_bases.append(base)

	return observed_bases

def confirm_parent_no_alt(alt_base, dad_bases, mom_bases):
	if int(dad_bases[alt_base]) < 1 and int(mom_bases[alt_base]) < 1:
		return True

	return False

def check_for_denovo(site_info, dad_bases, mom_bases, child_bases):
	base_idx = ref_base_to_idx(site_info[2])
	alt_alleles = get_child_alt_alleles(base_idx, child_bases)
	if len(alt_alleles) != 1:
		return None

	if not confirm_parent_no_alt(alt_alleles[0], dad_bases, mom_bases):
		return None

	alt_allele = idx_to_base(alt_alleles[0])
	id = '_'.join(site_info[:3] + [alt_allele])
	out_data = [site_info[0], site_info[1], id, site_info[2], alt_allele]

	return out_data


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-r", "--ref", required=True, help="filtered haplotype caller output")
	# ap.add_argument("-f", "--father", required=True, nargs = "+", help="father files")
	# ap.add_argument("-m", "--mother", required=True, nargs = "+", help="mother files")
	# ap.add_argument("-c", "--child", required=True, nargs = "+", help="child files")
	# ap.add_argument("-o", "--output", required=True, type=str, help="output file")
	# ap.add_argument("-d", "--dnms", required=True, type=str, help="dnm file")
	#
	# args = ap.parse_args()

	file_list = get_file_list(snakemake.input.dad_intervals, snakemake.input.mom_intervals, snakemake.input.child_intervals)
	# file_list = get_file_list(args.father, args.mother, args.child)

	out_header = [
		'chr', 'pos', 'ref_base',
		'dad_A', 'dad_C', 'dad_G', 'dad_T',
		'mom_A', 'mom_C', 'mom_G', 'mom_T',
		'kid_A', 'kid_C', 'kid_G', 'kid_T'
	]

	dnm_header = ['chr', 'pos', 'id', 'ref', 'alt']

	# inititalize by opening the first interval file
	int_idx = 0
	read_pos = file_list[int_idx][0]
	opened_files = open_files(file_list[int_idx][2:])

	print(f"starting with interval 1 of {len(file_list)}")
	with gzip.open(snakemake.input.genotyped_sites, 'rt') as ref_file, gzip.open(snakemake.output.callable_sites, 'wt') as out_file, open(snakemake.output.dnms, 'w') as dnm_file:
	# with gzip.open(args.ref, 'rt') as ref_file, gzip.open(args.output, 'wt') as out_file, open(args.dnms, 'w') as dnm_file:
		print('\t'.join(out_header), file = out_file)
		print('\t'.join(dnm_header), file = dnm_file)

		next(ref_file)
		for line in ref_file:
			site = line.rstrip().split('\t')
			pos = int(site[1])

			new_idx = advance_interval(pos, int_idx, file_list)
			if new_idx != int_idx:
				print(f"moving on to interval {new_idx + 1} of {len(file_list)}")
				int_idx = new_idx
				opened_files = open_files(file_list[int_idx][2:])
				read_pos = file_list[int_idx][0]

			while read_pos < pos:
				for file in opened_files:
					next(file)
				read_pos += 1

			read_pos += 1

			dad, mom, child = [next(f).rstrip().split('\t') for f in opened_files]

			# confirm that we're at the right position
			assert len({dad[0], mom[0], child[0]}) == 1, "read files not in the same position"
			assert int(dad[0]) == pos, "reads file are out of sync with ref file"

			# only consider sites with reads in all three samples
			if depth_filter(dad, mom, child):
				out_data = site[:3] + [*dad[1:], *mom[1:], *child[1:]]
				print('\t'.join(out_data), file = out_file)

				# check if there's an alternate allele in the child and not the parents
				denovo_data = check_for_denovo(site, dad, mom, child)
				if denovo_data:
					print('\t'.join(denovo_data), file = dnm_file)
