import argparse
import numpy as np

def get_var_info(var):
	var_info = var[0:5]
	var_info.append(var[7].split('=')[-1])
	return var_info

def import_unsourced_reads(readfile, child_id):
	read_dict = {}
	with open(readfile, 'r') as reads:
		header = reads.readline().rstrip().split('\t')[2:]
		for line in reads:
			read = line.rstrip().split('\t')
			var_id, sample = read[0:2]
			if (var_id, sample) not in read_dict:
				read_dict[(var_id, sample)] = []
			read_dict[(var_id, sample)].append((read[2:]))

	return header, read_dict

def filter_data(read_list, bq_idx, mq_idx):
	filtered_data = []

	for read in read_list:
		try:
			bq = int(read[bq_idx])
			mq = int(read[mq_idx])
		except ValueError:
			continue
		if bq >= 20 and mq >= snakemake.params.mapq:
			filtered_data.append(read)

	return filtered_data

def get_filtered_reads(site_id, sample_id, child_id, read_data, header):
	base_qual = header.index('base_quality')
	mapq = header.index('mapq')

	try:
		sample_data = read_data[(site_id, sample_id)]
	except KeyError:
		return []

	# if sample_id == child_id:
	filtered_reads = filter_data(sample_data, base_qual, mapq)
	return filtered_reads

	return sample_data

def filter_child_reads(site_id, sample_id, read_data, header):
	base_qual = header.index('base_quality')
	mapq = header.index('mapq')

	try:
		blood_data = read_data[(site_id, sample_id, 'blood')]
		filtered_blood = filter_data(blood_data, base_qual, mapq)
	except KeyError:
		filtered_blood = []

	try:
		cell_data = read_data[(site_id, sample_id, 'cells')]
		filtered_cells = filter_data(cell_data, base_qual, mapq)
	except KeyError:
		filtered_cells = []

	return filtered_blood, filtered_cells

def allele_count(read_data, allele_base, header):
	allele = header.index('allele')

	return [x[allele] for x in read_data].count(allele_base)

def allele_balance(read_data, alt_alleles):
	total_alleles = len(read_data)
	return round(alt_alleles / total_alleles, 5)

def get_allele_stats(reads, ref_allele, alt_allele, header):
	if len(reads) > 0:
		alt_allele_count = allele_count(reads, alt_allele, header)
		ref_allele_count = allele_count(reads, ref_allele, header)
		ab = allele_balance(reads, alt_allele_count)
		dp = len(reads)
	else:
		alt_allele_count = 0
		ref_allele_count = 0
		ab = np.nan
		dp = 0

	return (ref_allele_count, alt_allele_count), ab, dp

def get_status(mom_count, child_count, child_dp, pf):
	if pf == 'hifi':
		inherited_threshold = 1
		fp_threshold = 2
	else:
		inherited_threshold = 2
		fp_threshold = 1

	if np.nan in mom_count:
		return "missing_parent_data"
	if child_dp == 0:
		return "missing_child_data"
	if mom_count[1] > inherited_threshold:
		return "inherited"
	if child_count[1] >= fp_threshold and child_count[0] >= 2:
		return "biallelic"

	if mom_count[1] <= inherited_threshold and child_count[1] >= fp_threshold:
		return "true_de_novo"

	return "false_positive"



if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	# ap.add_argument("-v", "--variants", required=True, help="Candidate de novo file")
	# ap.add_argument("-r", "--reads", required=True, help="File of read data")
	# ap.add_argument("-f", "--father", required=True, help="father ID", type=str)
	# ap.add_argument("-c", "--child", required=True, help="sample ID", type=str)
	# ap.add_argument("-p", "--platform", required=True, help="Platform")
	# ap.add_argument("-o", "--output", required=True, help="Output file")
	# ap.add_argument("-mq", "--mapq_threshold", required=False, type=int, default=59)
	# args = ap.parse_args()

	sample_ids = {'mom': snakemake.params.mother, 'child': snakemake.params.child}
	# sample_ids = {'dad': args.father, 'child': args.child}

	read_header, reads = import_unsourced_reads(snakemake.input.read_data, snakemake.params.child)
	# read_header, reads = import_unsourced_reads(args.reads, args.child)
	all_variants = []

	# with open(args.input, 'r') as variants, open(args.output, 'w') as outfile:
	with open(snakemake.input.snvs, 'r') as variants, open(snakemake.output.validated, 'w') as outfile:
		out_columns = ['chr', 'pos', 'id', 'ref', 'alt', 'callers', 'mom_ref_count', 'mom_alt_count', 'child_ref_count', 'child_alt_count', 'child_allele_balance', 'child_depth', 'filtered_status']
		print('\t'.join(out_columns), file=outfile)
		for line in variants:
			var_data = line.rstrip().split('\t')
			site = var_data[2]
			ref = var_data[3]
			alt = var_data[4]

			variant_info = get_var_info(var_data)

			sample_counts = {}

			for sample_name, sample_id in sample_ids.items() :
				# filtered_reads = get_filtered_reads(site, sample_id, args.child, reads, read_header)
				filtered_reads = get_filtered_reads(site, sample_id, snakemake.params.child, reads, read_header)

				if sample_name == 'mom':
					if len(filtered_reads) > 0:
						num_alt_alleles = allele_count(filtered_reads, alt, read_header)
						num_ref_alleles = allele_count(filtered_reads, ref, read_header)
						mom_count = (num_ref_alleles, num_alt_alleles)
					else:
						mom_count = (np.nan, np.nan)
				else:
					child_count, child_ab, child_dp = get_allele_stats(filtered_reads, ref, alt, read_header)

			# validation = get_status(mom_count, child_count, child_dp, args.platform)
			validation = get_status(mom_count, child_count, child_dp, snakemake.params.platform)
			variant_info += [str(x) for x in [mom_count[0], mom_count[1], child_count[0], child_count[1], child_ab, child_dp, validation]]
			print('\t'.join(variant_info), file=outfile)
