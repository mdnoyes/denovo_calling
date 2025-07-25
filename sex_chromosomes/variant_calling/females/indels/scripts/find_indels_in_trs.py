import pandas as pd
from pybedtools import BedTool
from tqdm import tqdm
import argparse

def find_subunit(s):
	i = (s+s).find(s, 1, -1)
	return None if i == -1 else s[:i]

def read_variant_file(variant_file):
	variant_dict = {}
	with open(variant_file, 'r') as var_file:
		for line in var_file:
			var_data = line.rstrip().split()
			relevant_data = var_data[:5] + [var_data[7].split(';')[-1].split('=')[-1]]
			pos = int(var_data[1])
			bed_info = (var_data[0], pos - 1, pos)
			variant_dict[bed_info] = relevant_data

	return variant_dict

def make_tr_beds(tr_interval_dict):
	bed_dict = {}
	for chromsome, interval_list in tr_interval_dict.items():
		bed_dict[chromsome] = BedTool(interval_list)

	return bed_dict

def read_trf_file(trf_file):
	trf_dict = {}
	tr_interval_dict = {}
	with open(trf_file, 'r') as trf_file:
		for line in trf_file:
			trf_data = line.rstrip().split()
			bed_info = (trf_data[0], int(trf_data[1]), int(trf_data[2]))
			trf_dict[bed_info] = (trf_data[-1], trf_data[5], trf_data[1], int(trf_data[7]))
			if trf_data[0] not in tr_interval_dict:
				tr_interval_dict[trf_data[0]] = []
			tr_interval_dict[trf_data[0]].append(bed_info)

	return trf_dict, tr_interval_dict

def get_variants_in_trs(variant_dict, trf_bed):
	variant_bed = BedTool(list(variant_dict.keys()))
	tr_variants = variant_bed.intersect(trf_bed)
	tr_df = tr_variants.to_dataframe().drop_duplicates()
	return list(tr_df.itertuples(index=False, name=None))

def get_motif_info(variant_interval, trf_bed, trf_dict):
	motifs = []
	variant_bed = BedTool([variant_interval])
	tr = trf_bed.intersect(variant_bed, u=True)
	tr_df = tr.to_dataframe().drop_duplicates()
	for idx, row in tr_df.iterrows():
		tr_interval = (row['chrom'], row['start'], row['end'])
		motifs.append(trf_dict[tr_interval])

	return motifs

def get_indel(variant_info):
	if len(variant_info[3]) > 1:
		return variant_info[3], 'deletion'
	else:
		return variant_info[4], 'insertion'

def get_motif_change(indel_sequence, indel_type, motif):
	motif_count_change = indel_sequence.count(motif)
	if indel_type == 'insertion':
		return f'+{motif_count_change}'
	elif indel_type == 'deletion':
		return f'-{motif_count_change}'

def get_motif(variant_info, motif_list):
	indel, indel_type = get_indel(variant_info)

	indel_motif = find_subunit(indel)
	if not indel_motif:
		indel_motif = indel

	for motif_sequence, consensus_len, start_pos, purity in motif_list:
		if purity != 100:
			continue
		if motif_sequence == indel_motif:
			motif_size = len(indel_motif)
			motif_change = get_motif_change(indel, indel_type, indel_motif)

			if motif_size < 7:
				str_status = 'True'
			else:
				str_status = 'False'

			return [str_status, indel_motif, str(motif_size), start_pos, consensus_len, motif_change]

	return [str(len(indel))]

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	# ap.add_argument("-t", "--trf", required=True, help="Manifest of bam files")
	# ap.add_argument("-o", "--output", required=True, type=str, help="family id")
	# ap.add_argument("-s", "--strs", required=True, type=str, help="family id")
	# args = ap.parse_args()

	variants = read_variant_file(snakemake.input.indel_calls)
	# variants = read_variant_file(args.input)
	out_header = ['chr', 'pos', 'id', 'ref', 'alt', 'caller', 'in_tr']
	tr_header = out_header + ['str', 'repeat_motif', 'motif_length', 'motif_reference_start', 'consensus_length','denovo_motif_change']
	non_tr_header = out_header + ['indel_size']

	trs, trs_by_chromosome = read_trf_file(snakemake.params.trf)
	# trs, trs_by_chromosome = read_trf_file(args.trf)
	tr_beds = make_tr_beds(trs_by_chromosome)

	with open(snakemake.output.non_trs, 'w') as outfile, open(snakemake.output.trs, 'w') as str_file:
	# with open(args.output, 'w') as outfile, open(args.strs, 'w') as str_file:
		print('\t'.join(non_tr_header), file = outfile)
		print('\t'.join(tr_header), file = str_file)
		for interval, var_info in tqdm(variants.items(), total=len(variants)):
			motif_options = get_motif_info(interval, tr_beds[var_info[0]], trs)

			if len(motif_options) > 0:
				out_info = var_info + ['True']
				out_info = out_info + get_motif(var_info, motif_options)
				if len(out_info) == 13:
					print('\t'.join(out_info), file = str_file)
				else:
					print('\t'.join(out_info), file = outfile)
			else:
				out_info = var_info + ['False', str(len(get_indel(var_info)))]
				print('\t'.join(out_info), file = outfile)
