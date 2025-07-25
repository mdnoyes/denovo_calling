import pysam
import argparse
import gzip
from tqdm import tqdm
from collections import defaultdict, OrderedDict

def import_sources(source_file):
	source_dict = {}
	with open(source_file, 'r') as sources:
		header = sources.readline().rstrip().split('\t')
		source = header.index('source')
		prefix = header.index('readname_prefix')

		for line in sources:
			bam_info = line.rstrip().split('\t')
			source_dict[bam_info[prefix]] = bam_info[source]

	return source_dict

def get_interval_info(interval_string):
	interval = f"{interval_string}".split("_")

	return interval[0], int(interval[1]), int(interval[2])

def get_bam(bam_path):
	return pysam.AlignmentFile(bam_path, "rc")

def get_number_reads_in_bam(bam_path, chrom):
	chromosome_counts = pysam.idxstats(bam_path).split('\n')
	chrom_reads = [x for x in chromosome_counts if x.startswith(f'{chrom}\t')][0]
	return int(chrom_reads.split('\t')[2])

def check_hifi_source(read, read_sources):
	readname = read.alignment.query_name
	prefix = readname.split('/')[0]

	try:
		source = read_sources[prefix]
	except KeyError:
		return True

	if source == 'cells':
		return False

	return True

def get_base_counts(column, sources):
	position_bases = OrderedDict(A=0, C=0, G=0, T=0)

	for pileupread in column.pileups:
		if not check_hifi_source(pileupread, sources):
			continue
		if not pileupread.is_del and not pileupread.is_refskip:
			position_bases[pileupread.alignment.query_sequence[pileupread.query_position]] += 1

	return position_bases

if __name__ == "__main__":
# 	ap = argparse.ArgumentParser()
# 	ap.add_argument("-i", "--interval", required=True, help="interval")
# 	ap.add_argument("-b", "--bam", required=True, help="input bam file")
# 	ap.add_argument("-s", "--sample", required=True, help="sample id")
# 	ap.add_argument("-hs", "--hifi_sources", required=True, help="hifi sources")
# 	ap.add_argument("-o", "--output", required=True, help="Output file")
#
# 	args = ap.parse_args()

	print("importing hifi sources...")

	hifi_sources = import_sources(snakemake.params.hifi_sources)
	# hifi_sources = import_sources(args.hifi_sources)

	print(f'determining interval...')
	chrom, start, end = get_interval_info(snakemake.params.interval)
	# chrom, start, end = get_interval_info(args.interval)

	print(f'importing bam for sample {snakemake.params.sample}...')
	sample_bam = get_bam(snakemake.input.subsetted_bam)
	# print(f'importing bam for sample {args.sample}...')
	# sample_bam = get_bam(args.bam)
	pileup_columns = len(list(sample_bam.pileup(chrom, start-1, end, truncate = True)))

	print(f'querying reads for sample {snakemake.params.sample}...')
	# print(f'querying reads for sample {args.sample}...')

	pos = start - 1
	expected_position = start

	with gzip.open(snakemake.output.read_data, 'wt') as outfile:
	# with gzip.open(args.output, 'wt') as outfile:
		print('\t'.join(['pos', 'A', 'C', 'G', 'T']), file = outfile)
		for pileupcolumn in tqdm(sample_bam.pileup(chrom, start-1, end, min_base_quality = 20, truncate = True), total = pileup_columns):
			pos = pileupcolumn.pos + 1
			while pos > expected_position:
				out_info = [str(expected_position), '0', '0', '0', '0']
				print('\t'.join(out_info), file = outfile)
				expected_position += 1

			if pos == expected_position:
				base_counts = get_base_counts(pileupcolumn, hifi_sources)
				out_info = [str(pos)] + [str(x) for x in base_counts.values()]
				print('\t'.join(out_info), file = outfile)
				expected_position += 1

		if pos != end:
			while pos < end:
				pos += 1
				out_info = [str(pos), '0', '0', '0', '0']
				print('\t'.join(out_info), file = outfile)
