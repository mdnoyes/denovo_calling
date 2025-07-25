import pysam
import argparse
from tqdm import tqdm
from collections import Counter
import pandas as pd
import math

def import_sources(source_file: str):
	"""
	Import file of hifi data sources to dictionary of read names mapped to 'blood' or 'cells'.
	:param source_file: Path to tab-separated file with read prefix and source info. Only child data is necessary.
	:return: Dictionary of read prefix string mapped to 'blood' or 'cells'.
	"""
	source_dict = {}
	with open(source_file, 'r') as sources:
		header = sources.readline().rstrip().split('\t')
		sample = header.index('sample')
		source = header.index('source')
		prefix = header.index('readname_prefix')

		for line in sources:
			bam_info = line.rstrip().split('\t')
			source_dict[(bam_info[sample], bam_info[prefix])] = bam_info[source]

	return source_dict

def get_bam(sample, bams, sequence_platform):
	"""
	Retrieves the path to the sample's bam file.
	:param sample: Sample to find bam file for.
	:param bams: Manifest dataframe of bam files.
	:param sequence platform: illumina, ont, or hifi
	:returns: PySam AlignmentFile object of the sample's read data.
	"""
	fm, sm = sample.split("_")
	bam_path = bams.loc[(bams['sample'] == sm) & (bams['family'] == fm)][sequence_platform].values[0]
	bam = pysam.AlignmentFile(bam_path, "rc")
	return bam

def calculate_indel_size(ref, alt):
	"""
	Compare ref and alt allele to calculate indel size.
	:param ref: Refernce allele.
	:param alt: Alternate allele.
	:return: Integer of number of base pairs inserted (>0) or deleted (<0).
	"""
	size = len(alt) - len(ref)
	return size

def get_flanking_sequences(chrom, pos, indel_size, reference):
	"""
	Get the reference base pairs before/after the indel.
	:param chrom: Chromosome of the indel.
	:param pos: Starting coordinate of the indel.
	:param indel_size: Number of inserted or deleted bases.
	:param reference: PySam object of reference genome.
	:return: Strings of the padded number of base pairs before and after the indel.
	"""
	if indel_size < 0:
		end = pos - indel_size
	else:
		end = pos

	prev_seq = reference.fetch(chrom, pos - 10, pos)
	next_seq = reference.fetch(chrom, end, end + 10)
	return prev_seq.upper(), next_seq.upper()

def filter_read(read, platform):
	"""
	Checks that reads pass mapq filter and is a primary alignment.
	:return: True if variant passes filter, False if not.
	"""
	if read.is_secondary or read.is_supplementary:
		return False

	if platform != 'illumina':
		if read.mapping_quality < 40:
			return False

	return True

def define_flanks(read, cigar_position, lead_seq, tail_seq, indel_length, padding):
	"""
	Determine the expected and actual sequence on either side of the variant start site.
	:param read: PySam aligned segment.
	:param cigar_position: Position of the variant start site in the read.
	:param lead_seq: Predicted 10bp before the variant based on the reference.
	:param tail_seq: Predicted 10bp after the variant based on the reference.
	:param indel_length: Integer size of indel.
	:param padding: Number of bp on either side to check
	:return: Expected and observed sequence before and after variant, and the starting position of the right interval.
	"""
	expected_left = lead_seq[(-1 * padding):]
	expected_right = tail_seq[:padding]
	left_flank = read.query_sequence[cigar_position - 4: cigar_position]
	right = cigar_position + max([0, indel_length])
	right_flank = read.query_sequence[right: right + 4]

	return expected_left, left_flank, expected_right, right_flank, right

def check_cigar_string(read, lead_seq, tail_seq, indel_length, indel_position):
	"""
	Iterate through cigar operations in read to find the variant start site and sequence.
	:param read: PySam aligned segment option.
	:param lead_seq: 10bp of sequence before the variant, based on the reference.
	:param tail_seq: 10bp of sequence after the variant, based on the refference.
	:param indel_length: Integer size of indel.
	:param indel_position: Reference coordinate of indel.
	:return: Sequence of inserted or deleted bases between flanking seequence, or None if no indel is observed at variant site.
	"""
	cigar_position = 0
	reached_start = False
	iterations = 0
	inferred_coord = read.reference_start

	if indel_length < 0:
		expected_operation = 2
	else:
		expected_operation = 1

	if not read.cigartuples:
		return None

	for operation, length in read.cigartuples:

		# this last tuple will tell us if there's inserted or deleted sequence that corresponds to our indel
		if reached_start == True:
			if indel_position == inferred_coord:
				expected_left, left_flank, expected_right, right_flank, right_start = define_flanks(read, cigar_position, lead_seq, tail_seq, indel_length, 4)
			else:
				expected_left, left_flank, expected_right, right_flank, right_start = define_flanks(read, cigar_position, lead_seq, tail_seq, indel_length, 10)

			if left_flank == expected_left and right_flank == expected_right:
				sequence = read.query_sequence[cigar_position : right_start]
				if operation == expected_operation and length == abs(indel_length):
					# print(read.qname, left_flank, right_flank, operation, length, inferred_coord)
					return sequence

			reached_start = False

		# Once we think we reach the start, we still want to advance to the next cigar tuple before checking the operation
		if abs(indel_position - inferred_coord) <= length :
			reached_start = True

		# If there is an insertion (1) or soft clip (4), we advance in the read but not in the reference
		if operation not in (1,4):
			inferred_coord += length

		# If there is a deletion (2) or reference skip (3), we advance in the reference but not in the read
		if operation not in (2,3):
			cigar_position += length

		iterations += 1

	return None

def check_hifi_source(read, platform, read_sources, sample_name):
	"""
	Determine whether read is from blood or cell line derived data.
	:param read: PySam aligned segment.
	:param platform: Sequencing platform (one of "hifi", "ont", "illumina").
	:param read_sources: Dict of sample and read names mapped to data source.
	:param sample_name: Name of de novo sample.
	:return: True if read is from blood, False if read is from cells.
	"""
	if platform == 'hifi':
		readname = read.query_name
		prefix = readname.split('/')[0]
		source = read_sources[(sample_name, prefix)]
		if source == 'cells':
			return False

	return True

def check_reads(bam_info, chrom, pos, alt, indel_size, lead_seq, tail_seq, platform, hifi_sources, sample_id):
	"""
	Iterate through all reads aligned to indel coordinate and find the inserted or deleted sequence.
	:param bam_info: Tuple of sample name and AlignmentFile object.
	:param chrom: Chromosome name.
	:param pos: Indel start coordinate.
	:param alt: Alternate allele.
	:param indel_size: Number of inserted or deleted bases.
	:param lead_seq: The of sequence prior to the indel.
	:param tail_seq: The of sequence after the indel.
	:param platform: Sequencing platform (one of "hifi", "ont", "illumina").
	:param hifi_sources: Dict of sample and read names mapped to data source.
	:param sample id: Name of de novo sample.
	:return: Number of inserted or deleted bases in each read.
	"""
	sample_name, bam = bam_info

	if indel_size < 0:
		expected_alt = ''
	elif indel_size > 0:
		expected_alt = alt[1:]
	else:
		return(0,0,0)

	try:
		total_reads = 0
		supporting_blood_reads = 0
		supporting_cell_reads = 0

		for aligned_segment in bam.fetch(chrom, pos - 1, pos):

			if not filter_read(aligned_segment, platform):
				continue

			total_reads += 1

			sequence = check_cigar_string(aligned_segment, lead_seq, tail_seq, indel_size, pos)
			# print(aligned_segment.qname, sequence, alt)

			if sequence == expected_alt:
				# print(aligned_segment.qname, sequence, expected_alt)

				if sample_name == sample_id:
					if not check_hifi_source(aligned_segment, platform, hifi_sources, sample_name):
						supporting_cell_reads += 1
						continue

				supporting_blood_reads += 1

		return (supporting_blood_reads, supporting_cell_reads, total_reads)

	# # if the BAM is corrupted, return NA
	except OSError:
		return ('NA', 'NA', 'NA')

def main(input, hifi_source_file, reference, manifest, family, sample, platform, parent, output):
	hifi_sources = import_sources(hifi_source_file)

	variant_df = pd.read_csv(input, sep = '\t')
	ref = pysam.FastaFile(reference)
	bam_df = pd.read_csv(manifest, sep='\t', dtype=str)
	samples = [f'{family}_{sample}', parent]

	sample_bams = []
	for sample_id in samples:
		sample_bams.append((sample_id, get_bam(sample_id, bam_df, platform)))

	samples_out = open(output, 'w')

	out_header = ['sample'] + list(variant_df.columns) + ['supporting_reads', 'cell_line_supporting_reads', 'total_reads']
	print('\t'.join(out_header), file = samples_out)

	for idx, row in tqdm(variant_df.iterrows(), total=len(variant_df)):
		indel = calculate_indel_size(row['ref'], row['alt'])
		prev_seq, next_seq = get_flanking_sequences(row['chr'], row['pos'], indel, ref)

		for bam in sample_bams:
			output_info = [bam[0]] + [str(x) for x in row.to_list()][:-1] + [str(indel)]

			allele_counts = check_reads(bam, row['chr'], row['pos'], row['alt'], indel, prev_seq, next_seq, platform, hifi_sources, f'{family}_{sample}')

			output_info += [str(x) for x in allele_counts]
			print('\t'.join(output_info), file = samples_out)

	samples_out.close()



if __name__ == "__main__":
	# TO RUN FROM COMMAND LINE INSTEAD, USE COMMENTED OUT CODE INSTEAD OF SNAKEMAKE CODE
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	# ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	# ap.add_argument("-r", "--reference", required=False, help="Reference genome fasta")
	# ap.add_argument("-f", "--family", required=True, type=str, help="family id")
	# ap.add_argument("-s", "--sample", required=True, help="sample id")
	# ap.add_argument("-p", "--platform", required=True, help="sequencing platform")
	# ap.add_argument("-hs", "--hifi_sources", required=True, help="hifi source data")
	# ap.add_argument("-o", "--output", required=True, help="path to output file")
	# args = ap.parse_args()

	# main(args.input, args.hifi_sources, args.reference, args.bams, args.family, args.sample, args.platform, args.output)
	main(snakemake.input.indels, snakemake.params.hifi_sources, snakemake.params.reference, snakemake.params.bam_manifest, snakemake.params.family, snakemake.params.sample, snakemake.params.platform, snakemake.params.parent, snakemake.output.read_summary)
