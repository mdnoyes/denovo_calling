import pysam
import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm


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

def get_sample_list(bams, included_samples, family_id):
	"""
	Generates list of samples to check reads for.
	:param bams: Manifest dataframe of bam files.
	:param included_samples: String 'family' if just related individuals, 'all' for every sample in the manifest.
	:param family_id: String of sample's family ID.
	:returns: List of sample names for analysis.
	"""
	sample_list = []
	for idx, row in bams.iterrows():
		sm = row["sample"]; fm = row["family"]
		sm_id = fm + "_" + sm
		if included_samples == 'family':
			if fm == family_id:
				sample_list.append(sm_id)
		elif included_samples == 'all':
			sample_list.append(sm_id)

	return sample_list

def get_mean_qual(qualities):
	"""
	Computes the mean base quality for reads.
	:param qualities: List of integer base qualities from read data.
	:return: Number of quality values in list and their mean.
	"""
	qual_count = len(qualities)

	if qual_count != 0:
		return qual_count, np.mean(qualities)

	return qual_count, "NA"

def check_snv_read(read, position, ref, alt, ref_quality, alt_quality, other_quality, map_quality):
	"""
	Evaluate a read to record its quality data and the base at the query position.
	:param read: AlignedSegment object of read aligned to query position.
	:param position: Integer query position of DNM in reference genome.
	:param ref: Reference allele, A/C/G/T.
	:param alt: de novo alternate allele, A/C/G/T.
	:param ref_quality: List of base quality values for reads with reference allele at query position.
	:param alt_quality: List of base quality values for reads with alternate allele at query position.
	:param other_quality: List of base quality values for reads with other allele at query position.
	:param map_quality: List of map quality values for reads at query position.
	:returns: Name of the read, base and its quality at query position, median base quality across the read, length of the read, and its map quality.
	"""
	try:
		median_quality = np.median(read.query_qualities)
	except:
		return read.query_name, "NA", "NA", "NA", 0, 0, 0

	read_length = read.query_length
	mapq = read.mapping_quality
	map_quality.append(mapq)

	try:
		ref_idx = read.get_reference_positions(full_length=True).index(position - 1)
	except:
		return read.query_name, "deletion", "NA", "NA", median_quality, read_length, mapq

	base = read.query_sequence[ref_idx]
	quality = read.query_qualities[ref_idx]
	read_pos = ref_idx / read_length

	if base == ref:
		correct_quality = ref_quality
	elif base == alt:
		correct_quality = alt_quality
	else:
		correct_quality = other_quality

	correct_quality.append(quality)

	return read.query_name, base, quality, read_pos, median_quality, read_length, mapq

def get_snv_data(bam, site_id, chrom, pos, ref, alt, reads_out, samples_out):
	"""
	Write individual and summary stats to files for all the reads aligned to a query position in one sample.
	:param bam: AlignmentFile object of all reads from a sample's BAM.
	:site_id: Unique ID for DNM site, with sample name, chromosome, position, ref, and alt alleles.
	:chrom: Chromosome name for DNM site.
	:pos: Integer reference position of DNM site.
	:param ref: Reference allele, A/C/G/T.
	:param alt: de novo alternate allele, A/C/G/T.
	:param reads_out: Opened file of read data.
	:param samples_out: Opened file of summary data for each sample.
	"""
	sample_name, bam = bam
	ref_qual = []
	alt_qual = []
	other_qual = []
	map_qual = []

	for aligned_segment in bam.fetch(chrom, pos - 1, pos):
		line = "\t".join((
			site_id,
			sample_name,
			*map(str, check_snv_read(aligned_segment, pos, ref, alt, ref_qual, alt_qual, other_qual, map_qual)),
		))
		reads_out.write(line)
		reads_out.write("\n")


	ref_count, avg_ref_q = get_mean_qual(ref_qual)
	alt_count, avg_alt_q = get_mean_qual(alt_qual)
	other_count, avg_other_q = get_mean_qual(other_qual)
	avg_mapq = np.mean(map_qual)

	sample_line = "\t".join(map(str, (site_id, sample_name,
							ref_count, avg_ref_q,
							alt_count, avg_alt_q,
							other_count, avg_other_q,
							avg_mapq)))
	samples_out.write(sample_line)
	samples_out.write("\n")


def main(variants, bams, family, platform, denovo_sample, who):
	variant_df = pd.read_csv(variants, sep='\t', header = None)

	bam_df = pd.read_csv(bams, sep='\t')
	bam_df = bam_df.astype(str)

	samples = get_sample_list(bam_df, who, family)

	if who == 'family':
		reads_out = open(f'read_data/{family}_{denovo_sample}_{platform}.read_data.tsv', 'w')
		samples_out = open(f'read_data/{family}_{denovo_sample}_{platform}.read_summary.tsv', 'w')
	elif who == 'all':
		reads_out = open(f'read_data/{family}_{denovo_sample}_{platform}.all_samples.read_data.tsv', 'w')
		samples_out = open(f'read_data/{family}_{denovo_sample}_{platform}.all_samples.read_summary.tsv', 'w')

	# write headers to output files

	print("\t".join(("id", "sample", "read_id", "allele", "base_quality", "read_position",
				"median_read_quality", "read_length", "mapq")), file = reads_out)
	print("\t".join(("id", "sample", "ref_count", "mean_ref_qual", "alt_count", "mean_alt_qual",
				"other_count", "mean_other_qual", "avg_mapq")), file = samples_out)

	for sample_id in samples:
		print(f"currently working on sample {sample_id}, now processing:")
		sample_bam = (sample_id, get_bam(sample_id, bam_df, platform))

		for idx, row in tqdm(variant_df.iterrows(), total=len(variant_df)):
			if row[0] == 'chr':
				continue
			get_snv_data(sample_bam, row[2], row[0], int(row[1]), row[3], row[4], reads_out, samples_out)

	reads_out.close()
	samples_out.close()

if __name__ == "__main__":
	# to run from the command line, use commented out lines
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="Candidate de novo file")
	# ap.add_argument("-b", "--bams", required=True, help="Manifest of bam files")
	# ap.add_argument("-f", "--family", required=True, type=str, help="family id")
	# ap.add_argument("-c", "--child", required=True, type=str, help="sample id")
	# ap.add_argument("-p", "--platform", required=True, help="sequencing platform")
	# ap.add_argument("-w", "--who", required=True, help="who to run this on")
	#
	#
	# args = ap.parse_args()
	#
	# main(args.input, args.bams, args.family, args.platform, args.child, args.who)
	main(snakemake.input.snv_calls, snakemake.params.bam_manifest, snakemake.params.family, snakemake.params.platform, snakemake.params.child, snakemake.params.who)
