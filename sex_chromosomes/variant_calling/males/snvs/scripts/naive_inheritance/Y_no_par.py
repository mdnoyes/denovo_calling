import sys
import pandas as pd
import argparse
import gzip
from tqdm import tqdm

def read_vcf(file):
	"""
	Reads in vcf.
	:param file: Vcf file.
	:returns: List of header rows and list of variant rows.
	"""
	header = []
	vcf = []

	with gzip.open(file, 'rt') as f:
		for line in f:
			if line.startswith("#"):
				header.append(line)
			else:
				vcf.append(line.rstrip().split('\t'))
	return header, vcf


def assign_variant(row):
	"""
	Assigns inheritance and transmission based on genotypes.
	:param row: Row of dataframe corresponding to variant.
	:returns: Transmission status ("TRANSMITTED=no/yes") and inheritance ("INH=mother_only/denovo_pro/etc."")
	"""
	dad = row['fa'].split(":")[0]
	kid = row['kid'].split(":")[0]

	if dad == '0' and kid == '0':
		return "TRANSMITTED=no", "INH=not_variant"
	elif dad == '1' and kid == '1':
		return "TRANSMITTED=yes", "INH=fa_to_kid"
	elif dad == '0' and kid == '1':
		return "TRANSMITTED=no", "INH=denovo_kid"
	elif dad == '1' and kid == '0':
		return "TRANSMITTED=no", "INH=denovo_kid"
	else:
		return "TRANSMITTED=unknown", "INH=unknown"


if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="vcf of all family variant calls")
	# ap.add_argument("-s", "--snvs", required=True, help="output file")
	# ap.add_argument("-d", "--indels", required=True, help="output file")
	# args = ap.parse_args()

	print("\nReading VCF...")

	header, vcf = read_vcf(snakemake.input.y_no_par_vcf)
	# header, vcf = read_vcf(args.input)

	column_names = header[-1].rstrip().split('\t')
	sample = column_names[-1]
	column_names[-2:] = ['fa', 'kid']
	df = pd.DataFrame(vcf)
	df.columns = column_names

	transmission = []
	inheritance = []

	print("\nAssigning inheritance...")

	for index, row in tqdm(df.iterrows(), total = len(df)):
		tra, inh = assign_variant(row)
		row["INFO"] = row["INFO"] + ";" + ';'.join([tra, inh, 'CALLER=GATK'])
		chrom = "chrY"
		row["ID"] = "_".join([sample, chrom, str(row["POS"]), row["REF"], row["ALT"]])
		if row['fa'].split(":")[0] == '1' and row['kid'].split(":")[0] == '0':
			ref = row["ALT"]
			alt = row["REF"]
			row["ALT"] = alt
			row["REF"] = ref
			row["ID"] = "_".join([sample, chrom, str(row["POS"]), ref, alt])

	denovo = df[df["INFO"].str.contains("denovo_kid")]

	snp_mask = (denovo['REF'].str.len() == 1) & (denovo['ALT'].str.len() == 1)
	indel_mask = (denovo['REF'].str.len() > 1) | (denovo['ALT'].str.len() > 1)

	print("\nWriting output file...")

	denovo[snp_mask].to_csv(snakemake.output.snvs, sep = '\t', index = None, header = False)
	denovo[indel_mask].to_csv(snakemake.output.indels, sep = '\t', index = None, header = False)
	# denovo[snp_mask].to_csv(args.snvs, sep = '\t', index = None, header = False)
	# denovo[indel_mask].to_csv(args.indels, sep = '\t', index = None, header = False)
