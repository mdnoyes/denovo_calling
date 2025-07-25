import pandas as pd
import argparse
from tqdm import tqdm
import gzip

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
	dad = row['fa'].split(":")[0].replace("|", "/")
	mom = row['mo'].split(":")[0].replace("|", "/")
	kid = row['kid'].split(":")[0].replace("|", "/")

	if dad == '0/0' and mom == '0/0' and kid == '0/0':
		return "TRANSMITTED=no", "INH=not_variant"
	elif dad == '0/0' and (mom == '0/1' or mom == '1/1') and kid == '0/0':
		return "TRANSMITTED=no", "INH=mother_only"
	elif (dad == '0/1' or dad == '1/1') and mom == '0/0' and kid == '0/0':
		return "TRANSMITTED=no", "INH=father_only"
	elif dad == '0/0' and mom == '0/0' and (kid == '0/1' or kid == '1/1'):
		return "TRANSMITTED=no", "INH=denovo_kid"
	elif dad == './.' and mom == '0/0' and (kid == '0/1' or kid == '1/1'):
		return "TRANSMITTED=no", "INH=denovo_kid"
	elif (dad == '0/1' or dad == '1/1') and mom == '0/0' and kid == '0/1':
		return "TRANSMITTED=yes", "INH=fa_to_kid"
	elif dad == '0/0' and (mom == '0/1' or mom == '1/1') and kid == '0/1':
		return "TRANSMITTED=yes", "INH=mo_to_kid"
	elif dad == '1/1' and mom == '1/1' and kid == '1/1':
		return "TRANSMITTED=yes", "INH=unknown"
	elif (dad == '0/1' or dad == '1/1') and (mom == '0/1' or mom == '1/1') and (kid == '0/1' or kid == '1/1'):
		return "TRANSMITTED=yes", "INH=unknown"
	elif (dad == '0/1' or dad == '1/1') and (mom == '0/1' or mom == '1/1') and kid == '0/0':
		return "TRANSMITTED=no", "INH=parents_only"
	else:
		return "TRANSMITTED=unknown", "INH=unknown"

def main(vcf_file, family_id, output_file):
	print("\nReading VCF...")
	header, vcf = read_vcf(vcf_file)
	sample  = header[-1].rstrip().split('\t')[-1]

	column_names = header[-1].rstrip().split('\t')
	column_names[-3:] = ['fa', 'mo', 'kid']
	df = pd.DataFrame(vcf)
	df.columns = column_names
	df.columns = df.columns.str.replace(f'{family_id}_', '')

	transmission = []
	inheritance = []

	print("\nAssigning inheritance...")

	for index, row in tqdm(df.iterrows(), total = len(df)):
		tra, inh = assign_variant(row)
		row["INFO"] = row["INFO"] + ";" + tra + ";" + inh
		row["ID"] = '_'.join([sample, row['#CHROM'], row['POS'], row['REF'], row['ALT']])

	denovo = df[df["INFO"].str.contains("denovo_kid")]

	print("\nWriting output file...")

	denovo.to_csv(output_file, sep = '\t', index = None, header = False)


if __name__ == "__main__":
	# to run from the command line, use commented out lines
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-i", "--input", required=True, help="vcf of all family variant calls")
	# ap.add_argument("-f", "--family", required=True, help="family ID")
	# ap.add_argument("-o", "--output", required=True, help="output file")
	# args = ap.parse_args()
	#
	# main(args.input, args.family, args.output)
	main(snakemake.input.vcf, snakemake.params.family, snakemake.output.naive_calls)
