import argparse
from collections import Counter

def open_files(f_list):
	files = []

	for file in f_list:
		f = open(file, 'r')
		next(f)
		files.append(f)

	return files

def final_validation(hifi, ont):
	if hifi[-1] == "cannot_determine" and ont[-1] == "cannot_determine" and 'NA' not in hifi:
		dad_ref, dad_alt, mom_ref, mom_alt = [int(x) for x in hifi[6:10]]
		if dad_alt > 1 and mom_alt > 1:
			return False
		if dad_alt == 0 and mom_alt == 0:
			return False

	return True

def resolve_variant_type(hifi, ont, ab):
	# prefer hifi or ont assignment, default to allele balance when neither is available
	if hifi == "cannot_determine":
		if ont != "cannot_determine":
			return ont
		return ab

	if ont ==  "cannot_determine":
		return hifi

	# if hifi and ont agree, this is the best case
	if hifi == ont:
		return hifi

	# if they disagree, we assume pzm if there's any evidence. if not, we assume dnm
	if hifi == "dnm":
		if ont == "pzm" and ab == "pzm":
			return "pzm"
		else:
			return "dnm"

	if hifi == "pzm":
		return "pzm"

	raise Exception(f"unhandled type case: hifi {hifi} ont {ont} ab {ab}")

# def resolve_inheritance(hifi, ont, whatshap):
def resolve_inheritance(hifi, ont):
	if hifi == "cannot_determine":
		# if ont != "cannot_determine":
		return ont
		# return whatshap

	if ont ==  "cannot_determine":
		return hifi

	if hifi == ont:
		return hifi

	if hifi != ont:
		# if hifi == whatshap:
		return hifi
		# if ont == whatshap:
		# 	return ont

		return "cannot_determine"

	raise Exception(f"unhandled inheritance case: hifi {hifi} ont {ont}")

def main(snv_file, hifi_file, ont_file, output_file):
	haplotype_files = open_files([hifi_file, ont_file])

	types = []
	parents = []

	with open(snv_file, 'r') as snv_file, open(output_file, 'w') as out_file:
		header = next(snv_file).rstrip().split('\t')
		out_header = [*header[:5], 'variant_type', 'inheritance']
		print('\t'.join(out_header), file = out_file)

		type_idx = header.index('variant_type')
		# inheritance_idx = header.index('inheritance')

		for line in snv_file:
			dnm = line.rstrip().split('\t')
			hifi_hap, ont_hap = [next(file).rstrip().split('\t') for file in haplotype_files]

			assert len({dnm[2], hifi_hap[5], ont_hap[5]}) == 1, "files are not in sorted order"

			if not final_validation(hifi_hap, ont_hap):
				continue

			final_variant_type = resolve_variant_type(hifi_hap[-1], ont_hap[-1], dnm[type_idx])
			types.append(final_variant_type)

			# final_inheritance = resolve_inheritance(hifi_hap[-2], ont_hap[-2], dnm[inheritance_idx])
			final_inheritance = resolve_inheritance(hifi_hap[-2], ont_hap[-2])
			parents.append(final_inheritance)

			out_info = [*dnm[:5], final_variant_type, final_inheritance]
			print('\t'.join(out_info), file = out_file)

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-s", "--snvs", required=True, help="all input snvs")
	# ap.add_argument("-hifi", "--hifi_haps", required=True, help="hifi haplotype assignment")
	# ap.add_argument("-ont", "--ont_haps", required=True, help="ont haplotype assignment")
	# ap.add_argument("-o", "--output", required=True, help="path to output file")
	# args = ap.parse_args()
	#
	# main(args.snvs, args.hifi_haps, args.ont_haps, args.output)
	main(snakemake.input.annotated_snvs, snakemake.input.hifi_haps, snakemake.input.ont_haps, snakemake.output.final_haps)
