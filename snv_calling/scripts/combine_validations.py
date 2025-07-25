import argparse

def get_platform(file_name):
	"""
	Infer platform name from file name.
	:param file_name: Path to platform validation file.
	:return platform: 'hifi' 'ont' or 'illumina'
	"""
	platform = file_name.split('_')[-1].split('.')[0]
	return platform

def add_platform(platform, var_file):
	"""
	Add platform to the start of every line in corresponding validation file.
	:param platform: 'hifi' 'ont' or 'illumina'
	:param var_fil: Opened validation file.
	:return: Lines from file wtiht platform at the start.
	"""
	next(var_file)
	for line in var_file:
		yield platform, line.rstrip().split('\t')

def cross_platform_validated(pf_data):
	"""
	Determines final validation based on each of the three sequencing platforms.
	:param pf_data: Dictionary of platform validations mapped to platform names.
	:returns: Final validation (eg 'true_de_novo', 'inherited').
	"""
	hifi_val = pf_data['hifi']
	ont_val = pf_data['ont']
	ilm_val = pf_data['illumina']

	if hifi_val == 'true_de_novo' and ont_val == 'true_de_novo' and ilm_val == 'true_de_novo':
		return 'true_de_novo'

	if 'inherited' in [hifi_val, ont_val, ilm_val]:
		return 'inherited'

	if hifi_val  == 'true_de_novo':
		if ont_val == 'true_de_novo':
			if ilm_val == 'false_positive':
				return 'true_de_novo_ilm_false_positive'
			if 'missing' in ilm_val:
				return 'true_de_novo_no_ilm_data'

		if ont_val == 'false_positive':
			if ilm_val == 'true_de_novo':
				return 'true_de_novo_no_ont_support'
			if ilm_val in ['false_positive', 'missing_parent_data', 'missing_child_data']:
				return 'false_positive'

		if 'missing' in ont_val:
			if ilm_val == 'true_de_novo':
				return 'true_de_novo_no_ont_data'
			if ilm_val == 'false_positive':
				return 'false_positive'
			if 'missing' in ilm_val :
				return 'cannot_validate'

	if hifi_val in ['false_positive', 'missing_parent_data', 'missing_child_data']:
		if ont_val == 'true_de_novo' and ilm_val == 'true_de_novo':
			return 'true_de_novo_no_hifi_support'

	if hifi_val == 'cell_line_only':
		if ont_val == 'true_de_novo' and ilm_val == 'true_de_novo':
			return 'true_de_novo'
		elif ont_val == 'true_de_novo' and ilm_val == 'false_positive':
			return 'cell_line_only'
		elif ont_val == 'false_positive':
			return 'false_positive'
		elif 'missing' in ont_val and ilm_val == 'true_de_novo':
			return 'true_de_novo_no_ont_support'
		elif 'missing' in ilm_val and ont_val == 'true_de_novo':
			return 'cell_line_only'
		elif 'missing' in ont_val and 'missing' in ilm_val:
			return 'cannot_validate'

	if hifi_val == 'false_positive':
		return 'false_positive'

	if 'missing' in hifi_val:
		if 'missing_data' in [ont_val, ilm_val]:
			return 'cannot_validate'
		else:
			return 'false_positive'

	return f'missing_case_hifi_{hifi_val}_ont_{ont_val}_ilm_{ilm_val}'

def check_cell_line(hifi_validation, variant_data):
	if hifi_validation != 'cell_line_only':
		return False

	platform_ab = {x[0]: x[1][-2] for x in variant_data}
	if float(platform_ab['illumina']) > 0.05:
		return False

	return True

def get_output_data(hifi_validation, variant_data, platforms):
	"""
	Formats allele count and read depth information from each platform.
	:param hifi_validation: Validation based on hifi data (eg 'true_de_novo', 'inherited').
	:param variant_data: Dictionary of platform validations mapped to platform names.
	:param platforms: List of platforms with sequence data.
	:return: List of strings to be printed in output file.
	"""
	platforms = ['hifi', 'ont', 'illumina']

	info_dict = {'allele_count':{x[0]: x[1][-4] for x in variant_data},
				 'depth': {x[0]: x[1][-3] for x in variant_data},
				 'ab':  {x[0]: x[1][-2] for x in variant_data}}

	if hifi_validation == 'cell_line_only':
		hifi_data = [x[1] for x in variant_data if x[0] == 'hifi'][0]
		info_dict['allele_count']['hifi'] = hifi_data[6]
		info_dict['depth']['hifi'] = hifi_data[7]
		info_dict['ab']['hifi'] = hifi_data[8]

	var_out = [info_dict['allele_count'][x] for x in platforms] + [info_dict['depth'][x] for x in platforms] + [info_dict['ab'][x] for x in platforms]

	return var_out


def validate_variant(variant_data, platforms):
	"""
	Makes dictionary of platform validations and runs final validation script.
	:param variant_data: List of lists corresponding to each validation file.
	:param platforms: List of platforms with sequence data.
	:return: Final validation (eg 'true_de_novo', 'inherited').
	:return: List of strings to be printed in output file.
	"""
	platform_validation = {x[0]: x[1][-1] for x in variant_data}
	cell_line_artifact = check_cell_line(platform_validation['hifi'], variant_data)
	if cell_line_artifact:
		# print(variant_data[0][1][2])
		return 'cell_line_only', []

	final_validation = cross_platform_validated(platform_validation)

	var_output = get_output_data(platform_validation['hifi'], variant_data, platforms)

	return final_validation, var_output

def main(validation_files, output_file):
	variant_files = {get_platform(path): path for path in validation_files}
	files_with_platform = [add_platform(pf, open(path, 'r')) for pf, path in variant_files.items()]

	custom_order = {c: i for i, c in enumerate(['hifi', 'ont', 'illumina'])}
	platforms = sorted(variant_files.keys(), key=lambda c:custom_order.get(c, len(custom_order)))

	out_columns = ['chr', 'pos', 'id', 'ref', 'alt', 'caller'] + \
		[f'{x}_alt_allele_count' for x in  platforms] + \
		[f'{x}_dp' for x in  platforms] + \
		[f'{x}_ab' for x in  platforms] + \
		['final_validation']

	with open(output_file, 'w') as outfile:
		print('\t'.join(out_columns), file=outfile)

		for variants in zip(*files_with_platform):
			variant = variants[0][1][:6]
			for platform, line in variants:
				if line[:6] != variant:
					raise ValueError("your files are not in sorted order!")
			val, variant_out = validate_variant(variants, platforms)
			if 'true_de_novo' in val:
				full_variant_out = variant + variant_out + [val]
				print('\t'.join(full_variant_out), file=outfile)


	for file in files_with_platform:
		try:
			next(file)
		except StopIteration:
			pass
		else:
			raise ValueError("one of your files is longer than the other!")

if __name__ == "__main__":
	# ap = argparse.ArgumentParser()
	# ap.add_argument("-f", "--files", required=True, help="validation files", nargs='+')
	# ap.add_argument("-o", "--output", required=True, help="output")
	#
	# args = ap.parse_args()
	#
	# main(args.files, args.output)
	main(snakemake.input.platform_validations, snakemake.output.validated_snvs)