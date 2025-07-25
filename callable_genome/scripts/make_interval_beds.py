import argparse

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="Input file")

	args = ap.parse_args()

	interval = 10000000

	with open(args.input, 'r') as infile:
		for line in infile:
			chrom, start, end = line.rstrip().split('\t')
			start = 1
			end = int(end)
			pos = interval
			while pos < end:
				with open(f'bed/{chrom}_{start}_{pos}.bed', 'w') as outfile:
					print('\t'.join([chrom, str(start), str(pos)]), file = outfile)
				start = pos + 1
				pos += interval
			with open(f'bed/{chrom}_{start}_{end}.bed', 'w') as outfile:
				print('\t'.join([chrom, str(start), str(end)]), file = outfile)
