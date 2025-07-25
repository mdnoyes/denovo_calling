import argparse

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-i", "--input", required=True, help="input")
	ap.add_argument("-o", "--output", required=True, help="Output file")

	args = ap.parse_args()

	with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
		for line in infile:
			file_name, size = line.rstrip().split(' ')

			interval = file_name.split('.')[1]

			start, end = [int(x) for x in interval.split('_')[1:]]

			expected_size = end - start + 1

			actual_size = int(size) - 1

			if expected_size != actual_size:
				print(file_name, file = outfile)
