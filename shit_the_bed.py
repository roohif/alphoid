#!/usr/bin/python3

# Convert the CSV file created by BLAST into a 'bucketed' BED file for Ensembl
from colour import Color
import math

# Mapping T2T to hg19/GRCh38
chromosome_map = {
	'gb|CP068255.2|': 'chrX', 'gb|CP068256.2|': 'chr22', 'gb|CP068257.2|': 'chr21',
	'gb|CP068258.2|': 'chr20', 'gb|CP068259.2|': 'chr19', 'gb|CP068260.2|': 'chr18',
	'gb|CP068261.2|': 'chr17', 'gb|CP068262.2|': 'chr16', 'gb|CP068263.2|': 'chr15',
	'gb|CP068264.2|': 'chr14', 'gb|CP068265.2|': 'chr13', 'gb|CP068266.2|': 'chr12',
	'gb|CP068267.2|': 'chr11', 'gb|CP068268.2|': 'chr10', 'gb|CP068269.2|': 'chr9',
	'gb|CP068270.2|': 'chr8', 'gb|CP068271.2|': 'chr7', 'gb|CP068272.2|': 'chr6',
	'gb|CP068273.2|': 'chr5', 'gb|CP068274.2|': 'chr4', 'gb|CP068275.2|': 'chr3',
	'gb|CP068276.2|': 'chr2', 'gb|CP068277.2|': 'chr1', 'gb|CP086569.2|': 'chrY' }

INFILE = 'vs_T2T.csv'
OUTFILE = 'alphoid.bed'

# Read the CSV file, and create a dictionary of buckets
bed_map = {}
with open(INFILE) as in_file:
	for line in in_file:
		(qseqid, qstart, qend, sseqid, sstart, send, nident, pident, length, qlen, evalue) = line.split(",")
		bucket = (int(sstart) + int(send)) // 2000000
		chr = chromosome_map[sseqid]

		# Add it to the map
		if chr not in bed_map:
			bed_map[chr] = { bucket: 0 }

		if bucket not in bed_map[chr]:
			bed_map[chr][bucket] = 0

		# Increment
		bed_map[chr][bucket] = bed_map[chr][bucket] + 1

# Create the colour gradient
green = Color("green")
gradient = list(green.range_to(Color("red"), 1001))

# Write the file out in BED format
out_file = open(OUTFILE, "w")

for chr in bed_map:
	for bucket in bed_map[chr]:
		count = bed_map[chr][bucket]
		bucketStart = int(bucket * 1e6) + 1
		bucketEnd = int((bucket + 1) * 1e6)

		nameStr = "\"" + str(bucket) + " MM bucket has " + str(count) + " repeat(s)\""

		# Work out the colour on the gradient
		gradient_index = int(math.log10(min(count, 1000)) * 1000 / 3)
		c = gradient[gradient_index]
		Rgb = str(int(min(c.red * 256, 255))) + \
			"," + str(int(min(c.green * 256, 255))) + \
			"," + str(int(min(c.blue * 256, 255)))

		line = "{chrom}\t{chromStart}\t{chromEnd}\t{name}\t{score}\t+\t{thickStart}\t{thickEnd}\t{itemRgb}\n"
		out_file.write(line.format(chrom = chr, chromStart = bucketStart, chromEnd = bucketEnd,
			name = nameStr, score = min(count, 1000), thickStart = bucketStart, thickEnd = bucketEnd,
			itemRgb = Rgb))

out_file.close()
