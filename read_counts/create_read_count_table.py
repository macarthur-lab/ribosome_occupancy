import sys
sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools')
from utilities import get_profile_for_region, parse_pileup, get_list_of_all_bams

def write_read_count_table(out, results, samples, features):
	with open(out, 'a') as o:
		# Write header
		o.write('\t'.join(samples)+'\n')
		zipped = zip(*results)
		for i in range(len(zipped)):
			o.write(features[i]+'\t')
			counts = [str(x) for x in zipped[i]]
			o.write('\t'.join(counts)+'\n')


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--bed', dest="bed", help='Bedfile containing regions to compute coverage over', required=False)
	parser.add_argument('--output', dest="output", help='Output file containing table of read counts', required=False)
	parser.add_argument('--features', dest="features", help='Input file containing names of features in bedfile', required=False)
	parser.add_argument('--bam', dest="bam", help='Input file containing list of bams. If want all bams, enter "all"', required=True)
	args = parser.parse_args()

	if args.bam == 'all':
		bams = get_list_of_all_bams()
	else:
		bams = parse_commandline_argument(args.bam)
		bams = filter(lambda bam: not bam.endswith('bai'), bams)

	samples = [x.split('/')[-1].split('_')[0] for x in bams]

	with open(args.features, 'r') as f:
		features = [line.strip() for line in f]

	results = []
	for bam in bams:
		profile = get_profile_for_region(args.bed, bam, d=False)
		counts = parse_pileup(profile, d=False)
		results.append(counts)
	
	write_read_count_table(args.output, results, samples, features)




