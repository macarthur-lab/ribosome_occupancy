from utilities import get_profile_for_region, listdir_fullpath, get_list_of_all_bams, parse_commandline_argument
from parallelize_coverageBed import Parallelizer
import os

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--bed', dest="bed", help='Bedfile/Bedfile directory containing regions/bedfiles to compute coverage over', required=True)
	parser.add_argument('-d', dest="d", help='Use `-d` parameter in coverageBed?', type=int, required=True)
	parser.add_argument('--output_dir', dest="output_dir", help='Path to the directory to dump pileup files', required=False)
	parser.add_argument('--parallelize', dest="parallelize", help='Parallelize over bams when running coverageBed?', type=int, required=True)	
	parser.add_argument('--bam', dest="bam", help='Input file containing list of bams. If you want all bams, enter "all"', required=True)
	args = parser.parse_args()

	if args.bam == 'all':
		bams = get_list_of_all_bams()
	else:
		bams = parse_commandline_argument(args.bam)
		bams = filter(lambda bam: not bam.endswith('bai'), bams)

	# Get list of bed filenames

	# if user inputs filename
	if os.path.isfile(args.bed):
		beds = [args.bed]*len(bams)
	# if user inputs directory
	elif os.path.isdir(args.bed):
		beds = [f for f in listdir_fullpath(args.bed)]
	

	# Sort bams and beds by a common key, their sample ID 
	bams = sorted(bams, key=lambda bam: bam.split('/')[-1].split('_')[0])
	beds = sorted(beds, key=lambda bed: bed.split('/')[-1].split('.')[0])

	samples = [x.split('/')[-1].split('_')[0] for x in bams]

	if args.parallelize:
		cvg_dir = '../lsf-dump/cvg'
		log_dir = '../lsf-dump/logs'
		group = '/birnbaum/runs/cvg'

		p = Parallelizer(beds, bams, cvg_dir, log_dir, group, args.d)

		# if args.d:
		# 	p = Parallelizer(beds, bams, cvg_dir, log_dir, group, True)
		# else:
		# 	p = Parallelizer(beds, bams, cvg_dir, log_dir, group, False)

		p.run()

	else:
		for i in range(len(beds)):
			bed = beds[i]
			bam = bams[i]
			print bed, bam
			sample = samples[i]
			get_profile_for_region(bed, bam, d=args.d, output=os.path.join(args.output_dir, '%s.coverage' % sample))

# examples:
# python run_coverageBed.py --bed ../metagene/50cpm/middle/ --bam all --parallelize 0 --output_dir ../metagene/50cpm/cvg/middle/ -d 1

# python ../tools/run_coverageBed.py --bed input_files/middle-CDS.bed --bam RNA_seq_bams/GTEX-N7MS-0426-SM-2YUN6.bam,RNA_seq_bams/GTEX-P44G-0526-SM-2XCD1.bam --output_dir RNA-control/middle/ --parallelize 0 -d 1



