import re
from collections import defaultdict
import time
import sys
import os
sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools')
from utilities import get_depth_dict, listdir_fullpath, parse_commandline_argument, parse_single_column_file
from parsing_pileups import parse_pileup
from numpy import linspace, average, median

def collapse(collapsed, profile, strand, flank_size, num_samples, tiled=True, window_size=50):

	gene_length = len(profile.keys())
	start = flank_size+15
	end = gene_length-flank_size-5

	if (end - start) < 80:
		return

	mean = average([profile[pos] for pos in profile if (pos > start and pos < end)])
	if mean == 0:
		return
		# mean += 0.01
	med = median([profile[pos] for pos in profile])

	if tiled:
		for pos in profile:
			snapshot_idx = list(linspace(flank_size+1, gene_length-flank_size, num_samples, dtype=int))		
			cvg = profile[pos]
			# cvg = cvg/depth
			cvg = cvg/mean

			# flip position index for intervals that have negative strand direction
			if strand == '-':
				pos = gene_length - pos + 1			
			# is pos in 5' flank?
			if pos <= flank_size:
				collapsed[pos].append(cvg)
			# is pos in 3' flank?
			elif pos > gene_length - flank_size:
				pos = pos - (gene_length - flank_size)
				pos = pos + flank_size + num_samples
				collapsed[pos].append(cvg)
			# is pos one of the sampled positions in the CDS?
			elif pos in snapshot_idx:
				sample_idx = snapshot_idx.index(pos)+1
				collapsed[flank_size+sample_idx].append(cvg)

	else:
		for pos in profile:
			cvg = profile[pos]
			cvg = cvg/mean

			# flip position index for intervals that have negative strand direction
			if strand == '-':
				pos = gene_length - pos + 1			
			# is pos in 5' flank?
			if pos <= flank_size:
				collapsed[pos].append(cvg)
			# is pos in 3' flank?
			elif pos > gene_length - flank_size:
				pos = pos - (gene_length - flank_size)
				pos = pos + flank_size + num_samples
				collapsed[pos].append(cvg)
			# is pos in first window
			elif pos <= flank_size + window_size:
				collapsed[pos].append(cvg)
			# is pos in second window
			elif pos > gene_length - (flank_size+window_size):
				pos_adj = flank_size + window_size + window_size - ((gene_length-flank_size)-pos)
				collapsed[pos_adj].append(cvg)


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--cvg_dir', dest='cvg_dir', help='Directory containing coverageBed output files.', required=True)
	parser.add_argument('--plot', dest='plot', help='Output file', required=False)
	parser.add_argument('--flank', dest='flank', type=int, required=True)
	parser.add_argument('--strand', dest='strand', required=False)
	parser.add_argument('--tile', dest='tile', required=True, type=int)
	args = parser.parse_args()

	d = get_depth_dict()
	# pileups = [p for p in listdir_fullpath(args.cvg_dir)]
	pileups = parse_commandline_argument(args.cvg_dir)
	pileups = sorted(pileups, key=lambda p: p.split('/')[-1].split('.')[0])

	if not args.strand in ['+', '-']:
		strand_files = parse_commandline_argument(args.strand)
		strand_files = sorted(strand_files, key=lambda f: f.split('/')[-1].split('.')[0])
		if len(strand_files) == 1:
			strand_files = strand_files*len(pileups)
		strand_lists =[parse_single_column_file(f) for f in strand_files]

	j = 0
	for f in pileups:
		collapsed_profile = defaultdict(list)
		print "Collapsing %s .." % f
		start = time.time()
		gene_profiles = parse_pileup(f)
		sample = f.split('/')[-1].split('.')[0]
		# depth = d[sample]
		if not args.strand in ['+', '-']:
			strands = strand_lists[j]
		j += 1
		
		for i in range(len(gene_profiles)):
			profile = gene_profiles[i]

			if not args.strand in ['+', '-']:
				strand = strands[i]
			else:
				strand = args.strand

			if strand == '-':
				continue
			collapse(collapsed_profile, profile, strand, 50, 100, tiled=args.tile)
		print "Took %s minutes." % ((time.time()-start)/60.0)

		if args.plot:
			with open(args.plot, 'a') as o:
				for pos in collapsed_profile:
					o.write('%s\t%s\n' % (pos, average(collapsed_profile[pos])))


# python collapsing.py --cvg_dir 50cpm/cvg/middle-test/ --flank 50 --strand 50cpm/middle/strands-test/ --plot collapsing_test.plot



