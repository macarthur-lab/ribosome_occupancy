from subprocess import Popen, PIPE, STDOUT
import time
from collections import defaultdict
import os

"""
Input: coverageBed output file name OR list of lines in coverage files, extracted by "get_profile_for_region" 
"""
def parse_pileup(pileup, d=True):	
	if type(pileup) == str:
		# read pileup file 
		with open(pileup, 'r') as p:
			lines = [line.strip().split('\t') for line in p]
	else:
		lines = [line.split('\t') for line in pileup if not line == '']

	gene = None
	if d:
		genes = []
		# create a dictionary for each region of coverage, mapping pos to cvg
		# aggregate dictionaries in a list
		for line in lines[:-1]:
			# print line
			pos = int(line[3])
			cvg = int(line[4])
			if pos == 1:
				if gene:
					genes.append(gene)
				gene = {}
			gene[pos] = cvg
		return genes
	else:
		# return read counts for each interval in pileup
		return [int(line[3]) for line in lines]

"""
Utility for running quicking coverage jobs.
Accepts either a list of bed tuples (e.g. ('X', 123435, 643155)), or the name of a bedfile
"""

def get_profile_for_region(regions, bam, d=True, output=None):
	if output and os.path.isfile(output):
		print '%s already exists. Not computing coverage.' % output
		return

	sample = bam.split('/')[-1].split('_')[0]
	# If given python list of bed tuples
	if type(regions) == list:
		regions = ['\t'.join([str(x) for x in region]) for region in regions]
		regions = '\n'.join(regions)

		# Call coverageBed via subprocess module, piping in bedfile-like-region-object
		args = ['coverageBed', '-d', '-split', '-abam', bam, '-b', 'stdin']
		if not d:
			del args[1]
		proc = Popen(args, stdin=PIPE, stdout=PIPE, stderr=STDOUT)
		pileup = proc.communicate(regions)[0].split('\n')

	# If given a bed file
	elif type(regions) == str:
		args = ['coverageBed', '-d', '-split', '-abam', bam, '-b', regions]
		if not d:
			del args[1]
		proc = Popen(args, stdout=PIPE, stderr=STDOUT)
		print 'Getting coverage for %s over %s ...' % (bam, regions)
		start = time.time()
		pileup = proc.communicate()[0].split('\n')

	if output:
		with open(output, 'a') as o:
			for line in pileup:
				o.write('%s\n' % line)
		print 'Took %s minutes.' % ((time.time()-start)/60.0)
	else:
		print 'Took %s minutes.' % ((time.time()-start)/60.0)
		return pileup


def get_list_of_all_bams(bam_names='/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools/bams.txt'):
	with open(bam_names, 'r') as f:
		bams = [line.strip() for line in f]
	return bams


def listdir_fullpath(direc):
	import os
	direcf = [os.path.join(direc, f) for f in os.listdir(direc)]
	for f in direcf:
		try:
			yield os.readlink(f)
		except OSError:
			yield f


def get_depth_dict(info='/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/depth_distributions/depth_by_sample.txt'):
	with open(info, 'r') as h:
		depth_dict = dict([(line.split()[0], int(line.split()[1])) for line in h])
	return depth_dict


def bed_to_tuple(bed):
	with open(bed, 'r') as f:
		beds = [line.strip().split() for line in f]
	beds = [(x[0], int(x[1]), int(x[2])) for x in beds]
	return beds


def parse_commandline_argument(arg):
	# if user inputs filename
	if os.path.isfile(arg):
		return [arg]
	# if user inputs directory
	elif os.path.isdir(arg):
		return [f for f in listdir_fullpath(arg)]
	elif len(arg.split(',')) > 1:
		return arg.split(',')

def parse_single_column_file(fi, typ='string'):
	with open(fi, 'r') as f:
		if not typ == 'string':
			if typ == 'int':
				lines = [int(line.strip()) for line in f]
			elif typ == 'float':
				lines = [float(line.strip()) for line in f]
		else:
			lines = [line.strip() for line in f]
	return lines


if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--bed', dest="bed", help='Bedfile containing regions to compute coverage over', required=False)
	parser.add_argument('--counts', dest="counts", help='Output file containing table of read counts', required=False)
	parser.add_argument('--features', dest="features", help='Input file containing names of features in bedfile', required=False)
	parser.add_argument('--bam', dest="bam", help='Input file containing list of bams. If want all bams, enter "all"', required=True)
	args = parser.parse_args()

	if args.bam == 'all':
		bams = get_list_of_all_bams()
	elif args.bam:
		bams = get_list_of_all_bams(args.bam)

	if args.features:
		features = get_list_of_all_bams(args.features)

	samples = [x.split('/')[-1].split('_')[0] for x in bams]	
	bed = bed_to_tuple(args.bed)
	# print samples
	results = []
	counts = parallelize_coverageBed(bed, bams)
	# print counts
	write_read_count_table(args.counts, counts, samples, features)












