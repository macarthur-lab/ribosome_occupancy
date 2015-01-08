import argparse
from collections import defaultdict
import sys
sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools')
from utilities import listdir_fullpath

""" also deduplicates list of intervals """
def remove_squished_exons(exon_list, no_singles=False):
	ok = []
	# removes single-exon transcripts
	if no_singles:
		return None
	elif len(exon_list) == 1:
		return exon_list

	current = exon_list[0]
	ok.append(current)
	for after in exon_list[1:]:
		if current not in ok:
			left = int(current[1]) - int(before[2])
			right = int(after[1]) - int(current[2])
			if left > 50 and right > 50:
				ok.append(current)
		before = current
		current = after

	left = int(current[1]) - int(before[2])
	# last line
	if left > 50:
		ok.append(current)
	return ok

def squished(exons, flank):
	current = exons[0]
	before = (1, -1000, -900)
	for after in exons[1:]:
		left = int(current[1]) - int(before[2])
		right = int(after[1]) - int(current[2])
		if left < 50 or right < 50:
			return True
		before = current
		current = after

	left = int(current[1]) - int(before[2])
	if left < 50:
		return True

	return False

def clean(dirty_file):
	transcripts = defaultdict(list)
	with open(dirty_file, 'r') as f:
		for line in f:
			line = line.split()
			interval = line[2:]
			transcript = line[0]
			transcripts[transcript].append(interval)

	cleaned = defaultdict(list)
	for trans in transcripts:
		cdss = transcripts[trans]
		cdss = sorted(cdss, key=lambda interval:interval[1])
		# cdss = remove_squished_exons(cdss)
		if not squished(cdss, 50):
			cleaned[trans] = cdss
	
	return cleaned


""" Sort CDS within each transcript. Remove CDS within 50 bp of eachother. """

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='Input file containing unsorted, duplicated CDS list', dest='dirty', required=False)
	parser.add_argument('-o', help='Output file containing CDS list with crunched CDS removed', dest='cleaned', required=False)
	parser.add_argument('--idir', help='Input directory containing CDS files', dest='idir', required=False)
	parser.add_argument('--odir', help='Output directory containing cleaned CDS files', dest='odir', required=False)
	args = parser.parse_args()

	if args.idir:
		for f in listdir_fullpath(args.idir):
			sample = f.split('/')[-1]
			transcripts = clean(f)
			print sample
			if args.odir:
				with open(args.odir+sample+'.cleaned', 'a') as o:
					for trans in transcripts:
						o.write('#%s\n' % trans)
						for interval in transcripts[trans]:
							o.write('%s\t%s\t%s\n' % (interval[0], interval[1], interval[2]))



# COMMAND LINE:
# python clean_CDS_list.py --idir 50cpm/CDS/ --odir 50cpm/cleaned/
# python clean_CDS_list.py -i high_rpkm/high_rpkm.CDS -o high_rpkm/high_rpkm.CDS.cleaned





