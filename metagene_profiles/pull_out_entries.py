import argparse
from collections import Counter

def pull_end_exons(lines, strand=None, strand_dict=None, positive_strand_only=False, ignore_single_exon_genes=False):
	firsts = []
	lasts = []
	new_strands = []

	following = [lines[i+1] for i in range(len(lines)) if lines[i].startswith('#')]
	preceding = [lines[i-1] for i in range(len(lines)) if lines[i].startswith('#')]
	last_line = preceding[0]
	del preceding[0]
	preceding.append(last_line)

	idx_following = [i+1 for i in range(len(lines)) if lines[i].startswith('#')]
	idx_preceding = [i-1 for i in range(len(lines)) if lines[i].startswith('#')]
	del idx_preceding[0]
	idx_preceding.append(len(lines)-1)
	
	if strand_dict:
		transcripts = [line.lstrip('#').strip() for line in lines if line.startswith('#')]
		strands = [strands_dict[t] for t in transcripts]
	elif strand:
		strands = [strand]*len(following)
	
	idx_firsts = []
	idx_lasts = []

	for i in range(len(strands)):
		strand = strands[i]
		if ignore_single_exon_genes and following[i] == preceding[i]:
			continue
		if strand == '+':
			firsts.append(following[i])
			lasts.append(preceding[i])
			idx_firsts.append(idx_following[i])
			idx_lasts.append(idx_preceding[i])
			new_strands.append(strand)
		elif positive_strand_only:
				continue
		else:
			firsts.append(preceding[i])
			lasts.append(following[i])
			idx_firsts.append(idx_following[i])
			idx_lasts.append(idx_preceding[i])
			new_strands.append(strand)

	return firsts, lasts, idx_firsts, idx_lasts, new_strands


def pull_middle_exons(first_idx, last_idx, lines, strands):
	idx = []
	strand_mapping = []
	for x in range(len(first_idx)):
		strand = strands[x]
		strand_mapping += [strand]*len(range(first_idx[x]+1, last_idx[x]))
		idx += range(first_idx[x]+1, last_idx[x])
	pulled = [lines[i] for i in idx]
	return pulled, strand_mapping


def add_flanks(lines, strands, flank_size, flank_direction):
	newlines = []
	for i in range(len(lines)):
		line = lines[i].strip().split()
		strand = strands[i]
		if flank_direction == 'up':
			if strand == '+':
				line = [line[0], int(line[1])-50, line[2]]
			else:
				line = [line[0], line[1], int(line[2])+50]
		if flank_direction == 'down':
			if strand == '+':
				line = [line[0], line[1], int(line[2])+50]
			else:
				line = [line[0], int(line[1])-50, line[2]]
		line = [str(l) for l in line]
		newline = '\t'.join(line)+'\n'
		newlines.append(newline)
	return newlines


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', help='Input file containing cleaned, sorted CDS list', dest='in_file', required=True)
	parser.add_argument('--istrand', help='Input file containing strand direction for CDS in input', dest='strands', required=True)
	parser.add_argument('--first', help='Output file containing first CDS for each transcript', dest='first', required=False)
	parser.add_argument('--last', help='Output file containing last CDS for each transcript', dest='last', required=False)
	parser.add_argument('--middle', help='Output file containing middle CDS for each transcript', dest='middle', required=False)
	parser.add_argument('--ostrand-first', help='Output file containing strand directions for bed intervals in args.first', dest='first_strands',required=False)
	parser.add_argument('--ostrand-middle', help='Output file containing strand directions for bed intervals in args.middle', dest='middle_strands',required=False)
	parser.add_argument('--ostrand-last', help='Output file containing strand directions for bed intervals in args.last', dest='last_strands',required=False)
	parser.add_argument('--positive-strand-only', help='Only pull out CDS from genes on positive strand', dest='positive_strand_only', required=False, type=int)
	parser.add_argument('--ignore-single-exon-genes', help='Don\'t pull out single exon genes', type=int, dest='ignore_single_exon_genes')
	parser.add_argument('-a', '--annotation', dest='annotation', required=False)
	args = parser.parse_args()
	
	with open(args.in_file, 'r') as f:
		lines = [line for line in f]

	strands_dict = None
	if not args.strands in ['+', '-']:
		with open(args.strands, 'r') as f:
			strands_dict = dict([(line.split()[0], line.split()[1].strip()) for line in f])

	# if args.positive_strand_only:
	# 	firsts, lasts, idx_following, idx_preceding, first_strands = pull_end_exons(lines, strand_dict=strands_dict, positive_strand_only=True)
	# else:
	# 	firsts, lasts, idx_following, idx_preceding, first_strands = pull_end_exons(lines, strand_dict=strands_dict)
	
	if strands_dict:
		firsts, lasts, idx_following, idx_preceding, first_strands = pull_end_exons(lines, strand_dict=strands_dict, positive_strand_only=args.positive_strand_only, ignore_single_exon_genes=args.ignore_single_exon_genes)
	else:
		firsts, lasts, idx_following, idx_preceding, first_strands = pull_end_exons(lines, strand=args.strands, positive_strand_only=args.positive_strand_only, ignore_single_exon_genes=args.ignore_single_exon_genes)
	
	middle, middle_strands = pull_middle_exons(idx_following, idx_preceding, lines, first_strands)

	if args.first:
		firsts = add_flanks(firsts, first_strands, 50, 'up')
		new_firsts = add_flanks(firsts, first_strands, 50, 'down')
	if args.last:
		lasts = add_flanks(lasts, first_strands, 50, 'up')
		new_lasts = add_flanks(lasts, first_strands, 50, 'down')
	if args.middle:
		middle = add_flanks(middle, middle_strands, 50, 'up')
		middle = add_flanks(middle, middle_strands, 50, 'down')

	if args.first:
		with open(args.first, 'a') as f:
			for cds in new_firsts:
				f.write(cds)

	if args.middle:
		with open(args.middle, 'a') as f:
			for cds in middle:
				f.write(cds)

	if args.last:
		with open(args.last, 'a') as f:
			for cds in new_lasts:
				f.write(cds)

	if args.first_strands:
		with open(args.first_strands, 'a') as f:
			for strand in first_strands:
				f.write('%s\n' % strand)

	if args.middle_strands:
		with open(args.middle_strands, 'a') as f:
			for strand in middle_strands:
				f.write('%s\n' % strand)

	if args.last_strands:
		with open(args.last_strands, 'a') as f:
			for strand in first_strands:
				f.write('%s\n' % strand)

# Command line:

# python pull_out_entries.py -i high_rpkm/high_rpkm.CDS.cleaned --istrand high_rpkm/high_rpkm.transcripts.strands \
# --first high_rpkm/bed_input/first_CDS.bed --middle high_rpkm/bed_input/middle_CDS.bed --last high_rpkm/bed_input/last_CDS.bed \
# --ostrand-end high_rpkm/strands/ends.strands --ostrand-middle high_rpkm/strands/middle.strands

# After re-cleaning:

# python pull_out_entries.py -i recleaned_high_rpkm/high_rpkm.CDS.recleaned --istrand high_rpkm/high_rpkm.transcripts.strands \
# --first recleaned_high_rpkm/first.recleaned.no_singles_pulled.bed --middle recleaned_high_rpkm/middle.recleaned.bed \
# --last recleaned_high_rpkm/last.recleaned.no_singles_pulled.bed \
# --ostrand-first recleaned_high_rpkm/strands/first.with_singles.strands \
# --ostrand-middle recleaned_high_rpkm/strands/middle.strands \
# --ostrand-last recleaned_high_rpkm/strands/last.no_singles.strands






