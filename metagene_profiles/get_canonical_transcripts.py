import gzip
import sys
import os
sys.path.insert(0, '/humgen/atgu1/fs03/birnbaum/ribosomal_occupancy/tools')
from map_genes_to_intervals import parse_gene_table

CANON = '/humgen/atgu1/fs03/konradk/exac/exac_browser/canonical_transcripts_v75.txt.gz'
def get_gene_transcript_map():
	d = {}
	with gzip.open(CANON, 'r') as f:
		for line in f:
			line = line.split()
			d[line[0]] = line[1]
	return d

if __name__ == '__main__':
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('--feature-table', dest='table', required=True)
	parser.add_argument('--output-dir', dest='output_dir', required=False)
	args = parser.parse_args()

	with open(args.table, 'r') as f:
		lines = [line.split('\t') for line in f]

	for line in lines:
		sample = line[0]
		with open(os.path.join(args.output_dir, '%s.transcript' % sample), 'a') as o:
			o.write('\n'.join(line[1:]))
