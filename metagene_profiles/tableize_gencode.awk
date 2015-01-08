{
	gsub(/[\";]/, "", $12)
	split($12, transcript, ".")
	chrom=substr($1, 4)
	printf "%s\t%s\t%s\t%s\t%s\n", transcript[1],$3,chrom,$4,$5

}
