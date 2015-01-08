data <- read.delim('appris-principal-transcripts_pos_only_no_readthrough.table', header=TRUE, sep="\t", check.names=FALSE)
# library(edgeR)
# isexpr <- cpm(data)>50
isexpr <- data>100
samples <- colnames(data)
genes <- rownames(data)
for (i in 1:length(samples)) {
sample <- samples[i]
expr <- genes[isexpr[,sample]]
line <- c(sample, expr)
# print(line)
# print(length(line))	
write(line, file="../metagene/appris_principal/appris_principal_expressed_by_sample.pos_only.no_readthrough.txt", append=TRUE, sep="\t", ncolumns=length(line))
}
