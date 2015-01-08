args <- commandArgs(trailingOnly=TRUE)
data <- read.delim(args[1], header=TRUE, sep="\t", check.names=FALSE)
isexpr <- data>100
samples <- colnames(data)
genes <- rownames(data)
for (i in 1:length(samples)) {
sample <- samples[i]
expr <- genes[isexpr[,sample]]
line <- c(sample, expr)
write(line, file="../metagene/appris_principal/appris_principal_expressed_by_sample.pos_only.no_readthrough.txt", append=TRUE, sep="\t", ncolumns=length(line))
}
