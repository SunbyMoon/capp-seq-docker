library("plyr")
library("magrittr")
library("dplyr")

args <- commandArgs(TRUE)
sample = args[1] #sampleID
path = args[2] #/data/storage/capp-seq/docker_out/cnv
#CCD333-T-DNA_realigned_sorted

cnvkit_filter <- function(sample) {
sample <- read.csv2(file = paste0(path, "/", sample, "-T-DNA_realigned_sorted.cns"), header = T, stringsAsFactors = F, sep = "\t")
sample$copies <- 2*(2^(as.numeric(sample$log2)))
sample_subset <- data.frame(filter(sample, gene != "-") %>% filter(copies <= 1.5 | copies >= 2.5)) %>% arrange(copies)
return(sample_subset)
}

outfile <- cnvkit_filter(sample)
write.table(outfile, file=paste0(path, "/", sample, ".cnvkit.out"), sep = "\t", row.names = F, quote = F)

