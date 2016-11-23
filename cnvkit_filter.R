library("plyr")
library("magrittr")
library("dplyr")

args <- commandArgs(TRUE)
sample = args[1] #sampleID
cnsfile = args[2] #cnvkitcns=/data/storage/capp-seq/patients/CCD347/cnv/CCD347-T1-DNA-TES1-REP1_realigned_sorted.cns
cnvout = args[3] #cnvkitout=/data/storage/capp-seq/patients/CCD347/cnv/CCD347.cnvkit.out

#path = args[2] #/data/storage/capp-seq/patients/CCD333/cnv
#CCD333-T-DNA_realigned_sorted

cnvkit_filter <- function(sample) {
#sample <- read.csv2(file = paste0(path, "/", sample, "-T-DNA_realigned_sorted.cns"), header = T, stringsAsFactors = F, sep = "\t")
sample <- read.csv2(file = cnsfile, header = T, stringsAsFactors = F, sep = "\t")
sample$copies <- 2*(2^(as.numeric(sample$log2)))
sample_subset <- data.frame(filter(sample, gene != "-") %>% filter(copies <= 1.5 | copies >= 2.5)) %>% arrange(copies)
return(sample_subset)
}

outfile <- cnvkit_filter(sample)
#write.table(outfile, file=paste0(path, "/", sample, ".cnvkit.out"), sep = "\t", row.names = F, quote = F)
write.table(outfile, file=cnvout, sep = "\t", row.names = F, quote = F)

