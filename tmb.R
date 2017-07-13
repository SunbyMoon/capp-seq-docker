PANCeqTMB <- function(nmuts) {
  bed <- read.csv2('PANCeq_C200X50bp_6.bed', sep='\t', header = F)
  panelSize <- sum(bed$V3 - bed$V2) / 1000000
  tmb <- nmuts / panelSize
}
