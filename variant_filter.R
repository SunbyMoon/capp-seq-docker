library(VariantAnnotation)
library(GenomicRanges)

args <- commandArgs(TRUE)
infile = args[1] # Annovar VCF
outfile = args[2] # Output VCF
outfile2 = args[3] #TMB outfile
scriptdir = args[4] # Repo filepath

setwd(scriptdir)
#setwd("/home")

vcfFile <- infile
vcf <- readVcf(vcfFile, 'hg19')
vcf <- vcf[fixed(vcf)$FILTER=='PASS']

# remove germline
vcf <- vcf[info(vcf)$STATUS!='Germline']

# only Somatic
vcf <- vcf[info(vcf)$STATUS=='StrongSomatic' | info(vcf)$STATUS=='LikelySomatic' | info(vcf)$STATUS=='Cluster75bp' | info(vcf)$STATUS=='Cluster0bp']

# remove strand bias
vcf <- vcf[geno(vcf)$SBF[,1] > 0.06]

# filter low freq with low depth
# based on http://bcb.io/2016/04/04/vardict-filtering/
lq <- sapply(vcf, function(x) {
  if (geno(x)$AF[,1] * geno(x)$DP[,1] < 6) {
    if ((geno(x)$MQ[,1] < 55 & geno(x)$NM[,1] > 1.0) | (geno(x)$MQ[,1] < 60 & geno(x)$NM[,1] > 2.0)) {
      res <- 'LowAlleleDepth'
    }
    if (geno(x)$DP[,1] < 10) {
      res <- 'LowAlleleDepth'
    }
    if (geno(x)$QUAL[,1] < 45) {
      res <- 'LowAlleleDepth'
    }
  } else if (geno(x)$AF[,1] < 0.2 & geno(x)$QUAL[,1] < 55 & info(x)$SSF > 0.06) {
    res <- 'LowFreqQuality'
  } else {
    res <- 'PASS'
  }
  res
})

vcf <- vcf[lq=='PASS']

# remove blacklisted regions
# bed to granges function from https://davetang.org/muse/2015/02/04/bed-granges/
bedToGranges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}
dac <- bedToGranges('hg19.dac.bed')
duke <- bedToGranges('hg19.duke.bed')

keep <- which(countOverlaps(granges(vcf), dac)==0 & countOverlaps(granges(vcf), duke)==0)
vcf <- vcf[keep]

PANCeqTMB <- function(nmuts) {
  bed <- read.csv2('PANCeq_C200X50bp_6.bed', sep='\t', header = F)
  panelSize <- sum(bed$V3 - bed$V2) / 1000000
  tmb <- nmuts / panelSize
}

# conpute tmb
#tmb <- round(PANCeqTMB(dim(vcf)[1]),1)

# annotation based filtering
# only exonic regions
vcf <- vcf[unlist(info(vcf)$Func.refGene)=='exonic' | unlist(info(vcf)$Func.refGene)=='splicing']

# only nonsynonoumus & frameshifts & stop_gains
vcf <- vcf[unlist(info(vcf)$ExonicFunc.refGene=='frameshift_deletion' |
                    info(vcf)$ExonicFunc.refGene=='frameshift_insertion' |  
                    info(vcf)$ExonicFunc.refGene=='frameshift_substitution' |
                    info(vcf)$ExonicFunc.refGene=='nonframeshift_deletion' |
                    info(vcf)$ExonicFunc.refGene=='nonframeshift_insertion' |  
                    info(vcf)$ExonicFunc.refGene=='nonframeshift_substitution' |
                    info(vcf)$ExonicFunc.refGene=='nonsynonymous_SNV' |
                    info(vcf)$ExonicFunc.refGene=='stopgain' |
                    info(vcf)$ExonicFunc.refGene=='stoploss')]

tmb <- round(PANCeqTMB(dim(vcf)[1]),1)

writeVcf(vcf, outfile)
write.csv(tmb, outfile2, row.names = F, quote = F)
