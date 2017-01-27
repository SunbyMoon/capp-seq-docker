library(VariantAnnotation)
#library(myvariant)

args <- commandArgs(TRUE)
infile = args[1] #/data/storage/capp-seq/patients/CCD010/vcf/CCD010_vardict_somatic.annovar_filt.hg19_multianno.vcf
outfile = args[2] #/data/storage/capp-seq/patients/CCD010/vcf/CCD010_somatic_postfilter.vcf

#vcfFile <- '~/AWS/storage/capp-seq/patients/CCD347/vcf/CCD347_vardict_somatic.annovar.hg19_multianno_refGene.vcf'
#param <- ScanVcfParam(fixed="ALT")
vcfFile <- infile
vcf <- readVcf(vcfFile, 'hg19')

# remove germline
vcf <- vcf[info(vcf)$STATUS!='Germline']

vcf <- vcf[fixed(vcf)$FILTER=='PASS']

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

#writeVcf(vcf, '~/AWS/storage/capp-seq/patients/CCD347/vcf/CCD347_vardict_somatic_postfilter_refseq.vcf')
writeVcf(vcf, outfile)
