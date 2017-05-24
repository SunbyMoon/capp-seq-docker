###############################################################################
## AREP PDF
###############################################################################
# Copyright (c) 2016 Tobias Meissner

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights 
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
# copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
# THE SOFTWARE.

#!/usr/bin/env Rscript

###############################################################################
# command line options
###############################################################################
packages <- function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x, character.only=TRUE, quietly=TRUE)){
    install.packages(pkgs=x, repos="http://cran.r-project.org", quiet=TRUE)
    require(x, character.only=TRUE, quietly=TRUE)
  }
}
packages(optparse)
option_list = list(
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="annovar annoted filtered somatic vardict vcf file", metavar="character"),
  make_option(c("-r", "--raw"), type="character", default=NULL, 
              help="annovar annoted vardict vcf file", metavar="character"),
  make_option(c("-C", "--cnv"), type="character", default=NULL, 
              help="cnfkit output", metavar="character"),
  make_option(c("-i", "--id"), type="character", default=NULL, 
              help="patient ID", metavar="character"),
  make_option(c("-sc", "--sample_config"), type="character", default=NULL, 
              help="sample config file in .yaml format", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="~", 
              help="output file path [default= %default]", metavar="character"),
  make_option(c("-d", "--scriptdir"), type="character", default="~",
              help="script file path [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$vcf)){
  print_help(opt_parser)
  stop("A annotation file needs to be provided", call.=FALSE)
}
if (is.null(opt$id)){
  print_help(opt_parser)
  stop("A patient ID needs to be provided", call.=FALSE)
}

##
packages(knitr)
packages(VariantAnnotation)
packages(plyr)
packages(rtracklayer)
packages(myvariant)

setwd(opt$scriptdir)

# cancer genes
cg <- read.csv2('cancer_genes.csv', sep='\t', stringsAsFactors = F)

#isoforms
li <- read.csv2('longest_resfseq.txt', sep = ' ', stringsAsFactors = F)

# rando string fucntion
# code from https://ryouready.wordpress.com/2008/12/18/generate-random-string-name/
MHmakeRandomString <- function(n=1, lenght=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

vcf <- readVcf(opt$vcf, 'hg19')
vcfraw <- readVcf(opt$raw, 'hg19')

tab <- data.frame(Variant=names(ranges(vcf)),
                  Gene=unlist(lapply(info(vcf)$Gene.refGene, function(x) paste(x, collapse = ' '))),
                  Amino.Acid=unlist(lapply(info(vcf)$AAChange.refGene, function(x) paste(x, collapse = ' '))),
                  Type=info(vcf)$TYPE,
                  Status=info(vcf)$STATUS,
                  AlleleFrequency=geno(vcf)$AF[,1],
                  RadialSVM=unlist(info(vcf)$RadialSVM_pred),
                  CADD=as.numeric(unlist(info(vcf)$CADD_phred)),
                  ExAC=as.numeric(as.vector(unlist(info(vcf)$ExAC_ALL))),
                  CLINSIG=unlist(lapply(as.vector(info(vcf)$CLINSIG), paste, collapse='|')),
                  COSMIC=unlist(lapply(info(vcf)$cosmic70, function(x) paste(x, collapse = ','))),
                  ICGC=unlist(info(vcf)$ICGC_Id),
                  stringsAsFactors = F
                  )

# select isoform based on longest transcript
tab$Amino.Acid.Uq <- NA
for(i in 1:length(tab$Amino.Acid)) {
  s <- unlist(strsplit(unlist(as.vector(unlist(tab$Amino.Acid[i]))), ' '))
  if(length(s)>1) {
    ts <- li$name[match(tab$Gene[i], li$name2)]
    if(length(grep(ts, s))==1) {
      tab$Amino.Acid.Uq[i] <- s[grep(ts, s)]  
    } else {
      tab$Amino.Acid.Uq[i] <- s[1]
    }
  } else {
    tab$Amino.Acid.Uq[i] <- s
  }
}

#tab <- tab[order(tab$Gene),]

# add drive / actionable status
inf <- sapply(tab$Gene, function(x) {
  if (x %in% cg$symbol) {
    if(cg[match(x, cg$symbol),]$intogenDriver==1 & cg[match(x, cg$symbol),]$actionable==1) {
      'driver,actionable'
    } else if(cg[match(x, cg$symbol),]$intogenDriver==1 & cg[match(x, cg$symbol),]$actionable==0) {
      'driver'
    } else if (cg[match(x, cg$symbol),]$intogenDriver==0 & cg[match(x, cg$symbol),]$actionable==1) {
      'actionable'
    } else {
      ''
    }
  } else {
    ''
  }
})

tab <- cbind(tab, Info=as.vector(inf))

# ADD SNPEFF impact prediction
hgvs <- rownames(tab)
hgvs <- gsub('/', '>', hgvs)
hgvs <- gsub('_', '', hgvs)
hgvs <- gsub(':', ':g.', hgvs)

ts <- unlist(lapply(strsplit(tab$Amino.Acid.Uq, ':'), '[[', 2))

x <- getVariants(hgvs, fields = c('snpeff.ann'))

z <- NULL
for (i in 1:length(x$snpeff.ann)) {
  if(is.null(x$snpeff.ann[[i]])) {
    z[i] <- NA
  } else {
    a <- grep(ts[i], x$snpeff.ann[[i]]$feature_id)
    z[i] <- x$snpeff.ann[[i]]$putative_impact[a]  
  }
}

tab$SNPEFF <- z

## annotate with gavin
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1141-7
################################
## GAVIN Applied Rule Guide r0.3
################################
## 
## Variant can be interpreted by following the columns from left to right.
## This classification scheme was implemented in Java, and used to
## benchmark GAVIN in the paper (Van der Velde et al., using r0.2), see:
## https://github.com/molgenis/molgenis and https://github.com/molgenis/gavin
## 
## Genome-wide rules are used if the gene-specific rules fail to classify.
## These rules are applied as follows:
## 1) If impact equals MODIFIER -> benign,
## 2) if MAF greater than 0.003456145 -> benign,
## 3) if CADD greater than 15 -> pathogenic,
## 4) if CADD less than 15 -> benign.
## 
## Explanation of the gene calibration categories:
## C1 = CADD scores highly significantly predictive for pathogenicity (pval < 0.01).
## C2 = CADD scores significantly predictive for pathogenicity (pval < 0.05).
## C3 = CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).
## C4 = CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).
## C5 = CADD scores less predictive for pathogenicity (population CADD > pathogenic CADD).
## I1 = HIGH impact unique for, thus predictive, for pathogenic variants.
## I2 = MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.
## I3 = LOW or MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.
## T1 = Too few exac variants after filtering with pathogenic 95th percentile MAF.
## T2 = Too few exac variants after filtering with impact distribution.
## N1 = Too few ClinVar variants for calibration at this time.
## N2 = Too few ExAC variants found for calibration.
## 
## For C1 and C2, CADD score thresholds are means of stratified benign and pathogenic variants.
## For C3, C4 and C5, CADD score thresholds are 95th sensitivity/specificity percentiles of stratified benign and pathogenic variants.
## 
gavinRules <- read.csv2('GAVIN_ruleguide_r0.3.tsv', skip=33, sep='\t', dec='.', stringsAsFactors = F, na.strings = 'n/a')

gavin <- function(gene, cadd, exac, effect) {
  pred <- NA
  cat <- NA
  g <- 0
  
  ## calibrations
  cal <- list(C1 = 'CADD scores highly significantly predictive for pathogenicity (pval < 0.01).',
              C2 = 'CADD scores significantly predictive for pathogenicity (pval < 0.05).',
              C3 = 'CADD scores may be predictive for pathogenicity (pval > 0.05 but with few samples).',
              C4 = 'CADD scores less predictive for pathogenicity (pval > 0.05 with enough samples).',
              C5 = 'CADD scores less predictive for pathogenicity (population CADD > pathogenic CADD).',
              I1 = 'HIGH impact unique for, thus predictive, for pathogenic variants.',
              I2 = 'MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.',
              I3 = 'LOW or MODERATE or HIGH impact unique, thus predictive, for pathogenic variants.',
              T1 = 'Too few exac variants after filtering with pathogenic 95th percentile MAF.',
              T2 = 'Too few exac variants after filtering with impact distribution.',
              N1 = 'Too few ClinVar variants for calibration at this time.',
              N2 = 'Too few ExAC variants found for calibration.'
  )
  
  if(gene %in% gavinRules$Gene) {
    r <- subset(gavinRules, Gene==gene)
    if (!is.na(cadd) & !is.na(r$PathogenicIfCADDScoreGreaterThan) & cadd >= r$PathogenicIfCADDScoreGreaterThan) {
      pred <- 'pathogenic'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      catC <- r$CalibrationCategory
      g <- 1
    } else if (!is.na(cadd) & !is.na(r$PathogenicIfCADDScoreGreaterThan) & cadd < r$BenignIfCADDScoreLessThan) {
      pred <- 'benign'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      catC <- r$CalibrationCategory
      g <- 1
    } else if (!is.na(exac) & !is.na(r$BenignIfMAFGreaterThan) & exac >= r$BenignIfMAFGreaterThan) {
      pred <- 'benign'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      catC <- r$CalibrationCategory
      g <- 1
    } else if (!is.na(effect) & !is.na(r$PathogenicIfImpactEquals) & effect==r$PathogenicIfImpactEquals) {
      pred <- 'pathogenic'
      cat <- as.vector(unlist(cal[r$CalibrationCategory]))
      catC <- r$CalibrationCategory
      g <- 1
    }
  }
  # apply global fallback
  if (g==0 & !is.na(effect) & effect=='MODIFIER') {
    pred <- 'benign'
    cat <- 'Genome wide rule'
    catC <- 'G'
  } else if (g==0 &  !is.na(exac) & exac >= 0.003456145) {
    pred <- 'benign'
    cat <- 'Genome wide rule'
    catC <- 'G'
  } else if (g==0 &  !is.na(cadd) & cadd >= 15) {
    pred <- 'pathogenic'
    cat <- 'Genome wide rule'
    catC <- 'G'
  } else if (g==0 & !is.na(cadd) & cadd < 15) {
    pred <- 'benign'
    cat <- 'Genome wide rule'
    catC <- 'G'
  } else if (g==0) {
    pred <- NA
    cat <- NA
    catC <- ''
  }
  return(list=c(pred,cat,catC))
}

g <- t(apply(tab, 1, function(x) {
  gavin(x['Gene'], x['CADD'], x['ExAC'], x['SNPEFF'])
}))
colnames(g) <- c('GARVIN_PRED', 'GARVIN_CAT', 'GARVING_CAT_C')

tab <- cbind(tab, g)

#pathogenecity
gav <- ifelse(tab$GARVIN_PRED=='pathogenic', 'PATHOGENIC', ifelse(tab$GARVIN_PRED=='benign', 'BENIGN', NA))
#gav[is.na(gav)] <- 'VUS'
clin <- ifelse(grepl('Pathogenic', tab$CLINSIG), 'PATHOGENIC', ifelse(grepl('Benign', tab$CLINSIG), 'BENIGN', NA))
path <- cbind(gav,as.vector(tab$GARVING_CAT_C),clin)

pathS <- apply(path, 1, function(x) {
  if(is.na(x[1]) & is.na(x[3])) {
    'VUS'
  } else if(is.na(x[1]) & !is.na(x[3])) {
    if(is.na(x[1]) & x[3]=='BENIGN') {
      'likely benign'
    } else if (is.na(x[1]) & x[3]=='PATHOGENIC') {
      'likely pathogenic'
    } 
  } else if (!is.na(x[1]) & is.na(x[3])) {
    if(x[1]=='PATHOGENIC' & is.na(x[3])) {
      ifelse(x[2]=='C1' | x[2]=='C2' | x[2]=='C3' | x[2]=='I1' | x[2]=='I2' | x[2]=='I3', 'pathogenic', 'likely pathogenic')
    } else if(x[1]=='BENIGN' & is.na(x[3])) {
       'likely benign'
    } 
  } else if (!is.na(x[1]) & !is.na(x[3])) {
    if(x[1]=='BENIGN' & x[3]=='BENIGN') {
      'benign'
    } else if(x[1]=='BENIGN' & x[3]=='PATHOGENIC') {
      'likely pathogenic'
    } else if(x[1]=='PATHOGENIC' & x[3]=='BENIGN') {
      'likely pathogenic'
    } else if(x[1]=='PATHOGENIC' & x[3]=='PATHOGENIC') {
      'pathogenic'
    } 
  }
  
  # if(x[1]=='BENIGN' & x[3]=='BENIGN') {
  #   'benign'
  # } else if(x[1]=='BENIGN' & is.na(x[3])) {
  #   'likely benign'
  # } else if(x[1]=='BENIGN' & x[3]=='PATHOGENIC') {
  #   'likely pathogenic'
  # } else if(x[1]=='PATHOGENIC' & x[3]=='BENIGN') {
  #   'likely pathogenic'
  # } else if(x[1]=='PATHOGENIC' & is.na(x[3])) {
  #   ifelse(x[2]=='C1' | x[2]=='C2' | x[2]=='C3' | x[2]=='I1' | x[2]=='I2' | x[2]=='I3', 'pathogenic', 'likely pathogenic')
  # } else if(x[1]=='PATHOGENIC' & x[3]=='PATHOGENIC') {
  #   'pathogenic'
  # } else if(is.na(x[1]) & x[3]=='BENIGN') {
  #   'likely benign'
  # } else if(is.na(x[1]) & is.na(x[3])) {
  #   'VUS'
  # } else if(is.na(x[1]) & x[3]=='PATHOGENIC') {
  #   'likely pathogenic'
  # } 
  
})
pathS <- factor(pathS, levels=c('benign', 'likely benign', 'VUS', 'likely pathogenic', 'pathogenic'))

# Variant Rank
pred <- ifelse(tab$RadialSVM=='D' & as.numeric(as.vector(tab$CADD))>=20, 2, ifelse(tab$RadialSVM=='T' & as.numeric(as.vector(tab$CADD))<20, 0, ifelse(tab$RadialSVM=='D' | as.numeric(as.vector(tab$CADD))>=20, 1, NA)))
info <- ifelse(tab$Info!='', ifelse(tab$Info=='driver,actionable', 2, 1), 0)
exac <- ifelse(tab$ExAC<=0.01, 1, ifelse(tab$ExAC>0.01, 0, 0))
known <- ifelse(tab$ICGC!='.' | tab$COSMIC!='.', 1, 0)

rank <- apply(cbind(pred,info,exac,known), 1, sum, na.rm=T)

tab <- cbind(tab, rank, pathS)
tab <- arrange(tab, pathS, rank, Gene, decreasing=T)

# read in CNV data
if(any(names(opt)=='cnv')) {
  cnv <- read.csv2(opt$cnv, sep='\t', stringsAsFactors = F)
  
  parseGene <- function(i) {
    x <- strsplit(i, ',')
    xx <- lapply(x, strsplit, ';')
    xxx <- cbind(sapply(xx[[1]], '[[', 2), sapply(xx[[1]], '[[', 1))
    xxx[,1] <- gsub('gene_name=', '', xxx[,1])
    xxx[,2] <- gsub('exon=', '', xxx[,2])
    colnames(xxx) <- c('Gene', "Exon")
    xxx <- as.data.frame(xxx)
    if(length(unique(xxx$Gene))==1) {
      data.frame(Gene=xxx[1,1],Exon=paste(xxx$Exon, collapse=','))
    } else {
      aggregate(Exon ~ Gene, xxx, paste, collapse=',')
    }
  }
  
  cnvls <- lapply(cnv$gene, parseGene)
  
  for(i in 1:length(cnvls)) {
    cnvls[[i]] <- cbind(cnvls[[i]], CopyNumber=round(as.numeric(cnv$copies),2)[i])
  }
  
  cnvtab <- data.frame(do.call(rbind,cnvls), stringsAsFactors = F)
  cnvtab$Gene <- as.vector(cnvtab$Gene)
  cnvtab <- cnvtab[order(cnvtab$Gene),]
}

# ---------------------------
# pharm
# data from https://github.com/AshleyLab/stmp
pharm <- read.csv2('clinical_ann_metadata-snvs_canc.csv', sep=',', stringsAsFactors = F)

vcf.pharm <- vcfraw[unlist(info(vcfraw)$snp142) %in% pharm$Location]
pharmTab <- pharm[unlist(match(info(vcf.pharm)$snp142, pharm$Location)), c(2,4,5,13,7)]

a1 <- vector()
a2 <- vector()
x1 <- granges(vcf.pharm)
x2 <- geno(vcf.pharm)$GT[,2]
x2 <- strsplit(x2, '/')
for(i in 1:length(x2)) {
  if(x2[[i]][1]==0 & x2[[i]][2]==0) {
      a1[i] <- as.vector(x1[i]$REF)
      a2[i] <- as.vector(x1[i]$REF)
    } else if(x2[[i]][1]==0 & x2[[i]][2]==1) {
      a1[i] <- as.vector(x1[i]$REF)
      a2[i] <- as.vector(unlist(x1[i]$ALT))
      } else if(x2[[i]][1]==1 & x2[[i]][2]==0) {
        a1[i] <- as.vector(unlist(x1[i]$ALT))
        a2[i] <- as.vector(x1[i]$REF)
      } else if(x2[[i]][1]==1 & x2[[i]][2]==1) {
        a1[i] <- as.vector(unlist(x1[i]$ALT))
        a2[i] <- as.vector(unlist(x1[i]$ALT))
        }
}
pharmTab$GenoType <- paste0(a1,a2)
pharmTab$Gene <- unlist(info(vcf.pharm)$Gene.refGene)[match(pharmTab$Location, unlist(info(vcf.pharm)$snp142))]
colnames(pharmTab) <- c('dbSNP', 'Evidence', 'Effect',
                        'Disease', 'Variants', 'Genotype', 'Gene')
pharmTab <- pharmTab[,c(7,1,3,2,6,5,4)]
pharmTab$Effect <- gsub('/', ' \\\\newline ',pharmTab$Effect)
pharmTab$Effect <- gsub(';', ' \\\\newline ',pharmTab$Effect)
pharmTab$Disease <- gsub(';', ' \\\\newline ',pharmTab$Disease)

# exclude evidence levels 3 & 4
pharmTab <- subset(pharmTab, Evidence!='3' & Evidence!='4')
    
# reduce annotation to genotype
redAnn <- vector()
for (i in 1:length(pharmTab$Genotype)) {
  gen <- pharmTab$Genotype[i]
  
  x <- unlist(strsplit(pharmTab$Variants[i], ';'))
  xx <- unlist(lapply(strsplit(x, ':'), '[[', 2))
  names(xx) <- unlist(lapply(strsplit(x, ':'), '[[', 1))
  names(xx) <- gsub(' ', '', names(xx))
  
  if(gen %in% names(xx)) {
    redAnn[i] <- xx[gen]
  } else if(reverse(gen) %in% names(xx)) {
    redAnn[i] <- xx[reverse(gen)]
  } else {
    redAnn[i] <- NA
  }
}
pharmTab$Variants <- redAnn

#---------------------------------------------------------
# pahtogenic germline variants &
# germline varaints associated with drug response
# based on clinvar info

germline <- vcfraw[info(vcfraw)$STATUS=='Germline']
keep <- ! (grepl('d7', fixed(germline)$FILTER) | grepl('MSI', fixed(germline)$FILTER))
germline <- germline[keep]
germline <- germline[unlist(info(germline)$ExonicFunc.refGene)!='synonymous_SNV']

# remove strand bias
germline <- germline[geno(germline)$SBF[,1] > 0.06]

# annotation based filtering
# only exonic regions
germline <- germline[unlist(info(germline)$Func.refGene)=='exonic' | unlist(info(germline)$Func.refGene)=='splicing']
# 
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4318970/
# #≥0.5% in the 1000 Genomes Project, Exome Variant Server, Geuvadis European Exome Variants Server, or the Centre Nacional d'Anàlisi Genòmica in-house database were filtered out. 
# #Variants present in >10 of the 43 individuals in our data set were discarded because they most likely corresponded to polymorphisms. 
# #Also, only variants predicted to have a strong effect on gene function (frameshift, spsplice-site canonical, nonsense, and missense) were chosen. 
# #Regarding missense variants, we used six bioinformatics tools to select for a deleterious amino acid change, namely, PhyloP (score >0.85), SIFT (score <0.05), 
# #PolyPhen (score >0.85), GERP (score >2), Mutation Taster (score >0.5), and LRT (score >0.9), and only those with four or more deleterious predictions were further considered.
# 
# #keep <- ( as.numeric(unlist(info(germline)$ExAC_ALL))<=0.01 | info(germline)$`1000g2015aug_all`<=0.01 )
# #keep[is.na(keep)] <- FALSE
# #germline <- germline[keep]
# 
# pred <- data.frame(
# as.numeric(unlist(info(germline)$phyloP46way_placental))>0.85,
# as.numeric(unlist(info(germline)$SIFT_score))<=0.05,
# as.numeric(unlist(info(germline)$Polyphen2_HDIV_score))>0.85,
# as.numeric(unlist(info(germline)$`GERP++_RS`))>2,
# as.numeric(unlist(info(germline)$MutationTaster_score))>0.5,
# as.numeric(unlist(info(germline)$LRT_score))>0.9
# )
# predS <- apply(pred, 1, function(x) length(which(x==TRUE)))
# keep <- predS>=4
# germline <- germline[keep]

keep <- sapply(info(germline)$CLINSIG, function(x) {
  any(grepl('Pathogenic', x)) #|
  #any(grepl('drug_response', x))
})
#keep2 <- sapply(info(germline)$CLNDBN, function(x) {
  #any(grepl('cancer', x)) 
#})
#keep <- keep | keep2


pathoTab <- data.frame(Gene=unlist(info(germline[keep])$Gene.refGene),
                       CLNDBN=sapply(info(germline[keep])$CLNDBN, function(x) paste(x, collapse = ',')),
                       dbSNP=unlist(info(germline[keep])$snp142),
                       CADD=as.numeric(info(germline[keep])$CADD_phred),
                       ExAC=as.numeric(info(germline[keep])$ExAC_ALL)
                       )

pathoTab$Amino.Acid.Uq <- NA
for(i in 1:length(info(germline[keep])$AAChange.refGene)) {
  s <- unlist(strsplit(unlist(as.vector(unlist(info(germline[keep])$AAChange.refGene[i]))), ' '))
  if(length(s)>1) {
    ts <- li$name[match(pathoTab$Gene[i], li$name2)]
    if(length(grep(ts, s))==1) {
      pathoTab$Amino.Acid.Uq[i] <- s[grep(ts, s)]  
    } else {
      pathoTab$Amino.Acid.Uq[i] <- s[1]
    }
  } else {
    pathoTab$Amino.Acid.Uq[i] <- s
  }
}

pathoTab$Amino.Acid <- sapply(strsplit(pathoTab$Amino.Acid.Uq, ':'), function(x) {
  x[length(x)]
})

a1 <- vector()
a2 <- vector()
x1 <- granges(germline[keep])
x2 <- geno(germline[keep])$GT[,2]
x2 <- strsplit(x2, '/')
for(i in 1:length(x2)) {
  if(x2[[i]][1]==0 & x2[[i]][2]==0) {
    a1[i] <- as.vector(x1[i]$REF)
    a2[i] <- as.vector(x1[i]$REF)
  } else if(x2[[i]][1]==0 & x2[[i]][2]==1) {
    a1[i] <- as.vector(x1[i]$REF)
    a2[i] <- as.vector(unlist(x1[i]$ALT))
  } else if(x2[[i]][1]==1 & x2[[i]][2]==0) {
    a1[i] <- as.vector(unlist(x1[i]$ALT))
    a2[i] <- as.vector(x1[i]$REF)
  } else if(x2[[i]][1]==1 & x2[[i]][2]==1) {
    a1[i] <- as.vector(unlist(x1[i]$ALT))
    a2[i] <- as.vector(unlist(x1[i]$ALT))
  }
}
pathoTab$GenoType <- paste0(a1,a2)
  
pathoTab <- pathoTab[,c(1,8,3,7,4,5,2)]

#patient <- yaml.load_file(opt$sample_config)[[sampleID]]

## patient info
patient <- list(
  name              = NA,#'Peppermint',
  surname           = NA,#'Patty',
  street            = NA,#'123 Cray Court',
  state             = NA,#'CA',
  city              = NA,#'San Diego',
  zip               = NA,#'92122',
  date_of_birth     = NA,#'01/01/1899',
  clinical_diagnosis= NA,#'Breast Cancer',
  tumor_type        = NA,#'Breast Cancer',
  tumor_site        = NA,#'left breast',
  specimen_type     = NA,#'fresh frozen',
  cancer_stage      = NA,#'III',
  date_of_diagnosis = NA,#'01/01/1990',
  receptor_type     = NA,#'HER2+ ER- PR-',
  molecular_type    = NA,#'HER2',
  specimen_id       = NA,#'14MS-10038 1A',
  sex               = NA#'female'
)
physician <- list(
  name       = NA,#'House',
  surname    = NA,#'Gregory',
  street     = NA,#'1 Princeton-Plainsboro Teaching Hospital',
  state      = NA,#'NJ',
  city       = NA,#'Princeton',
  zip        = NA,#'12345',
  facility   = NA,#'Princeton-Plainsboro',
  facility_id= NA#'12345'
)
pathologist <- list(
  name       = NA,#'Murphy',
  surname    = NA#'Katie'
)
medical_director <- list(
  name       = 'Spinosa',
  surname    = 'John',
  title      = 'M.D., Medical Director'
)
lab <- list(
  name       = 'Avera Genomics Laboratory',
  street     = '11099 N Torrey Pines Rd, Suite 160',
  state      = 'CA',
  city       = 'La Jolla',
  zip        = '92037',
  country    = 'USA',
  clia_number= '05DXXXXX',
  cap_number = '12345',
  cal_number = 'CLFXXXXXX',
  phone      = '858 450 2805',
  web        = 'www.avera.org'
)
sample_details <- list(
  received       = NA,#'09/17/1014',
  biopsy_date    = NA,#'09/15/1014',
  sample_purity  = NA,#'85',
  amount_rna_used= NA,#'100ng',
  seq_type       = NA,#'RNA-Seq',
  seq_protocoll  = NA#'KAPPA'
)
run_info <- list(
  run_id                     = NA,#'1234',
  kit_len                    = NA#'150'
)

sample <- list(patID = opt$id,
               patient = patient,
               physician = physician,
               pathologist = pathologist,
               medical_director = medical_director,
               lab = lab,
               sample_details = sample_details,
               run_info = run_info
)

## create report
filename <- paste0(opt$id,'_',gsub(' ', '-',Sys.time()),'.pdf')

out <- paste0(tempdir(),'/',MHmakeRandomString())
dir.create(out)
knit('patientReport.Rnw', output = paste0(out,'/patientReport.tex'))
system(paste0('pdflatex -interaction=nonstopmode -output-directory=', out ,' ', out, '/patientReport.tex'))
file.copy(paste0(out, '/patientReport.pdf'), paste(opt$out,filename,sep='/')) # move pdf to file for downloading



