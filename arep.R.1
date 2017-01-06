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

tab <- data.frame(Variant=names(ranges(vcf)),
                  Gene=unlist(lapply(info(vcf)$Gene.refGene, function(x) paste(x, collapse = ' '))),
                  Amino.Acid=unlist(lapply(info(vcf)$AAChange.refGene, function(x) paste(x, collapse = ' '))),
                  Type=info(vcf)$TYPE,
                  Status=info(vcf)$STATUS,
                  AlleleFrequency=geno(vcf)$AF[,1],
                  RadialSVM=unlist(info(vcf)$RadialSVM_pred),
                  CADD=unlist(info(vcf)$CADD_phred),
                  ExAC=as.numeric(as.vector(unlist(info(vcf)$ExAC_ALL))),
                  CLINSIG=unlist(lapply(as.vector(info(vcf)$CLINSIG), paste, collapse='|')),
                  COSMIC=unlist(lapply(info(vcf)$cosmic70, function(x) paste(x, collapse = ','))),
                  ICGC=unlist(info(vcf)$ICGC_Id)
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

tab <- cbind(tab, Info=inf)

# Variant Rank
# path <- ifelse(tab$RadialSVM=='D' & as.numeric(as.vector(tab$CADD))>=20, 2, ifelse(tab$RadialSVM=='D' | as.numeric(as.vector(tab$CADD))>=20, 1, NA))
# info <- ifelse(tab$Info!='', ifelse(tab$Info=='driver,actionable', 2, 1), 0)
# rank <- apply(cbind(path, info), 1, sum, na.rm=T)
# tab <-  cbind(tab, path, info, rank)

pred <- ifelse(tab$RadialSVM=='D' & as.numeric(as.vector(tab$CADD))>=20, 1, ifelse(tab$RadialSVM=='T' & as.numeric(as.vector(tab$CADD))<20, -1, ifelse(tab$RadialSVM=='D' | as.numeric(as.vector(tab$CADD))>=20, 0.5, NA)))
info <- ifelse(tab$Info!='', ifelse(tab$Info=='driver,actionable', 1, 0.5), 0)
exac <- ifelse(tab$ExAC<=0.01, 1, ifelse(tab$ExAC>0.01, -1, NA))
known <- ifelse(tab$ICGC!='.' | tab$COSMIC!='.', 0.5, 0)
clin <- ifelse(grepl('Pathogenic', tab$CLINSIG), 1, ifelse(grepl('Benign', tab$CLINSIG), -1, 0))

rank <- 3+apply(cbind(pred,info,exac,known,clin), 1, sum, na.rm=T)

rankS <- sapply(rank, function(x) {
  if(x<1.5) {
    'benign'
  } else if(x>=1.5 & x<=2) {
    'likely benign'
  } else if(x>2 & x<=3) {
    'VUS, benign'
  } else if(x>3 & x<=4) {
    'VUS'
  } else if(x>4 & x<=5) {
    'VUS, pathogenic'
  } else if(x>5 & x <=6) {
    'likely pathogenic'
  } else if(x>6) {
    'pathogenic'
  }
})

tab <- cbind(tab, rank, rankS)
tab <- arrange(tab, rank, Gene, decreasing=T)

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



