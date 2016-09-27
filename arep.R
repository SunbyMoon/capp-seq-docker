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
  make_option(c("-a", "--anno"), type="character", default=NULL, 
              help="annotation file", metavar="character"),
  make_option(c("-i", "--id"), type="character", default=NULL, 
              help="patient ID", metavar="character"),
  make_option(c("-sc", "--sample_config"), type="character", default=NULL, 
              help="sample config file in .yaml format", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default="~", 
              help="output file path [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$anno)){
  print_help(opt_parser)
  stop("A annotation file needs to be provided", call.=FALSE)
}
if (is.null(opt$id)){
  print_help(opt_parser)
  stop("A patient ID needs to be provided", call.=FALSE)
}


##
packages(knitr)

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

## parse annotation file
tab <- read.csv2(opt$anno, 
                 stringsAsFactors = F, 
                 sep = '\t')

tabVarscan <- subset(tab, VariantCaller=='Varscan2 Somatic' | is.na(VariantCaller))[,c('Gene','Variant','Amino.Acid')]
tabMuTect <- subset(tab, VariantCaller=='MuTect')[,c('Gene','Variant','Amino.Acid')]

tab <- merge(tabVarscan, tabMuTect, by=c('Variant', 'Gene', 'Amino.Acid'), all=T)
# sort by gene nanme
tab <- tab[order(tab$Gene),]

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



