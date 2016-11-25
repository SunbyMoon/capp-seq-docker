########################################
# Next Generation Sequencing QC Report
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
  make_option(c("--fastqc"), type="character", default=NULL, 
              help="path to fastqc results", metavar="character"),
  make_option(c("--qualimap"), type="character", default=NULL, 
              help="path to qualimap results", metavar="character"),
  make_option(c("--id"), type="character", default=NULL, 
              help="patient ID", metavar="character"),
  make_option(c("--id_tumor"), type="character", default=NULL, 
              help="patient tumor extension", metavar="character"),
  make_option(c("--id_blood"), type="character", default=NULL, 
              help="patient germline extension", metavar="character"),
  make_option(c("--sample_config"), type="character", default=NULL, 
              help="sample config file in .yaml format", metavar="character"),
  make_option(c("--baseSpace_token"), type="character", default=NULL, 
              help="BaseSpace token", metavar="character"),
  make_option(c("--sample_tracking"), type="character", default=NULL, 
              help="sample tracking output from contest script", metavar="character"),
  make_option(c("--outdir"), type="character", default="~", 
              help="output file path [default= %default]", metavar="character"),
  make_option(c("--scriptdir"), type="character", default="~",
              help="script file path [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

setwd(opt$scriptdir)
##
packages(knitr)
packages(yaml)

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

# sample id
sampleID <- opt$id
sample <- yaml.load_file(opt$sample_config)[[sampleID]]
sampleT <- paste0(opt$id, '-', opt$id_tumor)
sampleB <- paste0(opt$id, '-', opt$id_blood)

fastqc1FileT <- paste0(opt$fastqc,'/', opt$id, '-', opt$id_tumor, '_1_fastqc.zip')
fastqc2FileT <- paste0(opt$fastqc,'/', opt$id, '-', opt$id_tumor, '_2_fastqc.zip')
fastqc1FileB <- paste0(opt$fastqc,'/', opt$id, '-', opt$id_blood, '_1_fastqc.zip')
fastqc2FileB <- paste0(opt$fastqc,'/', opt$id, '-', opt$id_blood, '_2_fastqc.zip')
  
qualimapTextT <- paste0(opt$qualimap,'/tumor/', 'genome_results.txt',sep='')
qualimapHTMLT <- paste0(opt$qualimap,'/tumor/', 'qualimapReport.html')
qualimapTextB <- paste0(opt$qualimap,'/normal/', 'genome_results.txt',sep='')
qualimapHTMLB <- paste0(opt$qualimap,'/normal/', 'qualimapReport.html')

# sample tracking data
st <- read.csv2(opt$sample_tracking, sep='\t', stringsAsFactors = F)

## create report
filename <- paste0(opt$id,'_',gsub(' ', '-',Sys.time()),'_QC.pdf')

out <- paste0(tempdir(),'/',MHmakeRandomString())
dir.create(out)
knit('QCReport.Rnw', output = paste0(out,'/QCReport.tex'))
system(paste0('pdflatex -interaction=nonstopmode -output-directory=', out ,' ', out, '/QCReport.tex'))
file.copy(paste0(out, '/QCReport.pdf'), paste(opt$outdir,filename,sep='/')) # move pdf to file for downloading


