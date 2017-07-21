#################################################################
# Dockerfile
#
# Description:      Docker image for Avera Cancer Institute's PANCeq pipeline to call variants and CNVs from targeted deep sequencing data. 
# Base Image:       ubuntu:16.04
# Pull Cmd:         docker pull anu9109/capp-seq
# Run Cmd:          docker run -it anu9109/capp-seq
#################################################################

# Set the base image to Ubuntu 16.04
FROM ubuntu:16.04

################## BEGIN INSTALLATION ###########################

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
  apt-get install -y \
  git \
  autoconf \
  curl \
  gcc \
  make \
  gawk \
  g++ \
  perl \
  pkg-config \ 
  zlib1g-dev \
  wget \
  libncurses5-dev \
  libcurl4-gnutls-dev \
  libgnutls-dev \
  libssl-dev \
  libexpat1-dev \
  libgd-gd2-perl \
  cpanminus \
  build-essential \
  libgd-dev \
  nettle-dev \
  bioperl \
  r-base \
  r-base-dev \
  r-cran-xml \
  libxml2-dev \
  python-pip \
  texlive-full \
  nano \
  default-jre \ 
  default-jdk && \
  apt-get clean && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# create /opt/software
RUN mkdir -p /opt/software

# install BWA-0.7.15
RUN cd /opt/software/ && \
  wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 && \
  tar -xvjf /opt/software/bwa-0.7.15.tar.bz2 && \
  cd /opt/software/bwa-0.7.15 && \
  make && \
  rm /opt/software/bwa-0.7.15.tar.bz2

# install samtools-1.3
RUN cd /opt/software/ && \

  wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 && \
  tar -xvjf /opt/software/samtools-1.3.tar.bz2 && \
  cd /opt/software/samtools-1.3 && \
  make && \
  make install
RUN cd /opt/software/samtools-1.3/htslib-1.3 && \
  ./configure && \
  make && \
  make install && \
  rm /opt/software/samtools-1.3.tar.bz2

# install sambamba-0.6.3
RUN cd /opt/software/ && \
  wget https://github.com/lomereiter/sambamba/releases/download/v0.6.3/sambamba_v0.6.3_linux.tar.bz2 && \
  tar -xvjf sambamba_v0.6.3_linux.tar.bz2 && \
  mv /opt/software/sambamba_v0.6.3 /opt/software/sambamba && \
  rm /opt/software/sambamba_v0.6.3_linux.tar.bz2

# install vcftools-0.1.14 (latest version)
RUN  cd /opt/software/ && \
  git clone https://github.com/vcftools/vcftools.git && \
  cd /opt/software/vcftools && \
  ./autogen.sh && \
  ./configure && \
  make && \
  make install

# install vcflib-1.0.0-rc1 (latest version)
RUN  cd /opt/software/ && \
  git clone --recursive https://github.com/vcflib/vcflib.git && \
  cd /opt/software/vcflib && \
  make openmp

# install abra-0.97
RUN  cd /opt/software/ && \
  wget https://github.com/mozack/abra/releases/download/v0.97/abra-0.97-SNAPSHOT-jar-with-dependencies.jar && \
  mv /opt/software/abra-0.97-SNAPSHOT-jar-with-dependencies.jar /opt/software/abra.jar

# install samblaster-0.1.22 (latest version)
RUN  cd /opt/software/ && \
  git clone git://github.com/GregoryFaust/samblaster.git && \
  cd /opt/software/samblaster && \
  make

# install R & R packages
RUN wget https://cran.r-project.org/src/base/R-3/R-3.3.0.tar.gz && tar -zxvf R-3.3.0.tar.gz && cd R-3.3.0 && ./configure && make && make install && rm -r /R-3.3.0.tar.gz
RUN echo "local({r <- getOption('repos'); r['CRAN'] <- 'https://cran.cnr.berkeley.edu'; options(repos=r)})" > ~/.Rprofile
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' \ 
  -e 'biocLite("Rsamtools")' \
  -e 'biocLite("VariantAnnotation")' \
  -e 'biocLite("GenomicRanges")' \
  -e 'biocLite("myvariant")'
RUN Rscript -e 'install.packages("plyr")'
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("magrittr")'
RUN Rscript -e 'install.packages("knitr")'
RUN Rscript -e 'install.packages("httr")'
RUN Rscript -e 'install.packages("jsonlite")'
RUN Rscript -e 'install.packages("yaml")'
RUN Rscript -e 'install.packages("xtable")'
RUN Rscript -e 'install.packages("XML")'

# install CNVkit-0.8.5
RUN pip install --upgrade pip && pip install cnvkit==0.8.5

# install VarDict-1.5.1
RUN  cd /opt/software/ && \
  git clone --recursive https://github.com/AstraZeneca-NGS/VarDictJava.git && \
  cd VarDictJava/ && \
  git checkout v1.5.1 -b v1.5.1 && \
  ./gradlew clean installApp

# install FastQC-0.11.5
RUN  cd /opt/software/ && \
  wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
  unzip fastqc_v0.11.5.zip && \
  cd /opt/software/FastQC/ && \
  chmod 755 fastqc  && \ 
  rm /opt/software/fastqc_v0.11.5.zip

# install Qualimap-2.2
RUN  cd /opt/software/ && \
  wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.zip && \
  unzip qualimap_v2.2.zip && \
  cd /opt/software/qualimap_v2.2/scripts/ && \
  Rscript installDependencies.r && \
  rm /opt/software/qualimap_v2.2.zip

# install BBMap-36.32
RUN cd /opt/software/ && \
  wget https://sourceforge.net/projects/bbmap/files/BBMap_36.32.tar.gz && \
  tar -xvzf BBMap_36.32.tar.gz && \
  rm /opt/software/BBMap_36.32.tar.gz

# install VarScan
#RUN cd /opt/software/ && \
#  git clone https://github.com/dkoboldt/varscan.git 

# install Docker
RUN apt-get update && apt-get install -y apt-transport-https ca-certificates software-properties-common
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -
RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu xenial stable"
RUN apt-get update && apt-get install -y docker-ce

# install Picard-2.10.3
RUN cd /opt/software/ && \ 
  git clone https://github.com/broadinstitute/picard.git && \
  cd /opt/software/picard/ && \
  git checkout 2.10.3 -b 2.10.3 && \
  ./gradlew shadowJar

# set environmental variables
ENV PATH=$PATH:/opt/software:/opt/software/vcflib/bin:/opt/software/samblaster:/opt/software/samtools-1.3/htslib-1.3:/opt/software/VarDict-1.5.1/bin:/opt/software/FastQC:/opt/software/qualimap_v2.2:/home

# Add scripts and files to image
ADD arep.R /home/
ADD patientReport.Rnw /home/
ADD cont.R /home/
ADD contPanel.csv /home/
ADD cnvkit_filter.R /home/
ADD PANCeq_CNV_capture_targets_6.bed /home/
ADD PANCeq_C200X50bp_6.bed /home/
#ADD table_annovar.pl /home/
#ADD annotate_variation.pl /home/
#ADD convert2annovar.pl /home/
ADD variant_filter.R /home/
ADD clinical_ann_metadata-snvs_canc.csv /home/
ADD createQCreport.R /home/
ADD QCReport.Rnw /home/
ADD longest_resfseq.txt /home/
ADD longest_isoform.R /home/
ADD BaseSpaceRunInfo.R /home/
ADD cancer_genes.csv /home/
ADD avera_letter.pdf /home/
ADD GAVIN_ruleguide_r0.3.tsv /home/
ADD hg19.dac.bed /home/
ADD hg19.duke.bed /home/

##################### INSTALLATION END ##########################

# File Author / Maintainer
MAINTAINER Anu Amallraja <anu9109@gmail.com>
