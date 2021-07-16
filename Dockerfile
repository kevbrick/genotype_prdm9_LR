FROM r-base

LABEL authors="Kevin Brick" \
      description="Docker image containing all software requirements for calling PRDM9 haplotypes from long read data"

## Install gcc
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y make perl cpanminus g++ tar zip libghc-zlib-dev autoconf pkg-config
RUN apt-get install -y wget curl libcurl4 libcurl4-openssl-dev libblas-dev
RUN apt-get install -y bioperl emboss ncbi-entrez-direct

RUN apt-get install -y libssl-dev 
RUN apt-get install -y gtk3-nocsd
RUN apt-get install -y ncbi-blast+

RUN cpanm -v Math::Round && \
    cpanm -v List::Util && \
    cpanm -v List::MoreUtils && \
    cpanm -v List::UtilsBy && \
    cpanm -v Getopt::Long

## GET BEDTOOLS
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && \
    mkdir -p /usr/local/bedtools2/bin && \
    mv bedtools.static.binary /usr/local/bedtools2/bin/bedtools && \ 
    echo '#!/bin/bash\n${0%/*}/bedtools slop "$@"' >/usr/local/bedtools2/bin/slopBed  && \
    echo '#!/bin/bash'                      >/usr/local/bedtools2/bin/intersectBed && \
    echo '${0%/*}/bedtools intersect "$@"' >>/usr/local/bedtools2/bin/intersectBed && \
    echo '#!/bin/bash'                      >/usr/local/bedtools2/bin/mergeBed && \
    echo '${0%/*}/bedtools merge "$@"'     >>/usr/local/bedtools2/bin/mergeBed && \
    echo '#!/bin/bash'                      >/usr/local/bedtools2/bin/mapBed && \
    echo '${0%/*}/bedtools map "$@"'       >>/usr/local/bedtools2/bin/mapBed && \
    chmod a+x /usr/local/bedtools2/bin/*

RUN R -e "if('optparse'     %in% rownames(installed.packages()) == FALSE){install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')}"

RUN apt-get install -y bzip2 libbz2-dev liblzma-dev && \
    cd /usr/local && \
    wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
    tar -xjvf htslib-1.12.tar.bz2 && \
    cd htslib-1.12 && \
    ./configure --prefix=/usr/local/htslib && \
    make && \
    make install

RUN cd /usr/local && \
    wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar -xjvf samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix=/usr/local/samtools && \
    make && \
    make install

RUN mkdir /usr/local/genotypePRDM9

ENV PATH=/usr/local/bedtools2/bin:/usr/local/samtools/bin:/usr/local/htslib/bin:$PATH
	
ADD Alleva_et_al_2021*txt      /usr/local/genotypePRDM9/
ADD extractPRDM9ZFsFromFA.pl   /usr/local/genotypePRDM9/
ADD getZFAhaplotypesFromFA.pl  /usr/local/genotypePRDM9/
ADD prdm9_test_*longreads.fa   /usr/local/genotypePRDM9/

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

