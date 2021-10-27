#FROM google/deepvariant:test-20201228-mod-gpu
FROM ubuntu:20.04
MAINTAINER Kishwar Shafin, kishwar.shafin@gmail.com

# update and install dependencies
RUN apt-get update && \
    apt-get -y install wget git cmake make gcc g++ autoconf bzip2 lzma-dev zlib1g-dev tabix libbz2-dev && \
    apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev liblzma-dev libhdf5-dev && \
    apt-get -y install python3-pip python3-virtualenv virtualenv && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && \
    apt-get -y install time git make wget autoconf gcc g++ && \
    apt-get -y install autoconf bzip2 lzma-dev zlib1g-dev && \
    apt-get -y install libcurl4-openssl-dev libpthread-stubs0-dev libbz2-dev && \
    apt-get -y install liblzma-dev libhdf5-dev libncurses5-dev libncursesw5-dev && \
    apt-get -y install python3-dev python3-pip && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install cython
RUN python3 -m pip install wheel
RUN python3 -m pip uninstall -y enum34
RUN python3 -m pip install pysam

# install samtools
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar -xvjf samtools-1.11.tar.bz2 && \
    cd samtools-1.11 && \
    make && \
    make install

# install bcftools
WORKDIR /opt
RUN wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && \
    tar -xvjf bcftools-1.11.tar.bz2 && \
    cd bcftools-1.11 && \
    make && \
    make install


RUN python3 -m pip install --upgrade pip

# setup scripts for variant calling
WORKDIR /opt
COPY filter.py /opt/
COPY vcf_merge_t2t.py /opt/

RUN rm -rf /opt/run_pepper_margin_deepvariant/

WORKDIR /opt
