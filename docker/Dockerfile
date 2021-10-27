FROM python:latest
MAINTAINER Kishwar Shafin, kishwar.shafin@gmail.com

RUN python3 --version
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install cython
RUN python3 -m pip install wheel
RUN python3 -m pip uninstall -y enum34
RUN python3 -m pip install pysam


RUN python3 -m pip install --upgrade pip

# setup scripts for variant calling
WORKDIR /opt
COPY filter.py /opt/
COPY vcf_merge_t2t.py /opt/

RUN rm -rf /opt/run_pepper_margin_deepvariant/

WORKDIR /opt
