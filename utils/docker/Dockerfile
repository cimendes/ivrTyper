FROM ubuntu:16.04

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y git wget

# -- ivrTyper General Dependencies ---
RUN apt-get install -y python3-dev python3-pip 
RUN pip3 install pysam==0.11.2.2

# --- ivrTyper ---
RUN git clone https://github.com/cimendes/ivrTyper.git
ENV PATH="/NGStools/ivrTyper:${PATH}"
ENV PATH="/NGStools/ivrTyper/utils:${PATH}"
#fixing permissions for ivrTyper executables
RUN chmod +x /NGStools/ivrTyper/src/bowtie2-2.2.9/* && chmod +x /NGStools/ivrTyper/src/samtools-1.3.1/bin/*

WORKDIR /data/