FROM ubuntu:16.04

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y git wget

# -- ivrTyper General Dependencies ---
RUN apt-get install -y python-dev default-jre
RUN pip install pip install pysam==0.11.2.2

# --- ReMatCh ---
RUN git clone https://github.com/B-UMMI/ReMatCh.git -b base_quality_filter
ENV PATH="/NGStools/ReMatCh:${PATH}"

# --- ivrTyper ---
RUN git https://github.com/cimendes/ivrTyper.git
ENV PATH="/NGStools/ivrTyper:${PATH}"

# Clean
RUN apt-get remove -y wget && apt-get autoclean -y

WORKDIR /data/