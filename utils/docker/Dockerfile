FROM sangerpathogens/seroba:latest

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y git wget
RUN apt-get install -y lbzip2 unzip

# -- ivrTyper General Dependencies ---
RUN apt-get install -y python-dev python-pip
RUN pip2 install pysam==0.11.2.2

# --- ReMatCh ---
RUN git clone https://github.com/B-UMMI/ReMatCh.git
ENV PATH="/NGStools/ReMatCh:${PATH}"

# --- ivrTyper ---
RUN git clone https://github.com/cimendes/ivrTyper.git
ENV PATH="/NGStools/ivrTyper:${PATH}"
ENV PATH="/NGStools/ivrTyper/utils:${PATH}"

# --- pneumoCaT Dependencies ---
RUN apt-get install -y make gcc libxml2-dev libxslt1-dev zlib1g-dev g++ build-essential libncurses5-dev libncursesw5-dev
#Bowtie2
RUN wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip && unzip bowtie2-2.2.9-linux-x86_64.zip && rm bowtie2-2.2.9-linux-x86_64.zip
ENV PATH="/NGStools/bowtie2-2.2.9:${PATH}"
#Samtools
RUN wget https://downloads.sourceforge.net/project/samtools/samtools/0.1.19/samtools-0.1.19.tar.bz2  && tar -xf samtools-0.1.19.tar.bz2 && rm samtools-0.1.19.tar.bz2
WORKDIR samtools-0.1.19
RUN make
ENV PATH="/NGStools/samtools-0.1.19:${PATH}"
WORKDIR /NGStools/
#Python dependencies
RUN pip2 install PyYaml numpy lxml biopython
# -- PneumoCaT ---
RUN git clone https://github.com/phe-bioinformatics/PneumoCaT.git
ENV PATH="/NGStools/PneumoCaT:${PATH}"

# Clean
RUN apt-get remove -y wget make gcc libxml2-dev libxslt1-dev zlib1g-dev g++ build-essential libncurses5-dev libncursesw5-dev && apt-get autoclean -y

WORKDIR /data/