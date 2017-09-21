ivrTyper - Docker
=================
-- ivr locus allele determination from genomic data --

<https://github.com/B-UMMI/INNUca>


This is a dockerfile for using INNUca, with all dependencies already installed.

Within this container you can find:
- ubuntu:16.04
- git
- Python v2.7
- pysam v0.11.2.2
- [ReMatCh](https://github.com/B-UMMI/ReMatCh) v3.2
- [ivrTyper](https://github.com/cimendes/ivrTyper.git) v0.4.6

### Using play-with-docker
[![Try in PWD](https://cdn.rawgit.com/play-with-docker/stacks/cff22438/assets/images/button.png)](http://labs.play-with-docker.com/)

Within [play-with-docker](http://labs.play-with-docker.com/) webpage click on **create session**. Then, another page
will open with a big counter on the upper left corner. Click on **+ add new instance** and a terminal like instance should be generated on the right. On
this terminal you can load this docker image as follows:

`docker pull cimendes/ivrTyper`

#### Build this docker on your local machine

For this, docker needs to be installed on your machine. Instructions for this can be found [here](https://docs.docker.com/engine/installation/).

##### Using DockerHub (automated build image)

`docker pull cimendes/ivrTyper`

##### Using GitHub (build docker image)

1) `git clone https://github.com/cimendes/ivrTyper.git -b v2.8`
2) `docker build -t ivrTyper ./ivrTyper/Docker/`

### Run (using automated build image)
    docker run --rm -u $(id -u):$(id -g) -it -v /local/folder/fastq_data:/data/ cimendes/ivrTyper ivrTyper.py --workdir /data/ --threads 8


Contact
-------
Catarina Mendes
<cimendes@medicina.ulisboa.pt>