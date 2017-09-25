# ivrTyper
-- ivr locus allele determination from genomic data --

The *ivr* locus, also known as *cod locus*, is a random six-phase switch in a Type I restriction-modification system.
It has been shown to be the responsible underlying mechanism for phase-variation in *streptococcus pneumonia* (Manso et al, 2014).

There are six possible alleles, identified from A to E:

    - Allele A (modules 1.1 and 2.1)
    - Allele B (modules 1.2 and 2.2)
    - Allele E (modules 1.2 and 2.3)
    - Allele D (modules 1.2 and 2.1)
    - Allele C (modules 1.2 and 2.2)
    - Allele F (modules 1.2 and 2.3)


Due to the high variation rate and structural rearrangement, this allele cannot reliably be determined using assembly or standard mapping approaches.
Therefore, we present **ivrTyper**, a standalone tool for the determination of the *ivr locus* allele from paired-end genomic data.
We've adapted the algorithm described by Lee et al (2017). The reads are mapped to the full *iver locus* and the mates of reads mapping to the reverse strand of the conserved 5' region are extracted for each sample and mapped with bowtie2 to the 1.1 and 1.2 modules. The representative module of the *hsdS* gene is chosen if the number of reads mapping is greater than the `--proportionCutOff` (default is 80%). The 3' allele is determined using the chosen 1.x module as target, and mapping the mates of reads to the module 2.x reference sequences.
For a module to be considered present it needs to have a number of reads mapped higher than `--minCoverage`.
The number of mates that map to the 1.x module is obtained though the conserved 5' region. The number of mates mapping to the 2.x module is obtained though the selected 1.x module. The calculation of the proportion takes into consideration the reads mapping to all 2.x possibilities only.


#### Dependencies

Required to download sequence data from ENA database:
- Aspera Connect 2 >= v3.6.1

Required to run analysis:
- Bowtie2 >= v2.2.9
- Samtools = v1.3.1

(these three executables are provided, but user's own executables can be used by providing --doNotUseProvidedSoftware option)
- pysam >= 0.11.2.2 (available though pip)
- python 2.7


#### Installation

`git clone https://github.com/cimendes/ivrTyper.git`


#### Usage

    usage: ivrTyper.py [-h] [--version] -w /path/to/workdir/directory/ [-j N] [-u]
                       [-k] [-m N] [-c N] [-gt N] [-rf /path/to/run.*.log]
                       [-a /path/to/asperaweb_id_dsa.openssh] [-kd] [-ip ILLUMINA]
                       [-pm HiSeq] [-ls GENOMIC]
                       [-l /path/to/list_IDs.txt | -t "Streptococcus pneumoniae"]

    Reads mapping against pneumococcal ivr locus for typing.

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information
      -l /path/to/list_IDs.txt, --listIDs /path/to/list_IDs.txt
                            Path to list containing the IDs to be downloaded (one
                            per line) (default: None)
      -t "Streptococcus pneumoniae", --taxon "Streptococcus pneumoniae"
                            Taxon name for which fastq files will be downloaded
                            (default: None)

    General options:
      -w /path/to/workdir/directory/, --workdir /path/to/workdir/directory/
                            Path to the directory containing the read files in
                            subdirectories, one per sample (fastq.gz or fq.gz and
                            pair-end direction coded as _R1_001 / _R2_001 or _1 /
                            _2) or to be downloaded (default: None)
      -j N, --threads N     Number of threads to use (default: 1)
      -u, --skipProvidedSoftware
                            Do not use provided software (default: False)
      -k, --keepFiles       Keep alignment files (default: False)

    ivrTyper module facultative options:
      -m N, --minCoverage N
                            minimumcoverage depth a module to be present in the
                            sample (default: 5)
      -c N, --proportionCutOff N
                            Proportion cut off for the module 1.x to be chosen as
                            target (default: 0.8)
      -gt N, --greaterThan N
                            proportion cut off for the "greater than" column in
                            the final report. (default: 0.5)
      -rf /path/to/run.*.log, --reportFile /path/to/run.*.log
                            Logfile to append run information to. (default: None)

    Download facultative options:
      -a /path/to/asperaweb_id_dsa.openssh, --asperaKey /path/to/asperaweb_id_dsa.openssh
                            Download fastq files from ENA using Aspera Connect.
                            With this option, the path to Private-key file
                            asperaweb_id_dsa.openssh must be provided (normaly
                            found in
                            ~/.aspera/connect/etc/asperaweb_id_dsa.openssh).
                            (default: None)
      -kd, --keepDownloadFiles
                            Keep downloaded read files (default: False)
      -ip ILLUMINA, --instrumentPlatform ILLUMINA
                            Download files with specific library layout (available
                            options: ILLUMINA, ALL) (default: ALL)
      -pm HiSeq, --platformModel HiSeq
                            Filter download by the model of the Illumina machine.
                            (default: None)
      -ls GENOMIC, --librarySource GENOMIC
                            Filter download by library source (available options:
                            GENOMIC, ALL) (default: None)

#### Running ivrTyper

######  - With Local Samples
To run ivrTyper in local fastq files, these need to be organized in sample folders.
It is advisable to use copied fastq files or symbolic links to the original ones. Then provide the directory containing sample folders to `--workdir`. The output report will be stored there.
```
  workir/
    sample_1/
      fastq_file_a_1.fq.gz
      fastq_file_a_2.fq.gz
    sample_2/
      fastq_file_b_R1_001.fastq.gz
      fastq_file_b_R2_001.fastq.gz
```

`ivrTyper.py --workdir /path/to/data/ --threads 8`

###### - With specific ENA sequencing data
To run ivrTyper in a specific set of ENA IDs, provide a file to `--listIDs` containing a list of ENA IDs that will be downloaded. This tool requires the genomic data to be paired-end and from Illumina technology.
ReMatCh will store the output files in the `--workdir`.
`ivrTyper.py --listIDs /path/to/list_IDs.txt --workdir /path/to/output/directory/  --threads 8`

###### - With ENA sequencing data of *streptococcus pneumoniae* taxon
To run ivrTyper in all ENA data of a given taxon, provide the taxon name to `--taxon`.
The ENA Run Accession numbers for the given taxon will be stored in IDs_list.seqFromWebTaxon.tab file.
Warning - This is option is in development.
`ivrTyper.py --taxon "Streptococcus pneumoniae" --workdir /path/to/output/directory/  --threads 8`



#### Output

**run.*.log**
ivrTyper running log file.

**report_*.csv**
Report for the samples that ran the ivrTyper successfully.
- **Sample** - Sample name
- **1.1** - Number of reads that mapped to the 1.1 module
- **1.2** - Number of reads that mapped to the 1.2 module
- **2.1** - Number of reads that mapped to the 2.1 module, using the propper 1.x module as target
- **2.2** - Number of reads that mapped to the 2.2 module, using the propper 1.x module as target
- **2.3** - Number of reads that mapped to the 2.3 module, using the propper 1.x module as target
- **pA** - proportion of reads that map to the allele A
- **pB** - proportion of reads that map to the allele B
- **pE** - proportion of reads that map to the allele E
- **pD** - proportion of reads that map to the allele D
- **pC** - proportion of reads that map to the allele C
- **pF** - proportion of reads that map to the allele F
- **Most Prevalent** - Most prevalent *ivr* type in the sample
- **gt0.5** - *ivr* type with a proportion greater than 0.5 (or whatever is set by `--greaterThan`)


#### References
Lees, J.A., Kremer, P.H.C., Manso, A.S., Croucher, N.J., Ferwerda, B., Serón, M.V., Oggioni, M.R., Parkhill, J., Brouwer, M.C., Ende, A. Van Der, Beek, D. Van De, Bentley, S.D., 2017. Large scale genomic analysis shows no evidence for pathogen adaptation between the blood and cerebrospinal fluid niches during bacterial meningitis. Microb. Genomics 3, 1–12. doi:10.1099/mgen.0.000103

Li, J., Li, J.W., Feng, Z., Wang, J., An, H., Liu, Y., Wang, Y., Wang, K., Zhang, X., Miao, Z., Liang, W., Sebra, R., Wang, G., Wang, W.C., Zhang, J.R., 2016. Epigenetic Switch Driven by DNA Inversions Dictates Phase Variation in Streptococcus pneumoniae. PLoS Pathog. 12, 1–36. doi:10.1371/journal.ppat.1005762

Manso, A.S., Chai, M.H., Atack, J.M., Furi, L., De Ste Croix, M., Haigh, R., Trappetti, C., Ogunniyi, A.D., Shewell, L.K., Boitano, M., Clark, T.A., Korlach, J., Blades, M., Mirkes, E., Gorban, A.N., Paton, J.C., Jennings, M.P., Oggioni, M.R., 2014. A random six-phase switch regulates pneumococcal virulence via global epigenetic changes. Nat. Commun. 5, 5055. doi:10.1038/ncomms6055


#### Contact
Catarina Mendes <cimendes@medicina.ulisboa.pt>
