# ivrTyper - Utils
-- ivr locus allele determination from genomic data --

Table of Contents
--

[Restart Typer](#restart-typer)
[Run ivrTyper and ReMatCh](#run-ivrtyper-rematch)

## Restart ivrTyper

Restart a ivrTyper run abruptly terminated.

**Dependencies:**
- Python 2.7.x

**Usage:**

    usage: restartTyper.py [-h] [--version] -i /path/to/initial/workdir/directory/
                           [-w /path/to/workdir/directory/] [-j N]

    Restart a ivr Typer run abruptly terminated

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information

    Required options:
      -i /path/to/initial/workdir/directory/, --initialWorkdir /path/to/initial/workdir/directory/
                            Path to the directory where ivr Typer was running
                            (default: None)

    General facultative options:
      -w /path/to/workdir/directory/, --workdir /path/to/workdir/directory/
                            Path to the directory where irv Typer will run again
                            (default: .)
      -j N, --threads N     New number of threads to use instead (default: None)


## Run ivrTyper ReMatCh

Reads mapping against pneumococcal *ivr* locus for typing and locus evaluation - ReMatCh and ivrTyper integration.
This script allows to obtain information on the *ivr* locus (type and conformation though ivrTyper and ReMatCh), as well as information on the strain (mlst though ReMatCh and serotype through PneumoCaT).

**Dependencies:**
- Python 2.7.x
- [ivrTyper](https://github.com/cimendes/ivrTyper)
- [ReMatCh](https://github.com/B-UMMI/ReMatCh) >= 3.2 (in path)
- [PneuCaT](https://github.com/phe-bioinformatics/PneumoCaT) 1.0 (in path)

**Facultative Dependencies**
- Aspera Connect 2 >= v3.6.1


**Usage:**

    usage: runTyperReMatCh.py [-h] [--version] -w /path/to/workdir/directory/
                              [-l /path/to/list_IDs.txt] [-j N]
                              [-a /path/to/asperaweb_id_dsa.openssh] [-kd]

    Reads mapping against pneumococcal ivr locus for typing and locus evaluation -
    ReMatCh and ivrTyper integration.

    optional arguments:
      -h, --help            show this help message and exit
      --version             Version information

    General options:
      -w /path/to/workdir/directory/, --workdir /path/to/workdir/directory/
                            Path to the output directory (default: None)
      -l /path/to/list_IDs.txt, --listIDs /path/to/list_IDs.txt
                            Path to list containing the IDs to be downloaded (one
                            per line) (default: None)
      -j N, --threads N     Number of threads to use (default: 1)

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