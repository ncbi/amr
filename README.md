# NCBI Antimicrobial Resistance Gene Finder (alpha)

## Overview

This software is designed to find antibiotic resistance genes in the
chromosomes and plasmids of bacterial genomes.

## Mechanism

## Installation

The AMR Finder only runs on Linux and depends upon two pieces of
software, Docker and CWL. We also support udocker, which is
installable by any stadard user account and does not need special
system access.

### Retrieving the AMR software

```shell
~$ svn co https://github.com/ncbi/pipelines/trunk/amr_finder
```

### Installing CWL + udocker

You will also need to install a CWL Runner to execute the
pipeline. The reference platform is suitable, and may be installed by
the user as a Python Package. Here is typical installion:

```shell
~$ virtualenv --python=python2 cwl
~$ source cwl/bin/activate
(cwl) ~$ pip install -U wheel setuptools
(cwl) ~$ pip install -U cwltool PyYAML cwlref-runner
(cwl) ~$ curl https://raw.githubusercontent.com/indigo-dc/udocker/master/udocker.py > cwl/bin/udocker
(cwl) ~$ chmod u+rx cwl/bin/udocker
(cwl) ~$ udocker install
```

Your milage may vary, as the required set of package dependencies will
be different from system to system, depending upon what is already
installed. Note that we support both Python 2 & 3, however, udocker
only works with Python 2.


### Installing Docker

If it is not restricted by your local security policy, it is
recommended that you install Docker instead of the more limited
udocker. Detailed instructions may be found on the docker
website. [Docker Install](https://docs.docker.com/install/). Please
install the latest version of docker, and not the one that comes with
your distribution.  Note that it requires root access to install, and
the user who will be running the software will need to be in the
docker group. The required docker containers images will download
automatically the first time the pipeline runs. Afterwards, they will
be cached and subsequent runs will execute much faster.

## Usage

### Initial test run.
```shell
(cwl) ~$ cd amr_finder
(cwl) ~/amr_finder$ ./amrfinder -p test_prot.fa
```

You may optionally add the amr_finder directory to the default path
and run it from anywhere.

### Typical options

Generally speaking, the only required argument is either `-p` for
proteins or `-n` for nucleotides. This is used to denote the sequence
type of the input FASTA.  You may also have the code update itself as
well by using `-u`. This will not only bring in the latest code
changes, it will also update the search databases.

### Complete set of options

```shell
$ ./amr_finder/amrfinder -h
usage: amrfinder [-h] [-v] [-u] [-d] [-o OUTFILE] [-i IDENT_MIN]
                 [-c COVERAGE_MIN] [-t TRANSLATION_TABLE] [-s] (-p | -n)
                 fasta

Run (and optionally update) the amr_finder pipeline.

positional arguments:
  fasta                 FASTA file containing the query sequence(s).

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version information
  -u, --update          Update this code from the source control system
                        (default: False)
  -d, --parse_deflines  Use -parse_deflines option for blast. This sometimes
                        fixes issues with format of the input FASTA file being
                        automatically parsed by BLAST. (default: False)
  -o OUTFILE, --output OUTFILE
                        tabfile output to this file instead of STDOUT
  -i IDENT_MIN, --ident_min IDENT_MIN
                        Minimum proportion identical translated AA residues
                        (default: 0.9).
  -c COVERAGE_MIN, --coverage_min COVERAGE_MIN
                        Minimum coverage of reference protein sequence
                        (default: 0.5).
  -t TRANSLATION_TABLE, --translation_table TRANSLATION_TABLE
                        Translation table for blastx (default: 11). More info
                        may be found at https://www.ncbi.nlm.nih.gov/Taxonomy/
                        Utils/wprintgc.cgi?mode=c
  -s, --show_output     Show the stdout and stderr output from the pipeline
                        execution.
  -p, --protein         Amino-acid sequences to search using BLASTP and HMMER
  -n, --nucleotide      genomic sequence to search using BLASTX
```