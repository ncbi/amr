# NCBI Antimicrobial Resistance Gene Finder (alpha)

## Overview

This software is designed to find antibiotic resistance genes in the chromosomes and plasmids of bacterial genomes.  

## Mechanism

## Installation

The AMR Finder only runs on Linux and depends upon two pieces of software, Docker and CWL.

### Retrieving the AMR software
```shell
~$ svn co https://github.com/ncbi/pipelines/trunk/amr_finder
```

### Installing Docker

Detailed instructions may be found on the docker website. [Docker Install](https://docs.docker.com/install/).
Note that it requires root access to install, and the user who will be running the software will need to be in the docker group. The required docker containers images will download automatically the first time the pipeline runs. Afterwards, they will be cached and subsequent runs will execute much faster.

### Installing CWL 
You will also need to install a CWL Runner to execute the pipeline. The reference platform is suitable, and may be installed by the user as a Python Package.

```shell
~$ python3 -m venv cwl
~$ source cwl/bin/activate
(cwl) ~$ pip install -U wheel setuptools
(cwl) ~$ pip install -U cwltool PyYAML
```

## Usage

### Initial test run.
```shell
(cwl) ~$ cd amr_finder
(cwl) ~/amr_finder$ ./amrfinder -p test_prot.fa
```

You may optionally add the amr_finder directory to the default path and run it from anywhere.

### Complete set of options

```shell
$ ./amr_finder/amrfinder -h
usage: amrfinder [-h] [-x] [-u] [-v] [-d] [-o OUTFILE] [-i IDENT_MIN]
                 [-c COVERAGE_MIN] [-t TRANSLATION_TABLE] [-s] (-p | -n)
                 fasta

Run (and optionally update) the amr_finder pipeline.

positional arguments:
  fasta                 FASTA file containing the query sequence(s).

optional arguments:
  -h, --help            show this help message and exit
  -x, --check-update    Check for any updates to this pipeline (default:
                        False)
  -u, --update          Update this code from the source control system
                        (default: False)
  -v, --version         Print version information
  -d, --no_parse_deflines
                        Do not use -parse_deflines option for blast (sometimes
                        fixes issues with format of the input FASTA file being
                        automatically parsed by BLAST)
  -o OUTFILE, --output OUTFILE
                        tabfile output to this file instead of STDOUT
  -i IDENT_MIN, --ident_min IDENT_MIN
                        Minimum proportion identical translated AA residues
                        (default: 0.9).
  -c COVERAGE_MIN, --coverage_min COVERAGE_MIN
                        Minimum coverage of reference protein sequence
                        (default: 0.5).
  -t TRANSLATION_TABLE, --translation_table TRANSLATION_TABLE
                        Translation table for blastx (default: 11).
  -s, --show_output     Show the stdout and stderr output from the pipeline
                        execution.
  -p, --protein         Amino-acid sequences to search using BLASTP and HMMER
  -n, --nucleotide      genomic sequence to search using BLASTX
```