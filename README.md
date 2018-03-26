# NCBI Antimicrobial Resistance Gene Finder (alpha)

## Overview

This software and the accompanying database are designed to find aquired
antimicrobial resistance genes in protein or nucleotide sequences.

## Mechanism

AMRFinder can be run in two modes with protein sequence as input or with DNA
sequence as input. When run with protein sequence it uses both BLASTP and HMMER
to search protein sequences for AMR genes along with a heirarchical tree of
gene families to classify and name novel sequences. With nucleotide sequences
it uses BLASTX translated searches and the heirarchical tree of gene families.
For details see <<CITATION>>

## Installation

### Quick start

These instructions assume that you have python installed and want to use the
userspace docker clone "udocker" in a python virtualenv. We suggest you read
the more detailed explanation of the install process that follows, but this is
included for the impatient.

```shell
~$ svn co https://github.com/ncbi/pipelines/trunk/amr_finder
~$ pip install virtualenv
~$ virtualenv --python=python2 cwl
~$ source cwl/bin/activate
(cwl) ~$ pip install -U wheel setuptools
(cwl) ~$ pip install -U cwltool PyYAML cwlref-runner
(cwl) ~$ curl https://raw.githubusercontent.com/indigo-dc/udocker/master/udocker.py > cwl/bin/udocker
(cwl) ~$ chmod u+rx cwl/bin/udocker
(cwl) ~$ udocker install
(cwl) ~$ cd amr_finder
(cwl) ~/amr_finder$ ./amrfinder -p test_prot.fa
```

### Installation summary

The AMR Finder only runs on Linux and depends upon two main pieces of
software, Docker and CWL (Common Workflow Language). 

We briefly show two possible routes to installation, one using docker and one
using udocker. Docker (http://docker.com) requires root to install, so we have
also made AMRFinder compatable with udocker
(https://github.com/indigo-dc/udocker).

Note that udocker requires python2. AMRFinder itself is compatible with 
either python2 or 3. 

There are two parts to installing AMRFinder. Installing the code itself and 
installing the prerequisites. 

Prerequisites:
   - python2 or python3 (udocker requires python2)
   - docker or udocker
   - python packages
        - wheel
        - setuptools
        - PyYAML 
        - cwlref-runner
        - cwltool

### Retrieving the AMR software

The AMRFinder software is available at github at https://github.com/ncbi/pipelines/tree/master/amr_finder,
and can be retrieved with svn like:

```shell
~$ svn co https://github.com/ncbi/pipelines/trunk/amr_finder
```

### Prerequisites

You will need to install the prerequisites if they're not already installed on
your system.

e.g.,

The instructions that follow use pip and virtualenv, which are usually 
included with most python installs, so try:

```shell
~$ pip --version
~$ virtualenv --version
```
If pip is not installed see https://pip.pypa.io/en/stable/installing/ for installation instructions.

Virtualenv can be easily installed with pip:

```shell 
~$ pip install virtualenv
```

To create a virtualenv for your installation of CWL and AMRFinder:

```shell
~$ virtualenv --python=python2 cwl
```
(Note that if you're running python2 by default you will skip the '--python=python2')

### Installing CWL

```shell
~$ source cwl/bin/activate
(cwl) ~$ pip install -U wheel setuptools
(cwl) ~$ pip install -U cwltool PyYAML cwlref-runner
```

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

### Installing udocker

A simple example of how you might install udocker given the virtualenv we
created above:

```shell
~$ source cwl/bin/activate   # to enter the virtualenv
(cwl) ~$ curl https://raw.githubusercontent.com/indigo-dc/udocker/master/udocker.py > cwl/bin/udocker
(cwl) ~$ chmod u+rx cwl/bin/udocker
(cwl) ~$ udocker install
```
Your milage may vary, as the required set of package dependencies will
be different from system to system, depending upon what is already
installed. Note that we support both Python 2 & 3, however, udocker
only works with Python 2.

### Initial test run.

```shell
(cwl) ~$ cd amr_finder
(cwl) ~/amr_finder$ ./amrfinder -p test_prot.fa
```

### Shell script to run virtualenv before running amrfinder.

Because we used a virtualenv we might want to use a tiny shell script to
invoke amrfinder:

```shell
#!/bin/sh

source $HOME/cwl/bin/activate   # activate virtualenv (not necessary if you are running outside of a venv)
$HOME/amr_finder/amrfinder $@
```

## Running AMRFinder

### Typical options

Generally speaking, the only required arguments are either 
`-p <protein_fasta>` for proteins or `-n <nucleotide_fasta>` for nucleotides. 
You may also have the code update itself as well by using `-u`. This will 
update to the latest AMR database, as well as any code changes in AMRFinder.

### Complete set of options

```shell
$ ./amr_finder/amrfinder -h
usage: amrfinder [-h] [-v] [-u] [-d] [-o OUTFILE] [-i IDENT_MIN]
                 [-c COVERAGE_MIN] [-t TRANSLATION_TABLE] [-s] 
				 (-p <fasta> | -n <fasta>)

Run (and optionally update) the amr_finder pipeline.

  -p FASTA, --protein    Amino-acid sequences to search using BLASTP and HMMER
  -n FASTA, --nucleotide genomic sequence to search using BLASTX

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
```

## License

### HMMER

This distribution includes HMMER (c) Sean Eddy and the Howard Hughes Medical
Institue and licensed under the GNU General Public License version 3 (GPLv3)
(https://www.gnu.org/licenses/) 

See http://hmmer.org for details.

### PUBLIC DOMAIN NOTICE                                       
                                                           
This software/database is "United States Government Work" under the terms of
the United States Copyright Act. It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted. This software/database is freely available to the public for use
without a copyright notice. Restrictions cannot be placed on its present or
future use.                                                       
                                                           
Although all reasonable efforts have been taken to ensure  the accuracy and
reliability of the software and data, the National Center for Biotechnology
Information (NCBI) and the U.S. Government do not and cannot warrant the
performance or results that may be obtained by using this  software or data.
NCBI and the U.S. Government disclaim all warranties as to performance,
merchantability or fitness for any particular purpose.                        
                                                           
In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using:
<<CITATION>> as the citation.
