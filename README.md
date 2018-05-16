# NCBI Antimicrobial Resistance Gene Finder (alpha)

## Overview

This software and the accompanying database are designed to find acquired
antimicrobial resistance genes in protein or nucleotide sequences.

## Mechanism

AMRFinder can be run in two modes with protein sequence as input or with DNA
sequence as input. When run with protein sequence it uses both BLASTP and HMMER
to search protein sequences for AMR genes along with a hierarchical tree of
gene families to classify and name novel sequences. With nucleotide sequences
it uses BLASTX translated searches and the hierarchical tree of gene families.

## Installation

To run AMRFinder you will need Linux, Docker and CWL (Common Workflow
Language). We provide instructions here for two installation modes. One
using Docker and the other using the docker emulator uDocker. We recommend the
Docker installation, but installing Docker requires root, so we have made
AMRFinder also compatible with uDocker and included instructions for uDocker
installation as well. AMRFinder runs considerably slower under uDocker than 
Docker.

### Quick start

These instructions assume that you have python, pip, and docker installed. We provide more details about how to install these prerequisites below. 

```shell
(cwl) ~$ pip install -U wheel setuptools
(cwl) ~$ pip install -U cwltool[deps] PyYAML cwlref-runner
(cwl) ~$ svn co https://github.com/ncbi/pipelines/trunk/amr_finder
(cwl) ~$ cd amr_finder
(cwl) ~/amr_finder$ ./amrfinder -p test_prot.fa
```

### Installation summary

The AMR Finder only runs on Linux and depends upon two main pieces of
software, Docker and CWL (Common Workflow Language). 

We briefly show two possible routes to installation, one using docker and one
using uDocker. Docker (http://docker.com) requires root to install, so we have
also made AMRFinder compatible with uDocker
(https://github.com/indigo-dc/udocker).

Note that uDocker requires python2. AMRFinder itself is compatible with 
either Python 2 or 3. 

There are two parts to installing AMRFinder. Installing the code itself and 
installing the prerequisites. 


Prerequisites:
   - python2 or python3 (uDocker requires python2)
   - docker or uDocker
   - python packages
        - wheel
        - setuptools
        - PyYAML 
        - cwlref-runner
        - cwltool

### Retrieving the AMR software

The AMRFinder software is available at GitHub at https://github.com/ncbi/pipelines/tree/master/amr_finder,
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
(cwl) ~$ pip install -U cwltool[deps] PyYAML cwlref-runner
```

### Installing Docker

It is recommended that you install Docker instead of the more limited uDocker.
Detailed instructions may be found on the docker website. [Docker
Install](https://docs.docker.com/install/). Please install the latest version
of docker, it is usually newer than the one that comes with your distribution.
Note that it requires root access to install, and the user who will be running
the software will need to be in the docker group. The required docker
containers images will download automatically the first time the pipeline runs.
Afterwards, they will be cached and subsequent runs will execute much faster.

### Installing uDocker

UDocker compatibility is provided only because some servers may have security
policies that make it difficult to install Docker. We recommend the docker
installation if possible.

A simple example of how you might install uDocker given the virtualenv we
created above:

```shell
~$ source cwl/bin/activate   # to enter the virtualenv
(cwl) ~$ curl https://raw.githubusercontent.com/indigo-dc/udocker/master/udocker.py > cwl/bin/udocker
(cwl) ~$ chmod u+rx cwl/bin/udocker
(cwl) ~$ udocker install
```
Your mileage may vary, as the required set of package dependencies will
be different from system to system, depending upon what is already
installed. Note that we support both Python 2 & 3, however, uDocker
only works with Python 2.

### Initial test run

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

### Testing AMRFinder

A small set of test data are included with AMRFinder just to make sure things
are working

```shell
~$ source cwl/bin/activate   # to enter the virtualenv
(cwl) ~$ cd amr_finder
(cwl) ~$ ./amrfinder -p test_prot.fa -g test_prot.gff
```
You should see something like:

    Target identifier  Contig id Start Stop Strand Gene symbol Protein name                         Method  Target lengthReference protein length % Coverage of reference protein % Identity to reference protein Alignment length Accession of closest protein Name of closest protein HMM id                                                              HMM description
    blaOXA-436_partial contig1    4001 4699      + blaOXA      OXA-48 family class D beta-lactamase PARTIAL                                   233                             265                           87.92           100.00                          233 WP_058842180.1          OXA-48 family carbapenem-hydrolyzing class D beta-lactamase OXA-436 NF000387.2      OXA-48 family class D beta-lactamase
    blaPDC-114_blast   contig1    2001 3191      + blaPDC      PDC family class C beta-lactamase    BLAST                                     397                             397                          100.00            99.75                          397 WP_061189306.1          class C beta-lactamase PDC-114                                      NF000422.2      PDC family class C beta-lactamase
    blaTEM-156         contig1       1  858      + blaTEM-156  class A beta-lactamase TEM-156       ALLELE                                    286                             286                          100.00           100.00                          286 WP_061158039.1          class A beta-lactamase TEM-156                                      NF000531.2      TEM family class A beta-lactamase
    vanG               contig1    5001 6047      + vanG        D-alanine--D-serine ligase VanG      EXACT                                     349                             349                          100.00           100.00                          349 WP_063856695.1          D-alanine--D-serine ligase VanG                                     NF000091.3      D-alanine--D-serine ligase VanG


```shell
~$ ./amrfinder -n test_dna.fa
```

You should see something like:

    Target identifier      Contig id              Start Stop Strand Gene symbol Protein name                                Method Target lengthReference protein length % Coverage of reference protein % Identity to reference protein Alignment length Accession of closest protein Name of closest protein                                             HMM id                          HMM description
    blaOXA-436_partial_cds blaOXA-436_partial_cds   101  802      + blaOXA      OXA-48 family class D beta-lactamasePARTIAL    234                                   265                           88.30                          100.00              234 WP_058842180.1               OXA-48 family carbapenem-hydrolyzing class D beta-lactamase OXA-436 NF000387.2                      OXA-48 family class D beta-lactamase
    blaPDC-114_blast       blaPDC-114_blast           1 1191      + blaPDC      PDC family class C beta-lactamase            BLAST                                   397                             397                          100.00            99.75 397                          WP_061189306.1                                                      class C beta-lactamase PDC-114  NF000422.2                           PDC family class C beta-lactamase
    blaTEM-156             blaTEM-156               101  958      + blaTEM-156  class A beta-lactamase TEM-156              ALLELE                                   286                             286                          100.00           100.00 286                          WP_061158039.1                                                      class A beta-lactamase TEM-156  NF000531.2                           TEM family class A beta-lactamase
    vanG                   vanG                     101 1147      + vanG        D-alanine--D-serine ligase VanG              EXACT                                   349                             349                          100.00           100.00 349                          WP_063856695.1                                                      D-alanine--D-serine ligase VanG NF000091.3                           D-alanine--D-serine ligase VanG

## Running AMRFinder

### Typical options

The only required arguments are either
`-p <protein_fasta>` for proteins or `-n <nucleotide_fasta>` for nucleotides.
We also provide an automatic update mechanism to update the code and database
by using `-u`. This will update to the latest AMR database, as well as any code
changes in AMRFinder. Use '--help' to see the complete set of options and
flags.

### Input file formats

`-p <protein_fasta>` and `-n <nucleotide_fasta>`:
FASTA files are in standard format. The identifiers reported in the output are
the first non-whitespace characters on the defline.

`-g <gff_file>`
GFF files are used to get sequence coordinates for AMRFinder hits from protein
sequence. The identifier from the identifier from the FASTA file is matched up
with the 'Name=' attribute from field 9 in the GFF file. See test_prot.gff for
a simple example. (e.g., `amrfinder -p test_prot.fa -g test_prot.gff` should
result in the sample output shown below)

### Output format

AMRFinder output is in tab-delimited format (.tsv). The output format depends
on the options `-p`, `-n`, and `-g`. Protein searches with gff files (`-p
<file.fa> -g <file.gff>` and translated dna searches (`-n <file.fa>`) will also 
include contig, start, and stop columns. 

A sample AMRFinder report:

    Target identifier  Contig id Start Stop Strand Gene symbol Protein name                                   Method  Target length Reference protein length % Coverage of reference protein % Identity to reference protein Alignment length Accession of closest protein Name of closest protein                                             HMM id     HMM description
    blaOXA-436_partial contig1    4001 4699      + blaOXA      OXA-48 family class D beta-lactamase           PARTIAL           233                      265                           87.92                          100.00              233 WP_058842180.1               OXA-48 family carbapenem-hydrolyzing class D beta-lactamase OXA-436 NF000387.2 OXA-48 family class D beta-lactamase
    blaPDC-114_blast   contig1    2001 3191      + blaPDC      PDC family class C beta-lactamase              BLAST             397                      397                          100.00                           99.75              397 WP_061189306.1               class C beta-lactamase PDC-114                                      NF000422.2 PDC family class C beta-lactamase
    blaTEM-156         contig1       1  858      + blaTEM-156  class A beta-lactamase TEM-156                 ALLELE            286                      286                          100.00                          100.00              286 WP_061158039.1               class A beta-lactamase TEM-156                                      NF000531.2 TEM family class A beta-lactamase
    nimIJ_hmm          contig1    1001 1495      + nimIJ       NimIJ family nitroimidazole resistance protein HMM               165                       NA                              NA                              NA               NA NA                           NA                                                                  NF000262.1 NimIJ family nitroimidazole resistance protein
    vanG               contig1    5001 6047      + vanG        D-alanine--D-serine ligase VanG                EXACT             349                      349                          100.00                          100.00              349 WP_063856695.1               D-alanine--D-serine ligase VanG                                     NF000091.3 D-alanine--D-serine ligase VanG


Fields:

- Target Identifier - This is from the FASTA defline for the protein or DNA sequence
- Contig id - (optional) Contig name
- Start - (optional) 1-based coordinate of first nucleotide coding for protein in DNA sequence on contig
- Stop - (optional) 1-based corrdinate of last nucleotide coding for protein in DNA sequence on contig
- Gene symbol - Gene or gene-family symbol for protein hit
- Protein name - Full-text name for the protein
- Method - Type of hit found by AMRFinder one of five options
  - ALLELE - 100% sequence match over 100% of length to a protein named at the allele level in the AMRFinder database
  - EXACT - 100% sequence match over 100% of length to a protein in the database that is not a named allele
  - BLAST - BLAST alignment is > 90% of length and > 90% identity to a protein in the AMRFinder database
  - PARTIAL - BLAST alignment is > 50% of length, but < 90% of length and > 90% identity
  - HMM - HMM was hit above the cutoff, but there was not a BLAST hit that met standards for BLAST or PARTIAL
- Target length - The length of the query protein. The length of the blast hit for translated-DNA searches
- Reference protein length - The length of the AMR Protein in the database (NA if HMM-only hit)
- % Coverage of reference protein - % covered by blast hit (NA if HMM-only hit)
- % Identity to reference protein - % amino-acid identity to reference protein (NA if HMM-only hit)
- Alignment length - Length of BLAST alignment in amino-acids (NA if HMM-only hit)
- Accession of closest protein - RefSeq accession for protin hit by BLAST (NA if HMM-only hit)
- Name of closest protein - Full name assigned to the AMRFinder database protein (NA if HMM-only hit)
- HMM id - Accession for the HMM
- HMM description - The family name associated with the HMM

## Known Issues

Handling of fusion genes is still under active development. Currently they are
reported as two lines, one for each portion of the fusion. Gene symbol, Protein
name, Name of closest protein, HMM id, and HMM description are with respect to
the individual elements of the fusion. This behavior is subject to change.

File format checking of input files is almost nonexistant. Software behavior
with incorrect input files is not defined.

If you find bugs not listed here or have other questions/comments please email
us at pd-help@ncbi.nlm.nih.gov.

## Methods

### Protein searches (-p)

AMRFinder-prot uses the database of AMR gene sequences, hidden Markov models
(HMMs), the hierarchical tree of AMR protein designations, and a custom
rule-set to generate names and coordinates for AMR genes, along with
descriptions of the evidence used to identify the sequence. Genes are reported
with the following procedure after both HMMER and BLASTP searches are run.

#### BLASTP matches

BLASTP is run with the -task blastp-fast -word_size 6 -threshold 21 -evalue
1e-20 -comp_based_stats 0 options against the AMR gene database described
above.  Exact BLAST matches over the full length of the reference protein are
reported. If there is no exact match, then the following rules are applied:
Matches with < 90% identity or with < 50% coverage of the protein are dropped.
If the hit is to a fusion protein than at least 90% of the protein must be
covered. A BLAST match to a reference protein is removed if it is covered by
another BLAST match which has more identical residues or the same number of
identical residues, but to a longer reference protein. A single match is chosen
as the best of what remains sorting by the following criteria in order (1) if
it is exact; (2) has more identical residues; (3) hits a shorter protein; or
(4) the gene symbol comes first in alphabetical order.

#### HMM matches

HMMER version 3.1b2 (http://hmmer.org/) is run using the --cut_tc -Z 10000
options with the HMM database described above. HMM matches with full_score
< TC1 or domain_score < TC2 are dropped. All HMM matches to HMMs for parent
nodes of other HMM matches in the hierarchy are removed. The match(es) with the
highest full score are kept. If there is an exact BLAST match or the family of
the BLAST match reference protein is descendant of the family of the HMM then
the information for the nearest HMM node to the BLAST match are returned.

### Translated DNA searches (-n)

Translated alignments using BLASTX of the assembly against the AMR protein
database can be used to help identify partial, split, or unannotated AMR
proteins using the -task tblastn-fast -word_size 3 -evalue 1e-20 -seg no
-comp_based_stats 0 options. The algorithm for selecting hits is as described
above for proteins, but note that HMM searches are not performed against the
unannotated assembly.

## Help

If you have questions about AMRFinder that aren't answered in this document you
can email us at pd-help@ncbi.nlm.nih.gov

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
NCBI, NLM, and the U.S. Government disclaim all warranties as to performance,
merchantability or fitness for any particular purpose.

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using:
https://ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/ as the 
citation.
 
 
 
 
