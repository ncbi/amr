# WARNING: This is an experimental version that includes StxTyper

Currently it will need to be compiled from source to install.

### Downloading and compiling AMRFinderPlus with StxTyper

    git clone https://github.com/ncbi/amr.git
    cd amr
    git checkout stxtype
    git submodule update --init
    make
    ./amrfinder -u
    make test

# NCBI Antimicrobial Resistance Gene Finder (AMRFinderPlus)

This software and the accompanying database are designed to find acquired antimicrobial resistance genes and point mutations in protein and/or assembled nucleotide sequences. We have also added "plus" stress, heat, and biocide resistance and virulence factors for [some organisms](https://github.com/evolarjun/amr/wiki/Curated-organisms).

## See [the wiki for documentation](https://github.com/ncbi/amr/wiki)
## [Citing AMRFinderPlus](https://github.com/ncbi/amr/wiki#how-to-cite)
## Please [subscribe to our announce list](https://www.ncbi.nlm.nih.gov/mailman/listinfo/amrfinder-announce) for announcements of database and software updates.

----
# Licenses

## PUBLIC DOMAIN NOTICE

### National Center for Biotechnology Information

This software/database is a "United States Government Work" under the
terms of the United States Copyright Act.  It was written as part of
the authors' official duties as a United States Government employee and
thus cannot be copyrighted.  This software/database is freely available
to the public for use. The National Library of Medicine and the U.S.
Government have not placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy
and reliability of the software and data, the NLM and the U.S.
Government do not and cannot warrant the performance or results that
may be obtained by using this software or data. The NLM and the U.S.
Government disclaim all warranties, express or implied, including
warranties of performance, merchantability or fitness for any particular
purpose.

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using the
following citation:

Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH,
Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W.
AMRFinderPlus and the Reference Gene Catalog facilitate examination of the
genomic links among antimicrobial resistance, stress response, and virulence.
Sci Rep. 2021 Jun 16;11(1):12728. doi: [10.1038/s41598-021-91456-0](https://doi.org/10.1038/s41598-021-91456-0). PMID: [34135355](https://pubmed.ncbi.nlm.nih.gov/34135355/); PMCID: [PMC8208984](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8208984/).


