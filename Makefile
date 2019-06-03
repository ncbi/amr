##############################################################################
# PUBLIC DOMAIN NOTICE This software/database is "United States Government
# Work" under the terms of the United States Copyright Act. It was written as
# part of the authors' official duties for the United States Government and
# thus cannot be copyrighted. This software/database is freely available to the
# public for use without a copyright notice. Restrictions cannot be placed on
# its present or future use.
#
# Although all reasonable efforts have been taken to ensure the accuracy and
# reliability of the software and data, the National Center for Biotechnology
# Information (NCBI) and the U.S. Government do not and cannot warrant the
# performance or results that may be obtained by using this software or data.
# NCBI, NLM, and the U.S. Government disclaim all warranties as to performance,
# merchantability or fitness for any particular purpose.
#
# In any work or product derived from this material, proper attribution of the
# authors as the source of the software or data should be made, using:
# https://ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/ as the
# citation.
###############################################################################

# the SVNREV is set automatically here for convenience,
# but when actually building we should override it like:
# make all SVNREV=-D\'SVN_REV=\"$VERSION\"\' or use
# a version.txt file
ifeq ($(wildcard version.txt),)
	VERSION_STRING := $(shell svnversion -n .)
else
	VERSION_STRING := $(shell cat version.txt)
endif
SVNREV := -D'SVN_REV="$(VERSION_STRING)"'

# Define default paths
prefix=/usr/local
bindir = $(prefix)/bin
datadir=$(prefix)/share
BLAST_BIN=$(BIN)
HMMER_BIN=$(BIN)
INSTALL_DIR = $(datadir)/amrfinder
INSTALL_PROGRAM=install

CPPFLAGS = -std=gnu++11 -pthread -malign-double -fno-math-errno -O3 $(SVNREV)

CXX=g++
COMPILE.cpp= $(CXX) $(CPPFLAGS)  -c 


.PHONY: all clean install release

BINARIES= amr_report amrfinder amrfinder_update fasta_check fasta2parts gff_check point_mut 

all:	$(BINARIES)

release: clean
	svnversion . > version.txt
	make all

common.o:	common.hpp common.inc
gff.o: gff.hpp common.hpp common.inc

amr_report.o:	common.hpp common.inc gff.hpp
amr_reportOBJS=amr_report.o common.o gff.o
amr_report:	$(amr_reportOBJS)
	$(CXX) -o $@ $(amr_reportOBJS)

amrfinder.o:  common.hpp common.inc
amrfinderOBJS=amrfinder.o common.o
amrfinder:	$(amrfinderOBJS)
	$(CXX) -o $@ $(amrfinderOBJS) -pthread 

amrfinder_update.o:  common.hpp common.inc
amrfinder_updateOBJS=amrfinder_update.o common.o
amrfinder_update:      $(amrfinder_updateOBJS)
	$(CXX) -o $@ $(amrfinder_updateOBJS) -lcurl

fasta_check.o:	common.hpp common.inc 
fasta_checkOBJS=fasta_check.o common.o 
fasta_check:	$(fasta_checkOBJS)
	$(CXX) -o $@ $(fasta_checkOBJS)

fasta2parts.o:	common.hpp common.inc
fasta2partsOBJS=fasta2parts.o common.o
fasta2parts:	$(fasta2partsOBJS)
	$(CXX) -o $@ $(fasta2partsOBJS)

gff_check.o:	common.hpp common.inc gff.hpp
gff_checkOBJS=gff_check.o common.o gff.o
gff_check:	$(gff_checkOBJS)
	$(CXX) -o $@ $(gff_checkOBJS)

point_mut.o:	common.hpp common.inc 
point_mutOBJS=point_mut.o common.o
point_mut:	$(point_mutOBJS)
	$(CXX) -o $@ $(point_mutOBJS)


clean:
	rm -f *.o
	rm -f $(BINARIES)

install:
	$(INSTALL_PROGRAM) --target-directory=$(INSTALL_DIR) $(BINARIES)

# amrfinder binaries for github binary release
GITHUB_FILE=amrfinder_binaries_v$(VERSION_STRING)
GITHUB_FILES = test_* $(BINARIES)
github_binaries:
	@if [ ! -e version.txt ]; \
	then \
		echo >&2 "version.txt required to make a distribution file"; \
		false; \
	fi
	mkdir $(GITHUB_FILE)
	echo $(VERSION_STRING) > $(GITHUB_FILE)/version.txt
	cp $(GITHUB_FILES) $(GITHUB_FILE)
	if [ -e $(GITHUB_FILE).tar.gz ]; then rm $(GITHUB_FILE).tar.gz; fi
	cd $(GITHUB_FILE); tar cvfz ../$(GITHUB_FILE).tar.gz *
#	tar cvfz $(GITHUB_FILE).tar.gz $(GITHUB_FILE)/*
	rm -r $(GITHUB_FILE)/*
	rmdir $(GITHUB_FILE)

DISTFILES=$(GITHUB_FILES) Makefile *.cpp *.hpp *.inc
dist:
	@if [ ! -e version.txt ]; \
	then \
		echo >&2 "version.txt required to make a distribution file"; \
		false; \
	fi
	mkdir $(DISTDIR)
	echo $(VERSION_STRING) > $(DISTDIR)/version.txt
	cp $(DISTFILES) $(DISTDIR)
	datalink=`basename $$(readlink -f ../releases/current)`; \
		mkdir -p $(DISTDIR)/data/$$datalink; \
		cd $(DISTDIR)/data; \
		ln -s $$datalink latest; \
		cd ../../../data; \
		$(MAKE) dist INSTALL_DIR=../src/$(DISTDIR)/data/$$datalink;
	if [ -e $(DISTDIR).tar.gz ]; then rm $(DISTDIR).tar.gz; fi
	tar cvfz $(DISTDIR).tar.gz $(DISTDIR)/*
	rm -r $(DISTDIR)/*
	rmdir $(DISTDIR)

