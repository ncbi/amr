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
	VERSION_STRING := $(shell git describe --tags)
else
	VERSION_STRING := $(shell cat version.txt)
endif
SVNREV := -D'SVN_REV="$(VERSION_STRING)"'

# Define default paths
PREFIX ?= /usr/local
INSTALL=install

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

amrfinder.o:  common.hpp common.inc amrfinder.inc
amrfinderOBJS=amrfinder.o common.o
amrfinder:	$(amrfinderOBJS)
	$(CXX) -o $@ $(amrfinderOBJS) -pthread 

amrfinder_update.o:  common.hpp common.inc amrfinder.inc
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
	$(INSTALL) --target-directory=$(PREFIX)/bin $(BINARIES)

# amrfinder binaries for github binary release
GITHUB_FILE=amrfinder_binaries_v$(VERSION_STRING)
GITHUB_FILES = test_* $(BINARIES)
github_binaries:
	@if [ ! -e version.txt ]; \
	then \
		echo >&2 "version.txt required to make a distribution file"; \
		false; \
	fi
#   first recompile amrfinder.o to pick up the new version info
	rm amrfinder.o amrfinder
	make
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


test:
	curl -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_dna.fa \
		-O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.fa \
		-O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.gff \
		-O https://raw.githubusercontent.com/ncbi/amr/v3b/test_both.expected \
		-O https://raw.githubusercontent.com/ncbi/amr/v3b/test_dna.expected \
		-O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.expected
	amrfinder --plus -p test_prot.fa -g test_prot.gff -O Campylobacter > test_prot.got
	diff test_prot.expected test_prot.got
	amrfinder --plus -n test_dna.fa -O Campylobacter > test_dna.got
	diff test_dna.expected test_dna.got
	amrfinder --plus -n test_dna.fa -p test_prot.fa -g test_prot.gff -O Campylobacter > test_both.got
	diff test_both.got test_both.expected


     -O https://raw.githubusercontent.com/ncbi/amr/v3b/test_prot.expected
