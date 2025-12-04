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

INSTALL=install

# make it possible to hard define a database directory
# Define default paths
# This is a little convoluted because I broke things and don't want
# to change two different ways of defining the paths. This could
# be simplified in a later release
PREFIX ?= /usr/local
ifneq '$(INSTALL_DIR)' ''
	bindir=$(INSTALL_DIR)
endif
bindir ?= $(PREFIX)/bin
ifneq '$(CONDA_DB_DIR)' ''
	DBDIR := -D'CONDA_DB_DIR="$(CONDA_DB_DIR)"'
endif
ifneq '$(DEFAULT_DB_DIR)' ''
	DBDIR := -D'CONDA_DB_DIR="$(DEFAULT_DB_DIR)"'
endif

# for testing database updates using 
ifdef TEST_UPDATE
	TEST_UPDATE_DB := '-D TEST_UPDATE'
endif

# detect system architecture and set appropriate flags
# this is probably not the best way (i.e. M1 Mac would be arm64)
# but it works for Nvidia Jetson boards (aarch64) 
ARCH := $(shell uname -m)
OS := $(shell uname -s)
# "hack": if amd64 we can set to aarch64
# as AArch64 and ARM64 refer to the same thing
# this should build for Mac M1 and other arm64 chips
ifeq ($(ARCH),arm64)
  ARCH := aarch64
endif
# report detected OS and arch in stdout
$(info Dectected architecture: $(OS) $(ARCH))
# set CFLAGS based on arch
ifeq ($(ARCH),aarch64)
  # set arm CFLAGS
  CPPFLAGS = -std=gnu++17 -pthread --signed-char -falign-jumps -fno-math-errno -O3 
else
  # set x86_x64 CFLAGS
  CPPFLAGS = -std=gnu++17 -pthread -malign-double -fno-math-errno -O3
endif
# was: -std=gnu++14 

CXX=g++
COMPILE.cpp= $(CXX) $(CPPFLAGS) $(SVNREV) $(DBDIR) $(TEST_UPDATE_DB) -c 


.PHONY: all clean install release stxtyper test

BINARIES= amr_report amrfinder amrfinder_index amrfinder_update fasta_check \
		  fasta_extract fasta2parts gff_check dna_mutation mutate disruption2genesymbol

all:	$(BINARIES) stxtyper

release: clean
	svnversion . > version.txt
	make all

common.o:	common.hpp common.inc
curl_easy.o: curl_easy.hpp common.hpp common.inc
gff.o: gff.hpp common.hpp common.inc
alignment.o:	alignment.hpp seq.hpp common.hpp common.inc
seq.o: seq.hpp graph.hpp common.hpp common.inc

amr_report.o:	common.hpp common.inc gff.hpp alignment.hpp tsv.hpp seq.hpp columns.hpp version.txt
amr_reportOBJS=amr_report.o common.o gff.o alignment.o seq.o graph.o
amr_report:	$(amr_reportOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(amr_reportOBJS)

amrfinder.o:  common.hpp common.inc gff.hpp seq.hpp tsv.hpp columns.hpp version.txt
amrfinderOBJS=amrfinder.o common.o gff.o tsv.o
amrfinder:	$(amrfinderOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(amrfinderOBJS) -pthread $(DBDIR)

amrfinder_update.o:  common.hpp common.inc curl_easy.hpp version.txt
amrfinder_updateOBJS=amrfinder_update.o common.o curl_easy.o
amrfinder_update:      $(amrfinder_updateOBJS) 
	@if [ "$(TEST_UPDATE)" != "" ]  ; \
	then  \
		touch amrfinder_update.cpp ;\
	fi # make sure the next make command rebuilds amrfinder_update
	$(CXX) $(LDFLAGS) -o $@ $(amrfinder_updateOBJS) -lcurl 

amrfinder_index.o:  common.hpp common.inc version.txt
amrfinder_indexOBJS=amrfinder_index.o common.o
amrfinder_index:      $(amrfinder_indexOBJS) 
	$(CXX) $(LDFLAGS) -o $@ $(amrfinder_indexOBJS) 

fasta_check.o:	common.hpp common.inc version.txt
fasta_checkOBJS=fasta_check.o common.o 
fasta_check:	$(fasta_checkOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(fasta_checkOBJS)

fasta_extract.o:	common.hpp common.inc version.txt
fasta_extractOBJS=fasta_extract.o common.o
fasta_extract:	$(fasta_extractOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(fasta_extractOBJS)

fasta2parts.o:	common.hpp common.inc version.txt
fasta2partsOBJS=fasta2parts.o common.o
fasta2parts:	$(fasta2partsOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(fasta2partsOBJS)

gff_check.o:	common.hpp common.inc gff.hpp version.txt
gff_checkOBJS=gff_check.o common.o gff.o
gff_check:	$(gff_checkOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(gff_checkOBJS)

dna_mutation.o:	common.hpp common.inc alignment.hpp seq.hpp tsv.hpp columns.hpp version.txt
dna_mutationOBJS=dna_mutation.o common.o alignment.o seq.o graph.o
dna_mutation:	$(dna_mutationOBJS)
	$(CXX) $(LDFLAGS) -o $@ $(dna_mutationOBJS)

mutate.o:	common.hpp common.inc alignment.hpp seq.hpp version.txt
mutateOBJS=mutate.o common.o alignment.o seq.o graph.o
mutate:	$(mutateOBJS)
	$(CXX) -o $@ $(mutateOBJS)

disruption2genesymbol.o:	common.hpp common.inc seq.hpp version.txt
disruption2genesymbolOBJS=disruption2genesymbol.o common.o alignment.o seq.o graph.o
disruption2genesymbol:	$(disruption2genesymbolOBJS)
	$(CXX) -o $@ $(disruption2genesymbolOBJS)

stxtyper:
		$(MAKE) -C stx

clean:
	rm -f *.o
	rm -f $(BINARIES)
	$(MAKE) -C stx clean

install:
	@if [ ! -e $(DESTDIR)$(bindir) ]; \
	then \
		mkdir -p $(DESTDIR)$(bindir); \
	fi
	$(INSTALL) $(BINARIES) $(DESTDIR)$(bindir)
	make -C stx install PREFIX=$(PREFIX) bindir=$(bindir)
	mkdir $(DESTDIR)$(bindir)/stx
	ln -s ../stxtyper $(DESTDIR)$(bindir)/stx/stxtyper

# amrfinder binaries for github binary release
GITHUB_FILE=amrfinder_binaries_v$(VERSION_STRING)
GITHUB_FILES = test_amrfinder.sh test_*.expected test_*.fa test_*.gff $(BINARIES)

github_binaries:
	@if [ ! -e version.txt ]; \
	then \
		echo >&2 "version.txt required to make a distribution file"; \
		false; \
	fi
#   first recompile amrfinder.o to pick up the new version info
#	and remove leaky NCBI paths
	make clean
#	make all CXX=/usr/bin/g++ LD_RUN_PATH= 
	make all LD_RUN_PATH= 
	mkdir $(GITHUB_FILE)
	echo $(VERSION_STRING) > $(GITHUB_FILE)/version.txt
	cp $(GITHUB_FILES) $(GITHUB_FILE)
#	make -C stx
#	make -C stx install INSTALL_DIR=../$(GITHUB_FILE)/stx CXX=/usr/bin/g++ LD_RUN_PATH=
	make -C stx install INSTALL_DIR=../$(GITHUB_FILE)/stx LD_RUN_PATH=
	cp stx/test_stxtyper.sh stx/version.txt $(GITHUB_FILE)/stx
	mkdir $(GITHUB_FILE)/stx/test
	cp -R stx/test/*.fa stx/test/*.expected $(GITHUB_FILE)/stx/test
	if [ -e $(GITHUB_FILE).tar.gz ]; then rm $(GITHUB_FILE).tar.gz; fi
	cd $(GITHUB_FILE); ln -s stx/stxtyper .; tar cvfz ../$(GITHUB_FILE).tar.gz *
	rm -r $(GITHUB_FILE)/*
	rmdir $(GITHUB_FILE)

test: $(DISTFILES) Makefile *.cpp *.hpp *.inc test_dna.fa test_prot.fa test_prot.gff test_dna.fa test_dna.expected test_prot.expected test_both.expected
	make -C stx test
	# test the amrfinder in the current directory 
	# with the data in the current directory
	./test_amrfinder.sh -n 
