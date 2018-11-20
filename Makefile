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

SHELL=/bin/sh
INSTALL = install

CPPFLAGS = -std=gnu++11 \
  -malign-double -fno-math-errno \
  -Wall -Wextra \
  -Wcast-align -Wconversion -Wdeprecated-declarations -Wformat -Winit-self -Wlogical-op \
  -Wmissing-braces -Wmissing-field-initializers -Wmissing-format-attribute -Wmissing-include-dirs \
  -Woverloaded-virtual -pedantic -Wparentheses -Wpointer-arith -Wsequence-point -Wshadow -Wunused \
  -Wsuggest-attribute=format -Wswitch -Wuninitialized -Wsign-conversion -Wuseless-cast \
  -O3 \
  $(SVNREV)
#  -Wno-error=misleading-indentation -Wno-nonnull-compare \

CXX=g++
COMPILE.cpp= $(CXX) $(CPPFLAGS)  -c

#prefix = /panfs/pan1.be-md.ncbi.nlm.nih.gov/bacterial_pathogens/backup
prefix=/usr/local

# soft links are created here to INSTALL_DIR
bindir="$(prefix)/bin"

datadir=$(prefix)/share

# This is where AMRFinder and the default database will be installed
INSTALL_DIR=$(datadir)/amrfinder


.PHONY: all clean install dist
DISTFILES = Makefile *.cpp *.hpp *.inc test_* amrfinder.pl AMRFinder-dna.sh AMRFinder-prot.sh

all:	amr_report fasta_check gff_check

release: clean
	svnversion . > version.txt
	make all

common.o:	common.hpp common.inc
gff.o: gff.hpp common.hpp common.inc

amr_report.o:	common.hpp common.inc gff.hpp
amr_reportOBJS=amr_report.o common.o gff.o
amr_report:	$(amr_reportOBJS)
	$(CXX) -o $@ $(amr_reportOBJS)

fasta_check.o:	common.hpp common.inc
fasta_checkOBJS=fasta_check.o common.o
fasta_check:	$(fasta_checkOBJS)
	$(CXX) -o $@ $(fasta_checkOBJS)

gff_check.o:	common.hpp common.inc gff.hpp
gff_checkOBJS=gff_check.o common.o gff.o
gff_check:	$(gff_checkOBJS)
	$(CXX) -o $@ $(gff_checkOBJS)

clean:
	rm -f *.o
	rm -f amr_report fasta_check gff_check
	#rm version.txt

install:
	$(INSTALL) -D -t $(INSTALL_DIR) amr_report fasta_check gff_check amrfinder.pl
#	@dest=$(INSTALL_DIR); \
#		if [ ! -e $(bindir)/amrfinder.pl ]; \
#		then \
#			ln -s "$$dest/amrfinder.pl" $(bindir); \
#			ln -s "$$dest/amr_report" $(bindir); \
#			ln -s "$$dest/fasta_check" $(bindir); \
#			ln -s "$$dest/gff_check" $(bindir); \
#		else \
#			echo "$$dest/amrfinder.pl already exists, so skipping link creation"; \
#		fi

DISTDIR=amrfinder.$(VERSION_STRING)
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
#	rm -r $(DISTDIR)/*
#	rmdir $(DISTDIR)
