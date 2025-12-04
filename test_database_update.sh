#!/bin/sh

echo "test_database_update.sh - test update from staging area"
git status -uno
echo ""
echo "Attempts to update from ftp://ftp.ncbi.nlm.nih.gov/pathogen/Technical/AMRFinder_technical/test_database/"
echo "WARNING: recompiles AMRFinderPlus to use a different update URL"
echo "To continue press <CR>, to abort ^C"
read

set -x

touch amrfinder_update.cpp
make TEST_UPDATE=1
./amrfinder -U
./test_amrfinder.sh -n
#touch amrfinder_update.cpp
#make

