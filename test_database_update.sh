#!/bin/sh

touch amrfinder_update.cpp
make TEST_UPDATE=1
./amrfinder -U
./test_amrfinder.sh
touch amrfinder_update.cpp
make

