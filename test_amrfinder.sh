#!/bin/bash

if [ "$1" == "path" ]
then
    echo "Testing amrfinder command in your \$PATH"
    which amrfinder
    AMRFINDER=amrfinder
else
    echo "Testing ./amrfinder"
    AMRFINDER=./amrfinder
fi

echo Downloading fresh test data...
# curl -s -f \
     # -O https://raw.githubusercontent.com/ncbi/amr/master/test_dna.fa \
     # -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.fa \
     # -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.gff \
     # -O https://raw.githubusercontent.com/ncbi/amr/master/test_both.expected \
     # -O https://raw.githubusercontent.com/ncbi/amr/master/test_dna.expected \
     # -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.expected

if [ $? != 0 ]
then
    echo "WARNING: Could not download new test data, test data included with installation may"
    echo "not match the latest database release"
fi

if ! $AMRFINDER --plus -p test_prot.fa -g test_prot.gff -O Escherichia --print_node > test_prot.got \
    || ! diff -q test_prot.expected test_prot.got
then
    echo "Test failed: "
    echo "  $AMRFINDER --plus -p test_prot.fa -g test_prot.gff -O Escherichia --print_node > test_prot.got"
    echo "  diff test_prot.expected test_prot.got "
    diff test_prot.expected test_prot.got
    exit 1
fi

if ! $AMRFINDER --plus -n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got --print_node > test_dna.got \
    || ! diff -q test_dna.expected test_dna.got
then 
    echo "Test failed: "
    echo "  $AMRFINDER --plus -n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got --print_node > test_dna.got"
    echo "  diff test_dna.expected test_dna.got"
    diff test_dna.expected test_dna.got
    exit 1
fi
#if ! diff -q test_dna_mut_all.expected test_dna_mut_all.got
#then
#        echo "Test failed: "
#    echo "  $AMRFINDER --plus -n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got > test_dna.got"
#    echo "  diff test_dna_mut_all.expected test_dna_mut_all.got"
#    diff test_dna_mut_all.expected test_dna_mut_all.got
#    exit 1
#fi

if ! $AMRFINDER --plus -n test_dna.fa -p test_prot.fa -g test_prot.gff -O Escherichia --print_node > test_both.got \
    || ! diff -q test_both.expected test_both.got
then
    echo "Test failed: "
    echo "  $AMRFINDER --plus -n test_dna.fa -p test_prot.fa -g test_prot.gff -O Escherichia --print_node > test_both.got" 
    echo "  diff test_both.expected test_both.got "
    diff test_both.expected test_both.got
    exit 1
fi

echo "Success!"
