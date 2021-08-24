#!/bin/bash -x

if [ "$1" == "path" ]
then
    echo "Testing amrfinder command in your \$PATH"
    which amrfinder
    AMRFINDER=amrfinder
else
    echo "Testing ./amrfinder"
    AMRFINDER=./amrfinder
fi

if ! $AMRFINDER --plus -p test_prot.fa -g test_prot.gff -O Escherichia > test_prot.got \
    || ! diff -q test_prot.expected test_prot.got
then
    echo "Test failed: "
    echo "  $AMRFINDER --plus -p test_prot.fa -g test_prot.gff -O Escherichia > test_prot.got"
    echo "  diff test_prot.expected test_prot.got "
    exit 1
fi

if ! $AMRFINDER --plus -n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got > test_dna.got \
    || ! diff -q test_dna.expected test_dna.got
then 
    echo "Test failed: "
    echo "  $AMRFINDER --plus -n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got > test_dna.got"
    echo "  diff test_dna.expected test_dna.got"
    exit 1
fi
if ! diff -q test_dna_mut_all.expected test_dna_mut_all.got
then
	echo "Test failed: "
    echo "  $AMRFINDER --plus -n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got > test_dna.got"
    echo "  diff test_dna_mut_all.expected test_dna_mut_all.got"
    exit 1
fi

if ! $AMRFINDER --plus -n test_dna.fa -p test_prot.fa -g test_prot.gff -O Escherichia > test_both.got \
    || ! diff -q test_both.expected test_both.got
then
    echo "Test failed: "
    echo "  $AMRFINDER --plus -n test_dna.fa -p test_prot.fa -g test_prot.gff -O Escherichia > test_both.got" 
    echo "  diff test_both.expected test_both.got "
    exit 1
fi

echo "Success!"
