#!/bin/bash

AMRFINDER_OPTS=" --plus --print_node --threads 6 "
path=0
no_download=0
print_help=0
while getopts "pnh" opt; do
    case $opt in
        p) path=1 ;;
        n) no_download=1 ;;
        h) print_help=1 ;;
    esac
done

if [ "$print_help" -gt 0 ]
then
    echo "test_amrfinder.sh - Run tests"
    echo "Options: "
    echo "    -p Test the amrfinder command in path instead of ./amrfinder"
    echo "    -n Don't attempt to download fresh test data, use the test data in $PWD"
    echo "    -h print this help message"
    exit 1
fi

# some color macros
if [ "$TERM" == "" ] || [ "$TERM" == "dumb" ] || [ ! -t 1 ]
then
    green='' # no colors
    red=''
    reset=''
else
    green=`tput setaf 2`  # Set green foreground color (code 2)
    red=`tput setaf 1`    # Set red foreground color (code 1)
    reset=`tput sgr0`     # Reset color to default
fi

if [ "$path" -gt 0 ]
then
    echo "Testing amrfinder command in your \$PATH"
    which amrfinder
    AMRFINDER=amrfinder
else
    echo "Testing ./amrfinder"
    AMRFINDER=./amrfinder
fi

if [ "$no_download" -gt 0 ]
then
    echo "-n option detected, skipping download of test data"
else
    echo Downloading fresh test data...
    curl -k -s -f \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_dna.fa \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.fa \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.gff \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_both.expected \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_dna.expected \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.expected

    if [ $? != 0 ]
    then
        echo "WARNING: Could not download new test data."
        echo "Will attempt to use test data in $PWD."
        echo "Test data included with installation not match the latest database release."
    fi
fi

TESTS=0
TEST_TEXT=""
FAILURES=0

function test_input_file {
    local test_base="$1"
    local options="$2"

    TESTS=$(( $TESTS + 1 ))

    if ! $AMRFINDER $options $AMRFINDER_OPTS > "$test_base.got"
    then
        echo "${red}not ok: $AMRFINDER returned a non-zero exit value indicating a failure of the software${reset}"
        echo "#   $AMRFINDER $options $AMRFINDER_OPTS > $test_base.got"
        return 1
    else
        if ! diff -q "$test_base.expected" "$test_base.got"
        then
            echo "${red}not ok: $AMRFINDER returned output different from expected.${reset}"
            echo "#   diff $test_base.expected $test_base.got"
            echo "#   To approve run: "
            echo "#         mv $test_base.got $test_base.expected"
            TEST_TEXT="$TEST_TEXT"$'\n'"${red}Failed $test_base${reset}";
            echo ""
            return 1
        else
            echo "${green}ok:${reset} $test_base"
            return 0
        fi
    fi
}

test_input_file "test_prot" "-p test_prot.fa -g test_prot.gff -O Escherichia"
FAILURES=$(( $? + $FAILURES ))

test_input_file "test_dna" "-n test_dna.fa -O Escherichia --mutation_all test_dna_mut_all.got"
FAILURES=$(( $? + $FAILURES ))

test_input_file "test_both" "-n test_dna.fa -g test_prot.gff -p test_prot.fa -O Escherichia"
FAILURES=$(( $? + $FAILURES ))

echo "Done."
echo "$TEST_TEXT"
echo ""
if [ "$FAILURES" -gt 0 ]
then
    PASSED=$(( $TESTS - $FAILURES ))
    echo "${red}not ok overall: $FAILURES out of $TESTS amrfinder tests failed${reset}"
    exit 1
else
    echo "${green}ok: all $TESTS amrfinder tests passed ${reset}"
    echo "Success!"
fi
