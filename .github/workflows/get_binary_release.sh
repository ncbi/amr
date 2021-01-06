#!/bin/bash
# found and lightly modified these functions
get_latest_release() {
    curl --silent "https://api.github.com/repos/$1/releases/latest" |
    grep '"tag_name":' |
    cut -d '"' -f 4
}

get_tarball_url() {
    curl --silent "https://api.github.com/repos/$1/releases/latest" |
        fgrep '"browser_download_url":' |
        cut -d '"' -f 4
}


release=$(get_latest_release ncbi/amr)
URL=$(get_tarball_url ncbi/amr)

>&2 echo "Downloading AMRFinderPlus version $release"
>&2 echo "Binaries from $URL"

# download and unpack AMRFinder binaries
    curl --silent -L -O $URL
    tarball_name=$(echo $URL | perl -pe 's#^.*/(.*)#\1#')
    tar xfz $tarball_name
    rm $tarball_name

# download and unpack test
    curl --silent \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_dna.fa \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.fa \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.gff \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_both.expected \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_dna.expected \
         -O https://raw.githubusercontent.com/ncbi/amr/master/test_prot.expected

# download database
    ./amrfinder --update
