#!/bin/bash

if [[ $1 ]]; then
    sed -i "s/amr:[0-9\.]\+/amr:$1/g" *.cwl
else
    echo "Usage: $0 <version>"
fi


