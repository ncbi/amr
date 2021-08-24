#!/bin/sh

cat <<END
Note this requires a ~/.github-token file that contains a github token
with appropriate permissions

Actions that can be triggered:

1. install-test
2. linux-compile-test
3. linux-binary-test
4. linux-bioconda-test
5. mac-bioconda-test
6. docker-test
END

echo -n "% "
read action

if [ "X$action" == "X" ]
then
    echo "Enter an action to take next time you run"
    exit 1;
fi

echo $action
case $action in 
    1|install-test) event=install-test;;
    2|linux-compile-test) event=linux-compile-test;;
    3|linux-binary-test) event=linux-binary-test;;
    4|linux-bioconda-test) event=linux-bioconda-test;;
    5|mac-bioconda-test) event=mac-bioconda-test;;
    6|docker-test) event=docker-test;;
    *) echo "Unknown action. Please try again"
        exit 1;;
esac

echo "action=$action event=$event"

cmd="
curl \
  -H 'Accept: application/vnd.github.everest-preview+json' \
  -H 'Authorization: token `cat ~/.github-token`' \
  --data '{\"event_type\": \"$event\"}' \
  --request POST \
  https://api.github.com/repos/ncbi/amr/dispatches
"
echo $cmd
echo $cmd | bash

