#!/bin/bash

printf "`date`\t$1"
cd $1
a=`ls -1`
k=".. OK"
for f in ${a[@]}
do
    cd "$f"
    if [ ! -z `md5sum -c --status *_checksum.txt` ]; then 
	printf "\nfail\t$f";
	k=""
    fi
    cd ..
done
printf "$k\n"
