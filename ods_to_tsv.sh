#!/bin/bash

if (($#>0)); then
    Name=${1%%.ods}
    ssconvert -S --export-options="quoting-mode=never separator='	'" ${1} "${Name}_%s.txt"
    for i in ${Name}_*.txt; do
	mv "${i}" "${i%%txt}tsv";
    done    
fi
