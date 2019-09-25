#!/bin/bash

if [[ `which ssconvert` ]]; then
    if (($#>0)); then
	Name=${1%%.ods}
	Name=${Name##*/}
	echo "Model name: «${Name}»"
	ssconvert -S --export-options="quoting-mode=never separator='	'" ${1} "${Name}_%s.txt"
	for i in ${Name}_*.txt; do
	    mv "${i}" "${i%%txt}tsv";
	done    
    fi
else
    echo "gnumeric's «ssconvert» not found"
fi
