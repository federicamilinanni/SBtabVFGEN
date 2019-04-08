#!/bin/bash
if (($#>0)); then
    echo  "ssconvert -S -T Gnumeric_stf:stf_csv ${1}.ods CaMKII_%s.tsv"
    ssconvert -S -T "Gnumeric_stf:stf_csv" ${1} 'CaMKII_%s.tsv'
    #ssconvert -S -T Gnumeric_stf:stf_assistant --export-options="quoting-mode=never separator=	" ${1} 'CaMKII_%s.tsv'
    
    for i in CaMKII_*.tsv; do
     	sed -i -E 's/,/\t/g' "${i}"
    done
fi
