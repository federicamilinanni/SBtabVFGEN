#!/bin/bash
j=0;
Document=`head -n 1 "${1}" | sed -E "s/^.*Document='([^']+)'.*\\$/\\1/"`
for i in $*; do
    TableName=`head -n 1 "${i}" | sed -E "s/^.*TableName='([^']+)'.*\\$/\\1/"`
    echo "$TableName"
    f[j]="/dev/shm/${TableName}"
    cp "${i}" "${f[j]}"
    ((j++));
done
HASH=`md5sum "$@" | md5sum`
FILE="./${Document}_`date +%Y-W%W_${HASH:0:6}`.ods"
ssconvert --merge-to="$FILE" "${f[@]}"
