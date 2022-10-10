#!/bin/bash

cache=$(pwd)/.cache/download/CBD-KEY-PROTEOMICS

mkdir -p .cache/proteomics4
cd .cache/proteomics4

ls -1 $cache/*.mgf | head -1 | while read input; do
    name=$(basename $input)
    name=${name%.mgf}
    echo $name
    [ ! -e ${name}.mzML ] || continue

    wine64_anyuser msconvert $input --outfile ./${name}.mzML
done

