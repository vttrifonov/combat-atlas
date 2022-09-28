#!/bin/bash

cache=$(pwd)/.cache/download/CBD-KEY-PROTEOMICS
remote=vtrifonov@witsc02dw380md6.wks.jnj.com
ident=~/.ssh/jnj_domino_laptop_ed25519
msconvert="/usr/local/bin/docker run -v /var/run/docker.sock:/var/run/docker.sock docker docker run -v /tmp:/data --rm proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert /data/input.mgf"
ssh="ssh -n -i $ident $remote"
rsync='rsync -avz -e "ssh -i '$ident'"'

mkdir -p .cache/proteomics2
cd .cache/proteomics2

ls -1 $cache/*.mgf | while read input; do
    name=$(basename $input)
    name=${name%.mgf}.mzML
    echo $name
    [ ! -e $name ] || continue

    eval $rsync $input $remote:/tmp/input.mgf
    eval $ssh $msconvert
    eval $rsync $remote:/tmp/input.mzML $name
    eval $ssh rm /tmp/input.*
done

