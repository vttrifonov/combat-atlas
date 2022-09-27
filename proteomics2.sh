#!/bin/bash

msconvert="docker run -it -v /tmp:/data --rm proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert"
cache=$(pwd)/.cache/download/CBD-KEY-PROTEOMICS

mkdir -p .cache/proteomics2
cd .cache/proteomics2

input=xxx.d.zip.raw

unzip $input

eval $msconvert --help