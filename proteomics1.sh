#!/bin/bash

phil=/philosopher/philosopher
cache=$(pwd)/.cache/download/CBD-KEY-PROTEOMICS

mkdir -p .cache/proteomics1
cd .cache/proteomics1

$phil workspace --init
$phil database --id UP000005640 --contam --reviewed

ls -1 $cache/*.pepXML | while read x; do
    y=$(basename $x)
    echo $y
    [ ! -e  $y ] || continue

    cp $x input.pepXML
    $phil peptideprophet --database 2022-09-22-decoys-reviewed-contam-UP000005640.fas --decoy rev_ --ppm --accmass \
        --expectscore --decoyprobs --nonparam input.pepXML
    $phil proteinprophet interact-input.pep.xml
    $phil filter --sequential --razor --picked --tag rev_ --pepxml interact-input.pep.xml --protxml interact.prot.xml
    $phil report

    mkdir -p $y
    mv *.tsv $y/
    mv *.xml $y/
    mv protein.fas $y/
    rm input.pepXML
done 
