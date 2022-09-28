#!/bin/bash

phil=/philosopher/philosopher
cache=$(pwd)/.cache/download/CBD-KEY-PROTEOMICS

mkdir -p .cache/proteomics3
cd .cache/proteomics3

$phil workspace --init
$phil database --id UP000005640 --contam --reviewed

x=$( ls -1 $cache/*.pepXML | head -1 )

ls -1 $cache/*.pepXML | while read x; do
    y=$(basename $x)
    echo $y
    y=${y%.pepXML}
    [ ! -e  $y ] || continue
    z=$( ls -1 ../proteomics2/${y}*)

    cp $x $y.pepXML
    cp $z $y.mzML

    $phil workspace --init
    $phil database --id UP000005640 --contam --reviewed
    $phil peptideprophet --database 2022-09-28-decoys-reviewed-contam-UP000005640.fas \
        --decoy rev_ --ppm --accmass \
        --expectscore --decoyprobs --nonparam $y.pepXML
    $phil proteinprophet interact-$y.pep.xml
    $phil filter --sequential --razor --picked --tag rev_ \
        --pepxml interact-$y.pep.xml \
        --protxml interact.prot.xml    
    $phil freequant --dir .
    $phil report

    mkdir -p $y
    mv *.tsv $y/
    mv *.xml $y/    
    mv protein.fas $y/
    rm $y.pepXML $y.mzML    
done 
