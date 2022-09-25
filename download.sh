#!/bin/bash

url=https://zenodo.org/record/6120249/files
cache=~/.cache/combat-atlas/download
mkdir -p $cache

function download() {
    file=$1
    [ -e $cache/$file ] || curl $url/$file?download=1 -o $cache/$file
}

download CBD-KEY-CLINVAR.tar.gz

download COMBAT-CITESeq-DATA.h5ad
download COMBAT-CITESeq-EXPRESSION-ATLAS.h5ad

download CBD-KEY-RNASEQ-WB.tar.gz

download CBD-KEY-FACS.tar.gz

download CBD-KEY-LUMINEX.tar.gz

download CBD-KEY-CYTOF-MYELOID.tar.gz
download CBD-KEY-CYTOF-WB-D.tar.gz
download CBD-KEY-CYTOF-WB.tar.gz

ftp=https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD023175/
proteomics=$cache/CBD-KEY-PROTEOMICS
mkdir -p $proteomics
function proteomics_download() {
    file=$1
    [ -e $proteomics/$file ] || ( curl $ftp/$file -o $proteomics/tmp.$file && mv $proteomics/tmp.$file $proteomics/$file )
}

proteomics_download checksum.txt
proteomics_download README.txt
proteomics_download sample_key.xlsx

cat $proteomics/README.txt | grep pepXML | tr '\t' '#' | cut -f2 -d'#' | while read x; do
    echo $x
    proteomics_download $x
done

cat $proteomics/README.txt | grep raw | tr '\t' '#' | cut -f2 -d'#' | while read x; do
    echo $x
    proteomics_download $x
done
