#!/binbash

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

