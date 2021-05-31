#!/bin/bash
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic
## assignation.

## Variables
NCBI="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
TAXDUMP="taxdump.tar.gz"
TAXIDPATH="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
TAXID="nucl_gb.accession2taxid.gz"
NAMES="names.dmp"
NODES="nodes.dmp"
DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
USELESS_FILES="${TAXDUMP} ${DMP} gc.prt readme.txt"

## Download taxdump
rm -rf ${USELESS_FILES} "${NODES}" "${NAMES}"
wget "${NCBI}${TAXDUMP}" && \
    tar zxvf "${TAXDUMP}" && \
    rm -rf ${USELESS_FILES}

## Limit search space to scientific names
grep "scientific name" "${NAMES}" > "${NAMES/.dmp/_reduced.dmp}" && \
    rm -f "${NAMES}" && \
    mv "${NAMES/.dmp/_reduced.dmp}" "${NAMES}"

## Download gi_taxid_nucl
rm -f "${TAXID/.gz/}*"
wget "${TAXIDPATH}${TAXID}" && \
    gunzip "${TAXID}"

exit 0
