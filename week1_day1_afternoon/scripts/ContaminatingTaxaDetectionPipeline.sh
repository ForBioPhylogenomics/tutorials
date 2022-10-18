#!/bin/sh

#  ContaminatingTaxaDetectionPipeline.sh
#  
#
#  Created by Torsten Hugo Struck on 09/02/2017 and modified on 17/10/2022.
#

# USAGE: ContaminatingTaxaDetectionPipeline.sh [Assembly file] [Barcode sequence] [e value of blast searches] [number of hits saved] [minimum percentage of identity]

# Set variables with the relevant paths
ASSEMBLY=$1
BARCODE=$2
EVALUE=$3
MAXHITS=$4
MINPERIDENT=$5

function RetrieveTaxonomyDatabase () {
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

## Download NCBI's taxonomic data and GI (GenBank ID) taxonomic
## assignation.

## Variables
local NCBI="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/"
local TAXDUMP="taxdump.tar.gz"
local TAXIDPATH="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
local TAXID="nucl_gb.accession2taxid.gz"
local NAMES="names.dmp"
local NODES="nodes.dmp"
local DMP=$(echo {citations,division,gencode,merged,delnodes}.dmp)
local USELESS_FILES="${TAXDUMP} ${DMP} gc.prt readme.txt"

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
}

function RetrieveTaxonPath () {
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

local NAMES="names.dmp"
local NODES="nodes.dmp"
local GI_TO_TAXID="nucl_gb.accession2taxid"
local TAXONOMY=""
local GI="${1}"

# Obtain the name corresponding to a taxid or the taxid of the parent taxa
get_name_or_taxid()
{
    grep --max-count=1 "^${1}"$'\t' "${2}" | cut -f"${3}"
}

# Get the taxid corresponding to the GI number
TAXID=$(get_name_or_taxid "${GI}" "${GI_TO_TAXID}" "3")

# Loop until you reach the root of the taxonomy (i.e. taxid = 1)
while [[ "${TAXID}" -gt 1 ]] ; do
    # Obtain the scientific name corresponding to a taxid
    NAME=$(get_name_or_taxid "${TAXID}" "${NAMES}" "3")
    # Obtain the parent taxa taxid
    PARENT=$(get_name_or_taxid "${TAXID}" "${NODES}" "3")
    # Build the taxonomy path
    TAXONOMY="${NAME};${TAXONOMY}"
    TAXID="${PARENT}"
done

echo -e "${GI}\t${TAXONOMY}"
}

function SpeciesTaxaDetectionPipeline () {
# Set variables with the relevant paths
local ASSEMBLY=$1

#Generate blast-searchable database
makeblastdb -in ${ASSEMBLY} -parse_seqids -dbtype nucl

#Blast the barcode sequence(s) against the assembly library
blastn -evalue ${EVALUE} -outfmt "6 qseqid sseqid pident score evalue" -out ${ASSEMBLY}_Species_BLASThits.txt -db ${ASSEMBLY} -query ${BARCODE}

#Retrieve the corresponding sequences to a hit and write it out to new fasta files
cut -f2 < ${ASSEMBLY}_Species_BLASThits.txt | sort | uniq | while read ID
do
echo ${ID}
sed -n "/^>${ID}/,/>/p" < ${ASSEMBLY} | sed '$d' >> ${ASSEMBLY}_Sequences_of_best_hit.fas
done

#Blast the retrieved sequences against NCBI nr database remotely (version 2.6 is needed for that)
echo "Blastn search started with the settings: -max_target_seqs ${MAXHITS} -evalue ${EVALUE} -perc_identity ${MINPERIDENT} -outfmt "6 sacc qseqid sseqid length pident score evalue" -out ${ASSEMBLY}_nr_BLAST_Matches.txt -db /cluster/shared/databases/blast/2022-02-21/nt -query ${ASSEMBLY}_Sequences_of_best_hit.fas"
blastn -max_target_seqs ${MAXHITS} -evalue ${EVALUE} -perc_identity ${MINPERIDENT} -outfmt "6 sacc qseqid sseqid length pident score evalue" -out ${ASSEMBLY}_nr_BLAST_Matches.txt -db /cluster/shared/databases/blast/2022-02-21/nt -query ${ASSEMBLY}_Sequences_of_best_hit.fas

#Determine the taxonomic position of the hit
echo "Retrieving taxa from blast results"
cut -f1 < ${ASSEMBLY}_nr_BLAST_Matches.txt | sort | uniq | while read GI
do
RetrieveTaxonPath ${GI} >> ${ASSEMBLY}_Taxa_found.txt
done
}

# run the species detection script retrieving all contigs and unassembled reads matching the provided query sequences as well as the corresponding taxon name
RetrieveTaxonomyDatabase
SpeciesTaxaDetectionPipeline ${ASSEMBLY}

mkdir Results
mv ${ASSEMBLY}_Species_BLASThits.txt ${ASSEMBLY}_Sequences_of_best_hit.fas ${ASSEMBLY}_nr_BLAST_Matches.txt ${ASSEMBLY}_Taxa_found.txt Results/
rm ${ASSEMBLY}.n* *.dmp nucl_gb.accession2taxid 
