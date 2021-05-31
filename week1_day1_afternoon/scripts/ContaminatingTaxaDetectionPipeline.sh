#!/bin/sh

#  ContaminatingTaxaDetectionPipeline.sh
#  
#
#  Created by Torsten Hugo Struck on 09/02/2017.
#

# USAGE: ContaminatingTaxaDetectionPipeline.sh [Assembly file] [Barcode sequence] [e value of blast searches] [number of hits saved] [minimum percentage of identity]

# Set variables with the relevant paths
ASSEMBLY=$1
BARCODE=$2
EVALUE=$3
MAXHITS=$4
MINPERIDENT=$5

#Generate blast-searchable database
makeblastdb -in ${ASSEMBLY} -parse_seqids -dbtype nucl

#Blast the barcode sequence(s) against the assembly library
blastn -evalue ${EVALUE} -outfmt "6 qseqid sseqid pident score evalue" -out ${ASSEMBLY}_library_RNA_Contamination_BLASThits.txt -db ${ASSEMBLY} -query ${BARCODE}


#Retrieve the corresponding sequences to a hit and write it out to new fasta files
cut -f2 < ${ASSEMBLY}_library_RNA_Contamination_BLASThits.txt | while read ID
do
echo ${ID}
sed -n "/^>${ID}$/,/>/p" < ${ASSEMBLY} | sed '$d' >> ${ASSEMBLY}_Sequences_of_best_hit.fas
done

#Blast the retrieved sequences against NCBI nr database remotely (version 2.6 is needed for that)
echo "Blastn search started with the settings: -max_target_seqs ${MAXHITS} -evalue ${EVALUE} -perc_identity ${MINPERIDENT} -outfmt "6 sacc qseqid sseqid length pident score evalue" -remote -out ${ASSEMBLY}_nr_RNA_BLAST_Matches.txt -db nr -query ${ASSEMBLY}_Sequences_of_best_hit.fas"
blastn -max_target_seqs ${MAXHITS} -evalue ${EVALUE} -perc_identity ${MINPERIDENT} -outfmt "6 sacc qseqid sseqid length pident score evalue" -remote -out ${ASSEMBLY}_nr_RNA_BLAST_Matches.txt -db nr -query ${ASSEMBLY}_Sequences_of_best_hit.fas

#Determine the taxonomic position of the hit
echo "Retrieving taxa from blast results"
sh RetrieveTaxonomyDatabase.sh
cut -f1 < ${ASSEMBLY}_nr_RNA_BLAST_Matches.txt | sort | uniq | while read GI
do
sh RetrieveTaxonPath.sh ${GI} >> ${ASSEMBLY}_Taxa_found.txt
done

mkdir Results
mv ${ASSEMBLY}_library_RNA_Contamination_BLASThits.txt ${ASSEMBLY}_Sequences_of_best_hit.fas ${ASSEMBLY}_nr_RNA_BLAST_Matches.txt ${ASSEMBLY}_Taxa_found.txt Results/
rm ${ASSEMBLY}.n* *.dmp nucl_gb.accession2taxid 
