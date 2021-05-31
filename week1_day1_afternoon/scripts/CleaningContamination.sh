#!/bin/sh

#  CleaningContamination.sh
#  
#
#  Created by Torsten Hugo Struck on 09/02/2017.
#

# USAGE: CleaningContamination.sh [Assembly file] [e value]

# Set variables with the relevant paths
ASSEMBLY=$1
EVALUE=$2

for POSFILE in pos_*.fasta
do
echo ${POSFILE}
sed "s/^>/>${POSFILE}_/"< ${POSFILE} >> ReferenceLibrary_positive.fas
done

for NEGFILE in neg_*.fasta
do
echo ${NEGFILE}
sed "s/^>/>${NEGFILE}_/"< ${NEGFILE} >> ReferenceLibrary_negative.fas
done

cat ReferenceLibrary_*.fas > ReferenceLibraries.fas


#Generate blast-searchable reference database
makeblastdb -in ReferenceLibraries.fas -parse_seqids -dbtype nucl

#Blast the assembly sequences against the reference library
tblastx -max_target_seqs 1 -evalue ${EVALUE} -outfmt "6 qseqid sseqid pident score evalue" -out ${ASSEMBLY}_ContaminationScreening_BLASThits.txt -db ReferenceLibraries.fas -query ${ASSEMBLY}


#Retrieve the corresponding sequences to a hit and write it out to new fasta files
grep -P "\tneg_" ${ASSEMBLY}_ContaminationScreening_BLASThits.txt | cut -f1 | sort | uniq > ${ASSEMBLY}_negativeSeqIDs.txt
cp ${ASSEMBLY} tmp.fas
echo ">" >> tmp.fas
grep "^>" ${ASSEMBLY} | sed "s/^>//" | sort | uniq | while read ID
do
 NEGTRUE=0
 while read -r LINE
 do
    if [ ${LINE} == ${ID} ]
    then NEGTRUE=1
    echo Negative ${ID} found
    fi
 done < ${ASSEMBLY}_negativeSeqIDs.txt
 if [ ${NEGTRUE} == 0 ]
 then echo ${ID}
 sed -n "/>${ID}$/,/>/p" < tmp.fas | sed '$d' >> ${ASSEMBLY}_pruned.fas
 fi
done

rm tmp.fas

