#!/bin/sh

#  CleaningTaxaSupermatrix.sh
#  
#
#  Created by Torsten Hugo Struck on 26/05/2021.
#

# USAGE: CleaningTaxaSupermatrix.sh [Supermatrix fasta file] [File with taxa to exclude]

# Set variables with the relevant paths
SUPERMATRIX=$1
EXCLUSION=$2

#Exclude the corresponding sequences above the RCFV threshold and write it out the others ones to new fasta files
cp ${SUPERMATRIX} tmp.fas
echo ">" >> tmp.fas
grep "^>" ${SUPERMATRIX} | sed "s/^>//" | sort | uniq | while read ID
do
 NEGTRUE=0
 while read -r LINE
 do
    if [ ${LINE} == ${ID} ]
    then NEGTRUE=1
    echo Negative ${ID} found
    fi
 done < ${EXCLUSION}
 if [ ${NEGTRUE} == 0 ]
 then echo ${ID}
 sed -n "/>${ID}$/,/>/p" < tmp.fas | sed '$d' >> ${SUPERMATRIX}_pruned.fas
 fi
done

rm tmp.fas

