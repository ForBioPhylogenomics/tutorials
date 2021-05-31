#!/bin/bash
# (c) Frédéric Mahé & modified by Torsten Hugo Struck @12.05.2021 adjusting to new databases in NCBI

NAMES="names.dmp"
NODES="nodes.dmp"
GI_TO_TAXID="nucl_gb.accession2taxid"
TAXONOMY=""
GI="${1}"

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

exit 0
