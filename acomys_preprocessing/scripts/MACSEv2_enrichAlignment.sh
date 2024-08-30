#!/bin/bash

# report piping errors
set -euo pipefail

# get arguments
added_acomys_projections="$1"
multiple_alignment_link="$2"
cactus_alignments="$3"
fasta_output="$4"
MACSEv2_URL="$5"

# download MACSEv2 version 2.07 in current working directory 
wget "${MACSEv2_URL}"

while IFS=$'\t' read -r -a p && [ -n "$p" ]
do
     # skip if first line 
    if [ "${p[0]}" = "MSA_filename" ]; then
        continue
    fi
    
    # store MSA_filename
    MSA_filename_gz="${p[0]}"
    # remove .gz extention from MSA_filename_gz
    read MSA_filename < <(echo ${MSA_filename_gz} | sed 's/.gz$//g')

    # store MSA_transcript
    MSA_transcript="${p[1]}"

    # store projection
    projection="${p[2]}" 

    # download and decompress the MSAs.
    wget --no-check-certificate -q -O - "${multiple_alignment_link}${MSA_filename_gz}" | gzip -d > "${cactus_alignments}${MSA_filename}"

    java -jar macse_v2.07.jar -prog enrichAlignment -align "${cactus_alignments}${MSA_filename}" -seq "${fasta_output}HLacoCah2_${MSA_transcript}.fa" -gc_def 1
    echo

done < "${added_acomys_projections}"
