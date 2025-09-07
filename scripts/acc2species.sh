#!/bin/bash

echo "Cleaning up taxid_to_name.tsv and joining files..."

join -t $'\t' -1 1 -2 2 -o '2.1,1.2' \
    <( \
        awk '{print $1 "\t" $2}' taxid_to_name.tsv | \
        LC_ALL=C sort -k1,1 \
    ) \
    <( \
        gzip -dc prot.accession2taxid.gz | \
        cut -f 2,3 | \
        LC_ALL=C sort -k2,2 \
    ) \
> acc_to_species.tsv
