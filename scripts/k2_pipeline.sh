#!/bin/bash
set -euo pipefail

DB="k2_db_nr"
THREADS=64
TMPDIR="${TMPDIR:-/tmp}"
export LC_ALL=C

# 1. Rebuild the sorted taxonomy map using the versioned accession ($2)
awk 'NR>1 {print $2 "\t" $3}' "${DB}/taxonomy/prot.accession2taxid" \
  | sort -t $'\t' -k1,1 -T "$TMPDIR" --parallel="$THREADS" \
  > all_acc2taxid.sorted

# 2. Extract the versioned ACCNUM list
awk -F $'\t' '/^ACCNUM/ {print $2}' "${DB}/taxonomy/prelim_map.txt" \
  | sort -t $'\t' -k1,1 -T "$TMPDIR" --parallel="$THREADS" \
  > accmap.sorted

# 3. Join on the versioned accession
join -t $'\t' -1 1 -2 1 accmap.sorted all_acc2taxid.sorted \
  > seqid2taxid_acc.tmp

# 4. Prepend any manual TAXID entries
awk -F $'\t' '/^TAXID/ { sub(/^TAXID\t/,""); print }' \
    "${DB}/taxonomy/prelim_map.txt" \
  > seqid2taxid_map_file.tmp

cat seqid2taxid_acc.tmp >> seqid2taxid_map_file.tmp

# 5. Place the final map so Kraken skips its mapping step
mv seqid2taxid_map_file.tmp "${DB}/seqid2taxid.map"
echo "finished seq2tax successfully. sample seqid2taxid.map:"
head -n 5 k2_db_nr/seqid2taxid.map

./kraken2_scripts/kraken2-build --build --db k2_db_nr/ --protein --kmer-len 8 --minimizer-len 5 --minimizer-spaces 1 --threads 64 --no-masking --fast-build --max-db-size 125000000000
