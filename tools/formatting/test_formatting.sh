#!/bin/bash

# Test if prettify was run
# author: Ole Schuett

set -e
rm -rf preprettify

find ./src/ -type f -not -path "*/preprettify/*" -not -path "*/.svn/*" -print0 | xargs -0 md5sum > checksums.md5
md5sum ./data/POTENTIAL >> checksums.md5

cd makefiles
make --jobs=20 pretty
make --jobs=20 pretty  # run twice to ensure consistency with doxify
cd ..

cd data
cat GTH_POTENTIALS HF_POTENTIALS NLCC_POTENTIALS ALL_POTENTIALS > POTENTIAL
cd ..

nfiles=`wc -l checksums.md5 | cut -f 1 -d " "`
summary="Checked $nfiles files."
status="OK"

echo "Searching for doxify warnings ..."
if grep -r -e "UNMATCHED_PROCEDURE_ARGUMENT" \
           -e "UNKNOWN_DOXYGEN_COMMENT" \
           -e "UNKNOWN_COMMENT" \
           --exclude-dir=".svn" \
           --exclude-dir="preprettify" \
           ./src/* ; then
  summary="Found doxify warnings"
  status="FAILED"
fi

echo "Comparing MD5-sums ..."
if ! md5sum --quiet --check checksums.md5 ; then
  summary='Code not invariant under "make pretty"'
  status="FAILED"
fi

rm checksums.md5

echo "Summary:" $summary
echo "Status:" $status

#EOF
