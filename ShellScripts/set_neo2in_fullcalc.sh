#! /bin/bash

for folder in es_*; do
  filename="$folder/neo2.in"
  echo "filename: $filename"
  cp $filename $filename.sav

  sed -i "s/.*lsw_read_precom.*/ LSW_READ_PRECOM=T/" $filename
  sed -i "s/.*LSW_READ_PRECOM.*/ LSW_READ_PRECOM=T/" $filename
  sed -i "s/.*lsw_write_precom.*/ LSW_WRITE_PRECOM=F/" $filename
  sed -i "s/.*LSW_WRITE_PRECOM.*/ LSW_WRITE_PRECOM=F/" $filename
done
