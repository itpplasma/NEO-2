#!/bin/bash


filename="$1"

foldercommonpart="$2"

# Print help message and exit.
if [ "x$filename" = "x" -o "x$filename" = "x-?" -o "x$filename" = "x-h" -o "x$filename" = "x--help" ] ; then
  echo ""
  echo "Collect data from hdf5 files of a scan."
  echo ""
  echo "Invokation from directory where the scan is located:"
  echo "  collect_hdf5_data.sh filename foldercommonpart"
  echo "The script will look for files 'filename' in the subdirectories"
  echo "'foldercommonpart*'."
  echo "The output will be named 'final_filename'."
  echo ""

  exit
fi

collectionname="final_$filename"

for file in `ls ${foldercommonpart}*/$filename` ; do
  # See https://confluence.hdfgroup.org/display/HDF5/h5copy or use 'info
  # h5copy' for more information about h5copy.
  path=`dirname "$file"`
  h5copy -v -p -i $file -o "$collectionname" -s "/" -d "/$path"
done
