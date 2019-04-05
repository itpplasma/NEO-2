#!/bin/bash


filename="$1"

foldercommonpart="$2"

blacklistfile="$3"

collectionname="final_$filename"

function print_help {
  echo ""
  echo "Collect data from hdf5 files of a scan."
  echo ""
  echo "Invokation from directory where the scan is located:"
  echo "  collect_hdf5_data.sh filename foldercommonpart [blacklistfile]"
  echo "The script will look for files 'filename' in the subdirectories"
  echo "'foldercommonpart*'."
  echo "The output will be named 'final_filename'."
  echo "The optional argument blacklistfile is a file which contains"
  echo "foldernames that should be ignored."
  echo ""
}

function copy_file_to_group {
  input_file_local="$1"
  output_file_local="$2"
  path_local="$3"
  # See https://confluence.hdfgroup.org/display/HDF5/h5copy or use 'info
  # h5copy' for more information about h5copy.
  h5copy -v -p -i $input_file_local -o "$output_file_local" -s "/" -d "/$path_local"
}

# Print help message and exit.
if [ "x$filename" = "x" -o "x$filename" = "x-?" -o "x$filename" = "x-h" -o "x$filename" = "x--help" ] ; then
  print_help # function call
  exit
fi

for file in `ls ${foldercommonpart}*/$filename` ; do
  path=`dirname "$file"`
  if [ -e "$blacklistfile" ] ; then
    count_found=`grep -c "$path" "$blacklistfile"`
    if [ "x$count_found" = "x0" ] ; then
      copy_file_to_group "$file" "$collectionname" "$path" # function call
    fi
  else
    copy_file_to_group "$file" "$collectionname" "$path" # function call
  fi

done
