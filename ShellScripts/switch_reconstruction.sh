#!/bin/bash
SCRIPTNAME=$(basename $0 .sh)
EXIT_SUCCESS=0
EXIT_FAILURE=1

function usage {
 echo "Usage: ./$SCRIPTNAME [-d (dry-run)] Value" >&2
 [[ $# -eq 1 ]] && exit $1 || exit $EXIT_FAILURE
}

DRY=0
while getopts 'dh' OPTION ; do
 case $OPTION in
 d) echo "Running in dry-mode, no directory is created"; DRY=1
 ;;
 h) usage $EXIT_SUCCESS
 ;;
 \?) echo "Unknown Option \"-$OPTARG\"." >&2
 usage $EXIT_ERROR
 ;;
 esac
done

if [ -f jobs_list.txt ];
then
for i in $(cat jobs_list.txt | grep -v '^#' | grep -v "^$"); do
  echo "Running batch..."
  echo "Switching $i..."
  value="prop_reconstruct = $1"
  if [ $DRY -eq 0 ] ; then
    cd $i
    sed -i "s#prop_reconstruct.*=.*#${value}#" neo2.in
    cd ..
  fi
done;
else
  echo "Modifying neo2.in..."
  value="prop_reconstruct = $1"
  sed -i "s#prop_reconstruct.*=.*#${value}#" neo2.in
fi
