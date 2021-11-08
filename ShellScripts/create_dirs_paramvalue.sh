#!/bin/bash
SCRIPTNAME=$(basename $0 .sh)
EXIT_SUCCESS=0
EXIT_FAILURE=1

function usage {
  echo "Usage: ./$SCRIPTNAME [-d (dry-run)]" >&2
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
    usage $EXIT_FAILURE
    ;;
  esac
done

# Do for every line(=file) that is not either a coment or empty.
for i in $(cat dirs_list.txt | grep -v '^#' | grep -v "^$"); do
  echo "Creating directory $i..."
  if [ $DRY -eq 0 ] ; then
    cp -r TEMPLATE_DIR/ $i
    cd $i
  fi
  string=$i
  array=(${string//-/ })
  for j in "${!array[@]}"
  do
    echo "Detected param-value attribute $j => ${array[j]}"
    string=${array[j]}
    paramvalue=(${string//=/ })
    for k in "${!paramvalue[@]}"
    do
      if [ $k == 0 ]
      then
        echo "Param: $k=>${paramvalue[k]}"
        param=${paramvalue[k]}
      else
        echo "Value: $k=>${paramvalue[k]}"
        value=${paramvalue[k]}
      fi
    done
    if [ $DRY -eq 0 ] ; then
      value=$(echo $value | sed "s#minus#-#" | sed "s#m#d-#" | sed "s#p#.#")
      #value="prop_reconstruct = $1"
      sed -i "s#${param}\s*=.*#${param} = ${value}#" neo2.in
      #sed -i "s#<${param}>#${value}#" neo2.in
    fi
  done
  if [ $DRY -eq 0 ] ; then
    cd ..
  fi

  #echo "Creating $i..."
  #value=$(echo $i | sed "s#m#d-#" | sed "s#p#.#")
  #echo "Transforming to $value"
  #if [ $DRY -eq 0 ] ; then
  #  cp -r TEMPLATE_DIR/ eta_part_$i
  #  cd eta_part_$i
  #  sed -i "s#<eta_part>#${value}#" neo2.in
  #  cd ..
  #fi
done;
