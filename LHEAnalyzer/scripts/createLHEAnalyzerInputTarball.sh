#!/bin/sh


echo "Running $0"

EXTRA_INCLUDES=""
for var in "$@"; do
  if [[ ! -z ${var} ]];then
    if [[ -d $var ]];then
      EXTRA_INCLUDES="${EXTRA_INCLUDES} -C ${var} ."
    else
      fname=$var
      dirname=$(pwd)
      if [[ "$var" == *"/"* ]];then
        fname="${var##*/}"
        dirname="${var%/*}"
      fi
      EXTRA_INCLUDES="${EXTRA_INCLUDES} -C ${dirname} ${fname}"
    fi
  fi
done

if [[ ! -z "$EXTRA_INCLUDES" ]];then
  TARFILE="lheanalyzer_inputs"
  tar Jcvf ${TARFILE} ${EXTRA_INCLUDES}
fi
