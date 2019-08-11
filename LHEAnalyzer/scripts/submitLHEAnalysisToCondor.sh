#!/bin/bash

INFILE=""
DATE=""
OUTPUTDIR=""
CONDOROUTDIR=""
EXTRA_INCLUDES=""

let print_help=0

for var in "$@"; do
  farg="${var//\"}"
  fargl="$(echo $farg | awk '{print tolower($0)}')"
  if [[ "$fargl" == "infile="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    INFILE="$fcnargname"
  elif [[ "$fargl" == "date="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    DATE="$fcnargname"
  elif [[ "$fargl" == "outputdir="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    OUTPUTDIR="$fcnargname"
  elif [[ "$fargl" == "condoroutdir="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    CONDOROUTDIR="$fcnargname"
  elif [[ "$fargl" == "includefile="* ]];then
    fcnargname="$farg"
    fcnargname="${fcnargname#*=}"
    EXTRA_INCLUDES="${EXTRA_INCLUDES} ${fcnargname}"
  elif [[ "$fargl" == *"help" ]] || [[ "$fargl" == "-h" ]];then
    let print_help=1
  fi
done

if [[ $print_help -eq 1 ]];then
  echo "$0 help:"
  echo "infile: Input file for the list of arguments"
  echo "date: Submission date or id"
  echo "outputdir: Main output directory"
  echo "condoroutdir: Main directory to transfer Condor output"
  echo "includefile: Extra file/directory to include in the submission tarball"
  exit 0
fi


QUEUE="default"

hname=$(hostname)

CONDORSITE="DUMMY"
if [[ "$hname" == *"lxplus"* ]];then
  echo "Setting default CONDORSITE to cern.ch"
  CONDORSITE="cern.ch"
elif [[ "$hname" == *"ucsd"* ]];then
  echo "Setting default CONDORSITE to t2.ucsd.edu"
  CONDORSITE="t2.ucsd.edu"
fi

if [[ "$OUTPUTDIR" == "" ]];then
  OUTPUTDIR="./output"
fi
if [[ "$DATE" == "" ]];then
  DATE=$(date +%y%m%d)
fi

OUTDIR="${OUTPUTDIR}/${DATE}"

mkdir -p $OUTDIR

TARFILE="lheanalyzer.tar"
if [ ! -e ${OUTDIR}/${TARFILE} ];then
  cd ${OUTDIR}
  createLHEAnalyzerTarball.sh $EXTRA_INCLUDES
  cd -
fi


while IFS='' read -r line || [[ -n "$line" ]]; do
  THECONDORSITE="${CONDORSITE+x}"
  THECONDOROUTDIR="${CONDOROUTDIR+x}"
  fcnarglist=($(echo $line))
  fcnargname=""
  fcnargnameextra=""
  fcnargnametmp=""
  for fargo in "${fcnarglist[@]}";do
    farg="${fargo//\"}"
    fargl="$(echo $farg | awk '{print tolower($0)}')"
    if [[ "$fargl" == "outfile="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      fcnargname="${fcnargname//.root}"
      fcnargname="${fcnargname##*/}"
    elif [[ "$fargl" == "condorsite="* ]];then
      line="${line//$fargo}"
      fcnargnametmp="$fargo"
      fcnargnametmp="${fcnargnametmp#*=}"
      THECONDORSITE="$fcnargnametmp"
    elif [[ "$fargl" == "condoroutdir="* ]];then
      line="${line//$fargo}"
      fcnargnametmp="$farg"
      fcnargnametmp="${fcnargnametmp#*=}"
      THECONDOROUTDIR="$fcnargnametmp"
    fi
  done
  line="${line//  / }" # replace double whitespace with single
  line="${line## }" # strip leading space
  line="${line%% }" # strip trailing space
  if [[ "${THECONDORSITE+x}" != "DUMMY" ]] && [[ -z "${THECONDOROUTDIR+x}" ]]; then
    echo "Need to set THECONDOROUTDIR"
    continue
  fi
  if [[ "$fcnargnameextra" != "" ]];then
    fcnargname=$fcnargname"_"$fcnargnameextra
  fi
  if [[ ! -z $fcnargname ]];then
    theOutdir="${OUTDIR}/${fcnargname}"
    mkdir -p $theOutdir
    ln -sf ${PWD}/${OUTDIR}/${TARFILE} ${PWD}/${theOutdir}/

    submitLHEAnalyzerGenericJob.sh "$line" "$QUEUE" "$theOutdir" "$THECONDORSITE" "$THECONDOROUTDIR"
  fi
done < "$INFILE"
