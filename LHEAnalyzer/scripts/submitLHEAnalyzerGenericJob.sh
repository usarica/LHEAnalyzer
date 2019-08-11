#!/bin/bash

SCRIPTCMD="$1"
QUEUE="$2"
OUTDIR="$3"
CONDORSITE="$4"
CONDOROUTDIR="$5"


echo "Calling the main submission script with the following arguments:"
echo "SCRIPTCMD: ${SCRIPTCMD}"
echo "QUEUE: ${QUEUE}"
echo "OUTDIR: ${OUTDIR}"
echo "CONDORSITE: ${CONDORSITE}"
echo "CONDOROUTDIR: ${CONDOROUTDIR}"


if [[ -z "$OUTDIR" ]];then
  echo "You must set the output directory!"
  exit 1
fi

CMSENVDIR=$CMSSW_BASE
if [[ -z "$CMSENVDIR" ]];then
  echo "Set up CMSSW first!"
  exit 1
fi


LOGSDIR=$OUTDIR"/Logs"
mkdir -p $LOGSDIR

extrainputcmd=""

extLog="LHEAnalysis"
if [[ ! -z "$SCRIPTCMD" ]];then
  fcnargname=""
  fcnargnameextra=""
  fcnargnametmp=""
  fcnarglist=($(echo $SCRIPTCMD))
  indir=""
  indircmd=""
  declare -a inputfiles
  for fargo in "${fcnarglist[@]}";do
    farg="${fargo//\"}"
    fargl="$(echo $farg | awk '{print tolower($0)}')"
    if [[ "$fargl" == "outfile="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      fcnargname="${fcnargname//.root}"
      fcnargname="${fcnargname##*/}"
    elif [[ "$fargl" == "indir="* ]];then
      fcnargname="$farg"
      fcnargname="${fcnargname#*=}"
      indir="${fcnargname}"
      indircmd="$farg"
    elif [[ "$fargl" != *"="* ]];then
      inputfiles+=($farg)
    fi
  done
  if [[ "$fcnargnameextra" != "" ]];then
    fcnargname=$fcnargname$fcnargnameextra
  fi
  if [[ "$fcnargname" == "" ]];then
    fcnargname=${SCRIPTCMD//" "/"_"}
  fi
  if [[ "$indir" != "" ]];then
    if [[ "$indir" != *"/" ]];then
      indir="${indir}/"
    fi
  fi
  fcnargname=${fcnargname//"="/"_"}
  fcnargname=${fcnargname//".root"}
  fcnargname=${fcnargname//".lhe"}
  fcnargname=${fcnargname//\"}
  fcnargname=${fcnargname//\!}
  fcnargname=${fcnargname//\\}
  fcnargname=${fcnargname//"("}
  fcnargname=${fcnargname//")"}
  fcnargname=${fcnargname//","/"_"}
  extLog=$extLog"_"$fcnargname

  # Assemble commands to add local input files into a tarball
  for efile in "${inputfiles[@]}";do
    if [[ -f "${indir}${efile}" ]] && [[ -s "${indir}${efile}" ]];then
      extrainputcmd="${extrainputcmd} --add_input=${indir}${efile}"
    fi
  done
  echo $extrainputcmd
  # If an input tarball is to be made, change the original indir to ./
  if [[ ! -z "$extrainputcmd" ]];then
    SCRIPTCMD="${SCRIPTCMD/$indircmd/indir=./}"
  fi
fi


hname=$(hostname)
echo $hname
TARFILE="lheanalyzer.tar"
if [[ "$hname" == *"lxplus"* ]] || [[ "$hname" == *"ucsd"* ]];then
  echo "Host is on LXPLUS or UCSD, so need to use HTCONDOR"
  THEQUEUE="vanilla"
  if [[ "$QUEUE" != "default" ]];then
    THEQUEUE=$QUEUE
  fi
  checkGridProxy.sh
  configureLHEAnalyzerCondorJob.py ${extrainputcmd} --configonly --tarfile="$TARFILE" --batchqueue="$THEQUEUE" --outdir="$OUTDIR" --outlog="Logs/log_$extLog" --errlog="Logs/err_$extLog" --batchscript="submitLHEAnalyzerGenericJob.condor.sh" --command="$SCRIPTCMD" --condorsite="$CONDORSITE" --condoroutdir="$CONDOROUTDIR"
elif [[ "$hname" == *"login-node"* ]] || [[ "$hname" == *"bc-login"* ]]; then
  echo "Host is on MARCC, so need to use SLURM batch."
  echo "MARCC submission is not supported."
  exit 1
fi
