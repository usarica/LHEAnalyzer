#!/bin/sh

getvarpaths(){
  for var in "$@";do
    tmppath=${var//:/ }
    for p in $(echo $tmppath);do
      if [[ -e $p ]];then
        echo $p
      fi
    done
  done  
}
searchfileinvar(){
  for d in $(getvarpaths $1);do
    for f in $(ls $d | grep $2);do
      echo "$d/$f"
    done
  done
}


CMSSWVERSION="$1"
SCRAMARCH="$2"
SUBMIT_DIR="$3" # Must be within $CMSSW_BASE/src/
TARFILE="$4"
FCNARGS="$5"
CONDORSITE="$6"
CONDOROUTDIR="$7"

export SCRAM_ARCH=${SCRAMARCH}

echo -e "\n--- begin header output ---\n" #                     <----- section division
echo "CMSSWVERSION: $CMSSWVERSION"
echo "SCRAMARCH: $SCRAMARCH"
echo "SUBMIT_DIR: $SUBMIT_DIR"
echo "TARFILE: $TARFILE"
echo "FCNARGS: $FCNARGS"
echo "CONDORSITE: $CONDORSITE"
echo "CONDOROUTDIR: $CONDOROUTDIR"

echo "GLIDEIN_CMSSite: $GLIDEIN_CMSSite"
echo "hostname: $(hostname)"
echo "uname -a: $(uname -a)"
echo "time: $(date +%s)"
echo "args: $@"
echo "tag: $(getjobad tag)"
echo "taskname: $(getjobad taskname)"
echo -e "\n--- end header output ---\n" #                       <----- section division

echo -e "\n--- begin memory specifications ---\n" #                     <----- section division
ulimit -a
echo -e "\n--- end memory specifications ---\n" #                     <----- section division


if [ -r "$OSGVO_CMSSW_Path"/cmsset_default.sh ]; then
  echo "sourcing environment: source $OSGVO_CMSSW_Path/cmsset_default.sh"
  source "$OSGVO_CMSSW_Path"/cmsset_default.sh
elif [ -r "$OSG_APP"/cmssoft/cms/cmsset_default.sh ]; then
  echo "sourcing environment: source $OSG_APP/cmssoft/cms/cmsset_default.sh"
  source "$OSG_APP"/cmssoft/cms/cmsset_default.sh
elif [ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
  echo "sourcing environment: source /cvmfs/cms.cern.ch/cmsset_default.sh"
  source /cvmfs/cms.cern.ch/cmsset_default.sh
else
  echo "ERROR! Couldn't find $OSGVO_CMSSW_Path/cmsset_default.sh or /cvmfs/cms.cern.ch/cmsset_default.sh or $OSG_APP/cmssoft/cms/cmsset_default.sh"
  exit 1
fi

INITIALDIR=$(pwd)

# If the first file in the tarball filelist starts with CMSSW, it is a
# tarball made outside of the full CMSSW directory and must be handled
# differently
if [ ! -z $(tar -tf ${TARFILE} | head -n 1 | grep "^CMSSW") ]; then
  echo "This is a full cmssw tar file."
  tar xf ${TARFILE}
  cd $CMSSWVERSION
  echo "Current directory ${PWD} =? ${CMSSWVERSION}"
  echo "Running ProjectRename"
  scramv1 b ProjectRename
else
  # Setup the CMSSW area
  echo "This is a selective CMSSW tar file."
  eval $(scramv1 project CMSSW $CMSSWVERSION)
  cd $CMSSWVERSION
fi


# Setup the CMSSW environment
eval $(scramv1 runtime -sh)
echo "CMSSW_BASE: ${CMSSW_BASE}"
echo "SCRAM_ARCH: ${SCRAM_ARCH}"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH-\'unset\'}"
echo "PYTHONPATH: ${PYTHONPATH-\'unset\'}"
echo "ROOT_INCLUDE_PATH: ${ROOT_INCLUDE_PATH-\'unset\'}"



# Ensure gfortran can be found
if [[ "$SCRAM_ARCH" == "slc7"* ]];then
  gcc --version
  g++ --version
  gfortran --version

  gfortran --print-file-name libgfortran.so
  gcc --print-file-name libgfortran.so

  echo $(ldconfig -p | grep gfortran)

  LIBGFORTRANDIR="$(ldconfig -p | grep gfortran | awk '{print($4)}')"
  LIBGFORTRANDIR="${LIBGFORTRANDIR%/*}"
  echo "gfortran installed in ${LIBGFORTRANDIR-\'unset\'}"
  if [[ ! -z "${LIBGFORTRANDIR+x}" ]];then
    beg=""
    if [[ -z "${LD_LIBRARY_PATH+x}" ]]; then
      beg=""
    else
      beg="${LD_LIBRARY_PATH}:"
    fi
    export LD_LIBRARY_PATH=$beg:$LIBGFORTRANDIR
  fi
fi


# Remove the tarfile
if [[ -e $INITIALDIR/${TARFILE} ]]; then
  echo "Moving the tarball from $INITIALDIR/${TARFILE} into "$(pwd)
  mv $INITIALDIR/${TARFILE} ./
  tar xf ${TARFILE}
  rm ${TARFILE}
else
  echo "The tarball does not exist in $INITIALDIR/${TARFILE}"
fi

if [[ -e $INITIALDIR/lheanalyzer_inputs ]]; then
  echo "Moving the tarball from $INITIALDIR/lheanalyzer_inputs into "$(pwd)
  mv $INITIALDIR/lheanalyzer_inputs ./
else
  echo "The tarball does not exist in $INITIALDIR/${TARFILE}"
fi

# Check the lib area as uploaded
echo "=============================="
echo "${CMSSW_BASE}/lib/${SCRAM_ARCH} as uploaded:"
ls ${CMSSW_BASE}/lib/${SCRAM_ARCH}
echo "=============================="


# Clean CMSSW-related compilation objects and print the lib area afterward
scramv1 b clean &>> compilation.log
echo "================================="
echo "${CMSSW_BASE}/lib/${SCRAM_ARCH} after cleaning:"
ls ${CMSSW_BASE}/lib/${SCRAM_ARCH}
echo "================================="

echo "================================="
echo "Contents of ${CMSSW_BASE}/src:"
ls -la ${CMSSW_BASE}/src
echo "================================="

echo "================================="
echo "Contents of /usr/lib/x86_64-linux-gnu:"
ls -la /usr/lib/x86_64-linux-gnu
echo "================================="


# MELA includes
LHEANALYZERDIR=${CMSSW_BASE}/src/LHEAnalyzer/LHEAnalyzer
if [[ ! -d $LHEANALYZERDIR ]];then
  echo "LHEANALYZERDIR=$LHEANALYZERDIR does not exist!"
  exit 1
fi

ZZMEDIR=${CMSSW_BASE}/src/ZZMatrixElement
if [[ ! -d $ZZMEDIR ]];then
  echo "ZZMEDIR=$ZZMEDIR does not exist!"
  exit 1
fi

CMSDATATOOLSDIR=${CMSSW_BASE}/src/CMSDataTools
if [[ ! -d $CMSDATATOOLSDIR ]];then
  echo "CMSDATATOOLSDIR=$CMSDATATOOLSDIR does not exist!"
  exit 1
fi

MELAANALYTICSDIR=${CMSSW_BASE}/src/MelaAnalytics
if [[ ! -d $MELAANALYTICSDIR ]];then
  echo "MELAANALYTICSDIR=$MELAANALYTICSDIR does not exist!"
  exit 1
fi

end=""

# Ensure CMSSW can find libmcfm
if [[ -z "${LD_LIBRARY_PATH+x}" ]]; then
  end=""
else
  end=":${LD_LIBRARY_PATH}"
fi
export LD_LIBRARY_PATH=${CMSSW_BASE}/src/ZZMatrixElement/MELA/data/${SCRAM_ARCH}$end
echo "LD_LIBRARY_PATH:";getvarpaths ${LD_LIBRARY_PATH}

# Ensure LHAPDF library path is also included in LD_LIBRARY_PATH
if [[ ! -z "${LHAPDF_DATA_PATH+x}" ]]; then
  export LD_LIBRARY_PATH=${LHAPDF_DATA_PATH}:${LD_LIBRARY_PATH}
else
  echo "CMSSW  configuration error: LHAPDF_DATA_PATH is undefined!"
  exit 1
fi

# Needed to locate the include directory of MELA classes. It can get lost.
if [[ -z "${ROOT_INCLUDE_PATH+x}" ]]; then
  end=""
else
  end=":${ROOT_INCLUDE_PATH}"
fi
export ROOT_INCLUDE_PATH=${CMSSW_BASE}/src/ZZMatrixElement/MELA/interface$end
echo "ROOT_INCLUDE_PATH:";getvarpaths ${ROOT_INCLUDE_PATH}


# Compile ZZMatrixElement
echo "Entering ZZMEDIR=${ZZMEDIR}"
cd $ZZMEDIR
./setup.sh clean
./setup.sh -j 12 &>> compilation.log
COMPILE_STATUS=$?
if [ $COMPILE_STATUS != 0 ];then
  echo "ZZMatrixElement compilation exited with error ${COMPILE_STATUS}. Printing the log:"
  cat compilation.log
fi
rm -f compilation.log
cd -

# Compile CMSSW-dependent packages
cd src
scramv1 b -j 12 &>> compilation.log
COMPILE_STATUS=$?
if [ $COMPILE_STATUS != 0 ];then
  echo "CMSSW compilation exited with error ${COMPILE_STATUS}. Printing the log:"
  cat compilation.log
fi
rm -f compilation.log
cd -


# Move inputs into running directory
if [[ -e lheanalyzer_inputs ]]; then
  mv lheanalyzer_inputs "src/${SUBMIT_DIR}/"
fi

# Go into the submission directory within $CMSSW_BASE/src
cd src/$SUBMIT_DIR

# Extract the inputs in the running directory
if [[ -e lheanalyzer_inputs ]]; then
  tar xvf lheanalyzer_inputs
fi

echo "Submission directory before running: ls -lrth"
ls -lrth


##############
# ACTUAL RUN #
##############
# Transfer needs to be done through the script.
# Script is actually run through bash to eliminate the extra processor consumption by python
echo -e "\n--- Begin RUN ---\n"
RUNFILE=LHEAnalyzer
RUN_CMD=$(runGenericExecutable.py --exe="$RUNFILE" --command="$FCNARGS" --dry) # Must run it dry
theTransferDir=""
theTransferFile=""
if [[ "$RUN_CMD" == "Running "* ]];then
  RUN_CMD="${RUN_CMD//Running }"
  echo "Analyzing command: $RUN_CMD"

  let fileLevel=0
  let pythiaStep=1
  indir="./"
  outdir="./"
  declare -a inputfiles
  outfile="tmp.root"
  declare -a restofcmds
  
  runargs=($(echo $RUN_CMD))
  for fargo in "${runargs[@]}";do
    farg="${fargo//\"}"
    fargl="$(echo $farg | awk '{print tolower($0)}')"
    if [[ "$fargl" == "lheanalyzer" ]];then      # Skip the executable name
      continue
    elif [[ "$fargl" == "filelevel="* ]];then    # Does not modify the command
      tmparg="$farg"
      tmparg="${tmparg#*=}"
      let fileLevel=$tmparg
      restofcmds+=($farg)
    elif [[ "$fargl" == "outdir="* ]];then       # Does not modify the command
      tmparg="$farg"
      tmparg="${tmparg#*=}"
      outdir="$tmparg"
      restofcmds+=($farg)
    elif [[ "$fargl" == "outfile="* ]];then      # Does not modify the command
      tmparg="$farg"
      tmparg="${tmparg#*=}"
      outfile="$tmparg"
      restofcmds+=($farg)
    elif [[ "$fargl" == "pythiastep="* ]];then
      tmparg="$farg"
      tmparg="${tmparg#*=}"
      let pythiaStep=$tmparg
    elif [[ "$fargl" == "indir="* ]];then
      tmparg="$farg"
      tmparg="${tmparg#*=}"
      indir="$tmparg"
    elif [[ "$fargl" != *"="* ]];then
      inputfiles+=($farg)
    else                                         # Does not modify the command
      restofcmds+=($farg)
    fi
  done

  if [[ "$indir" != *"/" ]];then
    indir=$indir"/"
  fi
  if [[ "$outdir" != *"/" ]];then
    outdir=$outdir"/"
  fi
  mkdir -p $outdir

  theTransferDir="${outdir}"
  theTransferFile="${outfile}"

  if [[ $fileLevel -eq 1 ]];then # Pythia processing
    pyindir="./"
    pythiaMainStr="pythiaTemp"
    declare -a pyinputfiles # These are supposed to replace inputfiles in running

    for index in "${!inputfiles[@]}"; do
      #FIXME: ak4 should come out of the options
      echo "Running trimPythia ${indir}${inputfiles[$index]} $pyindir $pythiaStep ak4"
      trimPythia "${indir}${inputfiles[$index]}" $pyindir $pythiaStep ak4
      RUN_STATUS=$?
      if [ $RUN_STATUS != 0 ]; then
        echo "Run has crashed with exit code ${RUN_STATUS}"
        exit 1
      fi

      mv "${pythiaMainStr}.root" "${pythiaMainStr}_${index}.root"
      pyinputfiles+=("${pythiaMainStr}_${index}.root")
    done

    echo "Running LHEAnalyzer" "${restofcmds[@]}" pythiaStep=-1 indir=$pyindir "${pyinputfiles[@]}"
    LHEAnalyzer "${restofcmds[@]}" pythiaStep=-1 indir=$pyindir "${pyinputfiles[@]}"
    RUN_STATUS=$?
    if [ $RUN_STATUS != 0 ]; then
      echo "Run has crashed with exit code ${RUN_STATUS}"
      exit 1
    fi
  else
    echo "Running LHEAnalyzer" "${restofcmds[@]}" indir=$indir "${inputfiles[@]}"
    LHEAnalyzer "${restofcmds[@]}" indir=$indir "${inputfiles[@]}"
    RUN_STATUS=$?
    if [ $RUN_STATUS != 0 ]; then
      echo "Run has crashed with exit code ${RUN_STATUS}"
      exit 1
    fi

  fi

else
  echo "Run command ${RUN_CMD} is invalid."
  exit 1
fi
echo -e "\n--- End RUN ---\n"
##############


##################
# TRANSFER FILES #
##################
if [[ "${theTransferDir}" == "./"* ]];then
  theTransferDir="${theTransferDir:1}"
fi
if [[ ! -z ${theTransferFile+x} ]];then
  current_transfer_dir=$(pwd)${theTransferDir}
  copyFromCondorToSite.sh ${current_transfer_dir} ${theTransferFile} ${CONDORSITE} ${CONDOROUTDIR}
  TRANSFER_STATUS=$?
  if [ $TRANSFER_STATUS != 0 ]; then
    echo " - Transfer crashed with exit code ${TRANSFER_STATUS}"
  fi
fi
##############


echo "Submission directory after running: ls -lrth"
ls -lrth

echo "time at end: $(date +%s)"
