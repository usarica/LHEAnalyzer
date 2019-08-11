#!/bin/env python

import sys
import imp
import copy
import os
import filecmp
import shutil
import pickle
import math
import pprint
import subprocess
from datetime import date
from argparse import ArgumentParser as OptionParser
from CMSDataTools.AnalysisTree.TranslateStringBetweenPythonAndShell import *


class BatchManager:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()

      self.parser.add_argument("--batchqueue", type=str, help="Batch queue")
      self.parser.add_argument("--batchscript", type=str, help="Name of the HTCondor executable")
      self.parser.add_argument("--tarfile", type=str, help="Name of the tar file to upload")
      self.parser.add_argument("--outdir", type=str, help="Name of the job output directory")
      self.parser.add_argument("--condorsite", type=str, help="Name of the HTCondor site")
      self.parser.add_argument("--condoroutdir", type=str, help="Name of the HTCondor output directory")
      self.parser.add_argument("--outlog", type=str, help="Name of the output log file")
      self.parser.add_argument("--errlog", type=str, help="Name of the output error file")

      self.parser.add_argument("--command", type=str, help="The arguments of the function")

      self.parser.add_argument("--add_input", type=str, action="append", help="Inputs to include other than the tar file")

      self.parser.add_argument("--dry", action="store_true", default=False, help="Do not submit jobs, just set up the files")
      self.parser.add_argument("--configonly",  action="store_true", default=False, help="Only configure the jobs, do not submit them")
      self.parser.add_argument("--interactive", action="store_true", default=False, help="Do not submit jobs; run them interactively")


      self.opt = self.parser.parse_args()

      optchecks=[
         "batchqueue",
         "batchscript",
         "tarfile",
         "outdir",
         "condorsite",
         "condoroutdir",
         "outlog",
         "errlog",
         "command"
      ]
      for theOpt in optchecks:
         if not hasattr(self.opt, theOpt) or getattr(self.opt, theOpt) is None:
            sys.exit("Need to set --{} option".format(theOpt))

      # Get abolute path
      self.opt.outdir = os.path.abspath(self.opt.outdir)

      if not os.path.isfile(self.opt.batchscript):
         print "Batch executable does not exist in current directory, will search for CMSSW_BASE/bin"
         if os.path.isfile(os.getenv("CMSSW_BASE")+"/bin/"+os.getenv("SCRAM_ARCH")+"/"+self.opt.batchscript):
            self.opt.batchscript = os.getenv("CMSSW_BASE")+"/bin/"+os.getenv("SCRAM_ARCH")+"/"+self.opt.batchscript
            print "\t- Found the batch executable"
         else:
            sys.exit("Batch executable {} does not exist. Exiting...".format(self.opt.batchscript))

      for theOpt in optchecks:
         print "Option {}={}".format(theOpt,getattr(self.opt, theOpt))

      if self.opt.add_input is None:
         self.opt.add_input = []

      self.submitJobs()


   def produceCondorScript(self):
      currentdir = os.getcwd()
      currentCMSSWBASESRC = os.getenv("CMSSW_BASE")+"/src/" # Need the trailing '/'
      currendir_noCMSSWsrc = currentdir.replace(currentCMSSWBASESRC,'')
      if self.opt.command is not None:
         self.opt.command = translateFromPythonToShell(self.opt.command)

      scramver = os.getenv("SCRAM_ARCH")
      singularityver = "cms:rhel6"
      if "slc7" in scramver:
         singularityver = "cms:rhel7"

      inputfiles = ""
      for infile in self.opt.add_input:
         inputfiles += " {}".format(infile)

      inputfilescmd=""
      if inputfiles:
         inputfilescmd=",lheanalyzer_inputs" # Comma is needed because it separates from TARFILE

      condorscriptargs = {
         "home" : os.path.expanduser("~"),
         "uid" : os.getuid(),
         "batchScript" : self.opt.batchscript,
         "CONDORSITE" : self.opt.condorsite,
         "CONDOROUTDIR" : self.opt.condoroutdir,
         "outDir" : self.opt.outdir,
         "outLog" : self.opt.outlog,
         "errLog" : self.opt.errlog,
         "QUEUE" : self.opt.batchqueue,
         "CMSSWVERSION" : os.getenv("CMSSW_VERSION"),
         "SCRAMARCH" : scramver,
         "SINGULARITYVERSION" : singularityver,
         "SUBMITDIR" : currendir_noCMSSWsrc,
         "TARFILE" : self.opt.tarfile,
         "INPUTFILES" : inputfilescmd,
         "RUNCMD" : self.opt.command
      }

      condorscriptcontents = """
universe={QUEUE}
+DESIRED_Sites="T2_US_UCSD,T2_US_Caltech,T3_US_UCR,T2_US_MIT,T2_US_Vanderbilt,T2_US_Wisconsin,T3_US_Baylor,T3_US_Colorado,T3_US_NotreDame,T3_US_Rice,T3_US_Rutgers,T3_US_UMD,T3_US_Vanderbilt_EC2,T3_US_OSU"
executable              = {batchScript}
arguments               = {CMSSWVERSION} {SCRAMARCH} {SUBMITDIR} {TARFILE} {RUNCMD} {CONDORSITE} {CONDOROUTDIR}
Initialdir              = {outDir}
output                  = {outLog}.$(ClusterId).$(ProcId).txt
error                   = {errLog}.$(ClusterId).$(ProcId).err
log                     = $(ClusterId).$(ProcId).log
request_memory          = 4000M
+JobFlavour             = "tomorrow"
x509userproxy           = {home}/x509up_u{uid}
#https://www-auth.cs.wisc.edu/lists/htcondor-users/2010-September/msg00009.shtml
periodic_remove         = JobStatus == 5
transfer_executable=True
transfer_input_files    = {TARFILE}{INPUTFILES}
transfer_output_files = ""
+Owner = undefined
+project_Name = "cmssurfandturf"
notification=Never
should_transfer_files = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
Requirements = ((HAS_SINGULARITY=?=True) && (HAS_CVMFS_cms_cern_ch =?= true)) || (regexp("(uaf-[0-9]{{1,2}}|uafino)\.", TARGET.Machine) && !(TARGET.SlotID>(TotalSlots<14 ? 3:7) && regexp("uaf-[0-9]", TARGET.machine)))
+SingularityImage = "/cvmfs/singularity.opensciencegrid.org/bbockelm/{SINGULARITYVERSION}"

queue

"""
      condorscriptcontents = condorscriptcontents.format(**condorscriptargs)

      self.condorScriptName = "condor.sub"
      condorScriptFile = open(self.opt.outdir+"/"+self.condorScriptName,'w')
      condorScriptFile.write(condorscriptcontents)
      condorScriptFile.close()

      if inputfiles:
         tarcmd="cd {}; createLHEAnalyzerInputTarball.sh {}; cd -".format(self.opt.outdir, inputfiles)
         os.system( tarcmd )


   def submitJobs(self):
      self.produceCondorScript()

      jobcmd = None
      if not self.opt.configonly:
         jobcmd = "cd {}; condor_submit {}; cd -".format(self.opt.outdir, self.condorScriptName)
         if self.opt.dry:
            jobcmd = "echo " + jobcmd
      else:
         jobcmd = "echo Configured {}/{}".format(self.opt.outdir, self.condorScriptName)

      ret = os.system( jobcmd )



if __name__ == '__main__':
   batchManager = BatchManager()
