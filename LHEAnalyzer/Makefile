HERE = $(PWD)
BASEINCLUDE=$(CMSSW_BASE)/src
MELADIR=$(BASEINCLUDE)/ZZMatrixElement/MELA
CANDLORECASTERDIR=$(BASEINCLUDE)/MelaAnalytics/CandidateLOCaster
EVENTCONTAINERDIR=$(BASEINCLUDE)/MelaAnalytics/EventContainer
GENERICMECOMPUTERDIR=$(BASEINCLUDE)/MelaAnalytics/GenericMEComputer

ROOTCFLAGS    = $(shell root-config --cflags) -Lrootlib -lEG
ROOTGLIBS     = $(shell root-config --glibs) -Lrootlib -lEG

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -g -O2

ARCH         := $(shell root-config --arch)

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
CXX           = c++ -lm
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  -j 4
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
UNDEFOPT      = suppress
endif

CXXFLAGS      += $(ROOTCFLAGS) -I$(ROOFITSYS)/include/ -I$(BASEINCLUDE) -I$(HERE)/interface/ -I$(MELADIR)/interface/ -I$(CANDLORECASTERDIR)/interface/ -I$(EVENTCONTAINERDIR)/interface/ -I$(GENERICMECOMPUTERDIR)/interface/

NGLIBS         = $(ROOTGLIBS)
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L $(ROOFITSYS)/lib/ -lMinuit -lRooFitCore -lRooFit #-lboost_regex
NGLIBS  += -L $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lZZMatrixElementMELA -lMelaAnalyticsCandidateLOCaster -lMelaAnalyticsEventContainer -lMelaAnalyticsGenericMEComputer
GLIBS          = $(filter-out -lNew, $(NGLIBS))


COMPILEPATH      = $(PWD)/
INCLUDEDIR       = $(COMPILEPATH)src/
OUTLIB           = $(COMPILEPATH)lib/
RUNDIR           = $(COMPILEPATH)


SOURCES = $(wildcard $(INCLUDEDIR)*.cc)
OBJECTS = $(SOURCES:.cc=.o)

.PHONY: all help compile clean

all: analyzer dolinks


# General rule for making object files
%.d:	$(INCLUDEDIR)%.cc $(INCLUDEDIR)%.cxx
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CXX) -MM -MT $@ -MT ${@:.d=.o} $(CXXFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	$(INCLUDEDIR)%.cc
	$(QUIET) echo "Compiling $<"; \
	$(CXX) $(CXXFLAGS) $< -c -o $@ $(GLIBS)



dolinks:
	cp $(COMPILEPATH)macros/loadLib.py $(RUNDIR)loadLib.py
	cp $(COMPILEPATH)macros/trimPythia.py $(RUNDIR)trimPythia.py
	ln -sf $(COMPILEPATH)macros/loadLib.cc $(RUNDIR)loadLib.cc
	ln -sf $(COMPILEPATH)macros/trimPythia.cc $(RUNDIR)trimPythia.cc

analyzer: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)analyzer $(INCLUDEDIR)*.o $(GLIBS)


clean:
	rm -rf $(OUTLIB)
	rm -rf $(RUNDIR)Pdfdata
	rm -f $(RUNDIR)input.DAT
	rm -f $(RUNDIR)process.DAT
	rm -f $(RUNDIR)ffwarn.dat
	rm -f $(RUNDIR)br.sm*
	rm -f $(INCLUDEDIR)*.o
	rm -f $(INCLUDEDIR)*.so
	rm -f $(INCLUDEDIR)*.d
	rm -f $(RUNDIR)*.cc
	rm -f $(RUNDIR)analyze*
	rm -f $(RUNDIR)*.o
	rm -f $(RUNDIR)*.so
	rm -f $(RUNDIR)*.d
	rm -f $(RUNDIR)*.pcm
	rm -f $(RUNDIR)*.py
	rm -f $(RUNDIR)*.pyc
	
