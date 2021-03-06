#
# common.mk contains the configuration for this build setup.
#
#
ifeq ($(origin PREFIX), undefined)
  ifeq ($(origin SEYMOUR_HOME), undefined)
    PREFIX = /opt/seymour
  else
    PREFIX = $(SEYMOUR_HOME)
  endif
endif
ANALYSIS_HOME=$(PREFIX)/analysis

PBCPP_DIR = ..

include ../common.mk

#
# Define the targets before including the rules since the rules contains a target itself.
#

EXECS = wordCounter printReadWordCount blasr sdpMatcher swMatcher kbandMatcher sawriter saquery samodify printTupleCountTable cmpPrintTupleCountTable malign removeAdapters tabulateAlignment samatcher saprinter buildQualityValueProfile guidedalign extendAlign sals pbmask

# DISABLE for now
#cmpMatcher


all: bin make.dep $(EXECS)

BUILTEXECS = $(addprefix bin/, $(EXECS))
DISTRIB_SET = blasr swMatcher kbandMatcher sawriter samodify printTupleCountTable cmpPrintTupleCountTable removeAdapters sdpMatcher pbmask
DISTRIB_EXECS = $(addprefix bin/, $(DISTRIB_SET))
INSTALL_EXECS = $(addprefix install-, $(DISTRIB_SET))


include ../make.rules

ifneq ($(wildcard make.dep),)
    include make.dep
endif

wordCounter:        bin/wordCounter
printReadWordCount: bin/printReadWordCount
blasr:        bin/blasr
cmpMatcher:         bin/cmpMatcher
sdpMatcher:         bin/sdpMatcher
samatcher:          bin/samatcher
saprinter:          bin/saprinter
sals:               bin/sals
swMatcher:          bin/swMatcher
kbandMatcher:       bin/kbandMatcher
sawriter:           bin/sawriter
saquery:            bin/saquery
samodify:           bin/samodify
printTupleCountTable:bin/printTupleCountTable
cmpPrintTupleCountTable:bin/cmpPrintTupleCountTable
malign:	            bin/malign
removeAdapters:     bin/removeAdapters
tabulateAlignment:  bin/tabulateAlignment
buildQualityValueProfile: bin/buildQualityValueProfile
guidedalign: bin/guidedalign
extendAlign: bin/extendAlign
pbmask: bin/pbmask

bin/sawriter: bin/SAWriter.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/pbmask: bin/Mask.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lz

bin/guidedalign: bin/GuidedAlign.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lz

bin/extendAlign: bin/ExtendAlign.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@ -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lz

bin/saquery: bin/SAQuery.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@

bin/samodify: bin/SAModify.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@

bin/wordCounter: bin/WordCounter.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@

bin/printReadWordCount: bin/PrintReadWordCount.o
	$(CPP) $(CPPOPTS) $^ $(STATIC) -o $@

#
# Add nonstandard compilation rules to build compressed dna aligners
#

bin/CmpReadMatcher.o: ReadMatcher.cpp
	$(CPP) -c $(CPPOPTS) -DCOMPRESSED $< $(STATIC) -o $@

bin/cmpMatcher: bin/CmpReadMatcher.o bin/PositionTable.o
	$(CPP) $(CPPOPTS) $< bin/PositionTable.o $(STATIC) -o $@ -lpthread -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lz -lpthread

bin/CmpPrintTupleCountTable.o: PrintTupleCountTable.cpp
	$(CPP) -c $(CPPOPTS) -DCOMPRESSED $^ $(STATIC) -o $@

bin/cmpPrintTupleCountTable: bin/CmpPrintTupleCountTable.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

ifneq ($(shell uname -s),Darwin)
    LRT = -lrt
endif
bin/blasr: bin/Blasr.o
	$(CPP) $(CPPOPTS) $< -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) $(LINK_PROFILER) -lpthread -lz $(LRT) -ldl $(STATIC) -o bin/blasr

bin/samatcher: bin/SAMatcher.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/saprinter: bin/SAPrinter.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/sals: bin/SALS.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/sdpMatcher: bin/SDPMatcher.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/malign: bin/MAlign.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/swMatcher: bin/SWMatcher.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/removeAdapters: bin/RemoveAdapters.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/kbandMatcher: bin/KBandMatcher.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/printTupleCountTable: bin/PrintTupleCountTable.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/tabulateAlignment: bin/TabulateAlignment.o
	$(CPP) $(CPPOPTS) $< $(STATIC) -o $@

bin/buildQualityValueProfile: bin/BuildQualityValueProfile.o
	$(CPP) -g $(CPPOPTS) $< $(STATIC) -o $@ -lpthread -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lz

#
# Set up a default value for the install dir if one does
# not exist.
#
INSTALL_DIR ?= $(ANALYSIS_HOME)/bin
BUILD_DIR ?= bin
install:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)

install-%:
	/usr/bin/install -d $(INSTALL_DIR)
	/usr/bin/install -m 555 $(DISTRIB_EXECS) $(INSTALL_DIR)
