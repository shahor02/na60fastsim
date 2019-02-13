# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linux
CXX           =
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------

GENMUONLMRO     = GenMUONLMR.o
GENMUONLMRS	= GenMUONLMR.cpp GenMUONLMRDict.cpp
GENMUONLMRSO	= GenMUONLMR.o GenMUONLMRDict.o
GENMUONLMRH	= GenMUONLMR.h GenMUONLMRDict.h

OBJS            = source_interface.o stuprf.o magfld.o mgdraw.o usrout.o

#------------------------------------------------------------------------------

all:            libGenMUONLMR.so FluLib.o rootfluka

libGenMUONLMR.so:  $(GENMUONLMRSO) $(GENMUONLMRSH)
		   $(LD) $(SOFLAGS) $(LDFLAGS) $(GENMUONLMRSO) $(GLIBS) -o libGenMUONLMR.so
		   @echo "$@ done"

GenMUONLMRDict.cpp:GenMUONLMR.h GenMUONLMRLinkDef.h
		@echo "Generating dictionary $@..."
		@$(ROOTSYS)/bin/rootcint -f GenMUONLMRDict.cpp -c -p GenMUONLMR.h GenMUONLMRLinkDef.h

FluLib.o:   GenMUONLMR.h
	    @echo "Generating Library $@..."
	    g++ -c FluLib.cpp $(ROOTCFLAGS)

.f.$(ObjSuf):
	       $(FLUPRO)/flutil/fff $?

rootfluka:     $(OBJS) FluLib.$(ObjSuf) libGenMUONLMR.so
	       $(FLUPRO)/flutil/lfluka_root -o $@ -m fluka $?

clean:
		@rm -f rootfluka GenMUONLMR.o GenMUONLMRDict.o FluLib.o core *.so GenMUONLMRDict.cpp $(OBJS)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
