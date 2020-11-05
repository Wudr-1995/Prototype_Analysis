# Makefile for the ROOT test progrrpc.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ARCH          = linuxegcs

CXX           =
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o 

EVENTLIB      = $(EVENTSO)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


ifeq ($(ARCH),linuxegcs)
CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC -g
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
endif

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#------------------------------------------------------------------------------
EXEFILO       	= main.$(ObjSuf)
EXEFILES       	= main.$(SrcSuf)
EXEFILE        	= getGain$(ExeSuf)
PROCESSEVENTO   = Process.$(ObjSuf)
PROCESSEVENTS   = Process.$(SrcSuf) 
GETPARSEVENTO	= GetSpecParameters.$(ObjSuf)
GETPARSEVENTS	= GetSpecParameters.$(SrcSuf)

OBJS          = $(EXEFILO) $(PROCESSEVENTO) $(GETPARSEVENTO)
PROGRRPC      = $(EXEFILE)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

all:            $(PROGRRPC)

$(EXEFILE):     $(EXEFILO) $(PROCESSEVENTO) $(GETPARSEVENTO)
		$(LD) $(LDFLAGS) $(EXEFILO) $(PROCESSEVENTO) $(GETPARSEVENTO)\
		$(GLIBS) -lMinuit $(OutPutOpt)$(EXEFILE)
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core $(EXEFILE) *~

distclean:      clean
		@rm -f $(PROGRRPC) *Dict.* *.def *.exp \
		   *.root *.ps *.so .def so_locations
		@rm -rf cxx_repository

.SUFFIXES: .$(SrcSuf)

###
$(GETPARSEVENTO): GetSpecParameters.h GetSpecParameters.cc
$(PROCESSEVENTO): Process.h Process.cc
$(EXEFILO): main.cc Process.h Process.cc GetSpecParameters.h GetSpecParameters.cc
.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
