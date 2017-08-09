DIR := ${CURDIR}

HDR = ./interface/
SRC = ./src/
PLG = ./plugins/
PRG = ./main/
OBJ = ./lib/
LIB = ./lib/
BIN = ./bin/

HDRSuf = .h
SRCSuf = .cc
PRGSuf = .cpp
OBJSuf = .o
LIBSuf = .so

HDRS     =  $(wildcard $(HDR)*$(HDRSuf))
SRCS     =  $(wildcard $(SRC)*$(SRCSuf))
PLGS     =  $(wildcard $(PLG)*$(SRCSuf))
_OBJS    =  $(patsubst %$(SRCSuf), %$(OBJSuf), $(SRCS))
OBJS     =  $(patsubst $(SRC)%, $(OBJ)%, $(_OBJS))
_PLGOBJS =  $(patsubst %$(SRCSuf), %$(OBJSuf), $(PLGS))
PLGOBJS  =  $(patsubst $(PLG)%, $(OBJ)%, $(_PLGOBJS))
_LIBS    =  $(patsubst %$(SRCSuf), %$(LIBSuf), $(PLGS))
LIBS     =  $(patsubst $(PLG)%, $(LIB)lib%, $(_LIBS))
PRGS     =  $(wildcard $(PRG)*$(PRGSuf))
_BINS    =  $(wildcard $(PRG)*$(PRGSuf))
__BINS   =  $(_BINS:$(PRGSuf)=$(BINSuf))
___BINS  =  $(notdir $(__BINS))
BINS     =  $(addprefix $(BIN),${___BINS})

LINKDEF   =  $(wildcard ${HDR}*LinkDef.h)
DICTHDRS  =  $(patsubst $(LINKDEF),,$(HDRS)) $(LINKDEF)


ARCH  =  $(shell root-config --arch)

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs) -lGenVector -lFoam -lMinuit -lTMVA -lMLP -lXMLIO  -lTreePlayer -lMathMore



CXX  =  g++
CXXFLAGS  = -Wall -Wno-sign-compare -Wno-overloaded-virtual -O0 -g -fPIC -I$(DIR) $(ROOTCFLAGS) 

CPP  =  g++
CPPFLAGS  = -Wall -Wno-sign-compare -Wno-overloaded-virtual -O0 -g -I$(DIR) $(ROOTCFLAGS)

LD       =  g++
LDFLAGS  =  -rdynamic -shared -O2
SONAME	 =  libH4Analysis.so
SOFLAGS  =  -Wl,-soname,

GLIBS   =  -lm -L./DynamicTTree/lib -L./CfgManager/lib -Wl,-rpath=lib/:DynamicTTree/lib/:CfgManager/lib/ -lDTT -lCFGMan $(ROOTGLIBS)



#################################################
#if mac 64
ifeq ($(ARCH),macosx64)
LIBSuf  =  .dylib

CPPFLAGS  =  -Wall -W -O2 -pipe -I$(HDR) $(ROOTCFLAGS)

CXXFLAGS  =  -Wall -W -O2 -pipe -I$(HDR) $(ROOTCFLAGS)

LDFLAGS  =  -dynamiclib -shared -single_module -undefined dynamic_lookup
SONAME	 =  libH4Analysis.dylib
SOFLAGS  =
endif
#################################################



.PHONY: all clean test


all: dynTTree cfgMan $(LIB)$(SONAME) $(LIBS) $(BINS)

test:
	@echo "HDRS = $(HDRS)"
	@echo "DICTHDRS = $(DICTHDRS)"
	@echo "SRCS = $(SRCS)"
	@echo "PLGS = $(PLGS)"
	@echo "PRGS = $(PRGS)"
	@echo "OBJS = $(OBJS)"
	@echo "PLGOBJS = $(PLGOBJS)"
	@echo "LIBS = $(LIBS)"
	@echo "BINS = $(BINS)"

$(BIN)%: $(PRG)%$(PRGSuf) $(HDRS) $(LIB)$(SONAME) Makefile
	@echo " CXX $<"
	$(CPP) $(CPPFLAGS) $(GLIBS) -L$(LIB) -lH4Analysis -o $@ $<

$(OBJ)%$(OBJSuf): $(SRC)%$(SRCSuf) Makefile
	@echo " CXX $<"
	@$ $(CXX) -c $(CXXFLAGS) -o $@ $< 

$(LIB)mydict.cc: $(DICTHDRS)
	@echo "Generating dictionary..."
	@$ rootcling -f $(LIB)mydict.cc -c -p ${CXXFLAGS} $(DICTHDRS)

$(LIB)mydict.o: $(LIB)mydict.cc
	@echo " CXX $<"	
	@$ $(CXX) -c $(CXXFLAGS) -o $@ $<

$(LIB)$(SONAME): $(OBJS) $(LIB)mydict.o
	@echo "Linking $(SONAME):"
	@$ $(LD) $(LDFLAGS) $(OBJS) $(LIB)mydict.o -o $(LIB)$(SONAME) $(SOFLAGS)$(SONAME)

$(LIB)lib%$(LIBSuf): $(PLG)%$(SRCSuf) $(PLG)%$(HDRSuf) $(OBJS)
	@echo "Creating plugin library " $@
	@echo " CXX $<"	
	@$ $(LD) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(SOFLAGS)lib$*.so $(GLIBS)

dynTTree:
	cd DynamicTTree && $(MAKE)

cfgMan:
	cd CfgManager && $(MAKE)

clean:
	@echo "cleaning..."
	rm -f lib/* $(OBJ)*$(OBJSuf) $(LIB)*$(LIBSuf) $(LIB)mydict* $(BIN)*$(BINSuf)
	cd DynamicTTree && $(MAKE) clean
	cd CfgManager && $(MAKE) clean
