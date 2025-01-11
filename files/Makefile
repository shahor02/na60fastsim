# for C++ define  CC = g++

#to compile w/o aliroot support call with `make DOPTS=-D_NOALIROOT_`
#additionally, to compile it under root6 use `make DOPTS=-D_NOALIROOT_ C17=1`
#or `make DOPTS=-D_NOALIROOT_ C20=1`

CC = g++

ifdef C20
$(info compiling with c++20 support )
CFLAGS = -g -O0 -Wall -fPIC -m64 -std=gnu++20 $(DOPTS)
else
ifdef C17
$(info compiling with c++17 support )
CFLAGS = -g -O0 -Wall -fPIC -m64 -std=gnu++17 $(DOPTS)
else
$(info compiling with c++14 support )
CFLAGS = -g -O0 -Wall -fPIC -m64 -std=gnu++14 $(DOPTS)
endif
endif

ifdef DOPTS
$(info additional options:  $(DOPTS) )
endif

LFLAGS = -L$(ROOTSYS)/lib 
INC =	-I$(ROOTSYS)/include -I$(ALICE_ROOT)/include -I./
TGT =	libKMCDetFwd.so
DICT=	KMCDetFwdDict.cxx
DICTO=	KMCDetFwdDict.o

SRC = 	GenMUONLMR.cxx  KMCClusterFwd.cxx  KMCDetectorFwd.cxx  KMCFlukaParser.cxx \
	KMCLayerFwd.cxx  KMCProbeFwd.cxx  KMCUtils.cxx NaMaterial.cxx \
        KMCLayerRect.cxx KMCPolyLayer.cxx KMCMagnetBuilder.cxx \
        KMCMSStation.cxx \
        KMCVTCoolingPlane.cxx KMCPixelPlane.cxx \
        LocLog.cxx TrackPar.cxx TLocTreeStream.cxx

HDR =	$(SRC:.cxx=.h) 

OBJ = 	$(SRC:.cxx=.o)


.PHONY: depend clean

all: 	$(TGT)
	@echo creating libKMCDetFwd.so

$(TGT):	$(OBJ) $(DICTO)
	$(CC) $(CFLAGS)  -shared -o $(TGT) $(OBJ) $(DICTO) `root-config --ldflags` $(LFLAGS)

# pull in dependency info for *existing* .o files
-include $(OBJ:.o=.d)

%.o : %.cxx
	$(CC) $(CFLAGS) $(INC) -c $<  -o $@
	$(CC) -MM $(CFLAGS) $(INC) -c $*.cxx > $*.d
	@cp -f $*.d $*.d.tmp
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

clean:
	rm -f *.o *~ *.so *.pcm *.d *Dict.{h,cxx}

$(DICT): $(HDR) KMCDetFwfLinkDef.h
	rootcling -f $@ $(DOPTS) $(INC) $(HDR) $^


depend: $(SRC)
	makedepend $(INC) $^

# DO NOT DELETE THIS LINE -- make depend needs it
