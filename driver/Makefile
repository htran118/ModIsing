# simple makefile

#a lot of flags
FC = gfortran
FCFLAGS = -g -fbounds-check
FCOPT = -O3
FCINCLUDE += -I/usr/include
#LDFLAGS = (someLibraries)

FCCOMPLETE = $(FC) $(FCFLAGS) $(FCOPT)

PROGRAMS = isingModel isingAnal dataAnal 

MODULE = ../module/
#compiling commands
const.o: const.f
	$(FCCOMPLETE) -c const.f
randomSeed.o: $(MODULE)randomSeed.f
	$(FCCOMPLETE) -c $(MODULE)randomSeed.f
crysProc.o: $(MODULE)crysProc.f const.o
	$(FCCOMPLETE) -c $(MODULE)crysProc.f const.f
statProc.o: $(MODULE)statProc.f crysProc.o randomSeed.o const.o
	$(FCCOMPLETE) -c $(MODULE)statProc.f $(MODULE)crysProc.f $(MODULE)randomSeed.f const.f
dataProc.o: $(MODULE)dataProc.f const.o
	$(FCCOMPLETE) -c $(MODULE)dataProc.f const.f
isingModel.o: isingModel.f crysProc.o statProc.o dataProc.o randomSeed.o const.o
	$(FCCOMPLETE) -c isingModel.f $(MODULE)crysProc.f $(MODULE)statProc.f $(MODULE)dataProc.f $(MODULE)randomSeed.f const.f
isingModel: isingModel.o crysProc.o statProc.o dataProc.o randomSeed.o const.o
	$(FCCOMPLETE) $? -o $@
isingAnal.o: isingAnal.f crysProc.o dataProc.o const.o
	$(FCCOMPLETE) -c isingAnal.f $(MODULE)crysProc.f $(MODULE)dataProc.f const.f
isingAnal: isingAnal.o crysProc.o dataProc.o const.o
	$(FCCOMPLETE) $? -o $@
dataAnal.o: dataAnal.f dataProc.o const.o
	$(FCCOMPLETE) -c dataAnal.f $(MODULE)dataProc.f const.f
dataAnal: dataAnal.o dataProc.o const.o
	$(FCCOMPLETE) $? -o $@

#general commands
.o:
	$(FCCOMPLETE) $@ -o $@

#generic commands
all: $(PROGRAMS)
clean:
	rm -f *.o *.mod *.MOD *.exe *.stackdump
clear:
	rm -f *.o *.mod *.MOD
