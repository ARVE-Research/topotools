#makefile

NC_LIB=/usr/local/netcdf/lib
NC_INC=/usr/local/netcdf/include

CC = icc
FC = ifort
LDFLAGS  = -L$(NC_LIB) -openmp
CPPFLAGS = -I$(NC_INC)

#FCFLAGS = -debug -fpe0 -traceback -g -check -check noarg_temp_created -warn all -openmp  -openmp-report
FCFLAGS  = -xHost -vec-report -openmp  -openmp-report -warn all
LIBS     = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lnetcdff -lnetcdf -loverprint


OBJS = parametersmod.o  \
       netcdfmod.o      \
       sortmod.o        \
       utilitiesmod.o   \
       statsmod.o       \
       calcstatsmod.o   \
       coordsmod.o      \
       outputmod.o      \
       topostats.o

.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	topostats

topostats:  $(OBJS)
	$(FC) -o topostats $(OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm topostats *.mod *.o




#FCFLAGS = -ftree-vectorize -msse2 -ftree-vectorizer-verbose=2 #-fbacktrace -ggdb -g -fcheck=all -ffpe-trap=invalid,zero