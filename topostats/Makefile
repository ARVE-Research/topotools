#makefile

NC_LIB=/share1/netcdf/4.7.1-impi/lib
NC_INC=/share1/netcdf/4.7.1-impi/include

CC = icc
FC = ifort
LDFLAGS  = -L$(NC_LIB) -qopenmp
CPPFLAGS = -I$(NC_INC) -Difort

#FCFLAGS = -debug -fpe0 -traceback -g -check -check noarg_temp_created -warn all -openmp  -openmp-report
FCFLAGS  = -xHost -qopenmp #-xHost -vec-report -openmp  -openmp-report -warn all
LIBS     = -lnetcdff -lnetcdf #-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -loverprint


OBJS = overprint.o			\
			 parametersmod.o  \
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
