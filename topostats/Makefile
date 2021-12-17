# makefile

NC_LIB=/home/public/easybuild/software/netCDF-Fortran/4.5.3-gompi-2020a/lib
NC_INC=/home/public/easybuild/software/netCDF-Fortran/4.5.3-gompi-2020a/include

CC = gcc
FC = gfortran
LDFLAGS  = -fopenmp
# CPPFLAGS = -Difort

LIBS     = -lnetcdff -lnetcdf #-loverprint #-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
FCFLAGS  = -O3 -fopenmp -lnetcdff -lnetcdf -L$(NC_LIB) -I$(NC_INC) -pg #-xHost -vec-report -openmp  -openmp-report -warn all

# FCFLAGS = -debug -fpe0 -traceback -g -check -check noarg_temp_created -warn all #-qopenmp
# FCFLAGS = -ftree-vectorize -msse2 -ftree-vectorizer-verbose=2 #-fbacktrace -ggdb -g -fcheck=all -ffpe-trap=invalid,zero

OBJS = overprint.o			\
			 parametersmod.o  \
       netcdfmod.o      \
       sortmod.o        \
       utilitiesmod.o   \
       statsmod.o       \
       calcstatsmod.o   \
       coordsmod.o      \
       outputmod.o      \
			 calcfracmod.o    \
       topostats.o


.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c $(*F).c $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	topostats

topostats:  $(OBJS)
	$(FC) -lc -o topostats $(OBJS) $(LDFLAGS) $(LIBS)

clean::
	-rm topostats *.mod *.o
