#makefile

NC_LIB=/share1/netcdf/4.7.1-impi/lib
NC_INC=/share1/netcdf/4.7.1-impi/include

CC = icc
FC = ifort
LDFLAGS  = -L$(NC_LIB) -qopenmp
CPPFLAGS = -I$(NC_INC) -Difort

# debugging flags gfortran
#FCFLAGS = -ftree-vectorize -msse2 -ftree-vectorizer-verbose=2 #-fbacktrace -ggdb -g -fcheck=all -ffpe-trap=invalid,zero

# optimized flags for gfortran
#FCFLAGS = -debug -fpe0 -traceback -g -check -check noarg_temp_created -warn all -openmp  -openmp-report

# optimized flags for intel ifort
FCFLAGS  = -xHost -qopenmp #-qopt-report

LIBS     =  -lnetcdff -lnetcdf

OBJS = overprint.o       \
       parametersmod.o   \
       netcdfmod.o       \
       utilitiesmod.o    \
       distmod.o         \
       coordsmod.o       \
       outputmod.o       \
       calcslope_world.o

.SUFFIXES: .o .f90 .f .mod

%.o : %.c
	$(CC) $(CFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

%.o : %.f90
	$(FC) $(FCFLAGS) -c -o $(*F).o $(CPPFLAGS) $<

all::	calcslope

calcslope:  $(OBJS)
	$(FC) -o calcslope $(OBJS) $(LDFLAGS) $(LIBS)

clean::	
	-rm calcslope *.mod *.o
