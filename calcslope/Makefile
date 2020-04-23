#makefile

NC_LIB=/usr/local/netcdf/lib
NC_INC=/usr/local/netcdf/include

CC = icc
FC = ifort
LDFLAGS  = -L$(NC_LIB) -openmp
CPPFLAGS = -I$(NC_INC)

#FCFLAGS = -debug -fpe0 -traceback -g -check -check noarg_temp_created -warn all -openmp  -openmp-report
FCFLAGS  = -xHost -vec-report -openmp  -openmp-report
#FCFLAGS = -ftree-vectorize -msse2 -ftree-vectorizer-verbose=2 #-fbacktrace -ggdb -g -fcheck=all -ffpe-trap=invalid,zero
LIBS     =  -lnetcdff -lnetcdf -loverprint

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
