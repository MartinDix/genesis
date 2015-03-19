######################
# Makefile for genesis 
#
# make [ all | genesis | clean | new ]
#
######################

FC = ifort
FCFLAGS =-O3 -ip -no-prec-div

LIBDIR = -L${NETCDF_ROOT}/lib/Intel
LIBS = -lnetcdff -lnetcdf
INCDIR = -I${NETCDF_ROOT}/include/Intel


F_SOURCES=establ.f  ngtopt.f
F90_SOURCES=charney_phillips.f90 genesis.f90 gradz.f90 \
            ncread_data.f90 radian.f90 slingo.f90      \
	    create_namelist.f90 geostr.f90 interp_linear.f90 \
	    ncread_dim.f90 relh2q.f90 surfdist.f90 find_latlon.f90 \
	    global.f90 q2relh.f90 search_grid.f90 netcdf_type.f90 \
	    netcdf_check.f90

OBJS=$(F_SOURCES:.f=.o) $(F90_SOURCES:.f90=.o)

default : genesis

# Clear out the suffixes
.SUFFIXES:

# List only suffixes we use
.SUFFIXES: .o .f .f90

ncread_dim.o ncread_data.o genesis.o : netcdf_type.o
ncread_dim.o ncread_data.o : netcdf_check.o

.f.o:
	$(FC) -c $(FCFLAGS) $(LIBDIR) $(INCDIR) $(LIBS) -o $@ $<

.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBDIR) $(INCDIR) $(LIBS) -o $@ $<

$(OBJS) : global.o

#this fails
#MODDIR = -M/nm/keep/local/intel_builds/netcdf-4.0.1-classic/include
#LIBDIR = -L/nm/keep/local/intel_builds/netcdf-4.0.1-classic/lib
#LIBS = -lnetcdf
#INCDIR = -I/nm/keep/local/intel_builds/netcdf-4.0.1-classic/include

#this fails
#FC = f95
#LIBDIR = -L/nm/keep/local/sun_builds/netcdf-4.0.1-classic/lib
#LIBS = -lnetcdf
#INCDIR = -I/nm/keep/local/sun_builds/netcdf-4.0.1-classic/include

#OBJS = ngtopt.o global.o ncread_era.o ncread_jra.o ncread_ncep.o \
#       ncread_um.o establ.o find_latlon.o interp_linear.o \
#       search_grid.o geostr.o surfdist.o \
#       radian.o charney_phillips.o create_namelist.o modulus.o \
#       ncread_dim.o ncread_data.o

genesis: $(OBJS)
	$(FC) ${FCFLAGS} $(LIBDIR) $(LIBS) $(INCDIR) -o $@ $^

geostr.o : radian.o

all: genesis

clean:
	rm -f $(OBJS) genesis.o *.o *.mod

new:
	make clean
	make genesis

.PHONY: all clean new default

###################
# END
