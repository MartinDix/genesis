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
	    global.f90 modulus.f90 q2relh.f90 search_grid.f90

OBJS=$(F_SOURCES:.f=.o) $(F90_SOURCES:.f90=.o)

default : genesis
all : genesis

# Clear out the suffixes
.SUFFIXES:

# List only suffixes we use
.SUFFIXES: .o .f .f90

.f.o:
	$(FC) -c $(FCFLAGS) $(LIBDIR) $(INCDIR) $(LIBS) -o $@ $<

.f90.o:
	$(FC) -c $(FCFLAGS) $(LIBDIR) $(INCDIR) $(LIBS) -o $@ $<

$(OBJS) : global.o

genesis: $(OBJS)
	$(FC) ${FCFLAGS} $(LIBDIR) $(LIBS) $(INCDIR) -o $@ $^

all: genesis

clean:
	rm -f $(OBJS) genesis.o *.o *.mod

new:
	make clean
	make genesis

.PHONY: all clean new default

###################
# END
