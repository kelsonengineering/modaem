# Makefile created by mkmf.pl $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

SRCROOT = ../../src

FC = gfortran-mp-4.8
LD = gfortran-mp-4.8
FFLAGS = -O3 -ftree-vectorize -g
LDFLAGS = -lgfortran -framework Accelerate

vpath %f ../../src
vpath %F90 ../../src

.DEFAULT:

all: modaem-1.8

%o: %f
	$(FC) -c $(FFLAGS) -o $@ $<

%o: %F90
	$(FC) -c $(FFLAGS) -o $@ $<
	
a_grid.o: a_grid.F90 u_constants.o u_io.o m_aem.o m_aqu.o

a_obs.o: a_obs.F90 u_constants.o u_io.o m_aem.o m_ls2.o

a_stdio.o: a_stdio.F90 u_constants.o u_io.o m_aem.o

a_tr0.o: a_tr0.F90 u_constants.o u_io.o m_wl0.o m_ls0.o m_ls1.o m_hb0.o m_wl1.o m_aem.o

f_dipole.o: f_dipole.F90 u_constants.o u_io.o i_dipole.o i_linesink.o

f_pond.o: f_pond.F90 u_constants.o u_io.o i_pond.o

f_well.o: f_well.F90 u_constants.o u_io.o i_well.o

f_bwl.o: f_bwl.F90 u_constants.o u_io.o i_well.o u_math.o

i_areasink.o: i_areasink.F90 u_constants.o

i_dipole.o: i_dipole.F90 u_constants.o

i_linesink.o: i_linesink.F90 u_constants.o

i_pond.o: i_pond.F90 u_constants.o

i_well.o: i_well.F90 u_constants.o

m_aem.o: m_aem.F90 u_constants.o u_io.o u_matrix.o f_well.o f_dipole.o f_pond.o m_aqu.o m_wl0.o m_wl1.o \
                              m_pd0.o m_ls0.o m_ls1.o m_ls2.o m_ls3.o m_hb0.o m_as0.o

m_aqu.o: m_aqu.F90 u_constants.o u_io.o u_matrix.o f_dipole.o i_linesink.o m_in0.o

m_as0.o: m_as0.F90 u_constants.o u_io.o u_polygon.o u_matrix.o f_well.o f_dipole.o i_areasink.o m_aqu.o

m_ext.o: m_ext.F90 u_constants.o u_io.o m_wl0.o m_wl1.o m_ls0.o m_ls1.o m_hb0.o m_aem.o

m_hb0.o: m_hb0.F90 u_constants.o u_io.o u_matrix.o f_dipole.o

m_in0.o: m_in0.F90 u_constants.o u_io.o u_matrix.o u_polygon.o f_dipole.o

m_inq.o: m_inq.F90 u_constants.o u_io.o m_wl0.o m_wl1.o m_ls0.o m_ls1.o m_hb0.o m_aem.o

m_ls0.o: m_ls0.F90 u_constants.o u_io.o u_matrix.o f_well.o f_dipole.o m_aqu.o

m_ls1.o: m_ls1.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o

m_ls2.o: m_ls2.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o

m_ls3.o: m_ls3.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o

m_pd0.o: m_pd0.F90 u_constants.o u_io.o f_pond.o

m_wl0.o: m_wl0.F90 u_constants.o u_io.o f_well.o f_bwl.o m_aqu.o 

m_wl1.o: m_wl1.F90 u_constants.o u_io.o f_well.o u_matrix.o m_aqu.o

modaem.o: modaem.F90 u_constants.o u_io.o m_aqu.o m_wl0.o m_wl1.o m_ls0.o m_ls1.o m_hb0.o m_aem.o m_inq.o m_ext.o a_grid.o a_tr0.o a_stdio.o a_obs.o

u_constants.o: u_constants.F90

u_grid.o: u_grid.F90 u_constants.o u_io.o

u_io.o: u_io.F90 u_constants.o

u_math.o: u_math.F90 u_constants.o

u_matrix.o: u_matrix.F90 u_constants.o u_io.o

u_polygon.o: u_polygon.F90 u_constants.o polypack.o

polypack.o: polypack.f

test.o: test.F90 $(INF)

INF = i_well.o i_linesink.o i_dipole.o i_pond.o i_areasink.o

OBJ = a_tr0.o m_as0.o m_wl1.o m_ls0.o u_matrix.o u_polygon.o u_math.o modaem.o m_wl0.o a_stdio.o \
      f_pond.o a_obs.o m_ls1.o a_grid.o m_in0.o m_aqu.o m_pd0.o \
      f_dipole.o m_ext.o m_aem.o u_io.o m_hb0.o u_constants.o  \
      f_well.o m_inq.o m_ls2.o u_grid.o m_ls3.o f_bwl.o polypack.o $(INF)

clean: neat
	-rm -f .cppdefs *.o *.mod

neat:
	-rm -f $(TMPFILES)
	
modaem-1.8: $(OBJ) 
	$(LD) $(OBJ) -o modaem-1.8	  $(LDFLAGS)

install: modaem-1.8
	cp -f modaem-1.8 ~/bin/modaem-1.8

