# Makefile created by mkmf.pl $Id: mkmf,v 14.0 2007/03/20 22:13:27 fms Exp $ 

SRCROOT = .
POLYPACK = ../dependencies/polypack

FC = gfortran-mp-4.8
LD = gfortran-mp-4.8
#FFLAGS = -O2
FFLAGS = -O3 -ftree-vectorize
LDFLAGS = -lgfortran -framework Accelerate

.DEFAULT:

all: modaem-1.8

a_grid.o: $(SRCROOT)/a_grid.F90 u_constants.o u_io.o m_aem.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/a_grid.F90
a_obs.o: $(SRCROOT)/a_obs.F90 u_constants.o u_io.o m_aem.o m_ls2.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/a_obs.F90
a_stdio.o: $(SRCROOT)/a_stdio.F90 u_constants.o u_io.o m_aem.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/a_stdio.F90
a_tr0.o: $(SRCROOT)/a_tr0.F90 u_constants.o u_io.o m_wl0.o m_ls0.o m_ls1.o m_hb0.o m_wl1.o m_aem.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/a_tr0.F90
f_dipole.o: $(SRCROOT)/f_dipole.F90 u_constants.o u_io.o i_dipole.o i_linesink.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/f_dipole.F90
f_pond.o: $(SRCROOT)/f_pond.F90 u_constants.o u_io.o i_pond.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/f_pond.F90
f_well.o: $(SRCROOT)/f_well.F90 u_constants.o u_io.o i_well.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/f_well.F90
f_bwl.o: $(SRCROOT)/f_bwl.F90 u_constants.o u_io.o i_well.o u_math.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/f_bwl.F90
i_areasink.o: $(SRCROOT)/i_areasink.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/i_areasink.F90
i_dipole.o: $(SRCROOT)/i_dipole.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/i_dipole.F90
i_linesink.o: $(SRCROOT)/i_linesink.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/i_linesink.F90
i_pond.o: $(SRCROOT)/i_pond.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/i_pond.F90
i_well.o: $(SRCROOT)/i_well.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/i_well.F90
m_aem.o: $(SRCROOT)/m_aem.F90 u_constants.o u_io.o u_matrix.o f_well.o f_dipole.o f_pond.o m_aqu.o m_wl0.o m_wl1.o \
                              m_pd0.o m_ls0.o m_ls1.o m_ls2.o m_ls3.o m_hb0.o m_as0.o m_cw0.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_aem.F90
m_aqu.o: $(SRCROOT)/m_aqu.F90 u_constants.o u_io.o u_matrix.o f_dipole.o i_linesink.o m_in0.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_aqu.F90
m_as0.o: $(SRCROOT)/m_as0.F90 u_constants.o u_io.o u_polygon.o u_matrix.o f_well.o f_dipole.o i_areasink.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_as0.F90
m_cw0.o: $(SRCROOT)/m_cw0.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_cw0.F90
m_ext.o: $(SRCROOT)/m_ext.F90 u_constants.o u_io.o m_wl0.o m_wl1.o m_ls0.o m_ls1.o m_hb0.o m_aem.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_ext.F90
m_hb0.o: $(SRCROOT)/m_hb0.F90 u_constants.o u_io.o u_matrix.o f_dipole.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_hb0.F90
m_in0.o: $(SRCROOT)/m_in0.F90 u_constants.o u_io.o u_matrix.o u_polygon.o f_dipole.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_in0.F90
m_inq.o: $(SRCROOT)/m_inq.F90 u_constants.o u_io.o m_wl0.o m_wl1.o m_ls0.o m_ls1.o m_hb0.o m_aem.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_inq.F90
m_ls0.o: $(SRCROOT)/m_ls0.F90 u_constants.o u_io.o u_matrix.o f_well.o f_dipole.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_ls0.F90
m_ls1.o: $(SRCROOT)/m_ls1.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_ls1.F90
m_ls2.o: $(SRCROOT)/m_ls2.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_ls2.F90
m_ls3.o: $(SRCROOT)/m_ls3.F90 u_constants.o u_io.o f_well.o f_dipole.o u_matrix.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_ls3.F90
m_pd0.o: $(SRCROOT)/m_pd0.F90 u_constants.o u_io.o f_pond.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_pd0.F90
m_wl0.o: $(SRCROOT)/m_wl0.F90 u_constants.o u_io.o f_well.o f_bwl.o m_aqu.o 
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_wl0.F90
m_wl1.o: $(SRCROOT)/m_wl1.F90 u_constants.o u_io.o f_well.o u_matrix.o m_aqu.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/m_wl1.F90
modaem.o: $(SRCROOT)/modaem.F90 u_constants.o u_io.o m_aqu.o m_wl0.o m_wl1.o m_ls0.o m_ls1.o m_hb0.o m_aem.o m_inq.o m_ext.o a_grid.o a_tr0.o a_stdio.o a_obs.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/modaem.F90
u_constants.o: $(SRCROOT)/u_constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/u_constants.F90
u_grid.o: $(SRCROOT)/u_grid.F90 u_constants.o u_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/u_grid.F90
u_io.o: $(SRCROOT)/u_io.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/u_io.F90
u_math.o: $(SRCROOT)/u_math.F90 u_constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/u_math.F90
u_matrix.o: $(SRCROOT)/u_matrix.F90 u_constants.o u_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/u_matrix.F90
u_polygon.o: $(SRCROOT)/u_polygon.F90 u_constants.o polypack.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(SRCROOT)/u_polygon.F90
polypack.o: $(POLYPACK)/polypack.f
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	$(POLYPACK)/polypack.f

SRC = $(SRCROOT)/a_tr0.F90 $(SRCROOT)/m_as0.F90 $(SRCROOT)/m_wl1.F90 $(SRCROOT)/m_ls0.F90 \
      $(SRCROOT)/u_matrix.F90 $(SRCROOT)/u_polygon.F90 $(SRCROOT)/u_math.F90 $(SRCROOT)/modaem.F90 \
      $(SRCROOT)/m_wl0.F90 $(SRCROOT)/a_stdio.F90 $(SRCROOT)/f_pond.F90 $(SRCROOT)/a_obs.F90 \
      $(SRCROOT)/m_cw0.F90 $(SRCROOT)/m_ls1.F90 $(SRCROOT)/i_areasink.F90 $(SRCROOT)/a_grid.F90 \
      $(SRCROOT)/m_in0.F90 $(SRCROOT)/m_aqu.F90 $(SRCROOT)/i_pond.F90 $(SRCROOT)/m_pd0.F90 \
      $(SRCROOT)/f_dipole.F90 $(SRCROOT)/i_dipole.F90 $(SRCROOT)/i_well.F90 $(SRCROOT)/m_ext.F90 \
      $(SRCROOT)/m_aem.F90 $(SRCROOT)/u_io.F90 $(SRCROOT)/m_hb0.F90 $(SRCROOT)/u_constants.F90 \
      $(SRCROOT)/i_linesink.F90 $(SRCROOT)/f_well.F90 $(SRCROOT)/m_inq.F90 $(SRCROOT)/m_ls2.F90 \
      $(SRCROOT)/m_ls3.F90 $(SRCROOT)/f_bwl.F90 $(SRCROOT)/u_grid.F90 $(POLYPACK)/u_polygon.F90
OBJ = a_tr0.o m_as0.o m_wl1.o m_ls0.o u_matrix.o u_polygon.o u_math.o modaem.o m_wl0.o a_stdio.o \
      f_pond.o a_obs.o m_cw0.o m_ls1.o i_areasink.o a_grid.o m_in0.o m_aqu.o i_pond.o m_pd0.o \
      f_dipole.o i_dipole.o i_well.o m_ext.o m_aem.o u_io.o m_hb0.o u_constants.o i_linesink.o \
      f_well.o m_inq.o m_ls2.o u_grid.o m_ls3.o f_bwl.o polypack.o
OFF = $(SRCROOT)/m_ext.F90 $(SRCROOT)/i_linesink.F90 $(SRCROOT)/modaem.F90 $(SRCROOT)/i_pond.F90 \
      $(SRCROOT)/m_wl0.F90 $(SRCROOT)/i_well.F90 $(SRCROOT)/m_aqu.F90 $(SRCROOT)/u_math.F90 \
      $(SRCROOT)/u_constants.F90 $(SRCROOT)/m_cw0.F90 $(SRCROOT)/f_well.F90 $(SRCROOT)/u_matrix.F90 \
      $(SRCROOT)/m_as0.F90 $(SRCROOT)/u_io.F90 $(SRCROOT)/u_grid.F90 $(SRCROOT)/u_polygon.F90 \
      $(SRCROOT)/m_ls0.F90 $(SRCROOT)/m_inq.F90 $(SRCROOT)/f_dipole.F90 $(SRCROOT)/a_obs.F90 \
      $(SRCROOT)/i_dipole.F90 $(SRCROOT)/i_areasink.F90 $(SRCROOT)/a_grid.F90 $(SRCROOT)/m_aem.F90 \
      $(SRCROOT)/m_ls1.F90 $(SRCROOT)/m_wl1.F90 $(SRCROOT)/a_tr0.F90 $(SRCROOT)/f_pond.F90 \
      $(SRCROOT)/m_ls2.F90 $(SRCROOT)/a_stdio.F90 $(SRCROOT)/m_hb0.F90 $(SRCROOT)/m_in0.F90 \
      $(SRCROOT)/m_pd0.F90 $(SRCROOT)/m_ls3.F90 $(SRCROOT)/f_bwl.F90 $(SRCROOT)/u_grid.F90 \
      $(POLYPACK)/u_polygon.F90

clean: neat
	-rm -f .cppdefs $(OBJ) *.mod modaem-1.8

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRCROOT)
	etags $(SRCROOT)

tags: $(SRCROOT)
	ctags $(SRCROOT)

modaem-1.8: $(OBJ) 
	$(LD) $(OBJ) -o modaem-1.8  $(LDFLAGS)

install: modaem-1.8
	cp -f modaem-1.8 ~/bin/modaem-1.8

