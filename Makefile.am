.SUFFIXES = .f
.f.o:
	$(F77) $(FFLAGS) -c $<

# ls *.f | awk '{ printf("%s ", $1) }' > tmp and paste the contents of tmp to XXX_SOURCES
# INCLUDES = -I../basic
# lib_LIBRARIES = util.a
# include_HEADERS = util.h
# util_a_SOURCES = util.f
# check_PROGRAMS = test1 test2 test3
# TESTS = $(check_PROGRAMS)
# LDADD = liblin.a ../blas1/libblas1.a ../mathutil/libmathutil.a $(FLIBS) -lm
# test1_SOURCES = test1.cc f2c-main.cc
# test2_SOURCES = test2.cc f2c-main.cc  
# test3_SOURCES = test3.cc f2c-main.cc
bin_PROGRAMS = topmon
topmon_SOURCES = Atom_energy.f Atom_traxyz.f Intra_energy.f alignplanes.f analysis.f anes.f atomtype.f boltz.f bondlength.f calcsolpar.f calctor.f charge.f chempt.f cleanup.f close.f cone.f coneangle.f config.f corp.f coru.f coruz.f ctrmas.f deriv.f dipole.f dump.f ee_energy.f ee_index_swap.f ee_recip.f eemove.f eesetup.f energy.f erfunc.f exgrph.f exmuir.f expand.f explct.f exsami.f exsix.f exzeo.f exzeof.f exzeot.f flucq.f forcefield.f garo3body.f garofalini.f genlj.f geometry.f grid.f inclus.f init_vars.f initia.f ipswsetup.f jtrim.f linkcell.f ljmuir.f ljpsur.f ljsami.f matops.f mbessel.f mimage.f mmff.f molsetup.f monola.f monper.f newton.f ninesix.f nrtab.f pgrid.f place.f pressure.f prvolume.f qqcheck.f random.f readRES.f readdat.f recip.f recip_atom.f recippress.f recipsum.f rigfix.f rigrot.f rosenbluth.f rotxyz.f safecbmc.f safeschedule.f schedule.f setnn.f setpbc.f slitpore.f sphere.f suijtab.f sumup.f susami.f suvibe.f suzeo.f swap.f swatch.f tabulated.f test1.f testupdate.f thermopress.f topmon.f torsion.f traxyz.f trim_swap.f update.f update1.f updnn.f v_elect_field.f virial.f volume.f vtorso.f writepdb.f zeocoord.f ztest.f 
EXTRA_DIST = reconf configure
# SUBDIRS = MANUAL TEST_SUITE MPI_SCRIPTS MPI_EXTRA

#########################
### Compiler definitions
#########################

### Intel Compiler for x86 architecture
#AM_FFLAGS      =  -w90  -check all -traceback 
#AM_FFLAGS     =     -w90 -w95 -cm -check all -traceback#-check_bounds# -g
#AM_FFLAGS     =  -O2 -w90 -w95 -cm -p#-check_bounds# -g
## only for debugging
#AM_FFLAGS     = -w90 -w95 -cm -CB -CU -g

#AM_FFLAGS     =  -O2 -w90 -w95 #-cm -p
#AM_FFLAGS     =  -ipo -O3 -no-prec-div -xHost 
#AM_FFLAGS     =  -O3 -no-prec-div -check bounds 
#LD            = mpif90 
#LDS           = mpif90 -c $(AM_FFLAGS)
#FLINK       =  $(AM_FFLAGS)

### NAG Compiler for OSX
#AM_FFLAGS     = -O0 -gline -f77 -mismatch_all -C=all -C=undefined -nan
## only for debugging
#AM_FFLAGS     = -O2 
#LD            = f95
#LDS           = f95 -c $(AM_FFLAGS)
#FLINK       = $(AM_FFLAGS)

### SGI Altix: Intel Compiler for ia64 architecture
#AM_FFLAGS     = -O2 -tpp2 -w90 -w95 -cm #-check all -traceback
# only for debugging
#AM_FFLAGS     =  -check all -traceback -w90 -w95 -cm -CB -CU -g
#LD            = ifort
#LDS           = ifort -c $(AM_FFLAGS)
#FLINK       = -static $(AM_FFLAGS)

### GNU Compiler
#AM_FFLAGS     =  -ftrace=full -pg -fbounds-check
AM_FFLAGS     = -O3 -fomit-frame-pointer -finline-functions -ffast-math 
#LD	      = g95
#LDS 	      = g95 -c $(AM_FFLAGS)
#FLINK	      = $(AM_FFLAGS)

### Portland Compiler for x86 architecture
#AM_FFLAGS     =  -O2 -tp k8-64 -Mvect=sse#-check_bounds# -g
#AM_FFLAGS     =  -tp k8-64 -Mbounds -g -Mchkptr#-check_bounds# -g
#LD	      = /opt/pgi/linux86-64/6.1/bin/pgf95
#LDS 	      = /opt/pgi/linux86-64/6.1/bin/pgf95 -c $(AM_FFLAGS)
#FLINK	      = $(AM_FFLAGS)

### SGI Compiler for Origin
#AM_FFLAGS     = -O2 #-march=opteron # -mips4 -r12000 #-check_bounds# -g
#LD            = f77 
#LDS           = f77 -c $(AM_FFLAGS)
#FLINK       = $(AM_FFLAGS)

### IBM AIX Compiler for SP
#AM_FFLAGS     = -O3 -qstrict  -bmaxdata:2000000000 #-check_bounds# -g
#LD            = xlf
#LDS           = xlf -c $(AM_FFLAGS)
#FLINK	      = $(AM_FFLAGS)

### GNU Compiler
#AM_FFLAGS     =  -fbounds-check -g
#AM_FFLAGS     = -O3 #-fomit-frame-pointer -finline-functions -ffast-math #-fbounds-check #-mcpu=i686 -march=i686# 
#LD	      = g77
#LDS 	      = g77 -c $(AM_FFLAGS)
#FLINK	      = $(AM_FFLAGS)

#############################
### End compiler definitions 
#############################
