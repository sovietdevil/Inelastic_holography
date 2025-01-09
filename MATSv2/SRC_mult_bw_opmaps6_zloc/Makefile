###############################################################################
#                                                                             #
#  make              ... generate executable                                  #
#  make clean        ... delete unnecessary files                             #
#                                                                             #
###############################################################################

#------------------------------------------------------------------------------
#
# FC ......... compiler name
# FOPT ....... optimization flags
# FGEN ....... code generation flags (flags not related to optimization)
# LDFLAGS .... linker flags
# LIBS ....... libraries needed to build the executable
# DESTDIR .... directory where the executable should go (incl. trailing '/'!)
# EXECNAME ... name of the resulting executable (without prefixed path!)
#
#------------------------------------------------------------------------------
FC = ifort -Nmkl -i-static -mkl=sequential 
#FC = ifort -Nmkl -mkl=sequential 
FOPT     = -fp-model consistent -traceback -O1
#FOPT     = -fp-model precise -traceback -xHost
#FOPT     = -check -traceback
LDFLAGS = $(FOPT)
LIBS = -I/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64 
DESTDIR  = ./
EXECNAME = dyndif_v15_fe_sf_multbw_opmaps6_zloc
EXESTATIC = SRC.a 

###############################################################################
###############################################################################

FFLAGS = $(FGEN) $(FOPT)
EXEC   = $(DESTDIR)$(EXECNAME)

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS1   = modules.o            \
          deffile.o            \
          theta2dir.o          \
          read_inputs.o        \
	  read_vhklwien.o      \
	  load_multislice.o    \
	  vector_ops.o         \
	  newbeam.o            \
	  generate_beams.o     \
          scafac.o             \
	  fscatt.o             \
	  vhklwien.o           \
          vhkl.o               \
          calc_potentials.o    \
	  calc_ws.o            \
	  select_beams.o       \
	  quicksort.o          \
	  sort_beams.o         \
	  calc_matrices.o      \
	  vabs.o               \
	  calc_etas.o          \
          calc_blochs.o        \
          hash.o               \
	  sort_blochs.o        \
          dyndif.o

#..............................................................................
#
OBJS = $(OBJS1)

#..............................................................................
#
#  Build executable
#
$(EXEC):	$(OBJS)
	$(FC) -o $(EXEC) $(OBJS) $(LDFLAGS) $(LIBS) #$(EXESTATIC)
#..............................................................................
#
#  define rules for generating files needed to check the static semantics
#  of the FORTRAN source code and for generating object files from source files
#
%.o : %.f90
	$(FC) $(FFLAGS) -c $<

%.o : %.f95
	$(FC) $(FFLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
$(OBJS):	Makefile

#..............................................................................
#
#  remove object files, preprocessed source files and files needed for
#  static semantics checking
#
clean:
	rm  -f $(OBJS) *.mod
