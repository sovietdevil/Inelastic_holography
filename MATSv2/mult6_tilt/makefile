F90 = ifort
FFLAGS = -O2 -xHost -ipo -qopenmp -static -fp-model precise -traceback
CFLAGS =
#INCLUDE = -I/software/intel/composer_xe_2013.2.146/mkl/include/fftw
INCLUDE = -I/opt/intel/mkl/include/fftw
#LIB = -L/software/apps/fftw/3.3.3/i1214/dp_omp/lib -lfftw3_omp -lfftw3
LIB = -mkl=parallel
RM = rm -f

# clear out all suffixes
.SUFFIXES:
# list only those we use
.SUFFIXES: .o .f90

# define a suffix rule for .f90 -> .o and .c -> .o
.c.o 	:
		$(CC)  $(CFLAGS) $(INCLUDE) -c $<
.f90.o 	:
		$(F90) $(FFLAGS) $(INCLUDE) -c $<

TAGT =  mult6
#In correct order of dependencies
SRCS = modules.f90 misc.f90 output.f90 wf.f90 angm.f90 multislice.f90 3dwf.f90 preprocess.f90 main.f90
OBJS = ${SRCS:.f90=.o} 

all :	$(TAGT)

$(TAGT) : $(OBJS)
	  $(F90) -o $(TAGT) $(OBJS) $(LIB)

clean :
	 $(RM) *.o *.mod $(TAGT)
tarfile : 
	tar cvfz MultiSlice`date "+%d%h%Y"`.tgz *.f90 *.in makefile	
