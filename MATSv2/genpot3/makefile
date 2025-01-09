F90 = ifort
FFLAGS = -O2 -xHost -ipo -qopenmp -static
CFLAGS =
INCLUDE = -I/opt/intel/mkl/include/fftw
#INCLUDE = -I/software/intel/composer_xe_2013.2.146/mkl/include/fftw
#LIB = -L/software/apps/fftw/3.3.3/i1214/dp_omp/lib -lfftw3_omp -lfftw3
LIB = -Nmkl -mkl=parallel
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

TAGT = genpot 
#In correct order of dependencies
SRCS = fftw3.f90 modules.f90 vector_ops.f90 fscatt.f scafac.f90 vhkl.f90 main.f90
OBJS = ${SRCS:.f90=.o} 

all :	$(TAGT)

$(TAGT) : $(OBJS)
	  $(F90) -o $(TAGT) $(OBJS) $(LIB)

clean :
	 $(RM) *.o *.mod
tarfile : 
	tar cvfz MultiSlice.tgz *.f *.f90 plot.gnu makefile	
