#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libutils.a

# complier and options
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

# directories
flib = ../

# compile links
MODD = -I$(flib)/mod
OPTION = $(MODD)

#//// files to be compiled ////#
OBJ = constants.o funcs.o general.o cosmofunc.o bstool.o 

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flib}/mod; mv $(target) ${flib}/lib/

clean:
	rm -f *.o *.mod *.a


