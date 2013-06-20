.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

include Makefile.deps

OBJS    = $(DEPS:,=.o).o

FC = ifort
CC = icc

FOPT = -O3
COPT = -O3

primme.a: $(OBJS)
	ar rc primme.a $(OBJS)




.f.o:
		$(FC) $(FOPT) -c $<
.f90.o:
		$(FC) $(FOPT) -c $<
.F90.o:
		$(FC) $(FOPT) -c $<
.c.o:
		$(CC) $(COPT) -c $<


clean:
		rm -f *.o *.mod primme.a
