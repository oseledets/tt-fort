.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

include Makefile.deps
include Makefile.in

OBJS    = $(DEPS:,=.o).o


primme.a: setints $(OBJS)
	ar rc primme.a $(OBJS)




.f.o:
		$(FC) $(FOPT) -c $<
.f90.o:
		$(FC) $(FOPT) -c $<
.F90.o:
		$(FC) $(FOPT) -c $<
.c.o:
		$(CC) $(COPT) -c $<

setints:
	echo "$(INTDEF)#define int long" > int_redefine.h


clean:
		rm -f *.o *.mod primme.a
