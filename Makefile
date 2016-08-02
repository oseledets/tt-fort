.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

#FOPT    = -m64 -fopenmp #-fdefault-integer-8 #-xSSE4.2
DEPS    =  nan, timef, say, rnd, ptype, sort, trans, ort, mat, check, lr, maxvol, svd, \
	   matrix_util, tt, ttaux, ttop, ttio,  tts, python_conv, putstrmodule, dispmodule, tt_linalg, \
           ttals, tt_eigb, default, ttlocsolve, ttnodeop, ttamen

include Makefile.in

OBJS    = $(DEPS:,=.o).o
MODS    = *.mod
OBJF	= $(OBJS)
OBJC	=
all: mytt.a primme.a

mytt.a : $(OBJS)
	ar rc mytt.a $(OBJS)

primme.a:
	cp Makefile.cpu primme/ && $(MAKE) -C primme


########### This is a test for different compilers, so I didn't put it into "all"
test_eigb_i: mytt.a primme.a
	ifort -O2 test_eigb.f90 mytt.a primme/primme.a -o test_eigb  $(LIB)

test_eigb_g: mytt.a primme.a
	gfortran -O3 test_eigb.f90 mytt.a primme/primme.a -o test_eigb  $(LIB)

test_svd_g: mytt.a 
	gfortran -O3 test_svd.f90 mytt.a -o test_svd  $(LIB)


.f.o:
		$(FC) $(FOPT) -c $<
.f90.o:
		$(FC) $(FOPT) -c $<
.F90.o:
		$(FC) $(FOPT) -c $<
.c.o:
		$(CC) $(COPT) -c $<


clean:
		rm -f $(OBJF) $(OBJC) $(MODS) mytt.a primme/*.o primme/*.mod primme/primme.a

