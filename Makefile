.SUFFIXES:
.SUFFIXES: .c .f .f90 .F90 .o

#FOPT    = -m64 -fopenmp #-fdefault-integer-8 #-xSSE4.2
DEPS    =  nan, timef, say, rnd, ptype, sort, trans, ort, mat, check, lr, maxvol, svd, \
	 matrix_util, tt, ttaux, ttop, ttio,  tts, python_conv, putstrmodule, dispmodule, tt_linalg, tt_adapt_als, gmres_3d
          # d3, mimic, d3als, d3kryl, d3op, d3test, d3elp,

include Makefile.in

OBJS    = $(DEPS:,=.o).o
MODS    = *.mod
OBJF	= $(OBJS)
OBJC	=


mytt.a : $(OBJS)
	ar rc mytt.a $(OBJS)




.f.o:
		$(FC) $(FOPT) -c $<
.f90.o:
		$(FC) $(FOPT) -c $<
.F90.o:
		$(FC) $(FOPT) -c $<
.c.o:
		$(CC) $(COPT) -c $<


clean:
		rm -f $(OBJF) $(OBJC) $(MODS) mytt.a
