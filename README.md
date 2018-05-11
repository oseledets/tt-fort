tt-fort
=======

Fortran computing core of the TT-Toolbox Right now, it is not easy to use as a
standalone project.

It also includes two slightly customized libraries:

1. [PRIMME](http://www.cs.wm.edu/~andreas/software/) code
2. [Expokit](http://www.maths.uq.edu.au/expokit/) code

for doing fast local solvers in eigenvalue and dynamical problems.

Compilation
===========

1. Create file Makefile.cpu in the currect directory with the following line

   ```makefile
       CPU = your_preset
   ```

   where `your_preset` is one from Makefile.in (e.g. `CPU = i4-intel`). You
   may start from the example file [Makefile.cpu.default](Makefile.cpu.default)
   provided.

       For x64 systems!
       Note that the standalone programs (e.g. test_eigb below) require
       four-bytes integers without fPIC (i4-gnu-nopic or i4-intel-nopic), while
       the MATLAB MEX libraries require 8-bytes integers with fPIC, i8-gnu or
       i8-intel.

2. Run make

   It will compile the static libraries `mytt.a` and `primme/primme.a`

3. Use the libs in your code =).

   To quick start, you may check the standalone test program for the block
   eigenvalue solver `test_eigb` (see [Makefile](Makefile)). It computes 5
   lowest eigenpairs of the Heisenberg chain with 20 spin-1/2 particles.
