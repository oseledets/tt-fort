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

There are two ways to assembly project. The obsolete one is based on naked
`make`. This approach has some drawbacks like lack of configurability, reusing
as submodule of other project, and not cross-platform. On other hand `cmake`
and its backend paricularly solve that issues.

Barely `make` assembly will be deprecated in the future.

Compilation with CMake
----------------------

In order to compile TT library one could run the following shell commands.

```bash
    mkdir build && cd build
    cmake ../ \
        -DCMAKE_BUILD_TYPE=RelWithDebInfo \
        -DCMAKE_C_COMPILER=gcc \
        -DCMAKE_Fortran_COMPILER=ifort \
        -DTT_DEFAULT_INTEGER=8
    make
```

It is reasonably tunable to support for different compilers and various build
types. Some options are described bellow.

- Option `CMAKE_BUILD_TYPE`
- Option `CMAKE_C_COMPILER`
- Option `CMAKE_Fortran_COMPILER`
- Option `TT_DEFAULT_INTEGER`
- Option `WITH_BLAS`
- Option `WITH_LAPACK`
- Option `WITH_MKL`

Linking issues are usually related to library search directories. See docs of
`LD\_LIBRARY\_PATH` and `ldconfig`.

Compilation with Make only
---------------------

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
