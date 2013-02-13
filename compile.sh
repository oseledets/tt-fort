gfortran -c ttals.f90 tt_ksl.f90  -I../build/temp.macosx-10.5-x86_64-2.7/
gfortran test_ksl.f90 -o test -I../build/temp.macosx-10.5-x86_64-2.7 -L../build/temp.macosx-10.5-x86_64-2.7 -llapack -lblas 
