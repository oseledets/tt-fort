#gfortran -c ttals.f90 tt_ksl.f90  -I../build/temp.macosx-10.5-x86_64-2.7/
gfortran -g ttals.f90 tt_ksl.f90 test_ksl.f90 -o test -I../build/temp.macosx-10.5-x86_64-2.7 -L../build/temp.macosx-10.5-x86_64-2.7 -lprint_lib -lexpokit -lmytt -llapack -lblas 
