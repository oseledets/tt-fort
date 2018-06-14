#   compiler-options.cmake
#   All compiler setting are located here.

# Validate input options.
if (WITH_BLAS AND WITH_LAPACK AND NOT WITH_MKL)
    message(STATUS "Use BLAS and LAPACK instead of Intel MKL.")
    set(GNU_LIBS ON)
elseif (NOT WITH_BLAS AND NOT WITH_LAPACK AND WITH_MKL)
    message(STATUS "Use Intel MKL instead of BLAS and LAPACK.")
    set(GNU_LIBS OFF)
else()
    message(FATAL_ERROR "BLAS and LAPACK could not be used together with MKL.")
endif()

# Set up linking libraries.
if (GNU_LIBS)
    set(LIBRARIES blas lapack gomp)
else()
    set(LIBRARIES mkl iomp5)
endif()

# Set up Fortran compiler flags
if (CMAKE_Fortran_COMPILER MATCHES "gfortran")
    set(INTEGER_FLAG "-fdefault-integer-${TT_INTEGER_SIZE}")
    set(CMAKE_Fortran_FLAGS "${INTEGER_FLAG} -ffree-line-length-none -O3")
elseif (CMAKE_Fortran_COMPILER MATCHES "ifort")
    set(INTEGER_FLAG "-i${TT_INTEGER_SIZE}")
    set(CMAKE_Fortran_FLAGS "${INTEGER_FLAG} -O3")
elseif()
    message(STATUS
            "Fortran Compiler ${CMAKE_Fortran_COMPILER} is not supported."
            CACHE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
endif()

# Set up C compiler Flags.
set(CMAKE_C_FLAGS "-O3 ${CMAKE_C_FLAGS}")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
