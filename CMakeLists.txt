cmake_minimum_required(VERSION 3.11)

project("TT for Fortran"
        VERSION 0.1.0
        LANGUAGES C CXX Fortran)

enable_testing()

set(TT_INTEGER_SIZE 8 CACHE STRING "Default integer size.")

option(WITH_BLAS "Use BLAS." OFF)
option(WITH_LAPACK "Use LAPACK." OFF)
option(WITH_MKL "Use MKL instead of BLAS and LAPACK." ON)

if (${TT_INTEGER_SIZE} STREQUAL 4 AND ${TT_INTEGER_SIZE} STREQUAL 8)
    message(FATAL_ERROR "Default integer size must be 4 bytes or 8 bytes.")
endif()

set(TT_CSOURCES gmres_3d.c)
set(TT_FSOURCES check.f90
                default.f90
                dispmodule.f90
                lr.f90
                mat.f90
                matrix_util.f90
                maxvol.f90
                nan.f90
                ort.f90
                ptype.f90
                putstrmodule.F90
                python_conv.f90
                rnd.f90
                say.f90
                sort.f90
                svd.f90
                timef.f90
                trans.f90
                tt.f90
                tt_adapt_als.f90
                ttals.f90
                ttamen.f90
                ttaux.f90
                tt_eigb.f90
                ttio.f90
                tt_linalg.f90
                ttlocsolve.f90
                ttnodeop.f90
                ttop.f90
                tts.f90)
set(TT_SOURCES ${TT_CSOURCES} ${TT_FSOURCES})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

add_subdirectory(primme)

add_library(tt-fort-obj OBJECT ${TT_SOURCES})
add_library(tt-fort-static STATIC $<TARGET_OBJECTS:tt-fort-obj>)
add_library(tt-fort-dynamic SHARED $<TARGET_OBJECTS:tt-fort-obj>)
add_library(mytt ALIAS tt-fort-static)  # backward compatibility

if (BUILD_TESTING)
    add_executable(test_eigb test_svd.f90)
    target_link_libraries(test_eigb tt-fort-static primme-static blas lapack gomp)

    add_executable(test_svd test_svd.f90)
    target_link_libraries(test_svd tt-fort-static blas lapack gomp)

    add_test(NAME test_eigb COMMAND test_eigb)
    add_test(NAME test_svd COMMAND test_eigb)
endif()
