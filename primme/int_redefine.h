/**
 *  \file int_redefine.h.in
 *  \brief This is pretty weird header that redegine int with long if
 *  necessary. The reason is that Fortran code could be defined with different
 *  int length. This fact should be considered in C code as well.
 */

#pragma once

#if INTEGER_LENGTH == 4
#define int int
#elif INTEGER_LENGTH == 8
#define int long
#else
#error "Integerl length is not specified."
#endif
