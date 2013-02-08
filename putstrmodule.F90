#ifndef MATLAB


MODULE PUTSTRMODULE ! DUMMY VERSION
  ! An auxilliary module that accompanies DISPMODULE. This module contains dummy versions of the
  ! subroutines putstr and putnl that do nothing. It is needed to avoid an "undefined symbol" link
  ! error for these. In addition it defines the named constant (or parameter) DEFAULT_UNIT = -3,
  ! which makes the asterisk unit (usually the screen) the default to display on.
  !
  ! The purpose of having this module is to make displaying possible in situations where ordinary
  ! print- and write-statements do not work. Then this module should be replaced by one defining
  ! functional versions of putstr and putnl. An example is given by the commented out PUTSTRMODULE
  ! for Matlab mex files below.
  !
  integer, parameter :: DEFAULT_UNIT = -3
  !
CONTAINS
  subroutine putstr(s)
    character(*), intent(in) :: s
    integer ldummy, ldummy1  ! these variables exist to avoid unused variable warnings
    ldummy = len(s)
    ldummy1 = ldummy
    ldummy = ldummy1
  end subroutine putstr
  subroutine putnl()
  end subroutine putnl
END MODULE PUTSTRMODULE

#else

MODULE PUTSTRMODULE  ! for Matlab mex files.
  ! This module contains functional versions of subroutines putstr and putnl. It also sets
  ! DEFAULT_UNIT = -2, which makes putstr/putnl the default to display with. Using this module,
  ! instead of the dummy module above allows DISPMODULE to be used with Matlab mex files.
  ! used (commented in) instead of the one above (which should then be commented out), then
  ! DISPMODULE can be used with Matlab mex files. A shorter version (given in the user manual)
  ! may be used with g95, but the one below works for both g95 and gfortran.
  !
  use, intrinsic :: ISO_C_BINDING
  integer, parameter :: default_unit = -2
  interface
    subroutine mexprintf(s) bind(C, name = 'mexPrintf')
      import c_char
      character(c_char) s(*)
    end subroutine mexprintf
  end interface
  interface
    subroutine mexEvalString(s) bind(C, name = 'mexEvalString')
      import c_char
      character(c_char) s(*)
    end subroutine mexEvalString
  end interface
CONTAINS
  subroutine putstr(s)
    character(*), intent(in) :: s
    call mexprintf(s//char(0))
    call mexEvalString('drawnow;'); ! to dump string.
  end subroutine putstr
  subroutine putnl()
    call mexprintf(char(10)//char(0))
    call mexEvalString('drawnow;'); ! to dump string.
  end subroutine putnl
END MODULE PUTSTRMODULE

#endif
