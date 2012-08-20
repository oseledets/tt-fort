module blopex
! See http://arxiv.org/pdf/0705.2626.pdf
  use dispmodule
  integer n0
  abstract interface 
     subroutine python_matvec_template(n,m,x,y)
       integer, intent(in) :: n,m
       double precision, intent(in) :: x(n,m)
       double precision, intent(out) ::  y(n,m)
     end subroutine python_matvec_template
  end interface
  procedure(python_matvec_template), pointer :: mv

contains
   subroutine matvec_python(x,y,k,primme)
     implicit none
     integer, intent(in) :: k,primme
     double precision, intent(in) :: x(*)
     double precision :: y(*)
     call mv(n0,k,x,y)
   end subroutine 
 


 subroutine lobpcg2(n,m, avec,bvec,tvec,X,maxiter,tol,lambda)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: maxiter
    double precision, intent(in) :: tol
    double precision, intent(inout) :: X(n,m)
    double precision, intent(out) :: lambda(m)
    double precision W(n,m)
    double precision  rnorms(m)
    !Memory
    integer primme, ierr
    external avec,bvec, tvec
    include 'primme_f77.h'
      
    mv => avec
     n0 = n
     call primme_initialize_f77(primme)
     call primme_set_member_f77(primme, PRIMMEF77_n, n)
     call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec, matvec_python)
     call primme_set_method_f77(primme, PRIMMEF77_JDQMR_ETol, ierr)
     call primme_set_member_f77(primme, PRIMMEF77_initSize, 2) 
     call primme_set_member_f77(primme, PRIMMEF77_numEvals, 2) 

     call primme_display_params_f77(primme)
     call dprimme_f77(lambda, X, W, primme, ierr)
     return 

  end subroutine lobpcg2
end module blopex
