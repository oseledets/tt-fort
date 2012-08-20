module tts_lib
 use tt_lib
 implicit none
contains
 
 subroutine dtt_sparse(arg,tol)
  type(dtt),intent(inout) :: arg
  double precision,intent(in),optional :: tol
  character(len=*),parameter :: subnam='dtt_sparce'
  double precision :: eps,err,val,nrm
  integer :: l,m,k,i,j,n(tt_size),r(0:tt_size),nnz(tt_size)

  eps=1.d-2; if(present(tol))eps=tol
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; nnz=0
  
  do k=l,m
   val=maxval(abs(arg%u(k)%p))
   where(abs(arg%u(k)%p) < eps*val) arg%u(k)%p=0.d0
  end do
 end subroutine
 
 subroutine ztt_sparse(arg,tol)
  type(ztt),intent(inout) :: arg
  double precision,intent(in),optional :: tol
  character(len=*),parameter :: subnam='ztt_sparce'
  double precision :: eps,err,val,nrm
  integer :: l,m,k,i,j,n(tt_size),r(0:tt_size),nnz(tt_size)

  call svd(arg,tol=1.d-2)
  eps=1.d-2; if(present(tol))eps=tol
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; nnz=0
  
  do k=l,m
   val=maxval(abs(arg%u(k)%p))
   where(abs(arg%u(k)%p) < eps*val) arg%u(k)%p=0.d0
  end do
 end subroutine

end module
