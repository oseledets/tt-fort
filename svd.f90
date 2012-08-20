module svd_lib 
 implicit none

contains
 subroutine svd2d(a,u,v,s,tol,rmax,err)
 implicit none
  double precision,intent(in) :: a(:,:)
  double precision,pointer :: u(:,:),v(:,:),s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  character(len=*),parameter :: subnam='svd2d_d2'
  double precision,allocatable :: b(:,:),work(:),ss(:),uu(:,:),vv(:,:)
  integer :: m,n,mn,lwork,info,r
  m=size(a,1); n=size(a,2); mn=min(m,n); lwork=256*max(m,n)
  allocate(uu(m,mn),vv(mn,n),ss(mn),b(m,n),work(lwork))
  b=a
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgesvd info: ',info; stop; end if
  r=chop(ss,tol,rmax,err)
  if(present(err))err=err/dsqrt(dot_product(ss,ss))
  deallocate(b,work)
  allocate(u(m,r),v(n,r),s(r))
  u=uu(:,1:r); v=transpose(vv(1:r,:)); s=ss(1:r)
  deallocate(uu,vv,ss)
 end subroutine 

 integer function chop(s,tol,rmax,err) result (r)
  implicit none
  double precision,intent(in) :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  double precision :: nrm,er,er2,bound
  double precision,external :: dnrm2
  integer :: n
  n=size(s); er2=0.d0
  if(present(rmax))then
   r=min(rmax,n)
   if(r.lt.n) er2=dot_product(s(r+1:n),s(r+1:n))
  else
   r=n
  end if 
  if(present(tol))then
   nrm=dnrm2(n,s,1)
   bound=tol*tol*nrm*nrm 
   er=er2+s(r)*s(r)
   do while(er.lt.bound)
    er2=er; r=r-1; er=er+s(r)*s(r)
   end do
  end if
  if(present(err))err=dsqrt(er2)
  return
 end function

end module
