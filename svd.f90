module svd_lib 
 use nan_lib
 implicit none
 interface svd
  module procedure d_svd, z_svd
 end interface

contains
 subroutine d_svd(a,u,v,s,tol,rmax,err,info)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,pointer :: u(:,:),v(:,:),s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='d_svd'
  double precision,allocatable :: b(:,:),work(:),ss(:),uu(:,:),vv(:,:)
  integer :: m,n,mn,lwork,ierr,r,i,j
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); lwork=256*max(m,n)
  allocate(uu(m,mn),vv(mn,n),ss(mn),b(m,n),work(lwork),stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dcopy(m*n,a,1,b,1)
  call dgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,ierr)
  if(present(info))info=ierr
  if(ierr.ne.0)then
   write(*,*)subnam,': dgesvd info: ',ierr
   if(ierr.lt.0)stop
   if(nan(a))then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   else
    write(*,*) subnam,': min/max element of input: ',minval(a),maxval(a)
   end if 
   u=>null(); v=>null(); s=>null()
  else 
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),v(n,r),s(r))
   call dcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   forall(i=1:n,j=1:r)v(i,j)=vv(j,i)
  end if 
  deallocate(uu,vv,ss,b,work)
 end subroutine 
 subroutine z_svd(a,u,v,s,tol,rmax,err,info)
  implicit none
  double complex,intent(in) :: a(:,:)
  double complex,pointer :: u(:,:),v(:,:)
  double precision,pointer :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='z_svd'
  double complex,allocatable :: b(:,:),work(:),uu(:,:),vv(:,:)
  double precision,allocatable :: ss(:),rwork(:)
  integer :: m,n,mn,lwork,ierr,r,i,j
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2); mn=min(m,n); lwork=256*max(m,n)
  allocate(uu(m,mn),vv(mn,n),ss(mn),b(m,n),work(lwork),rwork(lwork),stat=ierr)
  if(ierr.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call zcopy(m*n,a,1,b,1)
  call zgesvd('s','s',m,n,b,m,ss,uu,m,vv,mn,work,lwork,rwork,ierr)
  if(present(info))info=ierr
  if(ierr.ne.0)then
   write(*,*)subnam,': zgesvd info: ',ierr
   if(ierr.lt.0)stop
   if (nan(a)) then
    write(*,*) subnam,': NaNs detected in the input array'
    stop
   end if 
   u=>null(); v=>null(); s=>null()
  else 
   r=chop(ss,tol,rmax,err)
   if(present(err))err=err/dnrm2(mn,ss,1)
   allocate(u(m,r),v(n,r),s(r))
   call zcopy(m*r,uu,1,u,1)
   call dcopy(r,ss,1,s,1)
   forall(i=1:n,j=1:r)v(i,j)=vv(j,i)
  end if
  deallocate(uu,vv,ss,b,work)
 end subroutine 

 integer function chop(s,tol,rmax,err) result (r)
  implicit none
  double precision,intent(in) :: s(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: rmax
  double precision,intent(out),optional :: err
  double precision :: nrm,er,er2,bound
  double precision,external :: dnrm2
  r=size(s); er2=0.d0
  if(present(rmax))then
   if(rmax.lt.r)then
    er2=dot_product(s(rmax+1:r),s(rmax+1:r))
    r=rmax
   end if 
  end if 
  if(present(tol))then
   nrm=dnrm2(size(s),s,1)
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
