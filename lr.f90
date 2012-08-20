module lr_lib
 use ort_lib
 implicit none

 interface lr
  module procedure lr_d2
 end interface

 contains
 subroutine lr_d2(a,u,b,tol,maxr,err)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,pointer :: u(:,:),b(:,:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: maxr
  double precision,intent(out),optional :: err
  character(len=*),parameter :: subnam='lr_d2'
  double precision,allocatable :: x(:,:),y(:,:),z(:,:),g(:,:),v(:)
  integer,allocatable :: q(:)
  integer :: m,n,mn,i,j,jj,ij(2),rmax,r,info
  double precision :: er,nrm,nr1,er1,zz,xx
  logical :: done
  double precision,external :: dnrm2
  integer,external :: idamax

  m=size(a,1); n=size(a,2); mn=min(m,n)
  if(present(maxr))then;if(maxr.ge.0)then;rmax=min(maxr,mn);else;rmax=mn;endif; else;rmax=mn;end if

  allocate(x(m,rmax),y(n,rmax),z(m,n),v(n),q(n), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory!';stop;endif
  
  call dcopy(m*n,a,1,z,1)
  nrm=dnrm2(m*n,a,1); er=nrm; nr1=-1.d0
  r=0; done=(r.ge.rmax)
  
  do while(.not.done)
!   ij=maxloc(abs(z)); i=ij(1); j=ij(2)
!...OMP parallel do shared(z,q,v,m,n)
   do jj=1,n
    q(jj)=idamax(m,z(1,jj),1); v(jj)=abs(z(q(jj),jj))
   end do
!...OMP END parallel do
   j=maxloc(v,1); i=q(j); zz=z(i,j); er1=dabs(zz)
   if(r.eq.0)nr1=er1
   r=r+1
   call dcopy(m,z(1,j),1,x(1,r),1); xx=dnrm2(m,x(1,r),1); call dscal(m,1.d0/xx,x(1,r),1)
   call dcopy(n,z(i,1),m,y(1,r),1); call dscal(n,xx/zz,y(1,r),1)
   call dger(m,n,-1.d0,x(1,r),1,y(1,r),1,z,m)
   er=dnrm2(m*n,z,1)
   done=r.eq.rmax
   if(present(tol))done=done.or.er.le.tol*nrm

   !write(*,'(i3,a,2i6,a,2e10.3)')r,' @ ',i,j,' er ',er/nrm,er1/nr1
  end do
  
  allocate(u(m,r),b(r,n),g(r,r), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory!';stop;endif
  
  call ort0(x(:,1:r),out=u,mat=g)
  call dgemm('n','t',r,n,r,1.d0,g,r,y,n,0.d0,b,r)
!-  
!  call dcopy(m*n,a,1,z,1)
!  call dgemm('n','n',m,n,r,-1.d0,u,m,b,r,1.d0,z,m)
!  er=dnrm2(m*n,z,1)
!  write(*,*) 'truerr : ',er/nrm
!- 
  deallocate(x,y,z,v,q,g)
  if(present(err))err=er/nrm
 end subroutine
end module
