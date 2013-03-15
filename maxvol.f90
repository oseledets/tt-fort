module maxvol_lib
use time_lib
use ort_lib, only : ort0
use mat_lib, only : matinv
use sort_lib
use say_lib

double precision,private,parameter :: tol_def=1.d-2
logical,private,parameter :: loud=.false.
interface maxvol
 module procedure maxvol_d,maxvol_z
end interface
interface maxvola
 module procedure maxvola_d,maxvola_z
end interface
interface volume
 module procedure volume_d,volume_z
end interface
interface domval
 module procedure domval_d,domval_z
end interface

private :: maxvol_d,volume_d,maxvol0_d
private :: maxvol_z,volume_z,maxvol0_z

contains
 subroutine maxvol0_d(a,p,z,tol)
  implicit none
  double precision,intent(in) :: a(:,:)
  integer,intent(inout) :: p(size(a,2))
  double precision,intent(inout) :: z(size(a,1),size(a,2))
  double precision,intent(in),optional  :: tol
  character(len=*),parameter :: subnam='maxvol'
  integer,parameter          :: maxiter=2**10
  integer          :: m,n,i,j,ip,ij(2),info,iter,lwork,ii,jj
  integer,allocatable         :: q(:),per(:)
  double precision :: zinv,vol,amax,err,nrm,eps,t0
  double precision,allocatable :: b(:,:),work(:),u(:),v(:)
  double precision,allocatable :: bb(:,:),tt(:,:) !debug
  character(len=4) :: verb
  double precision,external :: dnrm2
  integer,external :: idamax
  
  m=size(a,1); n=size(a,2); t0=timef()
  
  if(present(tol))then; eps=tol;else;eps=tol_def;endif
  lwork=256*n
  allocate(work(lwork),b(n,n),q(n),per(m),u(m),v(n))
  
  if(isperm(p,m))then
   verb="from"
   vol=volume(a(p(1:n),:))
  else
   verb="init"
   call dcopy(m*n,a,1,z,1)   ! z=a
   call dgetrf(m,n,z,m,p,info)
   if(info.ne.0)then; write(*,*)subnam,': dgetrf(mxn) info: ',info; stop; end if
   forall(i=1:n)v(i)=dlog(dabs(z(i,i))); vol=sum(v)
   q=(/ (i,i=1,n) /)
   call dgetri(n,z,m,q,work,lwork,info)
   if(info.ne.0)then; write(*,*)subnam,': dgetri(mxn) info: ',info; stop; end if
   b=z(1:n,:)
   per(1:m)=(/ (i,i=1,m) /)
   do i=1,n
    ip=per(i); per(i)=per(p(i)); per(p(i))=ip ! per(i) <-> per(p(i))
   end do
   p=per(1:n)
   call dgemm('n','n',m,n,n,1.d0,a,m,b,n,0.d0,z,m)
   
   !write(*,'(2a,i8,1x,i3,3a,f9.3,a,f7.1,a)') subnam, '[',shape(a),'] ',verb,' e',vol,' in',timef()-t0,' sec'
  end if
  
  do iter=1,maxiter
!--
!   allocate(tt(m,n)); call maxvola_d(a,p,tt); nrm=dnrm2(m*n,tt,1); tt=tt-z; err=dnrm2(m*n,tt,1); deallocate(tt)
!   write(*,'(a,i4,3(a,e9.3))')'iter: ',iter,' err: ',err,' nrm: ',nrm,' relerr: ',err/nrm
!--

   do jj=1,n
    q(jj)=idamax(m,z(1,jj),1)
    v(jj)=abs(z(q(jj),jj))
   end do
   j=maxloc(v,1); i=q(j)
   
   amax=dabs(z(i,j))
   if(amax.lt.1.d0+eps .or. p(j).eq.i)goto 100
   ip=p(j); p(j)=i
   !write(*,'(a,i4,a,e11.5,a,i7,a,i7)') 'iter: ',iter,' amax: ',amax,' swap: ',ip,' <-> ',i
   u=z(:,j); v=z(ip,:)-z(i,:)
   zinv=1.d0/z(i,j)

!!  call dger(m,n,zinv,u,1,v,1,z,m)

   do jj=1,n
    call daxpy(m,zinv*v(jj),u,1,z(1,jj),1)
   end do
   z(i,:)=0.d0; z(i,j)=1.d0   ! FORCE
  end do 
100 continue
  deallocate(work,b,q,per,u,v)
  

  if(loud)write(*,'(2a,i8,1x,i3,3a,f9.3,a,f9.3,a,i5,a,f7.1,a)') subnam, '[',shape(a),&
       '] ',verb,' e',vol,' upto e',volume(a(p,:)),' in ',iter,' iters in',timef()-t0,' sec'
  
!  call say(z)
!  call sort(p)
!   write(*,'(a,64(i7,1x))') 'after: ',p
!  do i=1,n
!   z(p(i),:)=0.d0; z(p(i),i)=1.d0
!  end do
!  call say(z)
!--
!allocate(tt(m,n)); call maxvola_d(a,p,tt); nrm=dnrm2(m*n,tt,1); tt=tt-z; err=dnrm2(m*n,tt,1); deallocate(tt)
!  write(*,'(2a,i9,1x,i4,3a,f10.5,a,f10.5,a,i5,a,f9.1,a,e9.3)') subnam, '[',shape(a),&
!       ']:',verb,' e',vol,' upto e',volume(a(p,:)),' in ',iter,' iters in',timef()-t0,' sec; relerr: ',err/nrm
!--
  return
 end subroutine
 
 subroutine maxvol0_z(a,p,z,tol)
  implicit none
  complex(8),intent(in) :: a(:,:)
  integer,intent(inout) :: p(size(a,2))
  complex(8),intent(inout) :: z(size(a,1),size(a,2))
  double precision,intent(in),optional  :: tol
  character(len=*),parameter :: subnam='maxvol'
  integer,parameter          :: maxiter=2**10
  integer          :: m,n,i,j,ip,ij(2),info,iter,lwork,ii,jj
  integer,allocatable         :: q(:),per(:)
  double precision :: vol,amax,err,nrm,eps,t0
  complex(8)   :: zinv
  complex(8),allocatable :: b(:,:),work(:),u(:)
  double precision,allocatable :: v(:)
  complex(8),allocatable :: bb(:,:),tt(:,:) !debug
  character(len=4) :: verb
  double precision,external :: dznrm2
  integer,external :: izamax
  real(8) abs
  m=size(a,1); n=size(a,2); t0=timef()
  
  if(present(tol))then; eps=tol;else;eps=tol_def;endif
  lwork=256*n
  allocate(work(lwork),b(n,n),q(n),per(m),u(m),v(n))
  
  if(isperm(p,m))then
   verb="from"
   vol=volume(a(p(1:n),:))
  else
   verb="init"
   call zcopy(m*n,a,1,z,1)   ! z=a
   call zgetrf(m,n,z,m,p,info)
   if(info.ne.0)then; write(*,*)subnam,': zgetrf(mxn) info: ',info; stop; end if
   forall(i=1:n)v(i)=dlog(abs(z(i,i))); vol=sum(v)
   q=(/ (i,i=1,n) /)
   call zgetri(n,z,m,q,work,lwork,info)
   if(info.ne.0)then; write(*,*)subnam,': zgetri(mxn) info: ',info; stop; end if
   b=z(1:n,:)
   per(1:m)=(/ (i,i=1,m) /)
   do i=1,n
    ip=per(i); per(i)=per(p(i)); per(p(i))=ip ! per(i) <-> per(p(i))
   end do
   p=per(1:n)
   call zgemm('n','n',m,n,n,(1.d0,0.d0),a,m,b,n,(0.d0,0.d0),z,m)
   
   !write(*,'(2a,i8,1x,i3,3a,f9.3,a,f7.1,a)') subnam, '[',shape(a),'] ',verb,' e',vol,' in',timef()-t0,' sec'
  end if
  
  do iter=1,maxiter
   do jj=1,n
    q(jj)=izamax(m,z(1,jj),1)
    v(jj)=abs(z(q(jj),jj))
   end do
   j=maxloc(v,1); i=q(j)
   
   amax=abs(z(i,j))
   if(amax.lt.1.d0+eps .or. p(j).eq.i)goto 100
   ip=p(j); p(j)=i
   !write(*,'(a,i4,a,e11.5,a,i7,a,i7)') 'iter: ',iter,' amax: ',amax,' swap: ',ip,' <-> ',i
   u=z(:,j); v=z(ip,:)-z(i,:)
   zinv=(1.d0,0.d0)/z(i,j)

   do jj=1,n
    call zaxpy(m,zinv*v(jj),u,1,z(1,jj),1)
   end do
   z(i,:)=(0.d0,0.d0); z(i,j)=(1.d0,0.d0)   ! FORCE
  end do 
100 continue
  deallocate(work,b,q,per,u,v)
  

  if(loud)write(*,'(2a,i8,1x,i3,3a,f9.3,a,f9.3,a,i5,a,f7.1,a)') subnam, '[',shape(a),&
       '] ',verb,' e',vol,' upto e',volume(a(p,:)),' in ',iter,' iters in',timef()-t0,' sec'
 end subroutine


 subroutine maxvol_d(a,p,z,tol)
  implicit none
  double precision, intent(in)            :: a(:,:)
  integer,intent(inout)                   :: p(size(a,2))
  double precision,pointer,optional       :: z(:,:)
  double precision,intent(in),optional    :: tol
  character(len=*),parameter :: subnam='maxvol'
  integer :: m,n,k
  double precision,allocatable :: d(:,:),y(:,:),f(:,:),g(:,:),dd(:,:)
  double precision,pointer     :: znew(:,:)
  double precision :: err,nrm
  double precision,external :: dnrm2
  
  if(size(a,1).lt.size(a,2))then
   write(*,'(2a,a,2i8,a)') subnam,': wide matrices not permited: ','a[',shape(a),']'
   stop
  end if
  
  if(.not.present(z))then
   allocate(y(size(a,1),size(a,2)))
   if(isperm(p,size(a,1)))call maxvola_d(a,p,y)
   call maxvol0_d(a,p,y,tol)
   deallocate(y)
   return
  end if
  if(.not.associated(z))then
   allocate(z(size(a,1),size(a,2)))
   if(isperm(p,size(a,1)))call maxvola_d(a,p,z)
   call maxvol0_d(a,p,z,tol)
   return
  end if
  
  if(size(z,1).ne.size(a,1))then
   write(*,'(2a,2(a,2i8),a)') subnam,': sizes mismatch: ','a[',shape(a),'], z[',shape(z),']'
   stop
  end if
  if(size(z,2).gt.size(a,2))then
   write(*,'(2a,2(a,2i8),a)') subnam,': z is wider than a:','a[',shape(a),'], z[',shape(z),']'
   stop
  end if

  m=size(a,1); k=size(a,2)-size(z,2); n=size(z,2)
  if(k.eq.0)then
   if(isperm(p,size(a,1)))call maxvola_d(a,p,z)
   call maxvol0_d(a,p,z,tol)
   return
  end if
  
  IF(isperm(p(1:n),m))THEN
   allocate(znew(m,n+k),d(m,k),y(m,k),f(n,k),g(k,n))
   d=a(:,n+1:n+k); f=a(p(1:n),n+1:n+k)
   call dgemm('n','n',m,k,n,-1.d0,z,m,f,n,1.d0,d,m)
   
   d(p(1:n),:)=0.d0     ! FORCE
   call ort0(d)
   p(n+1:n+k)=0; call maxvol0_d(d,p(n+1:n+k),y,tol)
   y(p(1:n),:)=0.d0     ! FORCE
   
   g=z(p(n+1:n+k),:)
   call dcopy(m*n,z,1,znew,1)
   call dcopy(m*k,y,1,znew(:,n+1:n+k),1)
   call dgemm('n','n',m,n,k,-1.d0,y,m,g,k,1.d0,znew,m)
   znew(p(n+1:n+k),1:n)=0.d0     ! FORCE
   deallocate(d,y,f,g,z)
!--
!   allocate(d(n+k,n+k),dd(n+k,n+k),y(m,n+k))
!   d=a(p(1:n+k),:); write(*,'(2(a,f10.5))')'before: vol: ',volume(d),' val: ',domval(a,p)
!   call matinv(d,dd); call dgemm('n','n',m,n+k,n+k,1.d0,a,m,dd,n+k,0.d0,y,m)
!   y=y-ab; err=dnrm2(m*(n+k),y,1); nrm=dnrm2(m*(n+k),ab,1)
!   write(*,'(3(a,e9.3))')'before: err: ',err,' nrm: ',nrm,' relerr: ',err/nrm
!   deallocate(d,dd,y)
!--
   call maxvol0_d(a,p,znew,tol)
   z=>znew
!--
!   allocate(d(n+k,n+k),dd(n+k,n+k),y(m,n+k))
!   d=a(p(1:n+k),:); write(*,'(2(a,f10.5))')'after: vol: ',volume(d),' val: ',domval(a,p)
!   call matinv(d,dd); call dgemm('n','n',m,n+k,n+k,1.d0,a,m,dd,n+k,0.d0,y,m)
!   y=y-ab; err=dnrm2(m*(n+k),y,1); nrm=dnrm2(m*(n+k),ab,1)
!   write(*,'(3(a,e9.3))')'after: err: ',err,' nrm: ',nrm,' relerr: ',err/nrm
!   deallocate(d,dd,y)
!--
  ELSE
   write(*,*)subnam,': reject invalid p on input: ',p(1:n)
   deallocate(z)
   allocate(z(m,n+k))
   p=0; call maxvol0_d(a,p,z,tol)
  ENDIF
  
  return
 end subroutine
 
 subroutine maxvol_z(a,p,z,tol)
  implicit none
  complex(8), intent(in)            :: a(:,:)
  integer,intent(inout)                   :: p(size(a,2))
  complex(8),pointer,optional       :: z(:,:)
  double precision,intent(in),optional    :: tol
  character(len=*),parameter :: subnam='maxvol'
  integer :: m,n,k
  complex(8),allocatable :: d(:,:),y(:,:),f(:,:),g(:,:),dd(:,:)
  complex(8),pointer     :: znew(:,:)
  double precision :: err,nrm
  double precision,external :: dznrm2
  
  if(size(a,1).lt.size(a,2))then
   write(*,'(2a,a,2i8,a)') subnam,': wide matrices not permited: ','a[',shape(a),']'
   stop
  end if
  
  if(.not.present(z))then
   allocate(y(size(a,1),size(a,2)))
   if(isperm(p,size(a,1)))call maxvola_z(a,p,y)
   call maxvol0_z(a,p,y,tol)
   deallocate(y)
   return
  end if
  if(.not.associated(z))then
   allocate(z(size(a,1),size(a,2)))
   if(isperm(p,size(a,1)))call maxvola_z(a,p,z)
   call maxvol0_z(a,p,z,tol)
   return
  end if
  
  if(size(z,1).ne.size(a,1))then
   write(*,'(2a,2(a,2i8),a)') subnam,': sizes mismatch: ','a[',shape(a),'], z[',shape(z),']'
   stop
  end if
  if(size(z,2).gt.size(a,2))then
   write(*,'(2a,2(a,2i8),a)') subnam,': z is wider than a:','a[',shape(a),'], z[',shape(z),']'
   stop
  end if

  m=size(a,1); k=size(a,2)-size(z,2); n=size(z,2)
  if(k.eq.0)then
   if(isperm(p,size(a,1)))call maxvola_z(a,p,z)
   call maxvol0_z(a,p,z,tol)
   return
  end if
  
  IF(isperm(p(1:n),m))THEN
   allocate(znew(m,n+k),d(m,k),y(m,k),f(n,k),g(k,n))
   d=a(:,n+1:n+k); f=a(p(1:n),n+1:n+k)
   call zgemm('n','n',m,k,n,-(1.d0,0.d0),z,m,f,n,(1.d0,0.d0),d,m)
   
   d(p(1:n),:)=(0.d0,0.d0)     ! FORCE
   call ort0(d)
   p(n+1:n+k)=0; call maxvol0_z(d,p(n+1:n+k),y,tol)
   y(p(1:n),:)=(0.d0,0.d0)     ! FORCE
   
   g=z(p(n+1:n+k),:)
   call zcopy(m*n,z,1,znew,1)
   call zcopy(m*k,y,1,znew(:,n+1:n+k),1)
   call zgemm('n','n',m,n,k,-(1.d0,0.d0),y,m,g,k,(1.d0,0.d0),znew,m)
   znew(p(n+1:n+k),1:n)=(0.d0,0.d0)     ! FORCE
   deallocate(d,y,f,g,z)
   
   call maxvol0_z(a,p,znew,tol)
   z=>znew
  
  ELSE
   write(*,*)subnam,': reject invalid p on input: ',p(1:n)
   deallocate(z)
   allocate(z(m,n+k))
   p=0; call maxvol0_z(a,p,z,tol)
  ENDIF
 end subroutine
 

 double precision function volume_d(a) result(vol)
  implicit none
  double precision,intent(in) :: a(:,:)
  integer :: m,n,info,i
  double precision,allocatable :: b(:,:),d(:)
  integer, allocatable         :: p(:)

  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   vol=-999.d0
   return
  end if
  allocate(b(m,n),d(n),p(n))
  b=a
  call dgetrf(n,n,b,n,p,info)
  if(info.ne.0)then
   vol=-888.d0
   return
  end if
  forall(i=1:n) d(i)=dlog(dabs(b(i,i)))
  vol=sum(d)
  deallocate(b,d,p)
  return
 end function
 
 double precision function volume_z(a) result(vol)
  implicit none
  complex(8),intent(in) :: a(:,:)
  integer :: m,n,info,i
  complex(8),allocatable :: b(:,:),d(:)
  integer, allocatable         :: p(:)

  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   vol=-999.d0
   return
  end if
  allocate(b(m,n),d(n),p(n))
  b=a
  call zgetrf(n,n,b,n,p,info)
  if(info.ne.0)then
   vol=-888.d0
   return
  end if
  forall(i=1:n) d(i)=dlog(abs(b(i,i)))
  vol=sum(d)
  deallocate(b,d,p)
  return
 end function
 
 double precision function domval_d(a,p) result(val)
  implicit none
  double precision,intent(in) :: a(:,:)
  integer,intent(in)          :: p(:)
  double precision,allocatable:: b(:,:)
  integer                     :: m,n
  m=size(a,1); n=size(a,2)
  if(size(p).ne.n)then; val=-999.d0; return; end if
  allocate(b(m,n))
  call maxvola(a,p,b)
  val=maxval(abs(b))
  deallocate(b)
  return
 end function
 double precision function domval_z(a,p) result(val)
  implicit none
  complex(8),intent(in) :: a(:,:)
  integer,intent(in)          :: p(:)
  complex(8),allocatable:: b(:,:)
  integer                     :: m,n
  m=size(a,1); n=size(a,2)
  if(size(p).ne.n)then; val=-999.d0; return; end if
  allocate(b(m,n))
  call maxvola(a,p,b)
  val=maxval(abs(b))
  deallocate(b)
  return
 end function

 subroutine maxvola_d(a,p,b)
  implicit none
  double precision,intent(in) :: a(:,:)
  integer,intent(in)          :: p(size(a,2))
  double precision,intent(out):: b(size(a,1),size(a,2))
  double precision,allocatable:: abox(:,:),ainv(:,:)
  integer                     :: m,n
  m=size(a,1); n=size(a,2)
  allocate(abox(n,n),ainv(n,n))
  abox=a(p,:)
  call matinv(abox,ainv,alg='trf')
  call dgemm('n','n',m,n,n,1.d0,a,m,ainv,n,0.d0,b,m)
  deallocate(abox,ainv)
 end subroutine
 subroutine maxvola_z(a,p,b)
  implicit none
  complex(8),intent(in) :: a(:,:)
  integer,intent(in)          :: p(size(a,2))
  complex(8),intent(out):: b(size(a,1),size(a,2))
  complex(8),allocatable:: abox(:,:),ainv(:,:)
  integer                     :: m,n
  m=size(a,1); n=size(a,2)
  allocate(abox(n,n),ainv(n,n))
  abox=a(p,:)
  call matinv(abox,ainv,alg='trf')
  call zgemm('n','n',m,n,n,(1.d0,0.d0),a,m,ainv,n,(0.d0,0.d0),b,m)
  deallocate(abox,ainv)
 end subroutine

end module
