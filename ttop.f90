module ttop_lib
 use tt_lib
 use ttaux_lib
 implicit none
 
 interface axpy
  module procedure dtt_axpy,ztt_axpy
 end interface
 

contains

! AXPY
 subroutine dtt_axpy(a,x,b,y)
  ! y=b*y+a*x
  implicit none
  type(dtt),intent(in),target :: x
  type(dtt),intent(inout),target :: y
  double precision,intent(in) :: a,b
  character(len=*),parameter :: subnam='dtt_axpy'
  type(dtt) :: z
  integer :: l,m,k,i,j,p
  integer,pointer :: r(:),q(:),n(:)
  logical :: zx,zy
  
  zx=.true.; if(x%m.ge.x%l)zx=.not.(all(x%n(x%l:x%m)>0).and.all(x%r(x%l-1:x%m)>0))
  zy=.true.; if(y%m.ge.y%l)zy=.not.(all(y%n(y%l:y%m)>0).and.all(x%r(x%l-1:x%m)>0))
  if(zx.and.zy)return
  if(zx)then;call dscal(y%r(l-1)*y%n(l)*y%r(l),b,y%u(l)%p,1);return;endif
  if(zy)then;y=a*x;return;endif
  
  if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=x%l; m=x%m
  if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>x%n;r=>x%r;q=>y%r
  if(r(l-1).ne.q(l-1) .or. r(m).ne.q(m))then;write(*,*)subnam,': border ranks mismatch';stop;endif

  z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1)=r(l-1); z%r(m)=r(m); z%r(l:m-1)=r(l:m-1)+q(l:m-1)
  if(l.eq.m)then;y%u(l)%p=a*x%u(l)%p + b*y%u(l)%p;return;endif
  call alloc(z)
  
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) z%u(l)%p(i,j,     k)=a*(x%u(l)%p(i,j,k))
  forall(i=1:r(l-1),j=1:n(l),k=1:q(l)) z%u(l)%p(i,j,r(l)+k)=b*(y%u(l)%p(i,j,k))
  do p=l+1,m-1
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(       i,j,     k)=x%u(p)%p(i,j,k)
   forall(i=1:r(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(       i,j,r(p)+k)=0.d0
   forall(i=1:q(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(r(p-1)+i,j,     k)=0.d0
   forall(i=1:q(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(r(p-1)+i,j,r(p)+k)=y%u(p)%p(i,j,k)
   !do i=1,r(p-1);do j=1,n(p);do k=1,r(p); z%u(p)%p(       i,j,     k)=x%u(p)%p(i,j,k); enddo;enddo;enddo
   !do i=1,r(p-1);do j=1,n(p);do k=1,q(p); z%u(p)%p(       i,j,r(p)+k)=0.d0;            enddo;enddo;enddo
   !do i=1,q(p-1);do j=1,n(p);do k=1,r(p); z%u(p)%p(r(p-1)+i,j,     k)=0.d0;            enddo;enddo;enddo
   !do i=1,q(p-1);do j=1,n(p);do k=1,q(p); z%u(p)%p(r(p-1)+i,j,r(p)+k)=y%u(p)%p(i,j,k); enddo;enddo;enddo
  end do
  forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(       i,j,k)=x%u(m)%p(i,j,k)
  forall(i=1:q(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(r(m-1)+i,j,k)=y%u(m)%p(i,j,k)
  
  y=z; call dealloc(z)
 end subroutine

 subroutine ztt_axpy(a,x,b,y)
  ! y=b*y+a*x
  implicit none
  type(ztt),intent(in),target :: x
  type(ztt),intent(inout),target :: y
  complex(8),intent(in) :: a,b
  character(len=*),parameter :: subnam='ztt_axpy'
  type(ztt) :: z
  integer :: l,m,k,i,j,p
  integer,pointer :: r(:),q(:),n(:)
  logical :: zx,zy
  
  zx=.true.; if(x%m.ge.x%l)zx=.not.(all(x%n(x%l:x%m)>0).and.all(x%r(x%l-1:x%m)>0))
  zy=.true.; if(y%m.ge.y%l)zy=.not.(all(y%n(y%l:y%m)>0).and.all(x%r(x%l-1:x%m)>0))
  if(zx.and.zy)return
  if(zx)then;call zscal(y%r(l-1)*y%n(l)*y%r(l),b,y%u(l)%p,1);return;endif
  if(zy)then;y=a*x;return;endif
  
  if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=x%l; m=x%m
  if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>x%n;r=>x%r;q=>y%r
  if(r(l-1).ne.q(l-1) .or. r(m).ne.q(m))then;write(*,*)subnam,': border ranks mismatch';stop;endif

  z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1)=r(l-1); z%r(m)=r(m); z%r(l:m-1)=r(l:m-1)+q(l:m-1)
  if(l.eq.m)then;y%u(l)%p=a*x%u(l)%p + b*y%u(l)%p;return;endif
  call alloc(z)
 
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l)) z%u(l)%p(i,j,     k)=a*(x%u(l)%p(i,j,k))
  forall(i=1:r(l-1),j=1:n(l),k=1:q(l)) z%u(l)%p(i,j,r(l)+k)=b*(y%u(l)%p(i,j,k))
  do p=l+1,m-1
   forall(i=1:r(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(       i,j,     k)=x%u(p)%p(i,j,k)
   forall(i=1:r(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(       i,j,r(p)+k)=(0.d0,0.d0)
   forall(i=1:q(p-1),j=1:n(p),k=1:r(p)) z%u(p)%p(r(p-1)+i,j,     k)=(0.d0,0.d0)
   forall(i=1:q(p-1),j=1:n(p),k=1:q(p)) z%u(p)%p(r(p-1)+i,j,r(p)+k)=y%u(p)%p(i,j,k)
  end do
  forall(i=1:r(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(       i,j,k)=x%u(m)%p(i,j,k)
  forall(i=1:q(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(r(m-1)+i,j,k)=y%u(m)%p(i,j,k)
  
  y=z; call dealloc(z)
 end subroutine



! HA
 subroutine dtt_ha1(x,y,z)
  ! z(i)=x(i).*y(i)
  implicit none
  type(dtt),intent(in),target :: x,y
  type(dtt),intent(inout),target :: z
  character(len=*),parameter :: subnam='dtt_ha1'
  integer :: l,m,k,i,j,ii,jj,kk,p
  integer,pointer :: r(:),q(:),n(:)
  logical :: zx,zy
  
  zx=.true.; if(x%m.ge.x%l)zx=.not.all(x%n(x%m:x%l)>0)
  zy=.true.; if(y%m.ge.y%l)zy=.not.all(y%n(y%m:y%l)>0)
  if(zx.or.zy)then;z%r=0;return;endif

  if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=x%l; m=x%m
  if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>x%n;r=>x%r;q=>y%r
  if(r(l-1).ne.q(l-1) .or. r(m).ne.q(m))then;write(*,*)subnam,': border ranks mismatch';stop;endif
  if(m.eq.l)then;write(*,*)subnam,': special code required for m=l!';stop;endif
  z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1)=r(l-1); z%r(m)=r(m); z%r(l:m-1)=r(l:m-1)*q(l:m-1)
  call alloc(z)
  
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l),kk=1:q(l)) z%u(l)%p(i,j,k+(kk-1)*r(l))=x%u(l)%p(i,j,k)*y%u(l)%p(i,j,kk)
  do p=l+1,m-1
   !forall(i=1:r(p-1),ii=1:q(p-1),j=1:n(p),k=1:r(p),kk=1:q(p)) z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   do i=1,r(p-1); do ii=1,q(p-1); do j=1,n(p); do k=1,r(p); do kk=1,q(p)
    z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   enddo; enddo; enddo; enddo; enddo
  end do
  forall(i=1:r(m-1),ii=1:q(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(i+(ii-1)*r(m-1),j,k)=x%u(m)%p(i,j,k)*y%u(m)%p(ii,j,k)
 end subroutine
 subroutine ztt_ha1(x,y,z)
  ! z(i)=x(i).*y(i)
  implicit none
  type(ztt),intent(in),target :: x,y
  type(ztt),intent(inout),target :: z
  character(len=*),parameter :: subnam='ztt_ha1'
  integer :: l,m,k,i,j,ii,jj,kk,p
  integer,pointer :: r(:),q(:),n(:)
  logical :: zx,zy
  
  zx=.true.; if(x%m.ge.x%l)zx=.not.all(x%n(x%m:x%l)>0)
  zy=.true.; if(y%m.ge.y%l)zy=.not.all(y%n(y%m:y%l)>0)
  if(zx.or.zy)then;z%r=0;return;endif

  if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=x%l; m=x%m
  if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>x%n;r=>x%r;q=>y%r
  if(r(l-1).ne.q(l-1) .or. r(m).ne.q(m))then;write(*,*)subnam,': border ranks mismatch';stop;endif
  if(m.eq.l)then;write(*,*)subnam,': special code required for m=l!';stop;endif
  z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1)=r(l-1); z%r(m)=r(m); z%r(l:m-1)=r(l:m-1)*q(l:m-1)
  call alloc(z)
  
  forall(i=1:r(l-1),j=1:n(l),k=1:r(l),kk=1:q(l)) z%u(l)%p(i,j,k+(kk-1)*r(l))=x%u(l)%p(i,j,k)*y%u(l)%p(i,j,kk)
  do p=l+1,m-1
   !forall(i=1:r(p-1),ii=1:q(p-1),j=1:n(p),k=1:r(p),kk=1:q(p)) z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   do i=1,r(p-1); do ii=1,q(p-1); do j=1,n(p); do k=1,r(p); do kk=1,q(p)
    z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   enddo; enddo; enddo; enddo; enddo
  end do
  forall(i=1:r(m-1),ii=1:q(m-1),j=1:n(m),k=1:r(m)) z%u(m)%p(i+(ii-1)*r(m-1),j,k)=x%u(m)%p(i,j,k)*y%u(m)%p(ii,j,k)
 end subroutine
 subroutine dtt_ha(x,y,z)
  ! z(ij)=x(i).*y(j)
  implicit none
  type(dtt),intent(in),target :: x,y
  type(dtt),intent(inout),target :: z
  character(len=*),parameter :: subnam='dtt_ha'
  integer :: l,m,k,i,j,ii,jj,kk,p
  integer,pointer :: r(:),q(:),n(:)
  logical :: zx,zy
  
  zx=.true.; if(x%m.ge.x%l)zx=.not.all(x%n(x%m:x%l)>0)
  zy=.true.; if(y%m.ge.y%l)zy=.not.all(y%n(y%m:y%l)>0)
  if(zx.or.zy)then;z%r=0;return;endif

  if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=x%l; m=x%m
  if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>x%n;r=>x%r;q=>y%r

  z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1:m)=r(l-1:m)*q(l-1:m)
  call alloc(z)
  do p=l,m
   !forall(i=1:r(p-1),ii=1:q(p-1),j=1:n(p),k=1:r(p),kk=1:q(p)) z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   do i=1,r(p-1); do ii=1,q(p-1); do j=1,n(p); do k=1,r(p); do kk=1,q(p)
    z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   enddo; enddo; enddo; enddo; enddo
  end do
 end subroutine
 
 subroutine ztt_ha(x,y,z)
  ! z(ij)=x(i).*y(j)
  implicit none
  type(ztt),intent(in),target :: x,y
  type(ztt),intent(inout),target :: z
  character(len=*),parameter :: subnam='ztt_ha'
  integer :: l,m,k,i,j,ii,jj,kk,p
  integer,pointer :: r(:),q(:),n(:)
  logical :: zx,zy
  
  zx=.true.; if(x%m.ge.x%l)zx=.not.all(x%n(x%m:x%l)>0)
  zy=.true.; if(y%m.ge.y%l)zy=.not.all(y%n(y%m:y%l)>0)
  if(zx.or.zy)then;z%r=0;return;endif

  if(x%l.ne.y%l .or. x%m.ne.y%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=x%l; m=x%m
  if(.not.all(x%n(l:m)==y%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>x%n;r=>x%r;q=>y%r

  z%l=l; z%m=m; z%n=n; z%r=0; z%r(l-1:m)=r(l-1:m)*q(l-1:m)
  call alloc(z)
  do p=l,m
   !forall(i=1:r(p-1),ii=1:q(p-1),j=1:n(p),k=1:r(p),kk=1:q(p)) z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   do i=1,r(p-1); do ii=1,q(p-1); do j=1,n(p); do k=1,r(p); do kk=1,q(p)
    z%u(p)%p(i+(ii-1)*r(p-1),j,k+(kk-1)*r(p))=x%u(p)%p(i,j,k)*y%u(p)%p(ii,j,kk)
   enddo; enddo; enddo; enddo; enddo
  end do
 end subroutine

! INV
 subroutine dtt_inv(arg,inv,tol,maxiter)
  implicit none
  type(dtt),intent(inout) :: arg
  type(dtt),intent(inout),optional :: inv
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: maxiter
  character(len=*),parameter :: subnam='dtt_inv'
  type(dtt) :: id,x,y,z,tmp
  integer :: l,m,n(tt_size),i,maxit,it
  double precision :: err,truerr,olderr,eps,nrm1,nrm2,nrm3, t0,t1,t2

  if(present(maxiter))then;maxit=maxiter;else;maxit=100;endif
  if(present(tol))then; eps=tol; else; eps=1.d-12; endif
  l=arg%l; m=arg%m; n=arg%n
  id%l=l; id%m=m; id%n=n; call ones(id)
  if(present(inv))then;x=inv;else;x=id;endif

  err=2.d0
  do while(err.gt.1.d0)
   call dtt_ha1(arg,x,y); call svd(y,tol=1.d-12) 
   z=y; call dtt_axpy(-1.d0,id,1.d0,z); call svd(z,tol=1.d-2)
   err=norm(z)/norm(id); truerr=err
   !write(*,'(2a,e10.3)')subnam,':start err: ',err
   if(err.gt.1.d0)call dscal(size(x%u(l)%p),1.d-1,x%u(l)%p,1)
  end do
  
  it=0; t0=timef(); olderr=1.d0
  do while(it.le.maxit .and. err.ge.eps .and.truerr.lt.olderr)
   it=it+1; olderr=truerr

!-[newton]-
   t1=timef()
   z=y; call dtt_axpy(2.d0,id,-1.d0,z)
   call dtt_ha1(x,z, y); call svd(y,tol=1.d-14)
   x=y
   t2=timef()
   call dtt_ha1(arg,x, y); call svd(y,tol=1.d-12)
!-

!-[modified newton]
!   t1=timef()
!   z=y;
!   call dtt_axpy(2.d0,id,-1.d0,z); call svd(z,tol=1.d-12); 
!   call dtt_ha1(x,z, tmp);         call svd(tmp,tol=1.d-12); x=tmp
!   call dtt_ha1(y,z, tmp);         call svd(tmp,tol=1.d-12); y=tmp
!   t2=timef()
!-

!-[another newton]
!  t1=timef()
!   z=y; call d3t_axpy(2.d0,id,-1.d0,z,rmax=(/2,2,2/))!tol=1.d-12)
!   call d3op('ha',x,z, tmp, tol=1.d-12); x=tmp
!   call d3op('ha',a,x, y, tol=1.d-12)
!   !tmp=y; call d3t_axpy(2.d0,id,-1.d0,tmp,tol=1.d-2); err=norm(tmp)/norm(id); truerr=err
!  t2=timef()
!-
   
   tmp=y; call dtt_axpy(1.d0,id,-1.d0,tmp); call svd(tmp,tol=1.d-2); err=norm(tmp)/norm(id)
   call dtt_ha1(arg,x, tmp); call svd(tmp,tol=1.d-12)
   call dtt_axpy(1.d0,id,-1.d0,tmp); call svd(tmp,tol=1.d-2); truerr=norm(tmp)/norm(id)
  
   !write(*,'(i3,3(a,f9.2),a,2e11.3,a,2f9.1)') it,' h[',rank(z),'] y[',rank(y),'] x[',rank(x),'] err ',err,truerr,' time ',t2-t1,t2-t0
  end do 
  if(present(inv))then;inv=x;else;arg=x;endif
  call dealloc(x); call dealloc(y); call dealloc(z); call dealloc(id); call dealloc(tmp)
 end subroutine 
 subroutine ztt_inv(arg,inv,tol,maxiter)
  implicit none
  type(ztt),intent(inout) :: arg
  type(ztt),intent(inout),optional :: inv
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: maxiter
  character(len=*),parameter :: subnam='ztt_inv'
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0),two=(2.d0,0.d0)
  type(ztt) :: id,x,y,z,tmp
  integer :: l,m,n(tt_size),i,maxit,it
  double precision :: err,truerr,olderr,eps,nrm1,nrm2,nrm3, t0,t1,t2

  if(present(maxiter))then;maxit=maxiter;else;maxit=100;endif
  if(present(tol))then; eps=tol; else; eps=1.d-12; endif
  l=arg%l; m=arg%m; n=arg%n
  id%l=l; id%m=m; id%n=n; call ones(id)
  if(present(inv))then;x=inv;else;x=id;endif

  err=2.d0
  do while(err.gt.1.d0)
   call ztt_ha1(arg,x,y); call svd(y,tol=1.d-12) 
   z=y; call ztt_axpy(-one,id,one,z); call svd(z,tol=1.d-2)
   err=norm(z)/norm(id); truerr=err
   !write(*,'(2a,e10.3)')subnam,':start err: ',err
   if(err.gt.1.d0)call zdscal(size(x%u(l)%p),1.d-1,x%u(l)%p,1)
  end do
  
  it=0; t0=timef(); olderr=1.d0
  do while(it.le.maxit .and. err.ge.eps .and.truerr.lt.olderr)
   it=it+1; olderr=truerr

!-[newton]-
   t1=timef()
   z=y; call ztt_axpy(two,id,-one,z)
   call ztt_ha1(x,z, y); call svd(y,tol=1.d-14)
   x=y
   t2=timef()
   call ztt_ha1(arg,x, y); call svd(y,tol=1.d-12)
!--   
   tmp=y; call ztt_axpy(one,id,-one,tmp); call svd(tmp,tol=1.d-2); err=norm(tmp)/norm(id)
   call ztt_ha1(arg,x, tmp); call svd(tmp,tol=1.d-12)
   call ztt_axpy(one,id,-one,tmp); call svd(tmp,tol=1.d-2); truerr=norm(tmp)/norm(id)
  
   !write(*,'(i3,3(a,f9.2),a,2e11.3,a,2f9.1)') it,' h[',rank(z),'] y[',rank(y),'] x[',rank(x),'] err ',err,truerr,' time ',t2-t1,t2-t0
  end do 
  if(present(inv))then;inv=x;else;arg=x;endif
  call dealloc(x); call dealloc(y); call dealloc(z); call dealloc(id); call dealloc(tmp)
 end subroutine 
 subroutine z1_inv(arg,inv,tol,maxiter)
  implicit none
  complex(8),intent(inout) :: arg(:)
  complex(8),intent(inout),optional :: inv(:)
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: maxiter
  character(len=*),parameter :: subnam='z1_inv'
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0),two=(2.d0,0.d0)
  complex(8),allocatable :: x(:),y(:),z(:),tmp(:)
  integer :: l,m,i,maxit,it,n,info
  double precision :: err,truerr,olderr,eps,nrm1,nrm2,nrm3, t0,t1,t2
  double precision,external :: dznrm2

  n=size(arg)
  allocate(x(n),y(n),z(n),tmp(n),stat=info)

  if(present(maxiter))then;maxit=maxiter;else;maxit=100;endif
  if(present(tol))then; eps=tol; else; eps=1.d-12; endif
  if(present(inv))then;x=inv;else;x=one;endif
  
  if(info.ne.0)then;write(*,*)subnam,': allocate fail: ',info;stop;endif

  err=2.d0
  do while(err.gt.1.d0)
   y=arg*x
   z=y-one
   err=dznrm2(n,z,1)/dsqrt(dble(n)); truerr=err
   write(*,'(2a,e10.3)')subnam,':start err: ',err
   if(err.gt.1.d0)x=x/10
  end do
  
  it=0; t0=timef(); olderr=1.d0
  do while(it.le.maxit .and. err.ge.eps .and.truerr.lt.olderr)
   it=it+1; olderr=truerr
!-[newton]-
   t1=timef()
   x=-(y-two)*x
   y=arg*x
!--   
   tmp=one-y; err=dznrm2(n,tmp,1)/dsqrt(dble(n))
   truerr=err
   write(*,'(i3,a,2e11.3)') it,' err ',err,truerr
  end do 
  if(present(inv))then;inv=x;else;arg=x;endif
  deallocate(x,y,z,tmp)
 end subroutine 


!PRM
 subroutine dtt_prm(arg,tol,low)
  implicit none
  type(dtt),intent(inout),target :: arg
  double precision,intent(in),optional :: tol
  integer,intent(in),optional :: low
  type(dtt) :: tmp
  integer :: k,m,l,ll,i,ii,j,n(tt_size)
  integer,pointer :: r(:)
  double precision,allocatable :: a(:,:),b(:,:)
  double precision :: eps
  l=arg%l;m=arg%m;n=arg%n;r=>arg%r
  if(present(tol))then;eps=tol;else;eps=1.d-12;endif
  if(present(low))then;ll=low;else;ll=l+1;endif

  do k=ll,m
   arg%l=l;arg%m=k
   tmp%l=l;tmp%m=k;tmp%n(l:k)=n(l:k);tmp%r(l-1:k)=2*r(l-1:k); call alloc(tmp)
   do j=l,k
    forall(i=1:r(j-1),ii=1:r(j))
     tmp%u(j)%p(       i,1,     ii)=arg%u(j)%p(i,1,ii); tmp%u(j)%p(       i,2,     ii)=0.d0
     tmp%u(j)%p(r(j-1)+i,2,     ii)=arg%u(j)%p(i,1,ii); tmp%u(j)%p(r(j-1)+i,1,     ii)=0.d0
     tmp%u(j)%p(       i,1,r(j)+ii)=arg%u(j)%p(i,2,ii); tmp%u(j)%p(       i,2,r(j)+ii)=0.d0
     tmp%u(j)%p(r(j-1)+i,2,r(j)+ii)=arg%u(j)%p(i,2,ii); tmp%u(j)%p(r(j-1)+i,1,r(j)+ii)=0.d0
    end forall
   end do
   allocate(a(r(l-1),2*n(l)*tmp%r(l)),b(tmp%r(k-1)*n(k)*r(k),2))
   call dcopy(size(a),tmp%u(l)%p,1,a,1); call dcopy(size(b),tmp%u(k)%p,1,b,1)
   tmp%r(l-1)=r(l-1); tmp%r(k)=r(k)
   deallocate(tmp%u(l)%p,tmp%u(k)%p)
   allocate(tmp%u(l)%p(r(l-1),n(l),tmp%r(l)),tmp%u(k)%p(tmp%r(k-1),n(k),r(k)))
   arg=tmp
   call d2submat(r(l-1),n(l)*r(l),a(1,1),r(l-1)*2,arg%u(l)%p)
   call d2submat(r(l-1),n(l)*r(l),a(1,2),r(l-1)*2,tmp%u(l)%p)
   call dcopy(r(k-1)*n(k)*r(k),b(1,1),1,arg%u(k)%p,1)
   call dcopy(r(k-1)*n(k)*r(k),b(1,2),1,tmp%u(k)%p,1)
   call dtt_axpy(1.d0,tmp,1.d0,arg)
   deallocate(a,b)
   
   call svd(arg,tol=eps)
   !call say(arg)
  end do 
  call dealloc(tmp)
 end subroutine
 

end module
