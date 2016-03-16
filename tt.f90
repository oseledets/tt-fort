module tt_lib
 use ptype_lib
 use svd_lib
 use mat_lib
 use trans_lib
 use time_lib
 use say_lib
 implicit none
 integer,parameter :: tt_size=1023
 double precision,parameter :: errval = -999.d0

 type,public:: dtt
  integer :: l=1
  integer :: m=0
  integer :: n(tt_size)=0
  integer :: r(0:tt_size)=0
  type(pointd3) :: u(tt_size)
 !contains
 ! final :: dtt_dealloc
 end type
 type,public:: ztt
  integer :: l=1
  integer :: m=0
  integer :: n(tt_size)
  integer :: r(0:tt_size)
  type(pointz3) :: u(tt_size)
 !contains
 ! final :: ztt_dealloc
 end type
 type,public:: ttind
  integer :: p(tt_size)=0
  integer :: n(tt_size)=0
  integer :: m=0
 end type

 

 interface nip
  module procedure dtt_nip, ztt_nip
 end interface
 interface full
  module procedure dtt_full, ztt_full
 end interface
 interface tijk
  module procedure dtt_ijk,ztt_ijk
 end interface
 interface seti
  module procedure dtt_seti,ztt_seti
 end interface
 interface sumall
  module procedure dtt_sumall,ztt_sumall
 end interface
 interface numel
  module procedure dtt_numel,ztt_numel
 end interface
 interface value
  module procedure dtt_value,dtt_value0,ztt_value,ztt_value0
 end interface
 
 interface mem          ! memory to keep all cores
  module procedure dtt_mem,ztt_mem
 end interface

 !interface elem
 ! module procedure dtt_elem
 !end interface

 interface alloc
  module procedure dtt_alloc,ztt_alloc
 end interface
 interface dealloc
  module procedure dtt_dealloc,ztt_dealloc
 end interface

 interface ort
  module procedure dtt_ort, ztt_ort
 end interface
 interface svd
  module procedure dtt_svd,dtt_svd0, ztt_svd,ztt_svd0
 end interface

 interface say
  module procedure dtt_say, ztt_say, ttind_say
 end interface
 interface sayfull
  module procedure dtt_sayfull, ztt_sayfull
 end interface

 interface norm
  module procedure dtt_norm,ztt_norm
 end interface
 interface rank
  module procedure dtt_rank,ztt_rank
 end interface

 interface dot
    module procedure dtt_dot, ztt_dot
 end interface


 interface copy
  module procedure dtt_copy,ztt_copy
 end interface

 interface assignment (=)
  module procedure dtt_assign, ztt_assign, ttind_assign
 end interface
 interface operator (*)
  module procedure dttmul_dt, zttmul_zt
 end interface
 interface operator (.eq.)
  module procedure ttind_eq
 end interface
 interface operator (.lt.)
  module procedure ttind_lt
 end interface
 interface operator (.le.)
  module procedure ttind_le
 end interface
 interface operator (.gt.)
  module procedure ttind_gt
 end interface
 interface operator (.ge.)
  module procedure ttind_ge
 end interface

 interface find
  module procedure find_ttind
 end interface
 interface push
  module procedure push_ttind
 end interface
 interface dble
  module procedure dble_ttind
 end interface
 interface int
  module procedure int_ttind
 end interface
contains


! ORT
 subroutine dtt_ort(arg)
  ![tt] ortogonalize from left
  implicit none
  type(dtt),intent(inout),target :: arg
  character(len=*),parameter :: subnam='dtt_ort'
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  double precision,allocatable :: work(:),tau(:),mat(:),u(:)
  double precision :: err,nrm
  double precision,external :: dnrm2
  !double precision :: t1,t2
  l=arg%l; m=arg%m
  if(m.lt.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  lwork=128*nn*rr
  allocate(work(lwork),tau(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  do k=l,m-1
   mm=r(k-1)*n(k); nn=r(k); mn=min(mm,nn); kk=n(k+1)*r(k+1)
   call dcopy(mm*nn, arg%u(k)%p,1,u,1)
   !t1=timef()
   call dgeqrf(mm,nn, u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': dgeqrf info: ',info; stop; end if
   do j=1,nn
    forall(i=1:min(j,mm))    mat(i+(j-1)*mn)=u(i+(j-1)*mm)
    forall(i=min(j,mm)+1:mn) mat(i+(j-1)*mn)=0.d0
   end do
   call dorgqr(mm,mn,mn,u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': dorgqr info: ',info; stop; end if
   call dcopy(mm*mn, u,1,arg%u(k)%p,1)
   call dgemm('n','n',mn,kk,nn,1.d0,mat,mn,arg%u(k+1)%p,nn,0.d0,u,mn)
   if(r(k).ne.mn)then
    call dcopy(mm*mn, arg%u(k)%p,1,mat,1)
    deallocate(arg%u(k)%p,arg%u(k+1)%p)
    r(k)=mn
    allocate(arg%u(k)%p(r(k-1),n(k),r(k)),arg%u(k+1)%p(r(k),n(k+1),r(k+1)))
    call dcopy(mm*mn, mat,1,arg%u(k)%p,1)
   end if
   call dcopy(mn*kk, u,1,arg%u(k+1)%p,1)
  end do
  deallocate(work,tau,mat,u)

 end subroutine
 subroutine ztt_ort(arg)
  ![tt] ortogonalize from left
  implicit none
  type(ztt),intent(inout),target :: arg
  character(len=*),parameter :: subnam='ztt_ort'
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  complex(8),allocatable :: work(:),tau(:),mat(:),u(:)
  double precision :: err,nrm
  double precision,external :: dznrm2

  l=arg%l; m=arg%m
  if(m.lt.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  lwork=128*nn*rr
  allocate(work(lwork),tau(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  do k=l,m-1
   mm=r(k-1)*n(k); nn=r(k); mn=min(mm,nn); kk=n(k+1)*r(k+1)
   call zcopy(mm*nn, arg%u(k)%p,1,u,1)
   call zgeqrf(mm,nn, u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': zgeqrf info: ',info; stop; end if
   do j=1,nn
    forall(i=1:min(j,mm))    mat(i+(j-1)*mn)=u(i+(j-1)*mm)
    forall(i=min(j,mm)+1:mn) mat(i+(j-1)*mn)=(0.d0,0.d0)
   end do

   call zungqr(mm,mn,mn,u,mm,tau,work,lwork,info)
   if(info.ne.0)then; write(*,*) subnam,': zungqr info: ',info; stop; end if
!
   nrm=dznrm2(mm*nn,arg%u(k)%p,1)
   call zgemm('n','n',mm,nn,mn,(1.d0,0.d0),u,mm,mat,mn,(-1.d0,0.d0),arg%u(k)%p,mm)
   err=dznrm2(mm*nn,arg%u(k)%p,1)
   if(err.gt.1.d-10*nrm)then
    write(*,*)subnam,': qr error: m,n: ',mm,nn
    write(*,*)subnam,': err: ',err,' nrm ',nrm
    stop
   endif
!
   call zcopy(mm*mn, u,1,arg%u(k)%p,1)
   call zgemm('n','n',mn,kk,nn,(1.d0,0.d0),mat,mn,arg%u(k+1)%p,nn,(0.d0,0.d0),u,mn)
   if(r(k).ne.mn)then
    call zcopy(mm*mn, arg%u(k)%p,1,mat,1)
    deallocate(arg%u(k)%p,arg%u(k+1)%p,stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot deallocate data';stop;endif
    arg%r(k)=mn
    allocate(arg%u(k)%p(r(k-1),n(k),r(k)),arg%u(k+1)%p(r(k),n(k+1),r(k+1)))
    call zcopy(mm*mn, mat,1,arg%u(k)%p,1)
   end if
   call zcopy(mn*kk, u,1,arg%u(k+1)%p,1)
  end do
  deallocate(work,tau,mat,u)
 end subroutine



! SVD
 subroutine dtt_svd(arg,tol,rmax)
  implicit none
  type(dtt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  character(len=*),parameter :: subnam='dtt_svd'
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  integer, optional :: rmax
  double precision,allocatable :: work(:),s(:),mat(:),u(:)
  double precision :: err,nrm
  double precision,external :: dnrm2

  l=arg%l; m=arg%m
  if(m.le.l)return
  r=>arg%r; n=>arg%n
  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  lwork=128*nn*rr
  allocate(work(lwork),s(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  call dtt_ort(arg)
  do k=m,l+1,-1
   mm=r(k-1); nn=n(k)*r(k); mn=min(mm,nn); kk=r(k-2)*n(k-1)
   call dgesvd('s','s',mm,nn,arg%u(k)%p,mm,s,mat,mm,u,mn,work,lwork,info)
   if(info.ne.0)then; write(*,*)subnam,': dgesvd info: ',info; stop; end if
   rr=chop(s(1:mn), tol/dsqrt(dble(m-l)))
   if (present(rmax)) then 
      rr = min(rr, rmax)
   end if
   forall(i=1:mm,j=1:rr)mat((j-1)*mm+i)=s(j)*mat((j-1)*mm+i)
   call d2submat(rr,nn,u,mn,arg%u(k)%p)

   call dgemm('n','n',kk,rr,mm,1.d0,arg%u(k-1)%p,kk,mat,mm,0.d0,u,kk)
   if(r(k-1).ne.rr)then
    call dcopy(rr*nn, arg%u(k)%p,1,mat,1)
    deallocate(arg%u(k-1)%p,arg%u(k)%p)
    r(k-1)=rr
    allocate(arg%u(k-1)%p(r(k-2),n(k-1),r(k-1)),arg%u(k)%p(r(k-1),n(k),r(k)))
    call dcopy(rr*nn, mat,1,arg%u(k)%p,1)
   end if
   call dcopy(kk*rr,u,1,arg%u(k-1)%p,1)
  end do

  deallocate(work,s,mat,u)
 end subroutine
 
 subroutine ztt_svd(arg, tol, rmax)
  implicit none
  type(ztt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  character(len=*),parameter :: subnam='ztt_svd'
  integer :: l,m,k,i,j,lwork,info,nn,rr,mn,mm,kk
  integer,pointer :: r(:),n(:)
  integer, intent(in), optional :: rmax
  double precision,allocatable :: s(:),rwork(:)
  complex(8),allocatable :: work(:),mat(:),u(:)
  double precision :: err,nrm
  double precision,external :: dznrm2

  l=arg%l; m=arg%m
  if(m.le.l)return
  r=>arg%r; n=>arg%n

  nn=maxval(n(l:m)); rr=maxval(r(l-1:m))
  lwork=128*nn*rr
  allocate(work(lwork),rwork(8*nn*rr),s(nn*rr), mat(rr*nn*rr),u(rr*nn*rr), stat=info)
  if(info.ne.0)then;write(*,*)subnam,': no memory';stop;endif

  call ztt_ort(arg)

  do k=m,l+1,-1
   mm=r(k-1); nn=n(k)*r(k); mn=min(mm,nn); kk=r(k-2)*n(k-1)
   call zgesvd('s','s',mm,nn,arg%u(k)%p,mm,s,mat,mm,u,mn,work,lwork,rwork,info)
   if(info.ne.0)then; write(*,*)subnam,': dgesvd info: ',info; stop; end if
   rr=chop(s(1:mn),tol/dsqrt(dble(m-l)))
   if (present(rmax)) then
      rr = min(rr, rmax)
   end if
   forall(i=1:mm,j=1:rr)mat((j-1)*mm+i)=s(j)*mat((j-1)*mm+i)
   call z2submat(rr,nn,u,mn,arg%u(k)%p)

   call zgemm('n','n',kk,rr,mm,(1.d0,0.d0),arg%u(k-1)%p,kk,mat,mm,(0.d0,0.d0),u,kk)
   if(r(k-1).ne.rr)then
    call zcopy(rr*nn, arg%u(k)%p,1,mat,1)
    deallocate(arg%u(k-1)%p,arg%u(k)%p)
    r(k-1)=rr
    allocate(arg%u(k-1)%p(r(k-2),n(k-1),r(k-1)),arg%u(k)%p(r(k-1),n(k),r(k)))
    call zcopy(rr*nn, mat,1,arg%u(k)%p,1)
   end if
   call zcopy(kk*rr,u,1,arg%u(k-1)%p,1)
  end do

  deallocate(work,rwork,s,mat,u)
 end subroutine

 subroutine dtt_svd0(n,a,tt,tol,rmax)
  implicit none
  integer,intent(in) :: n(:)
  double precision,intent(in) :: a(*)
  type(dtt),intent(inout),target :: tt
  double precision,intent(in) :: tol
  character(len=*),parameter :: subnam='dtt_svd0'
  double precision,dimension(:),allocatable :: s,b,u,work
  integer :: l,m,k,nn,mm,mn,info,i,lwork
  integer, intent(in), optional :: rmax
  integer,pointer :: r(:)
  double precision,external :: dnrm2

  l=tt%l; m=l+size(n)-1;tt%m=m;tt%n(l:m)=n; r=>tt%r; r(l-1)=1;r(m)=1
  nn=product(n)
  lwork=64*nn
  if(.not.lwork.gt.0)lwork=16*nn
  if(.not.lwork.gt.0)lwork=4*nn
  if(.not.lwork.gt.0)then;write(*,*)subnam,': nn, lwork: ',nn,lwork;stop;endif
  allocate(b(nn),u(nn),s(nn),work(lwork),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif
  if(dnrm2(nn,a,1).eq.0.d0)then;write(*,*)subnam,': zero norm input';tt%r=0;call dealloc(tt);return;endif

  call dcopy(nn,a,1,b,1)

  do k=m,l+1,-1
   mm=r(l-1)*product(tt%n(l:k-1)); nn=tt%n(k)*r(k); mn=min(mm,nn)
   call dgesvd('o','s',mm,nn,b,mm,s,b,1,u,mn,work,lwork,info)
   if(info.ne.0)then;write(*,*)subnam,': info: ',info;stop;endif
   !r(k-1)=chop(s(1:mn),tol/dsqrt(dble(m-l+1)))
   if(s(1).ne.0.d0)then 
       r(k-1)=chop(s(1:mn), tol/dsqrt(dble(m-l)))
   else
       r(k-1)=0
   endif
   if (present(rmax)) then
       r(k-1) = min(r(k-1), rmax)
   end if
   do i=1,r(k-1); call dscal(mm,s(i),b(1+(i-1)*mm),1);enddo
   if(associated(tt%u(k)%p))deallocate(tt%u(k)%p)
   allocate(tt%u(k)%p(r(k-1),tt%n(k),r(k)))
   call d2submat(r(k-1),nn,u,mn,tt%u(k)%p)
  end do
  if(associated(tt%u(l)%p))deallocate(tt%u(l)%p)
  allocate(tt%u(l)%p(r(l-1),tt%n(l),r(l)))
  call dcopy(r(l-1)*tt%n(l)*r(l),b,1,tt%u(l)%p,1)
  deallocate(work,b,u,s)
 end subroutine
 subroutine ztt_svd0(n,a,tt,tol,rmax)
  implicit none
  integer,intent(in) :: n(:)
  complex(8),intent(in) :: a(*)
  type(ztt),intent(inout),target :: tt
  double precision,intent(in) :: tol
  character(len=*),parameter :: subnam='ztt_svd0'
  double precision,dimension(:),allocatable :: s,rwork
  complex(8),dimension(:),allocatable :: b,u,work
  integer :: l,m,k,nn,mm,mn,info,i,lwork
  integer, intent(in), optional :: rmax
  integer,pointer :: r(:)
  double precision,external :: dznrm2

  l=tt%l; m=l+size(n)-1;tt%m=m;tt%n(l:m)=n; r=>tt%r; r(l-1)=1;r(m)=1
  nn=product(n)
  lwork=64*nn
  if(.not.lwork.gt.0)lwork=16*nn
  if(.not.lwork.gt.0)then;write(*,*)subnam,': nn, lwork: ',nn,lwork;stop;endif
  allocate(b(nn),u(nn),s(nn),work(lwork),rwork(8*nn),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif
  if(dznrm2(nn,a,1).eq.0.d0)then;write(*,*)subnam,': zero norm input';tt%r=0;call dealloc(tt);return;endif

  call zcopy(nn,a,1,b,1)
  do k=m,l+1,-1
   mm=r(l-1)*product(tt%n(l:k-1)); nn=tt%n(k)*r(k); mn=min(mm,nn)
   call zgesvd('o','s',mm,nn,b,mm,s,b,1,u,mn,work,lwork,rwork,info)
   if(info.ne.0)then;write(*,*)subnam,': info: ',info;stop;endif
   if(s(1).ne.0.d0)then 
       r(k-1)=chop(s(1:mn), tol/dsqrt(dble(m-l)))
   else
       r(k-1)=0
   endif
   if (present(rmax)) then
       r(k-1) = min(r(k-1), rmax)
   end if
   do i=1,r(k-1); call zdscal(mm,s(i),b(1+(i-1)*mm),1);enddo
   if(associated(tt%u(k)%p))deallocate(tt%u(k)%p)
   allocate(tt%u(k)%p(r(k-1),tt%n(k),r(k)))
   call z2submat(r(k-1),nn,u,mn,tt%u(k)%p)
  end do
  if(associated(tt%u(l)%p))deallocate(tt%u(l)%p)
  allocate(tt%u(l)%p(r(l-1),tt%n(l),r(l)))
  call zcopy(r(l-1)*tt%n(l)*r(l),b,1,tt%u(l)%p,1)
  deallocate(work,b,u,rwork,s)
 end subroutine




! GROUP
 subroutine dtt_group(arg,grp,side)
  !grp=[grp arg]
  implicit none
  type(dtt),intent(in),target :: arg
  type(dtt),intent(inout),target :: grp
  integer,intent(in),optional :: side
  character(len=*),parameter :: subnam='dtt_group'
  type(dtt) :: z
  integer,pointer :: r(:),q(:),n(:)
  integer :: l,m,k,sid,mm,ll,i,j

  if(arg%l.ne.grp%l .or. arg%m.ne.grp%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=arg%l; m=arg%m
  if(.not.all(arg%n(l:m)==grp%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>arg%n;r=>grp%r;q=>arg%r
  if(present(side))then;sid=side;else;if(r(l-1).ge.r(m))then;sid=0;else;sid=1;endif;endif

  z%l=l; z%m=m; z%n=n; z%r=0
  select case(sid)
   case(0)
    if(r(m).ne.q(m))then;write(*,*)subnam,': right border ranks mismatch:',r(m),q(m);stop;endif
    z%r(l-1:m-1)=r(l-1:m-1)+q(l-1:m-1); z%r(m)=r(m)
   case(1)
    if(r(l-1).ne.q(l-1))then;write(*,*)subnam,': left border ranks mismatch:',r(l-1),q(l-1);stop;endif
    z%r(l-1)=r(l-1); z%r(l:m)=r(l:m)+q(l:m)
   case default
    write(*,*)subnam,': illegal side:',sid; stop
  end select
  call alloc(z)

  if(sid.eq.1)then
   forall(j=1:r(l)) z%u(l)%p(:,:,     j)=grp%u(l)%p(:,:,j)
   forall(j=1:q(l)) z%u(l)%p(:,:,r(l)+j)=arg%u(l)%p(:,:,j)
   ll=l+1;mm=m
  endif
  if(sid.eq.0)then
   forall(i=1:r(m-1)) z%u(m)%p(       i,:,:)=grp%u(m)%p(i,:,:)
   forall(i=1:q(m-1)) z%u(m)%p(r(m-1)+i,:,:)=arg%u(m)%p(i,:,:)
   ll=l;mm=m-1
  endif

  do k=ll,mm
   z%u(k)%p=0.d0
   forall(i=1:r(k-1),j=1:r(k)) z%u(k)%p(       i,:,     j)=grp%u(k)%p(i,:,j)
   forall(i=1:q(k-1),j=1:q(k)) z%u(k)%p(r(k-1)+i,:,r(k)+j)=arg%u(k)%p(i,:,j)
  end do
  grp=z
  call dealloc(z)
 end subroutine

 subroutine ztt_group(arg,grp,side)
  !grp=[grp arg]
  implicit none
  type(ztt),intent(in),target :: arg
  type(ztt),intent(inout),target :: grp
  integer,intent(in),optional :: side
  character(len=*),parameter :: subnam='ztt_group'
  type(ztt) :: z
  integer,pointer :: r(:),q(:),n(:)
  integer :: l,m,k,sid,mm,ll,i,j

  if(arg%l.ne.grp%l .or. arg%m.ne.grp%m)then;write(*,*)subnam,': length mismatch';stop;endif
  l=arg%l; m=arg%m
  if(.not.all(arg%n(l:m)==grp%n(l:m)))then;write(*,*)subnam,': size mismatch';stop;endif
  n=>arg%n;r=>grp%r;q=>arg%r
  if(present(side))then;sid=side;else;if(r(l-1).ge.r(m))then;sid=0;else;sid=1;endif;endif

  z%l=l; z%m=m; z%n=n; z%r=0
  select case(sid)
   case(0)
    if(r(m).ne.q(m))then;write(*,*)subnam,': right border ranks mismatch:',r(m),q(m);stop;endif
    z%r(l-1:m-1)=r(l-1:m-1)+q(l-1:m-1); z%r(m)=r(m)
   case(1)
    if(r(l-1).ne.q(l-1))then;write(*,*)subnam,': left border ranks mismatch:',r(l-1),q(l-1);stop;endif
    z%r(l-1)=r(l-1); z%r(l:m)=r(l:m)+q(l:m)
   case default
    write(*,*)subnam,': illegal side:',sid; stop
  end select
  call alloc(z)

  if(sid.eq.1)then
   forall(j=1:r(l)) z%u(l)%p(:,:,     j)=grp%u(l)%p(:,:,j)
   forall(j=1:q(l)) z%u(l)%p(:,:,r(l)+j)=arg%u(l)%p(:,:,j)
   ll=l+1;mm=m
  endif
  if(sid.eq.0)then
   forall(i=1:r(m-1)) z%u(m)%p(       i,:,:)=grp%u(m)%p(i,:,:)
   forall(i=1:q(m-1)) z%u(m)%p(r(m-1)+i,:,:)=arg%u(m)%p(i,:,:)
   ll=l;mm=m-1
  endif

  do k=ll,mm
   z%u(k)%p=0.d0
   forall(i=1:r(k-1),j=1:r(k)) z%u(k)%p(       i,:,     j)=grp%u(k)%p(i,:,j)
   forall(i=1:q(k-1),j=1:q(k)) z%u(k)%p(r(k-1)+i,:,r(k)+j)=arg%u(k)%p(i,:,j)
  end do
  grp=z
  call dealloc(z)
 end subroutine



! ELEM TIJK
 pure double precision function dtt_ijk(arg,ind) result (a)
  implicit none
  type(dtt),intent(in) :: arg
  integer,intent(in) :: ind(:)
  character(len=*),parameter :: subnam='dtt_ijk'
  integer :: info,i,l,m,n(tt_size),r(0:tt_size)
  double precision,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;a=-3.d0;return;endif
  if(r(l-1).ne.1 .or. r(m).ne.1)then;a=-4.d0;return;endif
  allocate(x(r(m-1),r(m)),stat=info)
  if(info.ne.0)then;a=-1.d0;return;endif
  x=arg%u(m)%p(:,ind(m-l+1),:)
  do i=m-1,l,-1
   allocate(y(r(i-1),r(i)),z(r(i-1),r(m)),stat=info)
   if(info.ne.0)then;a=-2.d0;return;endif
   y=arg%u(i)%p(:,ind(i-l+1),:)
   z=matmul(y,x)
   deallocate(x,y); x=>z; nullify(z)
  end do
  a=x(1,1)
  deallocate(x)
 end function
 pure complex(8) function ztt_ijk(arg,ind) result (a)
  implicit none
  type(ztt),intent(in) :: arg
  integer,intent(in) :: ind(:)
  character(len=*),parameter :: subnam='ztt_ijk'
  integer :: info,i,l,m,n(tt_size),r(0:tt_size)
  complex(8),pointer :: x(:,:),y(:,:),z(:,:)
  complex(8),parameter :: one=(1.d0,0.d0),zero=(0.d0,0.d0),im1=(0.d0,1.d0)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;a=-3*one;return;endif
  if(r(l-1).ne.1 .or. r(m).ne.1)then;a=-4*one;return;endif
  allocate(x(r(m-1),r(m)),stat=info)
  if(info.ne.0)then;a=-one;return;endif
  x=arg%u(m)%p(:,ind(m-l+1),:)
  do i=m-1,l,-1
   allocate(y(r(i-1),r(i)),z(r(i-1),r(m)),stat=info)
   if(info.ne.0)then;a=-2*one;return;endif
   y=arg%u(i)%p(:,ind(i-l+1),:)
   z=matmul(y,x)
   deallocate(x,y); x=>z; nullify(z)
  end do
  a=x(1,1)
  deallocate(x)
 end function


 subroutine dtt_elem(arg,ind,a)
  implicit none
  type(dtt),intent(in):: arg
  integer,intent(in) :: ind(:)
  double precision,intent(out) :: a(*)
  character(len=*),parameter :: subnam='dtt_elem'
  integer :: info,i,l,m,n(tt_size),r(0:tt_size)
  double precision,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r
  if(any(ind(1:m-l+1)<=0).or.any(ind(1:m-l+1)>n(l:m)))then;write(*,*)subnam,': wrong index: ',ind;stop;endif
  allocate(x(r(l-1),r(l)),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
  x=arg%u(l)%p(:,ind(1),:)
  do i=l+1,m
   allocate(y(r(i-1),r(i)),z(r(l-1),r(i)),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
   y=arg%u(i)%p(:,ind(i-l+1),:)
   z=matmul(x,y)
   deallocate(x,y); x=>z; nullify(z)
  end do
  call dcopy(r(l-1)*r(m),x,1,a,1)
  deallocate(x)
 end subroutine

! VAL
 double precision function dtt_value(arg,x) result (val)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in) :: x(:)
  integer :: l,m,r(0:tt_size),n(tt_size), id,dd,i,j,pos,mm,ind(tt_size)=0
  double precision :: xx
  val=0.d0
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; dd=size(x)
  if(l.gt.m)return
  mm=(m-l+1)/dd
  do id=1,dd
   xx=x(id)
   if(xx.lt.0.d0)return
   if(xx.gt.1.d0)xx=xx-int(xx)
   do j=1,mm
    pos=l+(id-1)*mm+mm-j
    i=int(n(pos)*xx)
    if(i.eq.n(pos))i=n(pos)-1
    ind(pos-l+1)=i+1
    xx=xx*n(pos)-i
   end do
  end do
  val=dtt_ijk(arg,ind)
  !write(*,'(3f20.12,1x,e20.12)') x,val
  !write(*,'(127i1)')(mod(i,10),i=1,127)
  !write(*,'(127i1)')ind
 end function
 complex(8) function ztt_value(arg,x) result (val)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in) :: x(:)
  integer :: l,m,r(0:tt_size),n(tt_size), id,dd,i,j,pos,mm,ind(tt_size)=0
  double precision :: xx
  val=0.d0
  l=arg%l; m=arg%m; r=arg%r; n=arg%n; dd=size(x)
  if(l.gt.m)return
  mm=(m-l+1)/dd
  do id=1,dd
   xx=x(id)
   if(xx.lt.0.d0)return
   if(xx.gt.1.d0)xx=xx-int(xx)
   do j=1,mm
    pos=l+(id-1)*mm+mm-j
    i=int(n(pos)*xx)
    if(i.eq.n(pos))i=n(pos)-1
    ind(pos-l+1)=i+1
    xx=xx*n(pos)-i
   end do
  end do
  val=ztt_ijk(arg,ind)
 end function

 double precision function dtt_value0(arg,x) result (val)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in) :: x
  double precision :: xx(1)
  xx=x; val=dtt_value(arg,xx)
 end function
 complex(8) function ztt_value0(arg,x) result (val)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in) :: x
  double precision :: xx(1)
  xx=x; val=ztt_value(arg,xx)
 end function

! SETI
 subroutine dtt_seti(arg,pos,val)
  implicit none
  type(dtt),intent(inout) :: arg
  integer,intent(in) :: pos,val
  character(len=*),parameter :: subnam='dtt_seti'
  integer :: l,m,n,p,q,i,j
  double precision,allocatable :: a(:,:)
  l=arg%l; m=arg%m
  if(.not.(l.le.pos .and. pos.le.m))then;write(*,*)subnam,': pos not between l and m: ',pos,l,m;return;endif
  n=arg%n(pos); p=arg%r(pos-1);q=arg%r(pos)
  if(.not.(1.le.val .and. val.le.n))then;write(*,*)subnam,': val not between 1 and n: ',val,n;return;endif
  allocate(a(p,q))
  forall(i=1:p,j=1:q)a(i,j)=arg%u(pos)%p(i,val,j)
  deallocate(arg%u(pos)%p)
  allocate(arg%u(pos)%p(p,1,q))
  call dcopy(p*q,a,1,arg%u(pos)%p,1)
  deallocate(a)
  arg%n(pos)=1
 end subroutine
 subroutine ztt_seti(arg,pos,val)
  implicit none
  type(ztt),intent(inout) :: arg
  integer,intent(in) :: pos,val
  character(len=*),parameter :: subnam='ztt_seti'
  integer :: l,m,n,p,q,i,j
  complex(8),allocatable :: a(:,:)
  l=arg%l; m=arg%m
  if(.not.(l.le.pos .and. pos.le.m))then;write(*,*)subnam,': pos not between l and m: ',pos,l,m;return;endif
  n=arg%n(pos); p=arg%r(pos-1);q=arg%r(pos)
  if(.not.(1.le.val .and. val.le.n))then;write(*,*)subnam,': val not between 1 and n: ',val,n;return;endif
  allocate(a(p,q))
  forall(i=1:p,j=1:q)a(i,j)=arg%u(pos)%p(i,val,j)
  deallocate(arg%u(pos)%p)
  allocate(arg%u(pos)%p(p,1,q))
  call zcopy(p*q,a,1,arg%u(pos)%p,1)
  deallocate(a)
  arg%n(pos)=1
 end subroutine

! SUMALL
 double precision function dtt_sumall(arg) result(val)
  implicit none
  type(dtt),intent(in):: arg
  character(len=*),parameter :: subnam='dtt_sumall'
  integer :: info,i,j,p,q,l,m,n(tt_size),r(0:tt_size)
  double precision,pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r; val=0.d0
  if(r(l-1).gt.1 .or. r(m).gt.1) then; write(*,*)subnam,': matrix-valued output does not fit!';return;endif
  allocate(x(r(l-1),r(l)),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
  forall(i=1:r(l-1),j=1:r(l))x(i,j)=0.d0
  do p=1,n(l); forall(i=1:r(l-1),j=1:r(l))x(i,j)=x(i,j)+arg%u(l)%p(i,p,j); enddo
  do q=l+1,m
   allocate(y(r(q-1),r(q)),z(r(l-1),r(q)),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
   forall(i=1:r(q-1),j=1:r(q))y(i,j)=0.d0
   do p=1,n(q); forall(i=1:r(q-1),j=1:r(q))y(i,j)=y(i,j)+arg%u(q)%p(i,p,j); enddo
   z=matmul(x,y)
   deallocate(x,y); x=>z; nullify(z)
  end do
  val=x(1,1)
  deallocate(x)
 end function
 complex(8) function ztt_sumall(arg) result(val)
  implicit none
  type(ztt),intent(in):: arg
  character(len=*),parameter :: subnam='ztt_sumall'
  integer :: info,i,j,p,q,l,m,n(tt_size),r(0:tt_size)
  complex(8),pointer :: x(:,:),y(:,:),z(:,:)
  l=arg%l;m=arg%m;n=arg%n;r=arg%r; val=(0.d0,0.d0)
  if(r(l-1).gt.1 .or. r(m).gt.1) then; write(*,*)subnam,': matrix-valued output does not fit!';return;endif
  allocate(x(r(l-1),r(l)),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
  forall(i=1:r(l-1),j=1:r(l))x(i,j)=(0.d0,0.d0)
  do p=1,n(l); forall(i=1:r(l-1),j=1:r(l))x(i,j)=x(i,j)+arg%u(l)%p(i,p,j); enddo
  do q=l+1,m
   allocate(y(r(q-1),r(q)),z(r(l-1),r(q)),stat=info)
   if(info.ne.0)then;write(*,*)subnam,': allocation error: ',info;stop;endif
   forall(i=1:r(q-1),j=1:r(q))y(i,j)=(0.d0,0.d0)
   do p=1,n(q); forall(i=1:r(q-1),j=1:r(q))y(i,j)=y(i,j)+arg%u(q)%p(i,p,j); enddo
   z=matmul(x,y)
   deallocate(x,y); x=>z; nullify(z)
  end do
  val=x(1,1)
 end function

! NUMEL
 double precision function dtt_numel(arg) result (s)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: l,m,i
  s=0.d0; l=arg%l; m=arg%m; if(l.gt.m)return
  s=1.d0; do i=l,m;s=s*arg%n(i); enddo
  return
 end function
 double precision function ztt_numel(arg) result (s)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: l,m,i
  s=0.d0; l=arg%l; m=arg%m; if(l.gt.m)return
  s=1.d0; do i=l,m;s=s*arg%n(i); enddo
  return
 end function

! NIP (pack)
 subroutine dtt_nip(arg,ind)
  implicit none
  type(dtt),intent(inout),target :: arg
  integer,intent(in),optional :: ind(:)
  character(len=*),parameter :: subnam='dtt_nip'
  integer :: l,m,p,k,maxr,maxn,info
  integer,pointer :: r(:),n(:)
  double precision,allocatable :: tmp(:)
  l=arg%l; m=arg%m; r=>arg%r; n=>arg%n
  maxn=maxval(n(l:m)); maxr=maxval(r(l-1:m))
  allocate(tmp(maxr*maxn*maxr),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate tmp: ',info;stop;endif

  if(present(ind))then
   if(any(ind(1:m-l+1)<0) .or. any(ind(1:m-l+1)>n(l:m)))then
    write(*,*)subnam,': illegal elements in ind'
    write(*,*)ind(1:m-l+1)
    call say(arg)
    stop
   end if
   do p=l,m
    if(ind(p-l+1).ne.0)then
     do k=1,r(p)
      call dcopy(r(p-1),arg%u(p)%p(1,ind(p-l+1),k),1,tmp(1+(k-1)*r(p-1)),1)
     end do
     deallocate(arg%u(p)%p)
     allocate(arg%u(p)%p(r(p-1),1,r(p)))
     call dcopy(r(p-1)*r(p),tmp,1,arg%u(p)%p,1)
     n(p)=1
    end if
   end do
  end if

  do p=m,l+1,-1
   if(n(p).eq.1 .or. n(p-1).eq.1)then
    call dgemm('n','n',r(p-2)*n(p-1),n(p)*r(p),r(p-1),1.d0,arg%u(p-1)%p,r(p-2)*n(p-1),arg%u(p)%p,r(p-1),0.d0,tmp,r(p-2)*n(p-1))
    deallocate(arg%u(p-1)%p, arg%u(p)%p)
    allocate(arg%u(p-1)%p(r(p-2),n(p-1)*n(p),r(p)),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot allocate carriage: ',p,info;stop;endif
    call dcopy(r(p-2)*n(p-1)*n(p)*r(p),tmp,1,arg%u(p-1)%p,1)
    n(p-1)=n(p-1)*n(p);r(p-1)=r(p); n(p)=0;r(p)=0
   end if
  end do
  deallocate(tmp)

  k=l
  do p=l,m
   if(n(p).gt.0)then
    if(p.ne.k)then
     if(associated(arg%u(k)%p))then;write(*,*)subnam,': position is associated, it shouldnot happen: ',p,k;stop;endif
     allocate(arg%u(k)%p(r(k-1),n(p),r(p)))
     if( size(arg%u(p)%p,1).ne.r(k-1) .or. size(arg%u(p)%p,2).ne.n(p) .or. size(arg%u(p)%p,3).ne.r(p) )then
      write(*,*)subnam,': size mismatch! '
      write(*,*)size(arg%u(p)%p,1),size(arg%u(p)%p,2),size(arg%u(p)%p,3)
      write(*,*)r(k-1),n(p),r(p)
      stop
     end if
     call dcopy(r(k-1)*n(p)*r(p), arg%u(p)%p,1,arg%u(k)%p,1)
     deallocate(arg%u(p)%p)
     n(k)=n(p);r(k)=r(p);n(p)=0;r(p)=0
    end if
    k=k+1
   end if
  end do
  arg%m=k-1
 end subroutine
 subroutine ztt_nip(arg,ind)
  implicit none
  type(ztt),intent(inout),target :: arg
  integer,intent(in),optional :: ind(:)
  character(len=*),parameter :: subnam='ztt_nip'
  complex(8),parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
  integer :: l,m,p,k,maxr,maxn,info
  integer,pointer :: r(:),n(:)
  complex(8),allocatable :: tmp(:)
  l=arg%l; m=arg%m; r=>arg%r; n=>arg%n
  maxn=maxval(n(l:m)); maxr=maxval(r(l-1:m))
  allocate(tmp(maxr*maxn*maxr),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate tmp: ',info;stop;endif

  if(present(ind))then
   if(any(ind(1:m-l+1)<0) .or. any(ind(1:m-l+1)>n(l:m)))then
    write(*,*)subnam,': illegal elements in ind'
    write(*,*)ind(1:m-l+1)
    call say(arg)
    stop
   end if
   do p=l,m
    if(ind(p-l+1).ne.0)then
     do k=1,r(p)
      call zcopy(r(p-1),arg%u(p)%p(1,ind(p-l+1),k),1,tmp(1+(k-1)*r(p-1)),1)
     end do
     deallocate(arg%u(p)%p)
     allocate(arg%u(p)%p(r(p-1),1,r(p)))
     call zcopy(r(p-1)*r(p),tmp,1,arg%u(p)%p,1)
     n(p)=1
    end if
   end do
  end if

  do p=m,l+1,-1
   if(n(p).eq.1 .or. n(p-1).eq.1)then
    call zgemm('n','n',r(p-2)*n(p-1),n(p)*r(p),r(p-1),one,arg%u(p-1)%p,r(p-2)*n(p-1),arg%u(p)%p,r(p-1),zero,tmp,r(p-2)*n(p-1))
    deallocate(arg%u(p-1)%p, arg%u(p)%p)
    allocate(arg%u(p-1)%p(r(p-2),n(p-1)*n(p),r(p)),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot allocate carriage: ',p,info;stop;endif
    call zcopy(r(p-2)*n(p-1)*n(p)*r(p),tmp,1,arg%u(p-1)%p,1)
    n(p-1)=n(p-1)*n(p);r(p-1)=r(p); n(p)=0;r(p)=0
   end if
  end do
  deallocate(tmp)

  k=l
  do p=l,m
   if(n(p).gt.0)then
    if(p.ne.k)then
     if(associated(arg%u(k)%p))then;write(*,*)subnam,': position is associated, it shouldnot happen: ',p,k;stop;endif
     allocate(arg%u(k)%p(r(k-1),n(p),r(p)))
     if( size(arg%u(p)%p,1).ne.r(k-1) .or. size(arg%u(p)%p,2).ne.n(p) .or. size(arg%u(p)%p,3).ne.r(p) )then
      write(*,*)subnam,': size mismatch! '
      write(*,*)size(arg%u(p)%p,1),size(arg%u(p)%p,2),size(arg%u(p)%p,3)
      write(*,*)r(k-1),n(p),r(p)
      stop
     end if
     call zcopy(r(k-1)*n(p)*r(p), arg%u(p)%p,1,arg%u(k)%p,1)
     deallocate(arg%u(p)%p)
     n(k)=n(p);r(k)=r(p);n(p)=0;r(p)=0
    end if
    k=k+1
   end if
  end do
  arg%m=k-1
 end subroutine


! FULL
 subroutine dtt_full(arg,a,alpha,beta,part,ind)
  ! A = [beta]*A + [alpha]*FULL(TT), TT = arg([l:m]), l,m=[part]
  ! A size r(l-1) *n(l)*...*n(m)* r(m)
  implicit none
  type(dtt),intent(in),target :: arg
  double precision,intent(inout) :: a(*)
  double precision,intent(in),optional :: alpha,beta
  integer,intent(in),optional :: part(2),ind(:)
  character(len=*),parameter :: subnam='dtt_full'
  double precision :: alp
  type(pointd) :: p(0:1)
  integer,pointer :: r(:),n(:)
  integer :: l,m,na,nb,mem,rr,info,i,j,pp
  integer,allocatable :: ii(:),nn(:)
  double precision,allocatable :: q(:)

  if(present(alpha))then;alp=alpha;else;alp=1.d0;endif
  if(present(part))then;l=part(1);m=part(2);else;l=arg%l;m=arg%m;endif
  r=>arg%r; n=>arg%n

  if(l.gt.m)then
   write(*,*)subnam,': empty input';
   do i=1,r(l-1) * product(nn(l:m)) * r(m); a(i)=0.d0; enddo
   return
  end if

  allocate(ii(l:m),nn(l:m)); ii=0; nn(l:m)=arg%n(l:m)
  if(present(ind))then
   do i=l,m
    ii(i)=ind(i-l+1)
    if(ii(i).lt.0 .or. ii(i).gt.n(i))then;write(*,*)subnam,': invalid ind:',ind(1:m-l+1);stop;endif
    if(ii(i).ne.0)nn(i)=1
   end do
  end if

  na=r(l-1) * product(nn(l:m)) * r(m)
  mem=na; rr=1
  do i=l,m
   if(.not.(r(l-1)*product(nn(l:i))*r(i).gt.0))then;write(*,*)subnam,': oversized mem';stop;endif
   mem=max(mem, r(l-1)*product(nn(l:i))*r(i) )
   rr=max(rr,r(i-1)*r(i))
  end do
  allocate(p(0)%p(mem),p(1)%p(mem),q(rr),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif

  if(ii(l).eq.0)then
   call dcopy(r(l-1)*n(l)*r(l),arg%u(l)%p,1,p(0)%p,1)
  else
   do j=1,r(l);call dcopy(r(l-1),arg%u(l)%p(1,ii(l),j),1,p(0)%p(1+(j-1)*r(l-1)),1); enddo
  end if
  pp=0
  do i=l+1,m
   nb=r(l-1) * product(nn(l:i-1))
   if(ii(i).eq.0)then
    if(nb*n(i)*r(i).gt.mem)then;write(*,*)subnam,': nb-by-n-by-r > mem: ',nb,n(i),r(i),mem;stop;endif
    call dgemm('n','n',nb,n(i)*r(i),r(i-1),1.d0,p(pp)%p,nb,arg%u(i)%p,r(i-1),0.d0,p(1-pp)%p,nb)
   else
    if(nb*r(i).gt.mem)then;write(*,*)subnam,': nb-by-r > mem: ',nb,r(i),mem;stop;endif
    do j=1,r(i);call dcopy(r(i-1),arg%u(i)%p(1,ii(i),j),1,q(1+(j-1)*r(i-1)),1);enddo
    call dgemm('n','n',nb,r(i),r(i-1),1.d0,p(pp)%p,nb,q,r(i-1),0.d0,p(1-pp)%p,nb)
   end if
   pp=1-pp
  end do

  if(present(beta))then
   call dscal(na,beta,a,1)
   call daxpy(na,alp,p(pp)%p,1,a,1)
  else
   call dcopy(na,p(pp)%p,1,a,1)
  endif
  deallocate(p(0)%p,p(1)%p,q,ii,nn)
 end subroutine
 subroutine ztt_full(arg,a,alpha,beta,part,ind)
  ! A = [beta]*A + [alpha]*FULL(TT), TT = arg([l:m]), l,m=[part]
  ! A sizes r(l-1) *n(l)*...*n(m)* r(m)
  implicit none
  type(ztt),intent(in),target :: arg
  complex(8),intent(inout) :: a(*)
  complex(8),intent(in),optional :: alpha,beta
  integer,intent(in),optional :: part(2),ind(:)
  character(len=*),parameter :: subnam='ztt_full'
  complex(8),parameter :: one=(1.d0,0.d0),zero=(0.d0,0.d0)
  complex(8) :: alp
  type(pointz) :: p(0:1)
  integer,pointer :: r(:),n(:)
  integer :: l,m,na,nb,mem,rr,info,i,j,pp
  integer,allocatable :: ii(:),nn(:)
  complex(8),allocatable :: q(:)

  if(present(alpha))then;alp=alpha;else;alp=one;endif
  if(present(part))then;l=part(1);m=part(2);else;l=arg%l;m=arg%m;endif
  r=>arg%r; n=>arg%n

  allocate(ii(l:m),nn(l:m)); ii=0; nn(l:m)=arg%n(l:m)
  if(present(ind))then
   do i=l,m
    ii(i)=ind(i-l+1)
    if(ii(i).lt.0 .or. ii(i).gt.n(i))then;write(*,*)subnam,': invalid ind:',ind(1:m-l+1);stop;endif
    if(ii(i).ne.0)nn(i)=1
   end do
  end if

  na=r(l-1) * product(nn(l:m)) * r(m)
  mem=na; rr=1
  do i=l,m
   if(.not.(r(l-1)*product(nn(l:i))*r(i).gt.0))then;write(*,*)subnam,': oversized mem';stop;endif
   mem=max(mem, r(l-1)*product(nn(l:i))*r(i) )
   rr=max(rr,r(i-1)*r(i))
  end do
  allocate(p(0)%p(mem),p(1)%p(mem),q(rr),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': not enough memory';stop;endif

  if(ii(l).eq.0)then
   call zcopy(r(l-1)*n(l)*r(l),arg%u(l)%p,1,p(0)%p,1)
  else
   do j=1,r(l);call zcopy(r(l-1),arg%u(l)%p(1,ii(l),j),1,p(0)%p(1+(j-1)*r(l-1)),1); enddo
  end if
  pp=0
  do i=l+1,m
   nb=r(l-1) * product(nn(l:i-1))
   if(ii(i).eq.0)then
    if(nb*n(i)*r(i).gt.mem)then;write(*,*)subnam,': nb-by-n-by-r > mem: ',nb,n(i),r(i),mem;stop;endif
    call zgemm('n','n',nb,n(i)*r(i),r(i-1),one,p(pp)%p,nb,arg%u(i)%p,r(i-1),zero,p(1-pp)%p,nb)
   else
    if(nb*r(i).gt.mem)then;write(*,*)subnam,': nb-by-r > mem: ',nb,r(i),mem;stop;endif
    do j=1,r(i);call zcopy(r(i-1),arg%u(i)%p(1,ii(i),j),1,q(1+(j-1)*r(i-1)),1);enddo
    call dgemm('n','n',nb,r(i),r(i-1),one,p(pp)%p,nb,q,r(i-1),zero,p(1-pp)%p,nb)
   end if
   pp=1-pp
  end do

  if(present(beta))then
   call zscal(na,beta,a,1)
   call zaxpy(na,alp,p(pp)%p,1,a,1)
  else
   call zcopy(na,p(pp)%p,1,a,1)
  endif
  deallocate(p(0)%p,p(1)%p,q,ii,nn)
 end subroutine


!ALLOC
 subroutine dtt_memchk(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  integer :: i
  character(len=tt_size) :: a
  if(arg%m<arg%l)return
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))then
    if(size(arg%u(i)%p).gt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='>'
    else if (size(arg%u(i)%p).lt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='<'
    else
     a(i:i)='='
    endif
   else
    a(i:i)='-'
   end if
  end do
  write(*,'(a,i2,a,i2,2a)') 'dtt[',arg%l,':', arg%m,']: ',a(arg%l:arg%m)
 end subroutine
 subroutine ztt_memchk(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  integer :: i
  character(len=tt_size) :: a
  if(arg%m<arg%l)return
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))then
    if(size(arg%u(i)%p).gt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='>'
    else if (size(arg%u(i)%p).lt.arg%r(i-1)*arg%n(i)*arg%r(i))then
     a(i:i)='<'
    else
     a(i:i)='='
    endif
   else
    a(i:i)='-'
   end if
  end do
  write(*,'(a,i2,a,i2,2a)') 'ztt[',arg%l,':', arg%m,']: ',a(arg%l:arg%m)
 end subroutine

 subroutine dtt_alloc(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  integer :: i,info
  if(arg%m<arg%l)return
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
   allocate(arg%u(i)%p(arg%r(i-1),arg%n(i),arg%r(i)), stat=info)
   if(info.ne.0)then;write(*,*)'TT allocate fail: no memory';stop;endif
  end do
 end subroutine
 subroutine ztt_alloc(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  integer :: i,info
  if(arg%m<arg%l)return
  do i=arg%l,arg%m
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
   allocate(arg%u(i)%p(arg%r(i-1),arg%n(i),arg%r(i)), stat=info)
   if(info.ne.0)then;write(*,*)'TT allocate fail: no memory';stop;endif
  end do
 end subroutine

 subroutine dtt_dealloc(arg)
  implicit none
  type(dtt),intent(inout) :: arg
  integer :: i
  do i=1,tt_size
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
  end do
 end subroutine
 subroutine ztt_dealloc(arg)
  implicit none
  type(ztt),intent(inout) :: arg
  integer :: i
  do i=1,tt_size
   if(associated(arg%u(i)%p))deallocate(arg%u(i)%p)
  end do
 end subroutine



! MUL
 type(dtt) function dttmul_dt(a,b) result(c)
  double precision,intent(in) :: a
  type(dtt),intent(in) :: b
  integer :: k,l,m
  l=b%l; m=b%m; c%l=l; c%m=m; c%n=b%n; c%r=b%r; call alloc(c)
  do k=l,m
   call dcopy(b%r(k-1)*b%n(k)*b%r(k),b%u(k)%p,1,c%u(k)%p,1)
  end do
  call dscal(b%r(l-1)*b%n(l)*b%r(l),a,c%u(l)%p,1)
 end function
 type(ztt) function zttmul_zt(a,b) result(c)
  complex(8),intent(in) :: a
  type(ztt),intent(in) :: b
  integer :: k,l,m
  l=b%l; m=b%m; c%l=l; c%m=m; c%n=b%n; c%r=b%r; call alloc(c)
  do k=l,m
   call zcopy(b%r(k-1)*b%n(k)*b%r(k),b%u(k)%p,1,c%u(k)%p,1)
  end do
  call zscal(b%r(l-1)*b%n(l)*b%r(l),a,c%u(l)%p,1)
 end function

! ASSIGN COPY
 subroutine dtt_assign(b,a)
  implicit none
  type(dtt),intent(inout) :: b
  type(dtt),intent(in) :: a
  integer :: k,l,m
  l=a%l;m=a%m
  b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
  do k=l,m; call dcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1); end do
 end subroutine
 subroutine ztt_assign(b,a)
  implicit none
  type(ztt),intent(inout) :: b
  type(ztt),intent(in) :: a
  integer :: k,l,m
  l=a%l;m=a%m
  b%l=l; b%m=m; b%n(l:m)=a%n(l:m); b%r(l-1:m)=a%r(l-1:m); call alloc(b)
  do k=l,m; call zcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(k)%p,1); end do
 end subroutine

 subroutine dtt_copy(a,b,low)
  implicit none
  type(dtt),intent(in) :: a
  type(dtt),intent(inout) :: b
  integer,intent(in),optional :: low
  integer :: k,l,m,ll,mm
  l=a%l;m=a%m; ll=b%l; if(present(low))ll=low; mm=ll-l+m
  b%l=ll; b%m=mm; b%n(ll:mm)=a%n(l:m); b%r(ll-1:mm)=a%r(l-1:m);
  if(.not.all(a%n(l:m)>0))return;if(.not.all(a%r(l-1:m)>0))return;call alloc(b)
  do k=l,m; call dcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(ll-l+k)%p,1); end do
 end subroutine
 subroutine ztt_copy(a,b,low)
  implicit none
  type(ztt),intent(in) :: a
  type(ztt),intent(inout) :: b
  integer,intent(in),optional :: low
  integer :: k,l,m,ll,mm
  l=a%l;m=a%m; ll=b%l; if(present(low))ll=low; mm=ll-l+m
  b%l=ll; b%m=mm; b%n(ll:mm)=a%n(l:m); b%r(ll-1:mm)=a%r(l-1:m);
  if(.not.all(a%n(l:m)>0))return;if(.not.all(a%r(l-1:m)>0))return;call alloc(b)
  do k=l,m; call zcopy(a%r(k-1)*a%n(k)*a%r(k),a%u(k)%p,1,b%u(ll-l+k)%p,1); end do
 end subroutine

 function dtt_dot(tt1,tt2) result(dt)
 use matrix_util, only: dtransp
 implicit none
 type(dtt), intent(in), target :: tt1,tt2
 integer i
 real(8) dt(tt1%r(tt1%l-1)*tt2%r(tt2%l-1)*tt1%r(tt1%m)*tt2%r(tt2%m))
 real(8), allocatable :: phi(:), res(:) 
 integer :: l,m,mem,C
 integer, pointer :: r1(:), r2(:), n(:)
 l = tt1 % l
 m = tt1 % m
 r1 => tt1 % r
 r2 => tt2 % r
 n  => tt1 % n
 mem = 0
 do i = l-1,m
    mem = max(mem,r1(i)*r2(i))
 end do
 C = r1(l-1)*r2(l-1)
 mem = mem * C
 allocate(phi(mem))
 mem = 0
 do i = l,m
    mem = max(mem,C*r1(i-1)*n(i)*r2(i))
 end do 
 allocate(res(mem))
 call eye(phi,C)
 do i = l,m
    !Multiply phi by smth
    !phi  phi(a0,b0,a1,b1)*conj(cr1(a1,i1,a2))*cr2(b1,i1,b2)
    !phi(C,a1,b1)*cr2(b1,i1,b2)->res(C,a1,i1,b2)->res(a1,i1,b2,C)*cr1(a1,i1,a2)->a2,b2,C->C
    call dgemm('n','n',C*r1(i-1),n(i)*r2(i),r2(i-1),1d0,phi,C*r1(i-1),tt2%u(i)%p,r2(i-1),0d0,res,C*r1(i-1))
    call dtransp(C,r1(i-1)*n(i)*r2(i),res)
    call dgemm('t','n',r1(i),r2(i)*C,r1(i-1)*n(i),1d0,tt1%u(i)%p,r1(i-1)*n(i),res,r1(i-1)*n(i),0d0,phi,r1(i))
    call dtransp(r1(i)*r2(i),C,phi)
 end do
 call dcopy(C*r1(m)*r2(m),phi,1,dt,1)
 deallocate(phi)
 deallocate(res)
 end function dtt_dot
 
 function ztt_dot(tt1,tt2) result(dt)
 use matrix_util, only: ztransp
 implicit none
 type(ztt), intent(in), target :: tt1,tt2
 integer i
 complex(8) dt(tt1%r(tt1%l-1)*tt2%r(tt2%l-1)*tt1%r(tt1%m)*tt2%r(tt2%m))
 complex(8), allocatable :: phi(:), res(:) 
 integer :: l,m,mem,C
 integer, pointer :: r1(:), r2(:), n(:)
 real(8) dznrm2
 l = tt1 % l
 m = tt1 % m
 r1 => tt1 % r
 r2 => tt2 % r
 n  => tt1 % n
 mem = 0
 do i = l-1,m
    mem = max(mem,r1(i)*r2(i))
 end do
 C = r1(l-1)*r2(l-1)
 mem = mem * C
 allocate(phi(mem))
 mem = 0
 do i = l,m
    mem = max(mem,C*r1(i-1)*n(i)*r2(i))
 end do 
 allocate(res(mem))
 call eye(phi,C)
 do i = l,m
    !Multiply phi by smth
    !phi  phi(a0,b0,a1,b1)*conj(cr1(a1,i1,a2))*cr2(b1,i1,b2)
    !phi(C,a1,b1)*cr2(b1,i1,b2)->res(C,a1,i1,b2)->res(a1,i1,b2,C)*cr1(a1,i1,a2)->a2,b2,C->C
    call zgemm('n','n',C*r1(i-1),n(i)*r2(i),r2(i-1),(1d0,0d0),phi,C*r1(i-1),tt2%u(i)%p,r2(i-1),(0d0,0d0),res,C*r1(i-1))
    call ztransp(C,r1(i-1)*n(i)*r2(i),res)
    call zgemm('c','n',r1(i),r2(i)*C,r1(i-1)*n(i),(1d0,0d0),tt1%u(i)%p,r1(i-1)*n(i),res,r1(i-1)*n(i),(0d0,0d0),phi,r1(i))
    call ztransp(r1(i)*r2(i),C,phi)
 end do
 call zcopy(C*r1(m)*r2(m),phi,1,dt,1)
 deallocate(phi)
 deallocate(res)
 end function ztt_dot

! NORM

double precision function dtt_norm(arg,tol) result (nrm)
  implicit none
  type(dtt),intent(in) :: arg
  double precision,intent(in),optional :: tol
  integer :: l,m
  type(dtt) :: tmp
  double precision,external :: dnrm2
  double precision :: t1,t2
  l=arg%l;m=arg%m; nrm=0.d0
  tmp=arg
  if(present(tol))then
   call svd(tmp,tol)
   nrm=dnrm2(size(tmp%u(l)%p),tmp%u(l)%p,1)
  else
   call ort(tmp)
   nrm=dnrm2(size(tmp%u(m)%p),tmp%u(m)%p,1)
  end if
  call dealloc(tmp)
 end function

 double precision function ztt_norm(arg,tol) result (nrm)
  implicit none
  type(ztt),intent(in) :: arg
  double precision,intent(in),optional :: tol
  integer :: l,m
  type(ztt) :: tmp
  double precision,external :: dznrm2
  l=arg%l;m=arg%m; nrm=0.d0
  tmp=arg
  if(present(tol))then
   call svd(tmp,tol)
   nrm=dznrm2(size(tmp%u(l)%p),tmp%u(l)%p,1)
  else
   call ort(tmp)
   nrm=dznrm2(size(tmp%u(m)%p),tmp%u(m)%p,1)
  end if
  call dealloc(tmp)
 end function

! SAY
 subroutine dtt_say(arg)
  implicit none
  type(dtt),intent(in) :: arg
  write(*,'(a,i2,a,i3,a,f6.2)') 'dtt[',arg%l,':', arg%m,']: rank ',rank(arg)
  write(*,'(a,1x,128i3)') 'n: ',arg%n(arg%l:arg%m)
  write(*,'(a,128i3)') 'r: ',arg%r(arg%l-1:arg%m)
 end subroutine
 subroutine ztt_say(arg)
  implicit none
  type(ztt),intent(in) :: arg
  write(*,'(a,i2,a,i3,a,f6.2)') 'ztt[',arg%l,':', arg%m,']: rank ',rank(arg)
  write(*,'(a,1x,128i3)') 'n: ',arg%n(arg%l:arg%m)
  write(*,'(a,128i3)') 'r: ',arg%r(arg%l-1:arg%m)
 end subroutine
 subroutine ttind_say(arg)
  implicit none
  type(ttind),intent(in) :: arg
  write(*,'(a,i2)') 'ttind, m=', arg%m
  write(*,'(a,1x,128i3)') 'n: ',arg%n(1:arg%m)
  write(*,'(a,128i3)') 'p: ',arg%p(1:arg%m)
 end subroutine
 subroutine dtt_sayfull(arg)
  implicit none
  type(dtt),intent(in) :: arg
  character(len=*),parameter :: subnam='dtt_sayfull'
  integer :: l,m,n,i,info
  double precision,allocatable :: a(:)
  l=arg%l; m=arg%m; if(l.gt.m)return
  n=arg%r(l-1)*product(arg%n(l:m))*arg%r(m); if(n.le.0)return
  if(n.gt.2**20)then;write(*,*)subnam,': too much to say';return;endif
  allocate(a(n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call dtt_full(arg,a)
  do i=1,n
   write(*,'(e22.14)') a(i)
  end do
  deallocate(a)
 end subroutine
 subroutine ztt_sayfull(arg)
  implicit none
  type(ztt),intent(in) :: arg
  character(len=*),parameter :: subnam='ztt_sayfull'
  integer :: l,m,n,i,info
  complex(8),allocatable :: a(:)
  l=arg%l; m=arg%m; if(l.gt.m)return
  n=arg%r(l-1)*product(arg%n(l:m))*arg%r(m); if(n.le.0)return
  if(n.gt.2**20)then;write(*,*)subnam,': too much to say';return;endif
  allocate(a(n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call ztt_full(arg,a)
  do i=1,n
   write(*,'(e22.14,1x,e22.14)') a(i)
  end do
  deallocate(a)
 end subroutine

! RANK
 double precision function dtt_rank(arg) result (r)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: l,m,i,a,b,d
  l=arg%l;m=arg%m;d=m-l+1
  if(d.le.0)then;r=-1.d0;return;endif
  if(d.eq.1)then;r= 0.d0;return;endif
  r=0.d0
  do i=arg%l,arg%m
   r=r+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
  if(r.eq.0.d0)return
  b=arg%r(l-1)*arg%n(l) + arg%n(m)*arg%r(m)
  if(d.eq.2)then;r=r/b;return;endif
  a=sum(arg%n(l+1:m-1))
  r=(dsqrt(b*b+4.d0*a*r)-b)/(2.d0*a)
  return
 end function
 double precision function ztt_rank(arg) result (r)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: l,m,i,a,b,d
  l=arg%l;m=arg%m;d=m-l+1
  if(d.le.0)then;r=-1.d0;return;endif
  if(d.eq.1)then;r= 0.d0;return;endif
  r=0.d0
  do i=arg%l,arg%m
   r=r+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
  if(r.eq.0.d0)return
  b=arg%r(l-1)*arg%n(l) + arg%n(m)*arg%r(m)
  if(d.eq.2)then;r=r/b;return;endif
  a=sum(arg%n(l+1:m-1))
  r=(dsqrt(b*b+4.d0*a*r)-b)/(2.d0*a)
  return
 end function

! INDEX
 pure type(ttind) function ttindex(i,n) result (ind)
  implicit none
  integer,intent(in) :: i,n(:)
  integer :: m
  ind%p=0; m=1; ind%p(m)=i-1;
  do while(m.lt.tt_size .and. m.le.size(n) .and. ind%p(m).ge.n(m))
   m=m+1
   ind%p(m)=ind%p(m-1)/n(m-1)
   ind%p(m-1)=ind%p(m-1)-ind%p(m)*n(m-1)
  enddo
  ind%m=m
  ind%p(1:m)=ind%p(1:m)+1
  m=min(size(n),tt_size)
  ind%n(1:m)=n(1:m)
 end function
 subroutine mindex(ind,n,i)
  implicit none
  integer,intent(in) :: ind, n(:)
  integer,intent(out) :: i(:)
  integer :: j
  if(size(i).gt.size(n)+1)then;i=-999;return;endif
  i(1)=ind-1
  do j=2,size(i)
   i(j)=i(j-1)/n(j-1)
   i(j-1)=i(j-1)-i(j)*n(j-1)
  enddo
  i=i+1
 end subroutine

 integer function find_ttind(n,x,y) result (pos)
  ! for sorted vector x(1) <= x(2) <= ... <= x(n) and value y find pos, s.t. x(pos) <= y < x(pos+1)
  implicit none
  integer,intent(in) :: n
  type(ttind),intent(in) :: x(n),y
  integer :: s,t,i
  logical :: key
  if(n.eq.0)then;pos=0;return;endif
  if(y.lt.x(1))then;pos=0;return;endif
  if(x(n).le.y)then;pos=n;return;endif
  s=1;t=n;pos=(t+s)/2
  do while(t-s.gt.1)
   if(y.lt.x(pos))then;t=pos;else;s=pos;end if
   pos=(s+t)/2
  enddo
  return
 end function
 subroutine push_ttind(n,x)
  implicit none
  integer,intent(in) :: n
  type(ttind),intent(inout) :: x(n+1)
  integer :: i
  if(n.le.0)return
  do i=n,1,-1
   x(i+1)=x(i)
  end do
  x(1)%p=0;x(1)%m=0;x(1)%n=0
 end subroutine

 subroutine ttind_assign(b,a)
  implicit none
  type(ttind),intent(inout) :: b
  type(ttind),intent(in) :: a
  b%p=a%p; b%n=a%n; b%m=a%m
 end subroutine

 pure logical function ttind_eq(a,b) result (key)
  type(ttind),intent(in) :: a,b
  key=(a%m.eq.b%m); if(.not.key)return
  key=all(a%p(1:a%m)==b%p(1:b%m))
 end function
 pure logical function ttind_lt(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.false.;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.true. ;return;endif;endif
  key=.false.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).lt.b%p(i))
 end function
 pure logical function ttind_le(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.false.;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.true. ;return;endif;endif
  key=.true.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).lt.b%p(i))
 end function
 pure logical function ttind_gt(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.true. ;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.false.;return;endif;endif
  key=.false.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).gt.b%p(i))
 end function
 pure logical function ttind_ge(a,b) result (key)
  type(ttind),intent(in) :: a,b
  integer :: i
  if(a%m.gt.b%m)then;if(any(a%p(b%m+1:a%m)>1))then;key=.true. ;return;endif;endif
  if(b%m.gt.a%m)then;if(any(b%p(a%m+1:b%m)>1))then;key=.false.;return;endif;endif
  key=.true.; i=min(a%m,b%m)
  do while(i.gt.0 .and. a%p(i).eq.b%p(i)); i=i-1; end do
  if(i.gt.0)key=(a%p(i).gt.b%p(i))
 end function

 pure double precision function dble_ttind(ind) result (x)
  implicit none
  type(ttind),intent(in) :: ind
  integer :: j
  x=dble(ind%p(ind%m))-1
  do j=ind%m-1,1,-1
   x=x*ind%n(j)
   x=x+ind%p(j)-1
  end do
  x=x+1
 end function
 pure integer function int_ttind(ind) result (i)
  implicit none
  type(ttind),intent(in) :: ind
  integer :: j
  i=ind%p(ind%m)-1
  do j=ind%m-1,1,-1
   i=i*ind%n(j)
   i=i+ind%p(j)-1
  end do
  i=i+1
 end function

! MEM
 integer function dtt_mem(arg) result (sz)
  implicit none
  type(dtt),intent(in) :: arg
  integer :: i
  sz=0
  do i=arg%l,arg%m
   sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
 end function
 integer function ztt_mem(arg) result (sz)
  implicit none
  type(ztt),intent(in) :: arg
  integer :: i
  sz=0
  do i=arg%l,arg%m
   sz=sz+arg%r(i-1)*arg%n(i)*arg%r(i)
  end do
 end function
end module
