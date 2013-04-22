module mat_lib
 interface matinv
  module procedure matinv_d,matinv_z
 end interface
 private :: matinv_d,matinv_z,matinv_svd_d,matinv_svd_z

 interface eye
  module procedure eye_d1, eye_z1, eye_d2,eye_z2
 end interface
 interface laplace
  module procedure laplace_d2,laplace_z2
 end interface

contains
 subroutine matinv_d(a,ainv,alg,tol)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out),optional :: ainv(size(a,2),size(a,1))
  character(len=1),intent(in),optional :: alg
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_d'
  if(.not.present(alg) .or. size(a,1).ne.size(a,2))then
   call matinv_svd_d(a,ainv,tol)
  else
   select case(alg)
    case('s','S')
     call matinv_svd_d(a,ainv,tol)
    case('t','T')
     call matinv_lu_d(a,ainv)
    case default
     write(*,'(3a)')subnam,': unknown alg: ',alg
     stop
   end select
  end if
  return
 end subroutine 
 subroutine matinv_z(a,ainv,alg,tol)
  implicit none
  complex(8),intent(inout)  :: a(:,:)
  complex(8),intent(out),optional :: ainv(size(a,2),size(a,1))
  character(len=1),intent(in),optional :: alg
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_z'
  if(.not.present(alg) .or. size(a,1).ne.size(a,2))then
   call matinv_svd_z(a,ainv,tol)
  else
   select case(alg)
    case('s','S')
     call matinv_svd_z(a,ainv,tol)
    case('t','T')
     call matinv_lu_z(a,ainv)
    case default
     write(*,'(3a)')subnam,': unknown alg: ',alg
     stop
   end select
  end if
  return
 end subroutine 

 subroutine matinv_svd_d(a,ainv,tol)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out),optional :: ainv(size(a,2),size(a,1))
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_svd_d'
  double precision,allocatable :: b(:,:),u(:,:),v(:,:),s(:),work(:)
  integer :: r,lwork,info,i,ntrunc,m,n
  double precision :: si,s1,small
  character(len=180) :: str
  
  m=size(a,1); n=size(a,2)
  if(present(tol))then;small=tol;else;small=1.d-14;endif
  write(str,'(a,a,i8,1x,i8)')subnam,': m,n:',m,n 
  !call plog(2,str)
  r=min(m,n); lwork=64*max(m,n)
  allocate(b(m,n),u(m,r),v(r,n),s(r),work(lwork))
  b=a
  call dgesvd('s','s',m,n,b,m,s,u,m,v,r, work,lwork,info)
  if(info.ne.0)then
   write(*,'(a,a,i10)')subnam,': dgesvd info: ',info
   stop
  end if
  
  s1=s(1)
  !write(str,'(a,80(e9.3,1x))')'sv:',(fmem(s+i)/s1,i=1,r)
  !call plog(1,str)

  ntrunc=0
  do i=1,r
   si=s(i)
    if(si.lt.small*s1)then
     s(i)=0.d0; ntrunc=ntrunc+1
    else
     s(i)=1.d0/s(i)
    end if
  end do
  
  forall(i=1:r)u(:,i)=u(:,i)*s(i)
  if(present(ainv))then
   call dgemm('t','t',n,m,r,1.d0,v,r,u,m,0.d0,ainv,n)
  else
   call dgemm('t','t',n,m,r,1.d0,v,r,u,m,0.d0,a,n)
  end if 
  if(ntrunc.gt.0)then
   write(*,'(a,a,i5,a,i5)')subnam,': truncate: ',ntrunc,' of: ',r
   !call plog(1,str)
  end if

  deallocate(b,u,v,s,work)
  return
 end subroutine
 subroutine matinv_svd_z(a,ainv,tol)
  implicit none
  complex(8),intent(in)  :: a(:,:)
  complex(8),intent(out) :: ainv(size(a,2),size(a,1))
  double precision,intent(in),optional :: tol
  character(len=*),parameter   :: subnam='matinv_svd_z'
  complex(8),allocatable   :: b(:,:),u(:,:),v(:,:),work(:)
  double precision,allocatable :: s(:),rwork(:)
  integer :: r,lwork,lrwork,info,i,ntrunc,m,n
  double precision :: small=1.d-14
  double precision :: si,s1
  character(len=180) :: str
  
  m=size(a,1); n=size(a,2)
  if(present(tol))small=tol
  write(str,'(a,a,i8,1x,i8)')subnam,': m,n:',m,n 
  !call plog(2,str)
  r=min(m,n); lwork=64*max(m,n); lrwork=5*min(m,n)
  allocate(b(m,n),u(m,r),v(n,r),s(r),work(lwork),rwork(lrwork))
  b=a
  call zgesvd('s','s',m,n,b,m,s,u,m,v,n,work,lwork,rwork,info)
  if(info.ne.0)then
   write(*,'(a,a,i10)')subnam,': zgesvd info: ',info
   stop
  end if
  
  s1=s(1)
  !write(str,'(a,80(e9.3,1x))')'sv:',(fmem(s+i)/s1,i=1,r)
  !call plog(1,str)

  ntrunc=0
  do i=1,r
   si=s(i)
    if(si.lt.small*s1)then
     s(i)=0.d0; ntrunc=ntrunc+1
    else
     s(i)=1.d0/s(i)
    end if
  end do
  
  forall(i=1:r)u(:,i)=u(:,i)*s(i)
  call zgemm('c','c',n,m,r,(1.d0,0.d0),v,r,u,m,(0.d0,0.d0),ainv,n)
  if(ntrunc.gt.0)then
   write(str,'(a,a,i5,a,i5)')subnam,': truncate: ',ntrunc,' of: ',r
   !call plog(1,str)
  end if

  deallocate(b,u,v,s,work,rwork)
  return
 end subroutine

 subroutine matinv_lu_d(a,ainv)
  implicit none
  double precision,intent(inout)  :: a(:,:)
  double precision,intent(out),optional :: ainv(size(a,2),size(a,1))
  character(len=*),parameter   :: subnam='matinv_lu_d'
  double precision,allocatable :: work(:)
  double precision,dimension(:,:),pointer :: aa
  target  :: a,ainv
  integer :: lwork,m,n,info
  integer,allocatable :: piv(:)
  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   write(*,*)subnam,': matrix not square: ',m,n
   stop
  end if 
  if(present(ainv))then
   aa=>ainv; ainv=a
  else
   aa=>a
  end if

  lwork=64*n
  allocate(work(lwork),piv(n))
  call dgetrf(n,n,aa,n,piv,info)
  if(info.ne.0)then; write(*,*)subnam,': dgertf: info: ',info; stop; end if
  call dgetri(n,aa,n,piv,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgerti: info: ',info; stop; end if

  nullify(aa)
  deallocate(work,piv)
  return
 end subroutine
 subroutine matinv_lu_z(a,ainv)
  implicit none
  complex(8),intent(inout)  :: a(:,:)
  complex(8),intent(out),optional :: ainv(size(a,2),size(a,1))
  character(len=*),parameter   :: subnam='matinv_lu_z'
  complex(8),allocatable :: work(:)
  complex(8),dimension(:,:),pointer :: aa
  target  :: a,ainv
  integer :: lwork,m,n,info
  integer,allocatable :: piv(:)
  m=size(a,1); n=size(a,2)
  if(m.ne.n)then
   write(*,*)subnam,': matrix not square: ',m,n
   stop
  end if 
  if(present(ainv))then
   aa=>ainv; ainv=a
  else
   aa=>a
  end if

  lwork=64*n
  allocate(work(lwork),piv(n))
  call zgetrf(n,n,aa,n,piv,info)
  if(info.ne.0)then; write(*,*)subnam,': dgertf: info: ',info; stop; end if
  call zgetri(n,aa,n,piv,work,lwork,info)
  if(info.ne.0)then; write(*,*)subnam,': dgerti: info: ',info; stop; end if

  deallocate(work,piv)
  return
 end subroutine

 subroutine eye_d2(a)
  implicit none
  double precision,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  a=0.d0
  forall(i=1:min(m,n))a(i,i)=1.d0
  return
 end subroutine 
 
 subroutine eye_d1(a,n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(out) :: a(:)
    integer :: i
    a(1:n*n) = 0d0
    forall(i = 1:n) a(i+(i-1)*n) = 1d0
 end subroutine eye_d1
 
 subroutine eye_z1(a,n)
    implicit none
    integer, intent(in) :: n
    complex(8), intent(out) :: a(:)
    integer :: i
    a(1:n*n) = (0d0,0d0)
    forall(i = 1:n) a(i+(i-1)*n) = (1d0,0d0)
 end subroutine eye_z1



 subroutine eye_z2(a)
  implicit none
  complex(8),intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  a=(0.d0,0.d0)
  forall(i=1:min(m,n))a(i,i)=(1.d0,0.d0)
  return
 end subroutine 
 


 subroutine d2eye(n,a)
  implicit none
  integer,intent(in) :: n
  double precision,intent(out) :: a(n,n)
  integer :: i
  call dscal(n*n,0.d0,a,1)
  forall(i=1:n)a(i,i)=1.d0
 end subroutine 
 
 subroutine z2eye(n,a)
  implicit none
  integer,intent(in) :: n
  complex(8),intent(out) :: a(n,n)
  integer :: i
  call zdscal(n*n,0.d0,a,1)
  forall(i=1:n)a(i,i)=(1.d0,0.d0)
 end subroutine 

 subroutine laplace_d2(a)
  implicit none
  double precision,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  call dscal(m*n,0.d0,a,1)
  forall(i=1:min(m,n))a(i,i)=2.d0
  forall(i=1:min(m-1,n))a(i+1,i)=-1.d0
  forall(i=1:min(m,n-1))a(i,i+1)=-1.d0
  return
 end subroutine 
 subroutine laplace_z2(a)
  implicit none
  double complex,intent(out) :: a(:,:)
  integer :: m,n,i
  m=size(a,1); n=size(a,2)
  call zdscal(m*n,0.d0,a,1)
  forall(i=1:min(m,n))a(i,i)=(2.d0,0.d0)
  forall(i=1:min(m-1,n))a(i+1,i)=-(1.d0,0.d0)
  forall(i=1:min(m,n-1))a(i,i+1)=-(1.d0,0.d0)
  return
 end subroutine 
 
 subroutine d2submat(m,n,a,lda,b)
  implicit none
  integer,intent(in) :: m,n,lda
  double precision,intent(in) :: a(lda,n)
  double precision,intent(out) :: b(m,n)
  integer :: i,j
  if(m.gt.lda)then;write(*,*)'d2submat: lda,m: ',lda,m;stop;endif
  if(m.eq.lda)then;call dcopy(m*n,a,1,b,1);return;endif
  forall(i=1:m,j=1:n)b(i,j)=a(i,j)
 end subroutine
 subroutine z2submat(m,n,a,lda,b)
  implicit none
  integer,intent(in) :: m,n,lda
  complex(8),intent(in) :: a(lda,n)
  complex(8),intent(out) :: b(m,n)
  integer :: i,j
  if(m.gt.lda)then;write(*,*)'d2submat: lda,m: ',lda,m;stop;endif
  if(m.eq.lda)then;call zcopy(m*n,a,1,b,1);return;endif
  forall(i=1:m,j=1:n)b(i,j)=a(i,j)
 end subroutine
end module      
