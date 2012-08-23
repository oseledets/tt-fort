module check_lib
 implicit none
 double precision,parameter :: errval=-999.d0

 interface chk_abmi
  module procedure check_abmi_d,check_abmi_z
 end interface
 interface chk_abmc
  module procedure check_abmc_d,check_abmc_z
 end interface
 interface chk_ort
  module procedure check_ort_d,check_ort_z
 end interface
 interface nnz
  module procedure nnz_d,nnz_d2,nnz_d3
 end interface
 
 private :: check_abmi_d,check_abmi_z
 private :: check_abmc_d,check_abmc_z
 private :: check_ort_d,check_ort_z

contains
 double precision function check_abmi_d(a,b) result (err)
  implicit none
  double precision,intent(in) :: a(:,:),b(:,:)
  double precision,allocatable :: c(:,:)
  integer :: i,m,n
  double precision,external :: dnrm2
  if(size(a,1).ne.size(b,2) .or. size(a,2).ne.size(b,1))then
   err=errval; return
  end if
  m=size(a,1); n=size(a,2)
  allocate(c(m,m))
  c=matmul(a,b)
  !call dgemm('n','n',m,m,n,1.d0,a,m,b,n,0.d0,c,m)
  forall(i=1:m)c(i,i)=c(i,i)-1.d0
  err=dnrm2(m*m,c,1)/dsqrt(dble(m))
  deallocate(c)
  return
 end function check_abmi_d

 double precision function check_abmi_z(a,b) result (err)
  implicit none
  complex(8),intent(in) :: a(:,:),b(:,:)
  complex(8),allocatable :: c(:,:)
  integer :: i,m,n
  double precision,external :: dznrm2
  if(size(a,1).ne.size(b,2) .or. size(a,2).ne.size(b,1))then
   err=errval; return
  end if
  m=size(a,1); n=size(a,2)
  allocate(c(m,m))
  c=matmul(a,b)
  forall(i=1:m)c(i,i)=c(i,i)-(1.d0,0.d0)
  err=dznrm2(m*m,c,1)/dsqrt(dble(m))
  deallocate(c)
  return
 end function check_abmi_z
 
 double precision function check_abmc_d(a,b,c) result (err)
  implicit none
  double precision,intent(in) :: a(:,:),b(:,:),c(:,:)
  double precision,allocatable :: ab(:,:)
  integer :: m,n,k
  double precision :: nrm
  double precision,external :: dnrm2
  if(size(a,2).ne.size(b,1) .or. size(a,1).ne.size(c,1) .or. size(b,2).ne.size(c,2))then
   err=errval; return
  end if
  m=size(a,1); k=size(a,2); n=size(b,2)
  nrm=dnrm2(m*n,c,1)
  if(nrm.le.0.d0)then
   err=errval; return
  end if
  allocate(ab(m,n))
  ab=matmul(a,b)
  ab=ab-c; err=dnrm2(m*n,ab,1)/nrm
  deallocate(ab)
  return
 end function check_abmc_d
 double precision function check_abmc_z(a,b,c) result (err)
  implicit none
  complex(8),intent(in) :: a(:,:),b(:,:),c(:,:)
  complex(8),allocatable :: ab(:,:)
  integer :: m,n,k
  double precision :: nrm
  double precision,external :: dznrm2
  if(size(a,2).ne.size(b,1) .or. size(a,1).ne.size(c,1) .or. size(b,2).ne.size(c,2))then
   err=errval; return
  end if
  m=size(a,1); k=size(a,2); n=size(b,2)
  nrm=dznrm2(m*n,c,1)
  if(nrm.le.0.d0)then
   err=errval; return
  end if
  allocate(ab(m,n))
  ab=matmul(a,b)
  ab=ab-c; err=dznrm2(m*n,ab,1)/nrm
  deallocate(ab)
  return
 end function check_abmc_z
 double precision function check_ort_d(a) result (err)
  implicit none
  double precision,intent(in) :: a(:,:)
  double precision,allocatable :: g(:,:)
  integer :: m,n,i
  double precision,external :: dnrm2
  m=size(a,1); n=size(a,2)
  if(m.lt.n)then
   err=errval; return
  end if
  allocate(g(n,n))
  call dgemm('t','n',n,n,m,1.d0,a,m,a,m,0.d0,g,n)
  forall(i=1:n)g(i,i)=g(i,i)-1.d0
  err=dnrm2(n*n,g,1)/dsqrt(dble(n))
  deallocate(g)
  return
 end function check_ort_d
 double precision function check_ort_z(a) result (err)
  implicit none
  complex(8),intent(in) :: a(:,:)
  complex(8),allocatable :: g(:,:)
  integer :: m,n,i
  double precision,external :: dznrm2
  m=size(a,1); n=size(a,2)
  if(m.lt.n)then
   err=errval; return
  end if
  allocate(g(n,n))
  call zgemm('c','n',n,n,m,(1.d0,0.d0),a,m,a,m,(0.d0,0.d0),g,n)
  forall(i=1:n)g(i,i)=g(i,i)-(1.d0,0.d0)
  err=dznrm2(n*n,g,1)/dsqrt(dble(n))
  deallocate(g)
  return
 end function check_ort_z 
 
 double precision function check_sym_d(n,a) result (err)
  implicit none
  integer,intent(in) :: n
  double precision,intent(in) :: a(n,n)
  integer :: i,j
  double precision,external :: dnrm2
  err=0.d0
  do j=1,n
   do i=1,n
    err=err+(a(i,j)-a(j,i))**2
   end do
  end do
  err=dsqrt(err)
  if(err.gt.0.d0)err=err/dnrm2(n*n,a,1)
 end function 

 integer function nnz_d(x,tol) result (nnz)
  implicit none
  double precision,intent(in) :: x(:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(abs(x));eps=eps*nrm
  nnz=0
  do i=1,size(x)
   if(dabs(x(i)).gt.eps)nnz=nnz+1
  enddo
 end function 
 integer function nnz_d2(x,tol) result (nnz)
  implicit none
  double precision,intent(in) :: x(:,:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i,j
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(abs(x));eps=eps*nrm
  nnz=0
  do j=1,size(x,2)
   do i=1,size(x,1)
    if(dabs(x(i,j)).gt.eps)nnz=nnz+1
   enddo
  enddo
 end function 
 integer function nnz_d3(x,tol) result (nnz)
  implicit none
  double precision,intent(in) :: x(:,:,:)
  double precision,intent(in),optional :: tol
  double precision :: eps,nrm
  integer :: i,j,k
  if(present(tol))then;eps=tol;else;eps=1.d-13;endif
  nrm=maxval(abs(x));eps=eps*nrm
  nnz=0
  do k=1,size(x,3)
   do j=1,size(x,2)
    do i=1,size(x,1)
     if(dabs(x(i,j,k)).gt.eps)nnz=nnz+1
    enddo
   enddo
  enddo 
 end function 


end module
