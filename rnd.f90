module rnd_lib
 implicit none
 interface random
  module procedure d1rnd,z1rnd,d2rnd,z2rnd,d3rnd,z3rnd
 end interface
contains
 
 subroutine arnd()
  implicit none
  integer s,clock,crate,cmax
  integer,allocatable :: seed(:)
  call random_seed(size=s)
  allocate(seed(s))
  call random_seed(get=seed)
  write(*,*) 'oldseed: ',seed
  call system_clock(clock,crate,cmax)
  !write(*,*)clock,crate,cmax
  seed(1)=clock 
  call random_seed(put=seed)
  write(*,*) 'newseed: ',seed
 end subroutine 

 double precision function drnd( )
  implicit none
  call random_number(drnd)
 return
 end function

 subroutine d1rnd(d)
  double precision,intent(out)  :: d(:)
  call random_number(d)
 end subroutine
 subroutine z1rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:)
  character(len=*),parameter :: subnam='z1rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 subroutine d2rnd(d)
  double precision,intent(out)  :: d(:,:)
  call random_number(d)
 end subroutine
 subroutine z2rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:)
  character(len=*),parameter :: subnam='z2rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine

 subroutine d3rnd(d)
  double precision,intent(out)  :: d(:,:,:)
  call random_number(d)
 end subroutine
 subroutine z3rnd(z)
  implicit none
  double complex,intent(inout)  :: z(:,:,:)
  character(len=*),parameter :: subnam='z3rnd'
  double precision,allocatable :: d(:)
  integer :: n,info
  n=size(z)
  allocate(d(2*n),stat=info)
  if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
  call random_number(d)
  call dcopy(2*n,d,1,z,1)
  deallocate(d)
 end subroutine
 

 integer function irnd( maxi )
  integer,intent(in) :: maxi
  double precision :: d
  call random_number(d)
  irnd=int(d*maxi)+1
 return
 end function
 subroutine irand(maxi,ix)
  integer,intent(in)  :: maxi
  integer,intent(out) :: ix(:)
  integer :: i,n
  double precision,allocatable :: d(:)
  n=size(ix)
  allocate(d(n))
  call random_number(d)
  ix=int(d*maxi)+1
  deallocate(d)
 return
 end subroutine

end module
