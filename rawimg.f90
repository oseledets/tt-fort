module rawimg_lib
 implicit none
 integer,parameter,private :: bytes=1,rec1=4,maxint=2**(8*bytes)-1

contains

 subroutine read_raw(n,fnam,img,info)
  implicit none
  integer,intent(in) :: n
  character(len=*),intent(in) :: fnam
  double precision,intent(out) :: img(n)
  integer,intent(out) :: info
  character(len=*),parameter :: subnam='read_raw'
  integer(kind=bytes),allocatable :: a(:)
  logical :: ex,op
  integer :: i
  inquire (file=fnam, exist=ex, opened=op)
  if(.not.ex)then; write(*,*)subnam,': file not exist: ',fnam; info=-1; return; endif
  !open(unit=10,access='direct',form='unformatted', recl=bytes*n/rec1, file=fnam,err=101,iostat=info)
  open(unit=10,access='stream',form='unformatted', status='unknown', file=fnam,err=101,iostat=info)

  allocate(a(n))
  read(unit=10,err=102,iostat=info) a
  close(10)
  info=0
  do i=1,n
   if(a(i).lt.0)then
    img(i)=dble(a(i)+maxint+1)/maxint
   else 
    img(i)=dble(a(i))/maxint
   endif
  end do 
  deallocate(a)
  return
101 continue
  write(*,*)subnam,': file open error: ',info
  return
102 continue
  write(*,*)subnam,': file read error: ',info
  return
 end subroutine 
 
 subroutine write_raw(n,fnam,img,info)
  implicit none
  integer,intent(in) :: n
  character(len=*),intent(in) :: fnam
  double precision,intent(in) :: img(n)
  integer,intent(out) :: info
  character(len=*),parameter :: subnam='write_raw'
  integer(kind=bytes),allocatable :: a(:)
  integer :: inf,i
  
  !open(unit=10,access='direct',form='unformatted', recl=bytes*n/rec1, file=fnam,err=101,iostat=info)
  open(unit=10,access='stream',form='unformatted', status='replace', file=fnam,err=101,iostat=info)

  allocate(a(n),stat=inf)
  if(inf.ne.0)then; write(*,*)subnam,': allocate error: ',inf; info=inf; return; endif
  
  forall(i=1:n)a(i)=int(img(i)*maxint)
  do i=1,n
   if(img(i).gt.0.5d0)then
    a(i)=int(img(i)*maxint)-maxint-1; 
   else 
    a(i)=int(img(i)*maxint)
   end if 
  end do
  !write(*,'(32i4)') a
  write(unit=10,err=102,iostat=info) a
  close(10)
  info=0
  deallocate(a)
  return
101 continue
  write(*,*)subnam,': file open error: ',info
  return
102 continue
  write(*,*)subnam,': file write error: ',info
  return
 end subroutine 
end module
