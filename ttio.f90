module ttio_lib
 use tt_lib
 implicit none

 integer,parameter,private :: un=51, rec1=1
 character(len=*),parameter,private :: frm='unformatted', acc='stream'

 integer(kind=4),private  :: ver(2)=(/ 1, 0 /)

 type,private :: tthead
  sequence
  character(len=8) :: txt='TT      '
  integer(kind=4)  :: ver(2)=(/ 1, 0 /)
  integer(kind=4)  :: inf(4)=(/tt_size, 0, 0, 0/)
  character(len=64) :: comment
  integer(kind=4) :: i(8)
 end type

 interface read
  module procedure dtt_read,ztt_read
 end interface
 interface write
  module procedure dtt_write,ztt_write
 end interface

contains

 subroutine dtt_write(arg,fnam,info)
  implicit none
  type(dtt),intent(in) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='dtt_write'
  type(tthead) :: head
  integer :: io,u,i,l,m
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='write',position='rewind',status='replace',err=101,iostat=io)

  head%i(1)=arg%l
  head%i(2)=arg%m
  l=arg%l; m=arg%m
  write(u,err=111,iostat=io) head
  write(u,err=112,iostat=io) arg%l,arg%m
  write(u,err=113,iostat=io) arg%n(l:m),arg%r(l-1:m)
  do i=l,m
   write(u,err=114,iostat=io) arg%u(i)%p
  end do
  close(u,err=121,iostat=io)
  if(present(info))info=0
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error writing header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error writing lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error writing nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine

 subroutine dtt_read(arg,fnam,info)
  implicit none
  type(dtt),intent(inout) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='dtt_read'
  type(tthead) :: head
  integer(4) :: io,u,i,l,m
  integer(4), allocatable :: n4(:), r4(:)
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(.not.ex)then
   write(*,*)subnam,': file not exist: ',fnam
   if(present(info))info=-1
   return
  endif
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='read',position='rewind',status='old',err=101,iostat=io)
  read(u,err=111,iostat=io) head

  if(head%txt(1:2).ne.'TT')then
   write(*,*)subnam,': not TT header in file: ',fnam
   if(present(info))info=-2
   return
  end if
  if(head%ver(1).ne.ver(1))then
   write(*,*)subnam,': not correct version of TT file: ',head%ver
   if(present(info))info=-3
   return
  end if

  read(u,err=112,iostat=io) l,m
  arg%l=l; arg%m=m
  if(l.lt.0)then
   write(*,*)subnam,': read strange l,m: ',l,m
  end if

  allocate(n4(m-l+1), r4(m-l+2))

!   read(u,err=113,iostat=io) arg%n(l:m),arg%r(l-1:m)
  read(u,err=113,iostat=io) n4(1:m-l+1), r4(1:m-l+2)
  arg%n(l:m) = n4(1:m-l+1)
  arg%r(l-1:m) = r4(1:m-l+2)

  call alloc(arg)
  do i=l,m
   read(u,err=114,iostat=io) arg%u(i)%p
  end do
  close(u,err=121,iostat=io)
  if(present(info))info=0
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error reading header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error reading lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error reading nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine




 subroutine ztt_write(arg,fnam,info)
  implicit none
  type(ztt),intent(in) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='ztt_write'
  type(tthead) :: head
  integer :: io,u,i,l,m
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='write',position='rewind',status='replace',err=101,iostat=io)

  head%i(1)=arg%l
  head%i(2)=arg%m
  l=arg%l; m=arg%m
  write(u,err=111,iostat=io) head
  write(u,err=112,iostat=io) arg%l,arg%m
  write(u,err=113,iostat=io) arg%n(l:m),arg%r(l-1:m)
  do i=l,m
   write(u,err=114,iostat=io) arg%u(i)%p
  end do
  close(u,err=121,iostat=io)
  if(present(info))info=0
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error writing header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error writing lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error writing nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine
!
!
!
 subroutine ztt_read(arg,fnam,info)
  implicit none
  type(ztt),intent(inout) :: arg
  character(len=*),intent(in) :: fnam
  integer,intent(out),optional :: info
  character(len=*),parameter :: subnam='ztt_read'
  type(tthead) :: head
  integer :: io,u,i,l,m
  logical :: ex,op

  if(present(info))info=-11
  inquire (file=fnam, exist=ex, opened=op)
  if(.not.ex)then
   write(*,*)subnam,': file not exist: ',fnam
   if(present(info))info=-1
   return
  endif
  if(op)then
   write(*,*)subnam,': file is open, trying to close: ',fnam
   inquire(file=fnam, number=u)
   write(*,*)subnam,': establish unit: ',u
   close(unit=u,status='keep')
   write(*,*)subnam,': closed ok'
  end if

  u=un; op=.true.
  do  while(op)
   inquire(unit=u,opened=op)
   if(op)then
    write(*,*)subnam,': unit ',u,' is busy, trying next '
    u=u+1
   end if
  end do

  open(unit=u,file=fnam,form=frm,access=acc,action='read',position='rewind',status='old',err=101,iostat=io)
  read(u,err=111,iostat=io) head

  if(head%txt(1:2).ne.'TT')then
   write(*,*)subnam,': not TT header in file: ',fnam
   if(present(info))info=-2
   return
  end if
  if(head%ver(1).ne.ver(1))then
   write(*,*)subnam,': not correct version of TT file: ',head%ver
   if(present(info))info=-3
   return
  end if

  read(u,err=112,iostat=io) l,m
  arg%l=l; arg%m=m
  if(l.lt.0)then
   write(*,*)subnam,': read strange l,m: ',l,m
  end if

  read(u,err=113,iostat=io) arg%n(l:m),arg%r(l-1:m)
  call alloc(arg)
  do i=l,m
   read(u,err=114,iostat=io) arg%u(i)%p
  end do
  close(u,err=121,iostat=io)
  if(present(info))info=0
  return

101 continue
  write(*,*) subnam,': error opening file: ',io
  if(present(info))info=io
  return
111 continue
  write(*,*) subnam,': error reading header: ',io
  if(present(info))info=io
  return
112 continue
  write(*,*) subnam,': error reading lm: ',io
  if(present(info))info=io
  return
113 continue
  write(*,*) subnam,': error reading nr: ',io
  if(present(info))info=io
  return
114 continue
  write(*,*) subnam,': error writing cores: ',io
  if(present(info))info=io
  return
121 continue
  write(*,*) subnam,': error closing file: ',io
  if(present(info))info=io
  return
 end subroutine

end module
