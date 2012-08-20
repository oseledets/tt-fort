module sort_lib
 implicit none
 interface sort
  ! sort one vector, make the same permutation in the other
  module procedure sort_i,sort_d
 end interface
 private :: sort_i,sort_d
 interface find
  ! for sorted vector x(1) <= x(2) <= ... <= x(n) and value y find pos, s.t. x(pos) <= y < x(pos+1)
  module procedure find_i,find_d
 end interface
 private :: find_i,find_d
 interface push
  module procedure push_i,push_d
 end interface
 private :: push_i,push_d
contains
 subroutine sort_d(x,y,flag)
  implicit none
  double precision,intent(inout) :: x(:)
  integer,intent(inout),optional :: y(size(x))
  integer,intent(in),optional :: flag
  integer :: kflag,n,m,i,j,ij,k,l,il(21),iu(21),ty,tty
  double precision :: tt,t,r
  n=size(x); if(n.le.0)return
  if(present(flag))then
   if(flag.ge.0)then
    kflag=1
   else
    kflag=-1
   end if
  else
   kflag=1
  end if 
  x=x*kflag
  if(.not.present(y))then
    m = 1; i = 1; j = n; r = 0.375e0
 10 if (i .eq. j) goto 50
    if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
    else
         r = r-0.21875e0
    endif
 20 k = i
    ij = i + int((j-i)*r)
    t = x(ij)
    if (x(i) .gt. t) then
     x(ij) = x(i); x(i) = t; t = x(ij)
    endif 
    l = j
    if (x(j) .lt. t) then
     x(ij) = x(j); x(j) = t; t = x(ij)
      if (x(i) .gt. t) then
        x(ij) = x(i); x(i) = t; t = x(ij)
      endif
    endif
 30 l = l-1
    if (x(l) .gt. t) goto  30
 40 k = k+1
    if (x(k) .lt. t) goto  40
    if (k .le. l) then
      tt = x(l); x(l) = x(k); x(k) = tt
      goto 30
    endif
    if (l-i .gt. j-k) then
      il(m) = i; iu(m) = l; i = k; m = m+1
    else
      il(m) = k; iu(m) = j; j = l; m = m+1
    endif
    goto 60
 50 m = m-1
    if (m .eq. 0) goto 200
    i = il(m); j = iu(m)
 60 if (j-i .ge. 1) goto 20
    if (i .eq. 1) goto 10
    i = i-1
 70 i = i+1
    if (i .eq. j) goto 50
    t = x(i+1)
    if (x(i) .le. t) goto 70
    k = i
 80 x(k+1) = x(k)
    k = k-1
    if (t .lt. x(k)) goto 80
    x(k+1) = t
    goto 70
  else
    m = 1; i = 1; j = n; r = 0.375e0
110 if (i .eq. j) goto 150
    if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
    else
         r = r-0.21875e0
    endif
120 k = i
    ij = i + int((j-i)*r)
    t = x(ij); ty = y(ij)
    if (x(i) .gt. t) then
     x(ij) = x(i); x(i) = t; t = x(ij)
     y(ij) = y(i); y(i) = ty; ty = y(ij)
    endif 
    l = j
    if (x(j) .lt. t) then
     x(ij) = x(j); x(j) = t; t = x(ij)
     y(ij) = y(j); y(j) = ty; ty = y(ij)
      if (x(i) .gt. t) then
        x(ij) = x(i); x(i) = t; t = x(ij)
        y(ij) = y(i); y(i) = ty; ty = y(ij)
      endif
    endif
130 l = l-1
    if (x(l) .gt. t) goto 130
140 k = k+1
    if (x(k) .lt. t) goto 140
    if (k .le. l) then
      tt = x(l); x(l) = x(k); x(k) = tt
      tty = y(l); y(l) = y(k); y(k) = tty
      goto 130
    endif
    if (l-i .gt. j-k) then
      il(m) = i; iu(m) = l; i = k; m = m+1
    else
      il(m) = k; iu(m) = j; j = l; m = m+1
    endif
    goto 160
150 m = m-1
    if (m .eq. 0) goto 200
    i = il(m); j = iu(m)
160 if (j-i .ge. 1) goto 120
    if (i .eq. 1) goto 110
    i = i-1
170 i = i+1
    if (i .eq. j) goto 150
    t = x(i+1); ty = y(i+1)
    if (x(i) .le. t) goto 170
    k = i
180 x(k+1) = x(k); y(k+1) = y(k)
    k = k-1
    if (t .lt. x(k)) goto 180
    x(k+1) = t; y(k+1) = ty
    goto 170
  end if
200 continue
  x=x*kflag
  return
 end subroutine 
 subroutine sort_i(x,y,flag)
  implicit none
  integer,intent(inout) :: x(:)
  integer,intent(inout),optional :: y(size(x))
  integer,intent(in),optional :: flag
  integer :: kflag,n,m,i,j,ij,k,l,il(21),iu(21),ty,tty,tt,t
  double precision :: r
  n=size(x); if(n.le.0)return
  if(present(flag))then
   if(flag.ge.0)then
    kflag=1
   else
    kflag=-1
   end if
  else
   kflag=1
  end if 
  x=x*kflag
  if(.not.present(y))then
    m = 1; i = 1; j = n; r = 0.375e0
 10 if (i .eq. j) goto 50
    if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
    else
         r = r-0.21875e0
    endif
 20 k = i
    ij = i + int((j-i)*r)
    t = x(ij)
    if (x(i) .gt. t) then
     x(ij) = x(i); x(i) = t; t = x(ij)
    endif 
    l = j
    if (x(j) .lt. t) then
     x(ij) = x(j); x(j) = t; t = x(ij)
      if (x(i) .gt. t) then
        x(ij) = x(i); x(i) = t; t = x(ij)
      endif
    endif
 30 l = l-1
    if (x(l) .gt. t) goto  30
 40 k = k+1
    if (x(k) .lt. t) goto  40
    if (k .le. l) then
      tt = x(l); x(l) = x(k); x(k) = tt
      goto 30
    endif
    if (l-i .gt. j-k) then
      il(m) = i; iu(m) = l; i = k; m = m+1
    else
      il(m) = k; iu(m) = j; j = l; m = m+1
    endif
    goto 60
 50 m = m-1
    if (m .eq. 0) goto 200
    i = il(m); j = iu(m)
 60 if (j-i .ge. 1) goto 20
    if (i .eq. 1) goto 10
    i = i-1
 70 i = i+1
    if (i .eq. j) goto 50
    t = x(i+1)
    if (x(i) .le. t) goto 70
    k = i
 80 x(k+1) = x(k)
    k = k-1
    if (t .lt. x(k)) goto 80
    x(k+1) = t
    goto 70
  else
    m = 1; i = 1; j = n; r = 0.375e0
110 if (i .eq. j) goto 150
    if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
    else
         r = r-0.21875e0
    endif
120 k = i
    ij = i + int((j-i)*r)
    t = x(ij); ty = y(ij)
    if (x(i) .gt. t) then
     x(ij) = x(i); x(i) = t; t = x(ij)
     y(ij) = y(i); y(i) = ty; ty = y(ij)
    endif 
    l = j
    if (x(j) .lt. t) then
     x(ij) = x(j); x(j) = t; t = x(ij)
     y(ij) = y(j); y(j) = ty; ty = y(ij)
      if (x(i) .gt. t) then
        x(ij) = x(i); x(i) = t; t = x(ij)
        y(ij) = y(i); y(i) = ty; ty = y(ij)
      endif
    endif
130 l = l-1
    if (x(l) .gt. t) goto 130
140 k = k+1
    if (x(k) .lt. t) goto 140
    if (k .le. l) then
      tt = x(l); x(l) = x(k); x(k) = tt
      tty = y(l); y(l) = y(k); y(k) = tty
      goto 130
    endif
    if (l-i .gt. j-k) then
      il(m) = i; iu(m) = l; i = k; m = m+1
    else
      il(m) = k; iu(m) = j; j = l; m = m+1
    endif
    goto 160
150 m = m-1
    if (m .eq. 0) goto 200
    i = il(m); j = iu(m)
160 if (j-i .ge. 1) goto 120
    if (i .eq. 1) goto 110
    i = i-1
170 i = i+1
    if (i .eq. j) goto 150
    t = x(i+1); ty = y(i+1)
    if (x(i) .le. t) goto 170
    k = i
180 x(k+1) = x(k); y(k+1) = y(k)
    k = k-1
    if (t .lt. x(k)) goto 180
    x(k+1) = t; y(k+1) = ty
    goto 170
  end if
200 continue
  x=x*kflag
  return
 end subroutine 

 integer function find_d(n,x,y) result (pos)
  implicit none
  integer,intent(in) :: n
  double precision,intent(in) :: x(n),y
  integer :: s,t,i
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
 integer function find_i(n,x,y) result (pos)
  implicit none
  integer,intent(in) :: n
  integer,intent(in) :: x(n),y
  integer :: s,t,i
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

 subroutine push_d(n,x)
  implicit none
  integer,intent(in) :: n
  double precision,intent(inout) :: x(n+1)
  integer :: i
  if(n.le.0)return
  do i=n,1,-1
   x(i+1)=x(i)
  end do
  x(1)=0.d0
 end subroutine
 subroutine push_ii(n,x)
  implicit none
  integer,intent(in) :: n
  integer(kind=8),intent(inout) :: x(n+1)
  integer :: i
  if(n.le.0)return
  do i=n,1,-1
   x(i+1)=x(i)
  end do
  x(1)=0.d0
 end subroutine
 subroutine push_i(n,x)
  implicit none
  integer,intent(in) :: n
  integer,intent(inout) :: x(n+1)
  integer :: i
  if(n.le.0)return
  do i=n,1,-1
   x(i+1)=x(i)
  end do
  x(1)=0
 end subroutine
 
 
 logical function isperm(p,nmax)
  implicit none
  integer,intent(in) :: p(:)
  integer,intent(in),optional :: nmax
  integer,allocatable :: q(:)
  integer :: i,m,n

  m=size(p)
  if(present(nmax))then;n=nmax; else;n=m; endif

  isperm=(p(1).gt.0)
  if(.not.isperm)return
  do i=1,m
   if(p(i).lt.1 .or. p(i).gt.n)then
    isperm=.false.
    return
   end if
  end do 
  
  allocate(q(m))
  IF(m.eq.n)THEN
   q=0
   do i=1,m
    q(p(i))=i
   end do 
  
   do i=1,n
    if(q(i).lt.1 .or. q(i).gt.n)then
     isperm=.false.
     deallocate(q)
     return
    end if
   end do 
  
   do i=1,n
    if(p(q(i)).ne.i)then
     isperm=.false.
     deallocate(q)
     return
    end if
   end do 
  ELSE
   q=p; call sort_i(q)
   do i=1,m-1
    if(q(i+1).le.q(i))then
     isperm=.false.
     deallocate(q)
     return
    end if
   end do
  END IF 
  isperm=.true.
  deallocate(q)
  return
 end function

end module
