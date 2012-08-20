module tt_adapt_als

implicit none
real(8), allocatable :: result_core(:)
character*120 :: matlab_ist_dumme_kuh

contains
subroutine deallocate_result
  if ( allocated(result_core) ) then
    deallocate(result_core)
  end if
end subroutine deallocate_result


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General full matrix and tensor routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eye(n,a)
 integer, intent(in) :: n
 double precision, intent(inout) :: a(n,n)
 integer i
 a(:,:) = 0
 do i = 1,n
   a(i,i) = 1d0
 end do
end subroutine eye

subroutine qr(n,m, A, R, work, lwork, tau)
integer, intent(in):: n,m,lwork
real(8), intent(inout) :: A(n,m), R(min(n,m),n)
real(8), intent(inout) :: work(*), tau(*)
integer info
integer rnew, k,j

call dgeqrf(n, m, A, n, tau, work,lwork,info)
if (info.ne.0) then
  print *, 'qr: dgeqrf failed'
end if
rnew = min(n,m)
R(:,:)=0d0
do j=1,m
  R(1:min(j,n),j)=A(1:min(j,n),j)
!   forall(k=min(j,n)+1:rnew) R(k,j)=0.d0
end do
call dorgqr(n,rnew,rnew,A,n,tau,work,lwork,info)
if (info.ne.0) then
  print *, 'qr: dorgqr failed'
end if

end subroutine


subroutine row_add(m,n,k,A,B)
integer, intent(in) :: m, n, k
real(8), intent(inout) :: A(*)
real(8), intent(in) :: B(*)
real(8) swp(m)
integer i

do i=n,1,-1
  call dcopy(m, A(1+(i-1)*m), 1, swp, 1)
  call dcopy(m, swp, 1, A(1+(i-1)*(m+k)), 1)
  call dcopy(k, B(1+(i-1)*k), 1, A(m+1+(i-1)*(m+k)), 1)
end do
end subroutine

subroutine row_cut(m,n,k,A)
integer, intent(in) :: m, n, k
real(8), intent(inout) :: A(*)
integer i

do i=2,n
  call dcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
end do
end subroutine


subroutine transp(n, m, A, B)
  integer, intent(in):: n,m
  real(8), intent(in):: A(n,m)
  real(8), intent(inout):: B(m,n)

  integer i,j

  do i=1,n
    call dcopy(m, A(i,1), n, B(1,i),1)
  end do
end subroutine

subroutine perm1324(n1,n2,n3,n4, A, B)
  integer, intent(in) :: n1,n2,n3,n4
  real(8), intent (in) :: A(n1,n2,n3,n4)
  real(8), intent (inout) :: B(n1,n3,n2,n4)

  integer i2,i3,i4

  do i4=1,n4
    do i3=1,n3
      do i2=1,n2
        call dcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
      end do
    end do
  end do

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALS-related procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bfun3(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y, res1, res2)
! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
  integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
  real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
  real(8), intent(inout) :: y(*), res1(*), res2(*)

  call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
!    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
  call transp(rx1, m*ra2*ry2, res1, res2)
  call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
!     j1, a2, c2, b1
  call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
!     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
  call transp(ra1*n*ry2, rx1, res2, res1)
!     b1,a1,i1,c2
!    phi1: c1, b1, a1 : ry1, rx1, ra1
  call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res1, rx1*ra1, 0d0, y, ry1)
!     y: c1,i1,c2

end subroutine


subroutine bfun3_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, A, phi2, x, y, res1, res2)
! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
  integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
  real(8), intent(in) :: phi2(*), A(*), x(*)
  real(8), intent(inout) :: y(*), res1(*), res2(*)


!   phi2: b2,a2,c2: rx2, ra2, ry2
!  x: b1,j1,b2: rx1,m,rx2

  call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
!    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
  call transp(rx1, m*ra2*ry2, res1, res2)
  call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
!     j1, a2, c2, b1
  call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
!     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
  call transp(ra1*n*ry2, rx1, res2, y)
!   output is rx1*ra1,n,ry2
end subroutine


subroutine phi_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi2_old, A, x, y, phi2, res1, res2)
! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2, rx1*ra1*ry1)
  integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
  real(8), intent(in) ::  A(*), phi2_old(*), x(*), y(*)
  real(8), intent(inout) :: phi2(*), res1(*), res2(*)


  !   phi2: b2,a2,c2: rx2, ra2, ry2
  !  x: b1,j1,b2: rx1,m,rx2
  call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2_old, rx2, 0d0, res1, rx1*m)
  !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
  call transp(rx1, m*ra2*ry2, res1, res2)
  call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
  !     j1, a2, c2, b1
  call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
  !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
  call transp(ra1, n*ry2*rx1, res2, res1)
  !   i1,c2,b1,a1
  call dgemm('N', 'N', ry1, rx1*ra1, n*ry2, 1d0, y, ry1, res1, n*ry2, 0d0, phi2, ry1)
  ! phi2: ry1, rx1, ra1
  call transp(ry1, rx1*ra1, phi2, res1)
  call dcopy(rx1*ra1*ry1, res1, 1, phi2, 1)
end subroutine


! y'Ax
subroutine phi_left(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1_old, A, x, y, phi1, res1, res2)
! sizes of res1, res2: max(rx1*n*ra1*ry2, rx1*ra2*m*ry2, ry2*rx2*ra2)
  integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
  real(8), intent(in) ::  A(*), phi1_old(*), x(*), y(*)
  real(8), intent(inout) :: phi1(*), res1(*), res2(*)

  !   phi1: c1, b1, a1 : ry1, rx1, ra1
  !  y: c1,i1,c2: ry1,n,ry2
  call dgemm('T', 'N', rx1*ra1, n*ry2, ry1, 1d0, phi1_old, ry1, y, ry1, 0d0, res1, rx1*ra1)
  !    res1: rx1,ra1,n,ry2: b1,a1,i1,c2
  call transp(rx1, ra1*n*ry2, res1, res2)
  call dcopy(ra1*n*ry2*rx1, res2, 1, res1, 1)
  !     a1, i1, c2, b1
  call dgemm('T', 'N', ry2*rx1, m*ra2, ra1*n, 1d0, res1, ra1*n, A, ra1*n, 0d0, res2, ry2*rx1)
  !     res2: ry2,rx1,m,ra2 : c2, b1, j1, a2
  call transp(ry2, rx1*m*ra2, res2, res1)
  !   b1,j1,a2,y2
  call dgemm('T', 'N', ra2*ry2, rx2, rx1*m, 1d0, res1, rx1*m, x, rx1*m, 0d0, phi1, ra2*ry2)
  ! phi1: ra2, ry2, rx2
  call transp(ra2, ry2*rx2, phi1, res1)
  call dcopy(ry2*rx2*ra2, res1, 1, phi1, 1)
end subroutine


subroutine phi2_right(rx1, rx2, ry1, n, ry2, phi2_old, x, y, phi2, res1)
! sizes of res1: rx1*n*ry2
integer, intent(in) :: rx1, rx2, ry1, n, ry2
real(8), intent(in) ::  phi2_old(*), x(*), y(*)
real(8), intent(inout) :: phi2(*), res1(*)


!   phi2: b2,c2: rx2, ry2
!  x: b1,j1,b2: rx1,m,rx2
call dgemm('N', 'N', rx1*n, ry2, rx2, 1d0, x, rx1*n, phi2_old, rx2, 0d0, res1, rx1*n)
!    res1: rx1,m,ry2: b1,j1,c2
call dgemm('N', 'T', rx1, ry1, n*ry2, 1d0, res1, rx1, y, ry1, 0d0, phi2, rx1)
end subroutine


! y'x
subroutine phi2_left(rx1, rx2, ry1, n, ry2, phi1_old, x, y, phi1, res1)
! sizes of res1: rx1*n*ry2
integer, intent(in) :: rx1, rx2, ry1, n, ry2
real(8), intent(in) ::  phi1_old(*), x(*), y(*)
real(8), intent(inout) :: phi1(*), res1(*)

!   phi1: c1, b1, a1 : ry1, rx1
!  y: c1,i1,c2: ry1,n,ry2
call dgemm('T', 'N', rx1, n*ry2, ry1, 1d0, phi1_old, ry1, y, ry1, 0d0, res1, rx1)
!    res1: rx1,n,ry2: b1,j1,c2
call dgemm('T', 'N', ry2, rx2, rx1*n, 1d0, res1, rx1*n, x, rx1*n, 0d0, phi1, ry2)
end subroutine


subroutine Bfull(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, B, res1, res2)
! sizes of res1, res2: max(ry1*n*rx1*m*ra2, rx2*ra2*ry2, ry1*n*ry2*rx1*m*rx2)
integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
real(8), intent(in) :: phi1(*), A(*), phi2(*)
real(8), intent(inout) :: B(*), res1(*), res2(*)

! phi1: ry1,rx1,ra1
  call dgemm('N', 'N', ry1*rx1, n*m*ra2, ra1, 1d0, phi1, ry1*rx1, A, ra1, 0d0, res1, ry1*rx1)
! res1: ry1,rx1,n,m,ra2
  call perm1324(ry1, rx1, n, m*ra2, res1, res2)
  call dcopy(ry1*n*rx1*m*ra2, res2, 1, res1, 1)
! phi2: rx2,ra2,ry2
  call transp(rx2, ra2*ry2, phi2, res2)
  call dgemm('N', 'N', ry1*n*rx1*m, ry2*rx2, ra2, 1d0, res1, ry1*n*rx1*m, res2, ra2, 0d0, B, ry1*n*rx1*m);
  call perm1324(ry1*n, rx1*m, ry2, rx2, B, res1)
! now B: ry1,n,ry2,rx1,m,rx2
  call dcopy(ry1*n*ry2*rx1*m*rx2, res1, 1, B, 1)
end subroutine


subroutine djac_gen(ptype, rx1, n, rx2, ra1, ra2, Phi1, A, Phi2, jacs)
! sizes of work1, work2:
!   center:
!       rx*ra, rx1*n*n*rx2, rx2*ra2*rx2
!   left:
!       rx*ra, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2, rx2*ra2*rx2
!   right:
!       rx*ra, rx1*n*n*rx2*rx2, rx2*ra2*rx2
! sizes of jacs:
!   center:
!       rx1*n*n*ra2,  rx1*n*n*rx2
!   left:
!       rx1*rx1*n*n*ra2,  rx1*rx1*n*n*rx2
!   right:
!       n*n*rx2*rx2*rx1, rx1*n*n*ra2
character, intent(in):: ptype
integer, intent(in) :: rx1, n, rx2, ra1, ra2
real(8), intent(in) :: Phi1(*), A(*), Phi2(*)
real(8), intent(inout) :: jacs(*) !, work1(*), work2(*)
! integer, intent (inout) :: ipiv(*)
real(8), allocatable :: work1(:), work2(:)
integer, allocatable :: ipiv(:)
integer :: i, info


  if ((ptype=='c').or.(ptype=='C')) then
!     print *, rx1, n, rx2, ra1, ra2, ptype

    i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2)
    allocate(work1(i))
    allocate(work2(i))
    allocate(ipiv(n))

    call transp(rx2*ra2, rx2, Phi2, work1)
    do i=0,ra2-1
      call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
    end do
    do i=0,ra1-1
      call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
    end do
!     print *, 'diags prepared'

    call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
!     print *, 'first dgemm'
    call dgemm('N', 'T', rx1*n*n, rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2, 0d0, work1, rx1*n*n)
!     print *, '2nd dgemm'
    call transp(rx1, n*n*rx2, work1, work2)
    ! inversion
    work1(1:n*n)=0d0
    do i=1,n
      work1(i+(i-1)*n)=1d0
    end do
!     print *, rx2*rx1
    do i=0,rx2*rx1-1
      call dcopy(n*n, work1, 1, jacs(i*n*n+1), 1)
      call dgesv(n, n, work2(i*n*n+1), n, ipiv, jacs(i*n*n+1), n, info)
!       print *, i, '-th dgesv'
    end do
  end if

  if ((ptype=='l').or.(ptype=='L')) then
    i = max(rx1*ra1, rx2*ra2*rx2, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2)
    allocate(work1(i))
    allocate(work2(i))
    i = rx1*n
    allocate(ipiv(i))

    call transp(rx2*ra2, rx2, Phi2, work1)
    do i=0,ra2-1
      call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
    end do

    call dgemm('N', 'N', rx1*rx1, n*n*ra2, ra1, 1d0, Phi1, rx1*rx1, A, ra1, 0d0, work1, rx1*rx1)
    call perm1324(rx1, rx1, n, n*ra2, work1, jacs)
    call dgemm('N', 'T', rx1*n*rx1*n, rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2, 0d0, work1, rx1*n*rx1*n)
!     call dcopy(rx1*n*rx1*n*rx2, work1, 1, work2, 1)
    ! inversion
    work2(1:rx1*n*rx1*n)=0d0
    do i=1,rx1*n
      work2(i+(i-1)*rx1*n)=1d0
    end do
    do i=0,rx2-1
      call dcopy(rx1*n*rx1*n, work2, 1, jacs(i*rx1*n*rx1*n+1), 1)
      call dgesv(rx1*n, rx1*n, work1(i*rx1*n*rx1*n+1), rx1*n, ipiv, jacs(i*rx1*n*rx1*n+1), rx1*n, info)
    end do
  end if

  if ((ptype=='r').or.(ptype=='R')) then
    i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2*rx2)
    allocate(work1(i))
    allocate(work2(i))
    i = rx2*n
    allocate(ipiv(i))

    call transp(rx2*ra2, rx2, Phi2, work2)
    do i=0,ra1-1
      call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
    end do

    call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
    call dgemm('N', 'T', rx1*n*n, rx2*rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2*rx2, 0d0, work1, rx1*n*n)
    call perm1324(rx1*n, n, rx2, rx2, work1, work2)
    call transp(rx1, n*rx2*n*rx2, work2, work1)
    ! inversion
    work2(1:rx2*n*rx2*n)=0d0
    do i=1,rx2*n
      work2(i+(i-1)*rx2*n)=1d0
    end do
    do i=0,rx1-1
      call dcopy(rx2*n*rx2*n, work2, 1, jacs(i*rx2*n*rx2*n+1), 1)
      call dgesv(rx2*n, rx2*n, work1(i*rx2*n*rx2*n+1), rx2*n, ipiv, jacs(i*rx2*n*rx2*n+1), rx2*n, info)
    end do
  end if

  deallocate(work1,work2, ipiv)
end subroutine

subroutine djac_apply(ptype, rx1, n, rx2, jacs, x, y, work1)
! sizes of work1: rx1*n*rx2
character, intent(in) :: ptype
integer, intent(in) :: rx1, n, rx2
real(8), intent(in) :: jacs(*), x(*)
real(8), intent(inout) :: y(*), work1(*)
integer i

  if ((ptype=='c').or.(ptype=='C')) then
  ! jacs is n,n,rx2,rx1
    call transp(rx1, n*rx2, x, work1)
    do i=0,(rx2*rx1-1)
      call dgemv('N', n, n, 1d0, jacs(i*n*n+1), n, work1(i*n+1), 1, 0d0, y(i*n+1), 1)
    end do
    call transp(n*rx2, rx1, y, work1)
    call dcopy(rx1*n*rx2, work1, 1, y, 1)
  end if

  if ((ptype=='l').or.(ptype=='L')) then
  ! jacs is rx1*n,rx1*n, rx2
    call dcopy(rx1*n*rx2, x, 1, work1, 1)
    do i=0,rx2-1
      call dgemv('N', rx1*n, rx1*n, 1d0, jacs(i*rx1*n*rx1*n+1), rx1*n, work1(i*rx1*n+1), 1, 0d0, y(i*rx1*n+1), 1)
    end do
  end if

  if ((ptype=='r').or.(ptype=='R')) then
  ! jacs is n*rx2, n*rx2, rx1
    call transp(rx1, n*rx2, x, work1)
    do i=0,rx1-1
        call dgemv('N', n*rx2, n*rx2, 1d0, jacs(i*n*rx2*n*rx2+1), n*rx2, work1(i*n*rx2+1), 1, 0d0, y(i*n*rx2+1), 1)
    end do
    call transp(n*rx2, rx1, y, work1)
    call dcopy(rx1*n*rx2, work1, 1, y, 1)
  end if
!   print *, 'djac_apply done'
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TT and svd stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine compute_ps(d,r,n,ps)
integer, intent(in) :: d,r(*),n(*)
integer, intent(out) :: ps(*)
integer i

ps(1)=1;
do i=1,d
  ps(i+1) = ps(i) + r(i)*n(i)*r(i+1)
end do
end subroutine

! RELATIVE accuracy
real(8) function my_chop3(n, s, eps)
real(8), intent(in) :: s(*), eps
integer, intent(in) :: n
real(8) cursum, nrm
integer i
real(8) dnrm2

  nrm = dnrm2(n,s,1)
  nrm = (nrm*eps)*(nrm*eps);

  cursum = 0d0
  i = n;
  do while (i>0)
    cursum = cursum+(s(i)*s(i))
    if (cursum>nrm) then
      exit
    end if
    i = i-1
  end do

  my_chop3 = min(i+1, n)
end function


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AMR MatVec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tt_mvk4(d,n,m,rx,ra,crA, crX, crY0, ry, eps, rmax, kickrank, nswp, verb)
integer,intent(in) :: d,n(*),m(*),rx(*),ra(*), rmax
integer,intent(inout) :: ry(*)
integer, intent(in), optional :: kickrank, nswp, verb
! verb: 0: off; 1: matlab; 2: console
integer :: kickrank0, nswp0, verb0
real(8), intent(in) :: crA(*), crX(*), eps, crY0(*)
real(8) eps2
! integer cr_size
! real(8),allocatable ::  curcr(:), work(:), tau(:), R(:)
! real(8),allocatable :: crl(:,d), crr(:,d), phil(:,d), phir(:,d)
real(8),allocatable :: cr(:,:)
! real(8) :: crr(1000)
real(8),allocatable :: curcr(:)
real(8),allocatable :: work(:)
real(8),allocatable :: tau(:)
real(8),allocatable :: R(:)
! real(8) :: phil(1000)
real(8),allocatable :: phi(:,:)
real(8),allocatable :: full_core(:)

integer,allocatable :: pa(:), px(:)
integer :: i,j,k, swp, dir, lwork, mm, nn, rnew
integer info
real(8) :: err, err_max

real(8) dnrm2
real(8) sqrt

INTEGER :: i_, n_, clock_
INTEGER, DIMENSION(:), ALLOCATABLE :: seed_

CALL RANDOM_SEED(size = n_)
ALLOCATE(seed_(n_))
CALL SYSTEM_CLOCK(COUNT=clock_)
seed_ = clock_ + 37 * (/ (i_ - 1, i_ = 1, n_) /)
CALL RANDOM_SEED(PUT = seed_)
DEALLOCATE(seed_)


eps2 = eps/sqrt(real(d,8))
! print *, eps2

! real(4) etime
! real(4), dimension(2) :: tarray
! real(4) :: t0,t1
!
!   print *, rmax
!   print *, d
!   print *, rmax*maxval(n(1:d))*rmax
!   print *, rx
!   print *, ra
!   print *, ry

  kickrank0 = 5;
  if (present(kickrank)) then
    kickrank0 = kickrank
  end if
!   if (.not.present(lwork_core)) then
!     lwork_phi = sum(ry(2:d+1)*ry(1:d)*n(1:d))*2
!   end if
!   if (.not.present(lwork_phi)) then
!     lwork_phi = sum(rx*ry*ra)*2
!   end if
  nswp0 = 20
  if (present(nswp)) then
    nswp0 = nswp
  end if
  verb0 = 1
  if (present(verb)) then
    verb0 = verb
  end if

  allocate(pa(d+1), px(d+1))
!   print *, nswp0

  call compute_ps(d,ra,n(1:d)*m(1:d),pa)
  call compute_ps(d,rx,m,px)
!   px(1)=1;
!   do i=1,d
!     px(i+1) = px(i) + rx(i)*m(i)*rx(i+1)
!   end do


!   print *, pa
!   print *, px

!   working array: position is at i+1
!   previous array: position is at i
!   phi: at i
!   pyr = sum(ry(2:d+1)*ry(1:d)*n(1:d))+1
! !   print *, pyr
!   pyl = pyr-ry(d)*n(d)*ry(d+1);
!   ppr = sum(rx*ry*ra) + 1
!   ppl = ppr;
  lwork = rmax*maxval(n(1:d))*rmax

  allocate(cr(lwork, d))
  nn = maxval(rx(1:d+1))*maxval(ra(1:d+1))*rmax
  allocate(phi(nn, d+1))
!   print *, lwork
!   allocate(crr(pyr-1))
!   allocate(phir(ppr-1))
!   allocate(crl(pyr-1))
!   allocate(phil(ppl-1))
  allocate(curcr(lwork))
!   allocate(work(lwork))
  nn = maxval(ra(1:d+1))*maxval(rx(1:d+1))*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
  allocate(work(nn))
  allocate(R(lwork))
  allocate(tau(lwork))
  allocate(full_core(nn))

  cr(:,:)=0d0
  mm = 1
  do i=1,d
    call dcopy(ry(i)*n(i)*ry(i+1), crY0(mm), 1, cr(1,i), 1)
    mm = mm + ry(i)*n(i)*ry(i+1)
  end do
! They should be initialized to the identity matrix
 call eye(ry(1),phi(1,1))
 call eye(ry(d+1),phi(1,d+1))
!  phi(1,1)=1d0
!  phi(1,d+1)=1d0

!   t0 =  etime(tarray)
!   QR, psi
  dir = 1
  i = 1

  do while (i<d)
!     current core ppr-r1*n*r2:ppr-1. Transpose it and QR
!     mm=n(i)*ry(i+1); nn=ry(i)
    rnew = min(ry(i)*n(i), ry(i+1))
!     call dcopy(nn*mm, crr(pyr-ry(i)*n(i)*ry(i+1)), 1, curcr, 1)
!     print *, curcr(1:nn*mm)
!     call transp(nn, mm, cr(1,i), curcr)
    call qr(ry(i)*n(i), ry(i+1), cr(1,i), R, work, lwork, tau)
!     print *, curcr(1:nn*mm)

!     call transp(mm, rnew, curcr, cr(1,i))

    call dgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, R, rnew, cr(1,i+1), ry(i+1), 0d0, curcr, rnew)
!     call dgemm('N', 'T', ry(i-1)*n(i-1), rnew, ry(i), 1d0, cr(1,i-1), ry(i-1)*n(i-1), R, rnew, 0d0, curcr, ry(i-1)*n(i-1))
    call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, cr(1,i+1),1)
    ry(i+1) = rnew;

!     Phir
!     call phi_right(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phi(1,i+1), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i))
    call phi_left(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phi(1,i), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i+1), work, full_core)

!     write(*,"(A,I0)"), 'als_fort:  qr:', i

    i = i+dir
  end do

!   t1 = etime(tarray)
!   print *, 'qr time: ', t1-t0

!   print *, phi(1:2,1)
!   print *, phi(1:2,2)
!   print *, phi(1:2,3)
!   print *, phi(1:2,4)

  if (verb0==2) then
    print *, ''
  end if
  err_max = 0d0
  i = d;
  dir = -1;
  swp = 1
!   t0 = etime(tarray)
  do while (swp .le. nswp0)
!     write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  bfun3 started:[', swp, ',', i, ',', dir, ']'
    call bfun3(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phi(1,i), crA(pa(i)), phi(1,i+1), crX(px(i)), curcr, work, full_core)
!     write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  bfun3 done:[', swp, ',', i, ',', dir, ']'

!     print *, 'cr: ', curcr(1:ry(i)*n(i)*ry(i+1)), 'ry1: ', ry(i), 'ry2: ', ry(i+1)

    ! curcr is of size: r1,n,r2. we need r1n,(r2+kick)

!     print *, curcr(1:ry(i)*n(i)*ry(i+1))
    err = dnrm2(ry(i)*n(i)*ry(i+1), curcr(1:ry(i)*n(i)*ry(i+1)) - cr(1:ry(i)*n(i)*ry(i+1),i), 1) / dnrm2(ry(i)*n(i)*ry(i+1), curcr, 1)
    if (verb0>=2) then ! verb on, carriage_return, no matlab. Otherwise - a lot of shit comes...
      write(*,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3,20X,A$)"), 'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max, 13
!       write(*,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3)"),  'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max
    end if
!     write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3)"),  'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max


    err_max = max(err_max,err)


    if ((dir>0) .and. (i<d)) then
      rnew = min(ry(i+1), ry(i)*n(i))
      if (kickrank0>0) then
        call dgesvd('O', 'S', ry(i)*n(i), ry(i+1), curcr, ry(i)*n(i), tau, 0, 1, R, rnew, work, lwork, info)
    ! curcr is U and has the size r1*n,rnew; R is of size rnew,r2
        nn = my_chop3(rnew, tau, eps2)
! 	nn = rnew
	call row_cut(rnew, ry(i+1), nn, R)
	do j=1,nn
	  call dscal(ry(i+1), tau(j), R(j), nn)
	end do

! 	write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  svd:[', swp, ',', i, ',', dir, ']'

	rnew = nn
! 	print *, 'nn: ', nn

! 	print *, curcr(1:ry(i)*n(i)*nn)
! 	print *, '^curcr, vR'
! 	print *, R(1:nn*ry(i+1))

! 	call random_number(curcr(ry(i)*n(i)*nn+1:ry(i)*n(i)*(nn+kickrank0)))
! 	call random_number(work(1:ry(i)*n(i)*kickrank0))
! 	call dcopy(ry(i)*n(i)*kickrank0, work, 1, curcr(ry(i)*n(i)*nn+1), 1)
! 	tau(1:kickrank0*ry(i+1))=0d0
! 	print *, 'R: ', R(1:nn*ry(i+1))
! 	call row_add(nn,ry(i+1),kickrank0,R,tau)
! 	print *, 'R2: ', R(1:(nn+kickrank0)*ry(i+1))
! 	print *, nn, kickrank0, ry(i+1)
! 	print *, R(1:(nn+kickrank0)*ry(i+1))
! 	print *, R(1:nn*ry(i+1))
! 	print *, ''
! 	call dcopy((nn+kickrank0)*ry(i+1), work, 1, R, 1)

! 	call dgemm('N', 'N', ry(i)*n(i), ry(i+1), nn+kickrank0, 1d0, curcr, ry(i)*n(i), R, nn+kickrank0, 0d0, work, ry(i)*n(i))
! 	print *, 'cr2: ', work(1:ry(i)*n(i)*ry(i+1))


! 	call qr(ry(i)*n(i), nn+kickrank0, curcr, cr(1,i), work, lwork, tau)
! 	rnew = min(ry(i)*n(i), nn+kickrank0)
! 	print *, '^R  vR2'
! 	print *, cr(1:rnew*(nn+kickrank0),i)
! 	print *, '^R2  v curcr'
! 	print *, curcr(1:ry(i)*n(i)*rnew)
! 	print *, 'rnew ', rnew
! 	call dgemm('N', 'N', rnew, ry(i+1), nn+kickrank0, 1d0, cr(1,i), rnew, R, nn+kickrank0, 0d0, work, rnew)
! 	call dcopy(rnew*ry(i+1), work, 1, R, 1)

! 	call dgemm('N', 'N', ry(i)*n(i), ry(i+1), rnew, 1d0, curcr, ry(i)*n(i), R, rnew, 0d0, work, ry(i)*n(i))
! ! 	call dgemm('T', 'N', rnew, rnew, ry(i)*n(i), 1d0, curcr, ry(i)*n(i), curcr, ry(i)*n(i), 0d0, work, rnew)
! 	print *, 'cr2: ', work(1:ry(i)*n(i)*ry(i+1))
! 	print *, R(1:rnew*ry(i+1))
! 	print *, ''

! 	return
!         call transp(rnew, ry(i+1), R, work)
!         do j=1,nn
!           call dscal(ry(i+1), tau(j), work(1+(j-1)*ry(i+1)), 1)
!         end do
!         call transp(ry(i+1), nn, work, R)
      else
        call qr(ry(i)*n(i), ry(i+1), curcr, R, work, lwork, tau)
      end if

      call dcopy(ry(i)*n(i)*rnew, curcr, 1, cr(1,i), 1)
      call dgemm('N','N', rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, R, rnew, cr(1,i+1), ry(i+1), 0d0, curcr, rnew)
      call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, cr(1,i+1), 1)

      ry(i+1) = rnew
      call phi_left(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phi(1,i), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i+1), work, full_core)

!       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  phi:[', swp, ',', i, ',', dir, ']'
    end if

!     print *, phi(1:2,1)
!     print *, phi(1:2,2)
!     print *, phi(1:2,3)
!     print *, phi(1:2,4)

    if ((dir<0) .and. (i>1)) then
      rnew = min(n(i)*ry(i+1), ry(i))
      if (kickrank0>0) then
        call dgesvd('S', 'O', ry(i), n(i)*ry(i+1), curcr, ry(i), tau, R, ry(i), 0, 1, work, lwork, info)
! 	print *, 'VT: ', curcr(1:ry(i)*n(i)*ry(i+1))
! 	print *, 'R: ', R(1:ry(i)*rnew)
! 	print *, tau(1:rnew)
	! curcr is V and has the size rnew,n*r2; R is U and of size r1,rnew
        nn = my_chop3(rnew, tau, eps2)
	call transp(ry(i), n(i)*ry(i+1), curcr, work)
	call dcopy(n(i)*ry(i+1)*nn, work, 1, curcr, 1) ! now, n1,ry2,rnew
	call transp(ry(i), nn, R, work)
	call dcopy(nn*ry(i), work,1, R,1) ! rnew,r1
        do j=1,nn
          call dscal(ry(i), tau(j), R(j), nn)
        end do
!         call row_cut(ry(i), n(i)*ry(i+1), nn, curcr) ! its size is correct

! 	write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  svd:[', swp, ',', i, ',', dir, ']'
!         print *, 'calling full_core...'
!         if (i==d) then
!           print *, crA(pa(i):pa(i+1)-1)
!           print *, crX(px(i):px(i+1)-1)
!         end if
        ! kick
        call bfun3_right(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), crA(pa(i)), phi(1,i+1), crX(px(i)), full_core, full_core, work)
!         if (i==d) then
!           print *, full_core(1:rx(i)*ra(i)*n(i)*ry(i+1))
!         end if
! 	! full_core is of size rx1*ra1, n*ry2
        call duchol(rx(i)*ra(i), n(i)*ry(i+1), full_core,  min(kickrank0, n(i)*ry(i+1)), work, tau, mm);
        ! 	if (i==d) then
!           print *, work(1:n(i)*ry(i+1)*kickrank0)
!         end if
! 	print *, mm
	call dcopy(n(i)*ry(i+1)*mm, work, 1, curcr(n(i)*ry(i+1)*nn+1), 1)
! 	! work=du: n(i)*ry(i+1),kickrank

!         mm = kickrank0
! 	call random_number(curcr(n(i)*ry(i+1)*nn+1:n(i)*ry(i+1)*(nn+kickrank0)))
	tau(1:mm*ry(i))=0d0
	call row_add(nn, ry(i), mm, R, tau)

!         print *, 'VT2: ', curcr(1:nn*n(i)*ry(i+1))
!         call transp(ry(i), nn, R, work) ! R has to be transposed
! 	call dcopy(nn*ry(i), work, 1, R, 1)

	call qr(n(i)*ry(i+1), nn+mm, curcr, cr(1,i), work, lwork, tau)
	rnew = min(n(i)*ry(i+1), nn+mm)
	call dgemm('N', 'N', rnew, ry(i), nn+mm, 1d0, cr(1,i), rnew, R, nn+mm, 0d0, work, rnew)
	call dcopy(rnew*ry(i), work, 1, R, 1)

! 	write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  kick:[', swp, ',', i, ',', dir, ']'

!         rnew = nn
!         print *, rnew
      else
        call transp(ry(i), n(i)*ry(i+1), curcr, cr(1,i))
        call qr(n(i)*ry(i+1), ry(i), cr(1,i), R, work, lwork, tau)
	call dcopy(n(i)*ry(i+1)*ry(i), cr(1,i), 1, curcr, 1)
!         call transp(n(i)*ry(i+1), rnew, cr(1,i), curcr)
      end if

      call transp(n(i)*ry(i+1), rnew, curcr, cr(1,i))
!       call dcopy(rnew*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
!       call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
      call dgemm('N','T', ry(i-1)*n(i-1), rnew, ry(i), 1d0, cr(1,i-1), ry(i-1)*n(i-1), R, rnew, 0d0, curcr, ry(i-1)*n(i-1))
      call dcopy(ry(i-1)*n(i-1)*rnew, curcr, 1, cr(1,i-1), 1)

      ry(i) = rnew
      call phi_right(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phi(1,i+1), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i), work, full_core)
!       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  phi:[', swp, ',', i, ',', dir, ']'
    end if


    if ((dir>0) .and. (i==d)) then
      call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
      if (verb0>=1) then
	write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,I0,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry(1:d+1)), '  err_max:', err_max
! 	write(*,"(A,I0,A,I0,A,I0,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry), '  err_max:', err_max
! 	call mexPrintf(matlab_ist_dumme_kuh//achar(10))
      end if
      dir = -1
      i = i+1
      swp = swp+1
      if (err_max<eps) then
	exit
      end if
      err_max = 0d0
!       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  sweep_reverse:[', swp, ',', i, ',', dir, ']'
    end if

    if ((dir<0) .and. (i==1)) then
      call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
      if (verb0>=1) then
	write(*,"(A,I0,A,I0,A,I0,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry(1:d+1)), '  err_max:', err_max
! 	write(*,"(A,I0,A,I0,A,I0,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry), '  err_max:', err_max
! 	call mexPrintf(matlab_ist_dumme_kuh//achar(10))
      end if
      dir = 1
      i = i-1
!       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  sweep_reverse:[', swp, ',', i, ',', dir, ']'
    end if

    i=i+dir

  end do

  if (verb0==2) then
    print *, ''
  end if

!   t1 = etime(tarray)
!   print *, 'sweep time: ', t1-t0


  nn = sum(ry(2:d+1)*ry(1:d)*n(1:d));
  allocate(result_core(nn))

  nn = 1;
  do i=1,d
    call dcopy(ry(i)*n(i)*ry(i+1), cr(1,i), 1, result_core(nn), 1)
    nn = nn+ry(i)*n(i)*ry(i+1)
  end do
!   write(*,"(A)"), 'als_fort:  result_core dcopy'
!   print *, ry
!   call dcopy(sz, crr(pyr), 1, result_core, 1)

   deallocate(cr)
!    deallocate(crl)
   deallocate(R)
   deallocate(work)
   deallocate(curcr)
   deallocate(phi)
!    deallocate(phil)
   deallocate(tau)
   deallocate(full_core)
   deallocate(pa,px)

!   print *, pyl

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AMR lin. solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tt_amr_solve(d,n,m,ry,ra,crA, crY, crX0, rx, eps, rmax, kickrank, nswp, verb)
integer,intent(in) :: d,n(*),m(*),ry(*),ra(*), rmax
integer,intent(inout) :: rx(*)
integer, intent(in), optional :: kickrank, nswp, verb
! verb: 0: off; 1: matlab; 2: console
integer :: kickrank0, nswp0, verb0
real(8), intent(in) :: crA(*), crY(*), eps, crX0(*)
real(8) eps2
real(8),allocatable :: cr(:,:)
real(8),allocatable :: curcr(:)
real(8),allocatable :: work(:)
real(8),allocatable :: tau(:)
real(8),allocatable :: R(:)
real(8),allocatable :: B(:), res1(:), res2(:)
real(8),allocatable :: phiA(:,:), phiy(:,:)
real(8),allocatable :: full_core(:)

integer,allocatable :: pa(:), py(:)
integer :: i,j,k, swp, dir, lwork, mm, nn, rnew
integer info
real(8) :: err, err_max

real(8) dnrm2
real(8) sqrt

  kickrank0 = 5;
  if (present(kickrank)) then
    kickrank0 = kickrank
  end if

  nswp0 = 20
  if (present(nswp)) then
    nswp0 = nswp
  end if
  verb0 = 1
  if (present(verb)) then
    verb0 = verb
  end if

  allocate(pa(d+1), py(d+1))

  call compute_ps(d,ra,n(1:d)*m(1:d),pa)
  call compute_ps(d,ry,n,py)

  lwork = rmax*maxval(m(1:d))*rmax

  allocate(cr(lwork, d))
  nn = rmax*maxval(ra(1:d+1))*rmax
  allocate(phiA(nn, d+1))
  nn = maxval(ry(1:d+1))*rmax
  allocate(phiy(nn, d+1))
  allocate(curcr(lwork))
!   nn = maxval(ra(1:d+1))*maxval(ry(1:d+1))*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
  allocate(work(lwork))
  allocate(R(lwork))
  allocate(tau(lwork))
!   allocate(full_core(nn))

  mm = 1
  do i=1,d
    call dcopy(rx(i)*m(i)*rx(i+1), crX0(mm), 1, cr(1,i), 1)
    mm = mm + rx(i)*m(i)*rx(i+1)
  end do
  phiy(1,1)=1
  phiy(1,d+1)=1

  dir = 1
  i = 1
  do while (i<d)
    rnew = min(rx(i)*n(i), rx(i+1))
    call qr(rx(i)*n(i), rx(i+1), cr(1,i), R, work, lwork, tau)
    call dgemm('N', 'N', rnew, n(i+1)*rx(i+2), rx(i+1), 1d0, R, rnew, cr(1,i+1), rx(i+1), 0d0, curcr, rnew)
    call dcopy(rnew*n(i+1)*rx(i+2), curcr, 1, cr(1,i+1),1)
    rx(i+1) = rnew;

    allocate( res1(max(rx(i),rx(i+1))*max(rx(i),rx(i+1))*m(i)*max(ra(i),ra(i+1))), res2(max(rx(i),rx(i+1))*max(rx(i),rx(i+1))*m(i)*max(ra(i),ra(i+1))) )
    call phi_left(rx(i), m(i), rx(i+1), rx(i), m(i), rx(i+1), ra(i), ra(i+1), phiA(1,i), crA(pa(i)), cr(1,i), cr(1,i), phiA(1,i+1), res1, res2)
    deallocate (res1,res2)
    call phi2_left(ry(i), ry(i+1), rx(i), m(i), rx(i+1), phiy(1,i), crY(py(i)), cr(1,i), phiy(1,i+1), work)
    i = i+dir
  end do




  deallocate(work, R, tau, curcr, phiA, phiy, cr)
  deallocate(pa,py)

end subroutine

end module
