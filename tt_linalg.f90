module tt_linalg

  use iso_c_binding

  implicit none

  ! interface
  ! subroutine dbfun3(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y, res1, res2) bind(c,name="dbfun3")
  !   integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
  !   real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
  !   real(8), intent(inout) :: y(*), res1(*), res2(*)
  ! end subroutine
  !
  ! subroutine djac_apply(ptype, rx1, n, rx2, jacs, x, y, work1) bind(c,name="djac_apply")
  ! ! sizes of work1: rx1*n*rx2
  ! character, intent(in) :: ptype
  ! integer, intent(in) :: rx1, n, rx2
  ! real(8), intent(in) :: jacs(*), x(*)
  ! real(8), intent(inout) :: y(*), work1(*)
  ! end subroutine
  !
  ! end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! General full matrix and tensor routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dqr(n,m, A, R, work, lwork, tau)  bind(c)
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

  end subroutine dqr


  subroutine drow_add(m,n,k,A,B)  bind(c)
    integer, intent(in) :: m, n, k
    real(8), intent(inout) :: A(*)
    real(8), intent(in) :: B(*)
    ! real(8), intent(inout) :: C(m+k,n)
    real(8) swp(m)
    integer i

    do i=n,1,-1
       call dcopy(m, A(1+(i-1)*m), 1, swp, 1)
       call dcopy(m, swp, 1, A(1+(i-1)*(m+k)), 1)
       call dcopy(k, B(1+(i-1)*k), 1, A(m+1+(i-1)*(m+k)), 1)
       !   C(1:m,i) = A(1:m,i)
       !   C(m+1:m+k,i) = B(1:k,i)
    end do
  end subroutine drow_add

  subroutine drow_cut(m,n,k,A)   bind(c)
    integer, intent(in) :: m, n, k
    real(8), intent(inout) :: A(*)
    integer i

    do i=2,n
       call dcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
    end do
  end subroutine drow_cut


  subroutine dtransp(n, m, A, B)  bind(c)
    integer, intent(in):: n,m
    !   real(8), intent(in) :: A(*)
    !   real(8), intent(inout) :: B(*)
    real(8), intent(in):: A(n,m)
    real(8), intent(inout):: B(m,n)

    integer i,j

    !   do i=1,n
    !     do j=1,m
    !       B(j+(i-1)*m) = A(i+(j-1)*n)
    !     end do
    !   end do
    do i=1,n
       call dcopy(m, A(i,1), n, B(1,i),1)
    end do
  end subroutine dtransp

  subroutine dperm1324(n1,n2,n3,n4, A, B)  bind(c)
    integer, intent(in) :: n1,n2,n3,n4
    real(8), intent (in) :: A(n1,n2,n3,n4)
    real(8), intent (inout) :: B(n1,n3,n2,n4)

    integer i2,i3,i4

    do i4=1,n4
       do i3=1,n3
          do i2=1,n2
             call dcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
             !         do i1=1,n1
             !           B(i1,i3,i2,i4) = A(i1,i2,i3,i4)
             !         end do
          end do
       end do
    end do

  end subroutine dperm1324



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ALS-related procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dbfun3(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y, res1, res2)  bind(c)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    real(8), intent(inout) :: y(*), res1(*), res2(*)

    !   real(8), allocatable :: res1(:), res2(:), tmp(:)

    !   allocate(res1(rx1*m*ra2*ry2), res2(rx1*ra1*n*ry2), tmp(max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)))

    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    ! print *, 0
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !   print *, 1
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1, res2)
    call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
    !   print *, 2
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !   print *, 3
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call dtransp(ra1*n*ry2, rx1, res2, res1)
    !   print *, 4
    !     b1,a1,i1,c2
    !    phi1: c1, b1, a1 : ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res1, rx1*ra1, 0d0, y, ry1)
    !   print *, 5
    !     y: c1,i1,c2

    !   deallocate(res1, res2, tmp)
  end subroutine dbfun3


  subroutine dtransp2(n, m, A, B)
    integer, intent(in):: n,m
    real(8), intent(in):: A(n,m)
    real(8), intent(inout), optional, target :: B(m,n)
!     double precision, pointer :: C(:,:)
    real(8) C(m,n)
    integer i,j
!     if ( present(B) ) then
!        C => B
!     else
!        allocate(C(m,n))
!     end if
    do i=1,n
       call dcopy(m, A(i,1), n, C(1,i),1)
    end do
    if ( .not. present(B) ) then
       call dcopy(n*m, C, 1, A, 1)
    else
       call dcopy(n*m, C, 1, B, 1)
!        deallocate(C)
    end if
  end subroutine dtransp2

  subroutine dbfun32(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    real(8), intent(inout) :: y(*)
    double precision :: res1(rx1,m,ra2,ry2)
    double precision :: res2(ra1,n,ry2,rx1)
    double precision :: dnrm2
    !phi2(rx2,ra2,ry2)
    !phi1(ry1,rx1,ra1)
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp2(rx1, m*ra2*ry2, res1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n) !Here it would be a difference
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1                                                !rx1
    call dtransp2(ra1*n*ry2,rx1,res2)
    !     b1,a1,i1,c2
    !    phi1: c1, b1, a1 : ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res2, rx1*ra1, 0d0, y, ry1)
    !     y: c1,i1,c2

  end subroutine dbfun32


  subroutine dbfun3_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, A, phi2, x, y, res1, res2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi2(*), A(*), x(*)
    real(8), intent(inout) :: y(*), res1(*), res2(*)

    !   real(8), allocatable :: res1(:), res2(:), tmp(:)

    !   allocate(res1(rx1*m*ra2*ry2), res2(rx1*ra1*n*ry2), tmp(rx1*m*ra2*ry2))
    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    !   print *, rx1,m,rx2,n,ry2,ra1,ra2

    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !   if (rx1*m*ra2*ry2<10) then
    !   print *, res1(1:rx1*m*ra2*ry2)
    !   end if
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1, res2)
    call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
    !   print *, 2
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !   if (ra1*n*ry2*rx1<10) then
    !   print *, res2(1:ra1*n*ry2*rx1)
    !   end if
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call dtransp(ra1*n*ry2, rx1, res2, y)
    !   output is rx1*ra1,n,ry2
    !   deallocate(res1, res2, tmp)
  end subroutine dbfun3_right

  subroutine dbfun3_left(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, x, y, res1, res2)
    ! sizes of res1, res2: ry1*ra1*rx1, max(ra1*m,n*ra2)*ry1*rx2, ra1*n*m*ra2
    ! size of y: ry1*n, rx2*ra2
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), x(*)
    real(8), intent(inout) :: y(*), res1(*), res2(*)

    call dtransp(ry1*rx1, ra1, phi1, res1)
    call dgemm('N','N',ra1*ry1, m*rx2, rx1, 1d0, phi1, ra1*ry1, x, rx1, 0d0, res2, ra1*ry1)
    call dperm1324(ra1, ry1, m, rx2, res2, res1)
    call dperm1324(ra1, n, m, ra2, A, res2)
    call dgemm('T','N', ry1*rx2, n*ra2, ra1*m, 1d0, res1, ra1*m, res2, ra1*m, 0d0, y, ry1*rx2)
    call dperm1324(ry1, rx2, n, ra2, y, res1)
    call dcopy(ry1*n*rx2*ra2, res1, 1, y, 1)
  end subroutine dbfun3_left


  subroutine dphi_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi2_old, A, x, y, phi2, res1, res2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2, rx1*ra1*ry1)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) ::  A(*), phi2_old(*), x(*), y(*)
    real(8), intent(inout) :: phi2(*), res1(*), res2(*)

    !   real(8),allocatable :: res1(:), res2(:), tmp(:)

    !   allocate(res1(rx1*m*ra2*ry2), res2(rx1*ra1*n*ry2), tmp(max(rx1*m*ra2*ry2, rx1*ra1*n*ry2, rx1*ra1*ry1)))
    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2_old, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1, res2)
    call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call dtransp(ra1, n*ry2*rx1, res2, res1)
    !   i1,c2,b1,a1
    call dgemm('N', 'N', ry1, rx1*ra1, n*ry2, 1d0, y, ry1, res1, n*ry2, 0d0, phi2, ry1)
    ! phi2: ry1, rx1, ra1
    call dtransp(ry1, rx1*ra1, phi2, res1)
    call dcopy(rx1*ra1*ry1, res1, 1, phi2, 1)
    !   deallocate(res1,res2,tmp)
  end subroutine dphi_right
  ! y'Ax
  
  subroutine dphi_left(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1_old, A, x, y, phi1, res1, res2)
    ! sizes of res1, res2: max(rx1*n*ra1*ry2, rx1*ra2*m*ry2, ry2*rx2*ra2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) ::  A(*), phi1_old(*), x(*), y(*)
    real(8), intent(inout) :: phi1(*), res1(*), res2(*)

    !   real(8),allocatable :: res1(:), res2(:), tmp(:)

    !   allocate(res1(rx1*n*ra1*ry2), res2(rx1*ra2*m*ry2), tmp(max(rx1*n*ra1*ry2, rx1*ra2*m*ry2, ry2*rx2*ra2)))
    !   phi1: c1, b1, a1 : ry1, rx1, ra1
    !  y: c1,i1,c2: ry1,n,ry2
    call dgemm('T', 'N', rx1*ra1, n*ry2, ry1, 1d0, phi1_old, ry1, y, ry1, 0d0, res1, rx1*ra1)
    !    res1: rx1,ra1,n,ry2: b1,a1,i1,c2
    call dtransp(rx1, ra1*n*ry2, res1, res2)
    call dcopy(ra1*n*ry2*rx1, res2, 1, res1, 1)
    !     a1, i1, c2, b1
    call dgemm('T', 'N', ry2*rx1, m*ra2, ra1*n, 1d0, res1, ra1*n, A, ra1*n, 0d0, res2, ry2*rx1)
    !     res2: ry2,rx1,m,ra2 : c2, b1, j1, a2
    call dtransp(ry2, rx1*m*ra2, res2, res1)
    !   b1,j1,a2,y2
    call dgemm('T', 'N', ra2*ry2, rx2, rx1*m, 1d0, res1, rx1*m, x, rx1*m, 0d0, phi1, ra2*ry2)
    ! phi1: ra2, ry2, rx2
    call dtransp(ra2, ry2*rx2, phi1, res1)
    call dcopy(ry2*rx2*ra2, res1, 1, phi1, 1)
    !   deallocate(res1,res2,tmp)
  end subroutine dphi_left


  subroutine dphi2_right(rx1, rx2, ry1, n, ry2, phi2_old, x, y, phi2, res1)
    ! sizes of res1: rx1*n*ry2
    integer, intent(in) :: rx1, rx2, ry1, n, ry2
    real(8), intent(in) ::  phi2_old(*), x(*), y(*)
    real(8), intent(inout) :: phi2(*), res1(*)


    !   phi2: b2,c2: rx2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call dgemm('N', 'N', rx1*n, ry2, rx2, 1d0, x, rx1*n, phi2_old, rx2, 0d0, res1, rx1*n)
    !    res1: rx1,m,ry2: b1,j1,c2
    call dgemm('N', 'T', rx1, ry1, n*ry2, 1d0, res1, rx1, y, ry1, 0d0, phi2, rx1)
  end subroutine dphi2_right
  ! y'x
  subroutine dphi2_left(rx1, rx2, ry1, n, ry2, phi1_old, x, y, phi1, res1)
    ! sizes of res1: rx1*n*ry2
    integer, intent(in) :: rx1, rx2, ry1, n, ry2
    real(8), intent(in) ::  phi1_old(*), x(*), y(*)
    real(8), intent(inout) :: phi1(*), res1(*)

    !   phi1: c1, b1, a1 : ry1, rx1
    !  y: c1,i1,c2: ry1,n,ry2
    call dgemm('T', 'N', rx1, n*ry2, ry1, 1d0, phi1_old, ry1, y, ry1, 0d0, res1, rx1)
    !    res1: rx1,n,ry2: b1,j1,c2
    call dgemm('T', 'N', ry2, rx2, rx1*n, 1d0, res1, rx1*n, x, rx1*n, 0d0, phi1, ry2)
  end subroutine dphi2_left


  subroutine dBfull(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, B, res1, res2)
    ! sizes of res1, res2: max(ry1*n*rx1*m*ra2, rx2*ra2*ry2, ry1*n*ry2*rx1*m*rx2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*)
    real(8), intent(inout) :: B(*), res1(*), res2(*)

    ! phi1: ry1,rx1,ra1
    call dgemm('N', 'N', ry1*rx1, n*m*ra2, ra1, 1d0, phi1, ry1*rx1, A, ra1, 0d0, res1, ry1*rx1)
    ! res1: ry1,rx1,n,m,ra2
    call dperm1324(ry1, rx1, n, m*ra2, res1, res2)
    call dcopy(ry1*n*rx1*m*ra2, res2, 1, res1, 1)
    ! phi2: rx2,ra2,ry2
    call dtransp(rx2, ra2*ry2, phi2, res2)
    call dgemm('N', 'N', ry1*n*rx1*m, ry2*rx2, ra2, 1d0, res1, ry1*n*rx1*m, res2, ra2, 0d0, B, ry1*n*rx1*m);
    call dperm1324(ry1*n, rx1*m, ry2, rx2, B, res1)
    ! now B: ry1,n,ry2,rx1,m,rx2
    call dcopy(ry1*n*ry2*rx1*m*rx2, res1, 1, B, 1)
  end subroutine dBfull


  subroutine djac_gen(ptype, rx1, n, rx2, ra1, ra2, Phi1, A, Phi2, jacs)  bind(c)
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
       !     i = 100000
       allocate(work1(i))
       allocate(work2(i))
       allocate(ipiv(n))

       call dtransp(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          !       do info=1,rx2
          !         work2(info+(i-1)*rx2) = work1(info+(info-1)*rx2+(i-1)*rx2*rx2)
          !       end do
          call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do
       do i=0,ra1-1
          !       do info=1,rx1
          !         work1(info+(i-1)*rx1) = Phi1(info+(info-1)*rx1+(i-1)*rx1*rx1)
          !       end do
          call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do
       !     print *, 'diags prepared'

       call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
       !     print *, 'first dgemm'
       call dgemm('N', 'T', rx1*n*n, rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2, 0d0, work1, rx1*n*n)
       !     print *, '2nd dgemm'
       call dtransp(rx1, n*n*rx2, work1, work2)
       ! inversion
       do i=1,n*n
          work1(i)=0d0
       end do
       do i=1,n
          work1(i+(i-1)*n)=1d0
       end do
       !     print *, rx2*rx1
       do i=0,rx2*rx1-1
          call dcopy(n*n, work1, 1, jacs(i*n*n+1), 1)
          call dgesv(n, n, work2(i*n*n+1), n, ipiv, jacs(i*n*n+1), n, info)
          !       print *, i, '-th dgesv'
       end do
       !     print *, jacs(1:n*n*rx2*rx1)
       !     print *, 'dgesvs done'
    end if

    if ((ptype=='l').or.(ptype=='L')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*rx1*n*n*ra2, rx1*rx1*n*n*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx1*n
       allocate(ipiv(i))

       call dtransp(rx2*ra2, rx2, Phi2, work1)
       do i=0,ra2-1
          call dcopy(rx2, work1(i*rx2*rx2+1), rx2+1, work2(i*rx2+1), 1)
       end do

       call dgemm('N', 'N', rx1*rx1, n*n*ra2, ra1, 1d0, Phi1, rx1*rx1, A, ra1, 0d0, work1, rx1*rx1)
       call dperm1324(rx1, rx1, n, n*ra2, work1, jacs)
       call dgemm('N', 'T', rx1*n*rx1*n, rx2, ra2, 1d0, jacs, rx1*n*rx1*n, work2, rx2, 0d0, work1, rx1*n*rx1*n)
       !     call dcopy(rx1*n*rx1*n*rx2, work1, 1, work2, 1)
       ! inversion
       do i=1,rx1*n*rx1*n
          work2(i)=0d0
       end do
       do i=1,rx1*n
          work2(i+(i-1)*rx1*n)=1d0
       end do
       do i=0,rx2-1
          call dcopy(rx1*n*rx1*n, work2, 1, jacs(i*rx1*n*rx1*n+1), 1)
          call dgesv(rx1*n, rx1*n, work1(i*rx1*n*rx1*n+1), rx1*n, ipiv, jacs(i*rx1*n*rx1*n+1), rx1*n, info)
       end do
       !    print *, jacs(1:rx1*n*rx1*n*rx2)
    end if

    if ((ptype=='r').or.(ptype=='R')) then
       i = max(rx1*ra1, rx2*ra2*rx2, rx1*n*n*rx2*rx2)
       allocate(work1(i))
       allocate(work2(i))
       i = rx2*n
       allocate(ipiv(i))

       call dtransp(rx2*ra2, rx2, Phi2, work2)
       do i=0,ra1-1
          call dcopy(rx1, Phi1(i*rx1*rx1+1), rx1+1, work1(i*rx1+1), 1)
       end do

       call dgemm('N', 'N', rx1, n*n*ra2, ra1, 1d0, work1, rx1, A, ra1, 0d0, jacs, rx1)
       call dgemm('N', 'T', rx1*n*n, rx2*rx2, ra2, 1d0, jacs, rx1*n*n, work2, rx2*rx2, 0d0, work1, rx1*n*n)
       call dperm1324(rx1*n, n, rx2, rx2, work1, work2)
       call dtransp(rx1, n*rx2*n*rx2, work2, work1)
       ! inversion
       do i=1,rx2*n*rx2*n
          work2(i)=0d0
       end do
       do i=1,rx2*n
          work2(i+(i-1)*rx2*n)=1d0
       end do
       do i=0,rx1-1
          call dcopy(rx2*n*rx2*n, work2, 1, jacs(i*rx2*n*rx2*n+1), 1)
          call dgesv(rx2*n, rx2*n, work1(i*rx2*n*rx2*n+1), rx2*n, ipiv, jacs(i*rx2*n*rx2*n+1), rx2*n, info)
       end do
    end if

    deallocate(work1,work2, ipiv)
    !   print *, 'djac_gen done'
  end subroutine djac_gen

  subroutine djac_apply(ptype, rx1, n, rx2, jacs, x, y, work1)  bind(c)
    ! sizes of work1: rx1*n*rx2
    character, intent(in) :: ptype
    integer, intent(in) :: rx1, n, rx2
    real(8), intent(in) :: jacs(*), x(*)
    real(8), intent(inout) :: y(*), work1(*)
    integer i

    if ((ptype=='c').or.(ptype=='C')) then
       ! jacs is n,n,rx2,rx1
       call dtransp(rx1, n*rx2, x, work1)
       do i=0,(rx2*rx1-1)
          call dgemv('N', n, n, 1d0, jacs(i*n*n+1), n, work1(i*n+1), 1, 0d0, y(i*n+1), 1)
       end do
       call dtransp(n*rx2, rx1, y, work1)
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
       call dtransp(rx1, n*rx2, x, work1)
       do i=0,rx1-1
          call dgemv('N', n*rx2, n*rx2, 1d0, jacs(i*n*rx2*n*rx2+1), n*rx2, work1(i*n*rx2+1), 1, 0d0, y(i*n*rx2+1), 1)
       end do
       call dtransp(n*rx2, rx1, y, work1)
       call dcopy(rx1*n*rx2, work1, 1, y, 1)
    end if
    !   print *, 'djac_apply done'
  end subroutine djac_apply


!!!!!!!!!!!!!!
!!!! Unfinished Choletskiy !!!!!!!!
!!!!!!!!!!!!!

  subroutine duchol_fort(trans, m, n, dA, numvec, du, lam, real_numvec) bind(c)
    use dispmodule
    character, intent(in) :: trans
    integer, intent(in) :: m,n,numvec
    real(8), intent(in) :: dA(m,n)
    real(8), intent(inout) :: du(n,numvec), lam(*)
    integer, intent(inout) :: real_numvec

    real(8),allocatable :: ddiag(:), du2(:,:), work(:), dx(:), dv(:,:)
    integer i,j,cnt, maxi, lwork
    real(8) maxa, maxx

    ! functions
    real(8) ddot
    integer idamax

    if (numvec>n) then
       !write (*, "(A,I0,A,I0)"), 'UCHOL: warning: numvec=', numvec, ' > n=', n
       call disp('UCHOL: warning: numvec='//tostring(numvec*1d0)//' >n='//tostring(n*1d0))
   end if

    !   print *, m,n,numvec

    lwork = 3*numvec
    allocate(dx(numvec), dv(numvec, numvec), work(lwork), du2(n,numvec), ddiag(n))

    do j=1,numvec
       lam(j)=0d0
    end do

    !/* Initial diagonal */
    if ((trans=='n').or.(trans=='N')) then
       do i=1,n
          ddiag(i) = ddot(m, dA(1,i), 1, dA(1,i), 1)
       end do
    else
       do i=1,n
          ddiag(i) = ddot(m, dA(i,1), n, dA(i,1), n)
       end do
    end if

    !/* Main cycle */
    do cnt=1,numvec
       !/* maxel */
       maxi = idamax(n, ddiag, 1)
       maxa = ddiag(maxi)
       !         maxi = 0; maxa = 0d0;
       !    do i=1,n
       !             if (dabs(ddiag(i))>maxa) then
       !                 maxi = i;
       !                 maxa = dabs(ddiag(i));
       !             end if
       !         end do
       !         print *, 'diag: ', maxi, maxa

       !/* residual */
       if ((trans=='n').or.(trans=='N')) then
      call dgemv('T', m,n, 1d0, dA, m, dA(1,maxi), 1, 0d0, du(1,cnt), 1)
       else
      call dgemv('N', n,m, 1d0, dA, n, dA(maxi,1), n, 0d0, du(1,cnt), 1)
       end if

       do j=1,cnt-1
          call daxpy(n, -lam(j)*du(maxi,j), du(1,j), 1, du(1,cnt),1)
          !du(i,cnt) = du(i,cnt) -du(i,j)*lam(j)*du(maxi,j);
       end do
       maxa = dsqrt(dabs(du(maxi,cnt)));
       !    print *, 'after subtracion: ', maxa
       if (maxa<1d-300) then ! // We've got the exact zero matrix
          exit;
       end if
       call dscal(n, 1d0/maxa, du(1,cnt),1)
       do i=1,n
          ddiag(i) = ddiag(i)-(du(i,cnt)*du(i,cnt));
       end do
       !/* reorth */
       call dgeqrf(n, cnt, du, n, du2, work, lwork, i);
       do j=1,cnt
      dx(j) = du(j,cnt)
       end do
       call dorgqr(n, cnt, cnt, du, n, du2, work, lwork, i);

       !/* new eig */
       do i=1,cnt
      do j=1,cnt
             dv(j,i)=dx(j)*dx(i)
      end do
      dv(i,i)=dv(i,i) + lam(i)
       end do
       j=cnt+1;
       call dsyev('V', 'U', cnt, dv, numvec, lam, work, lwork, i);

       !/* update u */
       call dgemm('N','N',n,cnt,cnt,1d0,du,n,dv,numvec,0d0,du2,n);
       call dcopy(n*cnt, du2, 1, du, 1)
    end do ! cnt

    if (cnt>numvec) then
       cnt=numvec
    end if

    deallocate(ddiag, dx, dv, work, du2)
    real_numvec = cnt
  end subroutine duchol_fort


!!!!!!!!!!!!!!!!!!!!!!!
!!!!! GMRES !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!


  subroutine dgmresr_hh_fort(Phi1, A, Phi2, rhs, rx1, n, rx2, ra1, ra2, nrestart, tol, niters, ptype, jacs, sol, verb) bind(c)
    ! right preconditioned - for residual tolerance
    ! This one with Householder tranforms
    use dispmodule
    real(8),intent(in) :: Phi1(*), A(*), Phi2(*), rhs(*), jacs(*)
    integer, intent(in) :: rx1,n,rx2, ra1,ra2, nrestart, niters, verb
    real(8), intent(in) :: tol
    character, intent(in) :: ptype
    real(8), intent(inout) :: sol(*)

    integer :: i,j,it, sz;
    real(8), allocatable :: U(:,:), w(:), R(:,:), JJ(:,:), tau(:), res1(:), res2(:)
    real(8)  nrmr, curres, nrmrhs, dbeta, dalpha
    ! character last_iter

    real(8) dnrm2
    real(8) ddot


    sz = rx1*n*ra2*rx2;
    if (rx1*ra1*n*rx2>sz) then
       sz = rx1*ra1*n*rx2;
    end if
    allocate(res1(sz), res2(sz))

    sz = rx1*n*rx2;

    allocate(U(sz, nrestart+1), w(sz), tau(nrestart+1), JJ(2,nrestart), R(nrestart, nrestart))

    do j=1,sz
       sol(j)=0d0
    end do

    do it=0,niters-1
       ! r0
       if (.not.((ptype=='n').or.(ptype=='N'))) then
          call djac_apply(ptype, rx1, n, rx2, jacs, sol, w, res1)
!           call dbfun3(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, res1, res2)
          call dbfun32(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w)
       end if
       if ((ptype=='n').or.(ptype=='N')) then
!           call dbfun3(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, sol, w, res1, res2);
          call dbfun32(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, sol, w);
       end if
       call daxpy(sz,-1d0,rhs,1,w,1)
       !         call dscal(sz,-1d0,w,1)

       nrmr = dnrm2(sz,w,1);
       if (verb>1) then
      print *, 'restart: ' , it, ', res: ', nrmr
       end if
       if (it==0) then
      nrmrhs = nrmr
       end if
       if (nrmr==0d0) then
      exit
       end if
       ! initial HHT
       dbeta = nrmr;
       if (w(1)<0d0) then
      dbeta = -dbeta
       end if
       w(1) = w(1)+dbeta
       tau(1) = -dbeta
       nrmr = dnrm2(sz,w,1)
       dbeta = 1d0/nrmr
       call dscal(sz,dbeta,w,1)
       call dcopy(sz,w,1,U(1,1),1)

       do j=1,nrestart
          ! HHT on last U
          call dcopy(sz,U(1,j),1,w,1)
          dbeta = -2d0*U(j,j)
          call dscal(sz,dbeta,w,1)
          w(j) = w(j) + 1d0
          do i=j-1,1,-1
             dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
             call daxpy(sz,dbeta,U(1,i),1,w,1);
          end do
          dbeta = dnrm2(sz,w,1);
          dbeta = 1d0/dbeta;
          call dscal(sz,dbeta,w,1); ! w=w/norm(w);

          ! precvec, matvec
          if (.not.((ptype=='n').or.(ptype=='N'))) then
             call djac_apply(ptype, rx1, n, rx2, jacs, w, w, res1);
!              call dbfun3(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, res1, res2);
             call dbfun32(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w);
          end if
          if ((ptype=='n').or.(ptype=='N')) then
!              call dbfun3(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w, res1, res2);
             call dbfun32(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, w, w);
          end if

          ! Orthog w to j projectors
          do i=1,j
             dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
             call daxpy(sz,dbeta,U(1,i),1,w,1);
          end do

          ! new P_{j+1}
          if (j<sz) then
             do i=1,j
                U(i,j+1)=0d0;
             end do
             i = sz-j;
             call dcopy(i, w(j+1), 1, U(j+1, j+1), 1);
             dalpha = dnrm2(i, U(j+1, j+1), 1);
             if (.not.(dalpha==0d0)) then
                if (w(j+1)<0d0) then
                   dalpha = -dalpha;
                end if
                U(j+1, j+1) = U(j+1,j+1) + dalpha;
                dbeta = dnrm2(i, U(j+1, j+1), 1);
                dbeta = 1d0/dbeta;
                call dscal(i,dbeta,U(j+1, j+1),1);

                w(j+1) = -dalpha;
                do i=j+2,sz
                   w(i)=0d0;
                end do
             end if
          end if

          ! Givens rotators to the top of w
          do i=1,j-1
             dbeta = w(i);
             w(i) = JJ(1,i)*w(i) + JJ(2,i)*w(i+1);
             w(i+1) = -JJ(2,i)*dbeta + JJ(1,i)*w(i+1);
          end do

          ! New rotator
          if (j<sz) then
             dalpha = dnrm2(2, w(j), 1)
             !                 dalpha = sqrt((w(j)*w(j))+(w(j+1)*w(j+1)));
             JJ(1,j) = w(j)/dalpha;
             JJ(2,j) = w(j+1)/dalpha;
             tau(j+1) = -JJ(2,j)*tau(j);
             tau(j) = JJ(1,j)*tau(j);
             w(j) = dalpha;
             w(j+1) = 0d0;
          end if

          call dcopy(j, w, 1, R(1,j), 1);

          ! residual
          curres = abs(tau(j+1))/nrmrhs;
          if (verb>1) then
         !write(*,"(A,I0,A,I0,A,ES10.3)"),  'iter [', it, ',', j, '], res: ', curres
        call disp('iter ['//tostring(1d0*it)//','//tostring(j*1d0)//'], res:',tostring(curres))
        end if

          if (curres<tol) then
             exit
          end if
       end do ! inner

       if (j>nrestart) then
      j = nrestart
   !             j=nrestart-1;
   !             i=nrestart;
   !    else
   !      i=j+1;
       end if

       call dtrsv('u','n','n',j,R,nrestart,tau,1);

       ! Correction
       call dcopy(sz, U(1,j), 1, w, 1);
       dbeta = -2d0*U(j,j)*tau(j);
       call dscal(sz, dbeta, w, 1);
       w(j) = w(j) + tau(j);
       do i=j-1,1,-1
          w(i) = w(i)+tau(i);
          dbeta = -2d0*ddot(sz,U(1,i),1,w,1);
          call daxpy(sz,dbeta,U(1,i),1,w,1);
       end do
       dalpha=-1d0;
       call daxpy(sz,dalpha,w,1,sol,1);
       if (curres<tol) then
      exit;
       end if
    end do ! iters

    !     if (verb>0) sprintf(&gmres_3d_printf_buf_[0], "gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);
    !     if (verb>0) printf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);
    if (verb>0) then
       !       call dbfun3(rx1,n,rx2, rx1,n,rx2, ra1, ra2, Phi1,A,Phi2, sol, w, res1, res2);
       !       call daxpy(sz,-1d0,rhs,1,w,1)
       !write(*,"(A,I0,A,I0,A,ES10.3)"), 'gmres conducted [', it, ',', j, '] iters to relres ', curres
       call disp('gmres conducted['//tostring(it*1d0)//','//tostring(j*1d0)//'] iters to relres'//tostring(curres))
    end if

    deallocate(U,w,tau,JJ,R,res1,res2)
  end subroutine dgmresr_hh_fort

end module tt_linalg
