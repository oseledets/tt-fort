!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALS-related procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ttals
use matrix_util
contains
 !The algorithm for the dtransposed mv is quite simple
!
  subroutine dbfun3(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    real(8), intent(inout) :: y(*)
!    real(8) :: res1(rx1,m,ra2,ry2)
!    real(8) :: res2(ra1,n,ry2,rx1)
    real(8),allocatable :: res1(:), res2(:)
    allocate(res1(rx1*m*ra2*ry2), res2(ra1*n*ry2*rx1))
    !phi2(rx2,ra2,ry2)
    !phi1(ry1,rx1,ra1) 
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n) !Here it would be a difference
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1                                                !rx1
    call dtransp(ra1*n*ry2,rx1,res2)
    !     b1,a1,i1,c2
    !    phi1: c1, b1, a1 : ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res2, rx1*ra1, 0d0, y, ry1)
    !     y: c1,i1,c2
    deallocate(res1,res2)
  end subroutine dbfun3

  subroutine zbfun3(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    complex(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    complex(8), intent(inout) :: y(*)
    complex(8), allocatable :: res1(:), res2(:)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    allocate(res1(rx1*m*ra2*ry2), res2(ra1*n*ry2*rx1))
    !phi2(rx2,ra2,ry2)
    !phi1(ry1,rx1,ra1)  
    call zgemm('N', 'N', rx1*m, ra2*ry2, rx2, ONE, x, rx1*m, phi2, rx2, ZERO, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call ztransp(rx1, m*ra2*ry2, res1)
    !     j1, a2, c2, b1
    call zgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, ONE, A, ra1*n, res1, m*ra2, ZERO, res2, ra1*n) !Here it would be a difference
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1                                                !rx1
    call ztransp(ra1*n*ry2,rx1,res2)
    !     b1,a1,i1,c2
    !    phi1: c1, b1, a1 : ry1, rx1, ra1
    call zgemm('N', 'N', ry1, n*ry2, rx1*ra1, ONE, phi1, ry1, res2, rx1*ra1, ZERO, y, ry1)
    !     y: c1,i1,c2
    deallocate(res1, res2)
  end subroutine zbfun3


 subroutine dbfun3_transp(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    real(8), intent(inout) :: y(*)
    real(8) :: res1(rx1,m,ra2,ry2)
    real(8) :: res2(ra1,n,ry2,rx1)
    real(8) At(ra1,m,n,ra2)
    call dperm1324(ra1, n, m, ra2, A, At) !The dtranspose, but it is cheap

    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1, res2)
    call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1                                               
    call dtransp(ra1*n*ry2, rx1, res2, res1)
    !     b1,a1,i1,c2
    !    phi1: c1, b1, a1 : ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res1, rx1*ra1, 0d0, y, ry1)
    !     y: c1,i1,c2

  end subroutine dbfun3_transp

 subroutine zbfun3_transp(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    complex(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    complex(8), intent(inout) :: y(*)
    complex(8) :: res1(ry2, rx1, ra1, n)
    complex(8) :: res2(ry2, rx1, m, ra2)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    !conjx(ry1, n, ry2) * phi1(ry1, rx1, ra1) * A(ra1, n, m, ra2) * phi2(rx2, ra2, ry2)
    call zgemm('c','n', n * ry2, rx1 * ra1, ry1, ONE, x, &
               ry1 * n, phi1, ry1, ZERO, res1, n * ry2)
     !res1 is now n * ry2 * rx1 * ra1
    call ztransp(n, ry2 * rx1 * ra1, res1) 
     !res1 is now ry2 * rx1 * ra1 * n
    call zgemm('n', 'n', ry2 * rx1, m * ra2, ra1 * n, ONE, res1, ry2 * rx1, &
               A, ra1 * n, ZERO, res2, ry2 * rx1)
     !res2 is now ry2 * rx1 * m * ra2
    call ztransp(ry2, rx1 * m * ra2, res2) 
     !res2 is now rx1 * m * ra2 * ry2
    call zgemm('n', 't', rx1 * m, rx2, ra2 * ry2, ONE, res2, rx1 * m, &
               phi2, rx2, ZERO, y, rx1 * m)
    y(1:rx1 * m * rx2) = conjg(y(1:rx1 * m * rx2))
  end subroutine zbfun3_transp



  subroutine dbfun3_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, A, phi2, x, y, res1, res2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi2(*), A(*), x(*)
    real(8), intent(inout) :: y(*), res1(*), res2(*)
    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1, res2)
    call dcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call dtransp(ra1*n*ry1, rx1, res2, y)
    !   output is rx1*ra1,n,ry2
  end subroutine dbfun3_right


  subroutine zbfun3_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, A, phi2, x, y, res1, res2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    complex(8), intent(in) :: phi2(*), A(*), x(*)
    complex(8), intent(inout) :: y(*), res1(*), res2(*)
    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call zgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call ztransp(rx1, m*ra2*ry2, res1, res2)
    call zcopy(m*ra2*ry2*rx1, res2, 1, res1, 1)
    !     j1, a2, c2, b1
    call zgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call ztransp(ra1*n*ry1, rx1, res2, y)
    !   output is rx1*ra1,n,ry2
  end subroutine zbfun3_right


  subroutine dphi_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi2_old, A, x, y, phi2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2, rx1*ra1*ry1)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) ::  A(*), phi2_old(*), x(*), y(*)
    real(8), intent(inout) :: phi2(*)
!    real(8) :: res1(rx1*m*ra2*ry2)
!    real(8) :: res2(ra1*n*ry2*rx1)
    real(8), allocatable :: res1(:), res2(:)
    allocate(res1(rx1*m*ra2*ry2),res2(rx1*ra1*n*ry2))
    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2_old, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call dtransp(rx1, m*ra2*ry2, res1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n)
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call dtransp(ra1, n*ry2*rx1, res2)
    !   i1,c2,b1,a1
    call dgemm('N', 'N', ry1, rx1*ra1, n*ry2, 1d0, y, ry1, res2, n*ry2, 0d0, phi2, ry1)
    ! phi2: ry1, rx1, ra1
    call dtransp(ry1, rx1*ra1, phi2)
    deallocate(res1, res2)
  end subroutine dphi_right


  subroutine zphi_right(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi2_old, A, x, y, phi2)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2, rx1*ra1*ry1)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    complex(8), intent(in) ::  A(*), phi2_old(*), x(*), y(*)
    complex(8), intent(inout) :: phi2(*)
    complex(8) :: ycopy(ry1*n*ry2)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    complex(8) :: res1(rx1,m,ra2,ry2)
    complex(8) :: res2(ra1,n,ry2,rx1)
    ycopy(1:ry1*n*ry2) = conjg(y(1:ry1*n*ry2))
    !   phi2: b2,a2,c2: rx2, ra2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call zgemm('N', 'N', rx1*m, ra2*ry2, rx2, ONE, x, rx1*m, phi2_old, rx2, ZERO, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call ztransp(rx1, m*ra2*ry2, res1)
    !     j1, a2, c2, b1
    call zgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, ONE, A, ra1*n, res1, m*ra2, ZERO, res2, ra1*n)
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1
    call ztransp(ra1, n*ry2*rx1, res2)
    !   i1,c2,b1,a1
    call zgemm('N', 'N', ry1, rx1*ra1, n*ry2, ONE, ycopy, ry1, res2, n*ry2, ZERO, phi2, ry1)
    ! phi2: ry1, rx1, ra1
    call ztransp(ry1, rx1*ra1, phi2)
  end subroutine zphi_right

  ! y'Ax
  subroutine dphi_left(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1_old, A, x, y, phi1)
    ! sizes of res1, res2: max(rx1*n*ra1*ry2, rx1*ra2*m*ry2, ry2*rx2*ra2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) ::  A(*), phi1_old(*), x(*), y(*)
    real(8), intent(inout) :: phi1(*)
!    real(8) :: res1(rx1,ra1,n,ry2)
!    real(8) :: res2(ry2,rx1,m,ra2)
    real(8),allocatable :: res1(:), res2(:)

    allocate(res1(rx1*ra1*n*ry2), res2(ry2*rx1*m*ra2))
    !   phi1: c1, b1, a1 : ry1, rx1, ra1
    !  y: c1,i1,c2: ry1,n,ry2
    call dgemm('T', 'N', rx1*ra1, n*ry2, ry1, 1d0, phi1_old, ry1, y, ry1, 0d0, res1, rx1*ra1)
    !    res1: rx1,ra1,n,ry2: b1,a1,i1,c2
    call dtransp(rx1, ra1*n*ry2, res1)
    !     a1, i1, c2, b1
    call dgemm('T', 'N', ry2*rx1, m*ra2, ra1*n, 1d0, res1, ra1*n, A, ra1*n, 0d0, res2, ry2*rx1)
    !     res2: ry2,rx1,m,ra2 : c2, b1, j1, a2
    call dtransp(ry2, rx1*m*ra2, res2)
    !   b1,j1,a2,y2
    call dgemm('T', 'N', ra2*ry2, rx2, rx1*m, 1d0, res2, rx1*m, x, rx1*m, 0d0, phi1, ra2*ry2)
    ! phi1: ra2, ry2, rx2
    call dtransp(ra2, ry2*rx2, phi1)
    deallocate(res1,res2)
  end subroutine dphi_left

  subroutine zphi_left(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1_old, A, x, y, phi1)
    ! sizes of res1, res2: max(rx1*n*ra1*ry2, rx1*ra2*m*ry2, ry2*rx2*ra2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    complex(8), intent(in) ::  A(*), phi1_old(*), x(*), y(*)
    complex(8), intent(inout) :: phi1(*)
    complex(8) :: res1(rx1*ra1*n*ry2)
    complex(8) :: res2(ry2*rx1*m*ra2)
    complex(8) :: ycopy(ry1*n*ry2)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    ycopy(1:ry1*n*ry2) = conjg(y(1:ry1*n*ry2))

    !We have to figure out, what should be conjugated here

    !   phi1: c1, b1, a1 : ry1, rx1, ra1
    !  y: c1,i1,c2: ry1,n,ry2
    call zgemm('T', 'N', rx1*ra1, n*ry2, ry1, ONE, phi1_old, ry1, ycopy, ry1, ZERO, res1, rx1*ra1)
    !    res1: rx1,ra1,n,ry2: b1,a1,i1,c2
    call ztransp(rx1, ra1*n*ry2, res1)
    !     a1, i1, c2, b1
    call zgemm('T', 'N', ry2*rx1, m*ra2, ra1*n, ONE, res1, ra1*n, A, ra1*n, ZERO, res2, ry2*rx1)
    !     res2: ry2,rx1,m,ra2 : c2, b1, j1, a2
    call ztransp(ry2, rx1*m*ra2, res2)
    !   b1,j1,a2,y2
    call zgemm('T', 'N', ra2*ry2, rx2, rx1*m, ONE, res2, rx1*m, x, rx1*m, ZERO, phi1, ra2*ry2)
    ! phi1: ra2, ry2, rx2
    call ztransp(ra2, ry2*rx2, phi1)
  end subroutine zphi_left



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

  subroutine zphi2_right(rx1, rx2, ry1, n, ry2, phi2_old, x, y, phi2, res1)
    ! sizes of res1: rx1*n*ry2
    integer, intent(in) :: rx1, rx2, ry1, n, ry2
    complex(8), intent(in) ::  phi2_old(*), x(*), y(*)
    complex(8), intent(inout) :: phi2(*), res1(*)


    !   phi2: b2,c2: rx2, ry2
    !  x: b1,j1,b2: rx1,m,rx2
    call dgemm('N', 'N', rx1*n, ry2, rx2, 1d0, x, rx1*n, phi2_old, rx2, 0d0, res1, rx1*n)
    !    res1: rx1,m,ry2: b1,j1,c2
    call dgemm('N', 'T', rx1, ry1, n*ry2, 1d0, res1, rx1, y, ry1, 0d0, phi2, rx1)
  end subroutine zphi2_right


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

  subroutine zphi2_left(rx1, rx2, ry1, n, ry2, phi1_old, x, y, phi1, res1)
    ! sizes of res1: rx1*n*ry2
    integer, intent(in) :: rx1, rx2, ry1, n, ry2
    complex(8), intent(in) ::  phi1_old(*), x(*), y(*)
    complex(8), intent(inout) :: phi1(*), res1(*)

    !   phi1: c1, b1, a1 : ry1, rx1
    !  y: c1,i1,c2: ry1,n,ry2
    call dgemm('T', 'N', rx1, n*ry2, ry1, 1d0, phi1_old, ry1, y, ry1, 0d0, res1, rx1)
    !    res1: rx1,n,ry2: b1,j1,c2
    call dgemm('T', 'N', ry2, rx2, rx1*n, 1d0, res1, rx1*n, x, rx1*n, 0d0, phi1, ry2)
  end subroutine zphi2_left


  subroutine dBfull(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, B, res1, res2)
    use matrix_util
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

  subroutine zBfull(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, B)
    use matrix_util
    ! sizes of res1, res2: max(ry1*n*rx1*m*ra2, rx2*ra2*ry2, ry1*n*ry2*rx1*m*rx2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    complex(8), intent(in) :: phi1(*), A(*), phi2(*)
    complex(8), intent(inout) :: B(*)
    complex(8) ::  ZERO, ONE
    complex(8), allocatable :: res1(:),res2(:)
    integer :: mx
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    mx = max(ry1*n*rx1*m*ra2,rx2*ra2*ry2,ry1*n*ry2*rx1*m*rx2)
    allocate(res1(mx))
    allocate(res2(mx))
    ! phi1: ry1,rx1,ra1
    call zgemm('N', 'N', ry1*rx1, n*m*ra2, ra1, ONE, phi1, ry1*rx1, A, ra1, ONE, res1, ry1*rx1)
    ! res1: ry1,rx1,n,m,ra2
    call zperm1324(ry1, rx1, n, m*ra2, res1, res2)
    call zcopy(ry1*n*rx1*m*ra2, res2, 1, res1, 1)
    ! phi2: rx2,ra2,ry2
    call ztransp(rx2, ra2*ry2, phi2, res2)
    call zgemm('N', 'N', ry1*n*rx1*m, ry2*rx2, ra2, ONE, res1, ry1*n*rx1*m, res2, ra2, ONE, B, ry1*n*rx1*m);
    call zperm1324(ry1*n, rx1*m, ry2, rx2, B, res1)
    ! now B: ry1,n,ry2,rx1,m,rx2
    call zcopy(ry1*n*ry2*rx1*m*rx2, res1, 1, B, 1)
  end subroutine zBfull

end module ttals
