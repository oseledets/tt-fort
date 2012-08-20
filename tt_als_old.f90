!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ALS-related procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module ttals
use matrix_util
contains
 !The algorithm for the transposed mv is quite simple
!
  subroutine bfun3(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    real(8), intent(inout) :: y(*)
    double precision :: res1(rx1,m,ra2,ry2)
    double precision :: res2(ra1,n,ry2,rx1)
    !phi2(rx2,ra2,ry2)
    !phi1(ry1,rx1,ra1) 
    call dgemm('N', 'N', rx1*m, ra2*ry2, rx2, 1d0, x, rx1*m, phi2, rx2, 0d0, res1, rx1*m)
    !    res1: rx1,m,ra2,ry2: b1,j1,a2,c2
    call transp(rx1, m*ra2*ry2, res1)
    !     j1, a2, c2, b1
    call dgemm('N', 'N', ra1*n, ry2*rx1, m*ra2, 1d0, A, ra1*n, res1, m*ra2, 0d0, res2, ra1*n) !Here it would be a difference
    !     res2: ra1,n,ry2,rx1 : a1, i1, c2, b1                                                !rx1
    call transp(ra1*n*ry2,rx1,res2)
    !     b1,a1,i1,c2
    !    phi1: c1, b1, a1 : ry1, rx1, ra1
    call dgemm('N', 'N', ry1, n*ry2, rx1*ra1, 1d0, phi1, ry1, res2, rx1*ra1, 0d0, y, ry1)
    !     y: c1,i1,c2

  end subroutine bfun3


 subroutine bfun3_transp(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, x, y)
    ! sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*), x(*)
    real(8), intent(inout) :: y(*)
    double precision :: res1(rx1,m,ra2,ry2)
    double precision :: res2(ra1,n,ry2,rx1)
    double precision At(ra1,m,n,ra2)
    call perm1324(ra1, n, m, ra2, A, At) !The transpose, but it is cheap

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

  end subroutine bfun3_transp


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
  end subroutine bfun3_right


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
  end subroutine phi_right


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
  end subroutine phi_left


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
  end subroutine phi2_right


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
  end subroutine phi2_left


  subroutine Bfull(rx1, m, rx2, ry1, n, ry2, ra1, ra2, phi1, A, phi2, B, res1, res2)
    use matrix_util
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
  end subroutine Bfull
end module ttals
