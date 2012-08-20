!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General full matrix and tensor routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module matrix_util
contains

  subroutine eye(n,a)
    integer, intent(in) :: n
    double precision, intent(inout) :: a(n,n)
    integer i
    a(:,:) = 0d0
    do i = 1,n
       a(i,i) = 1d0
    end do
  end subroutine eye


  subroutine qr(n,m, A, R, work, lwork, tau)
    integer, intent(in):: n,m,lwork
    real(8), intent(inout) :: A(n,m), R(min(n,m),m)
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
    end do
    call dorgqr(n,rnew,rnew,A,n,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'qr: dorgqr failed'
    end if

  end subroutine qr


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
  end subroutine row_add

  subroutine row_cut(m,n,k,A)
    integer, intent(in) :: m, n, k
    real(8), intent(inout) :: A(*)
    integer i

    do i=2,n
       call dcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
    end do
  end subroutine row_cut


  subroutine transp(n, m, A, B)
    integer, intent(in):: n,m
    real(8), intent(in):: A(n,m)
    real(8), intent(inout), optional, target :: B(m,n)
    double precision, pointer :: C(:,:)
    integer i,j
    if ( present(B) ) then
       C => B
    else 
       allocate(C(m,n))
    end if
    do i=1,n
       call dcopy(m, A(i,1), n, C(1,i),1)
    end do
    if ( .not. present(B) ) then
       call dcopy(n*m, C, 1, A, 1)
       deallocate(C)
    end if
  end subroutine transp

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

  end subroutine perm1324

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
  end subroutine compute_ps

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
  end function my_chop3

end module matrix_util


