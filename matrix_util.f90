!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! General full matrix and tensor routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module matrix_util
!!$  interface qr
!!$     module procedure dqr
!!$  end interface qr
!!$
!!$  interface row_add
!!$     module procedure drow_add
!!$  end interface row_add
!!$
!!$  interface row_cut
!!$     module procedure drow_cut
!!$  end interface row_cut
!!$
  interface transp
     module procedure dtransp, ztransp, d1_transp, z1_transp
  end interface transp
!!$
!!$  interface perm1324
!!$     module procedure dperm1324
!!$  end interface perm1324

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


  subroutine dqr(n,m, A, R)
    integer, intent(in):: n,m
    real(8), intent(inout) :: A(n,m), R(min(n,m),m)
    integer, parameter :: nb = 256
    integer :: lwork 
    double precision :: tau(min(n,m))
    double precision :: work(nb*n)
    integer info
    integer rnew, k,j
    lwork = nb*n
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

  end subroutine dqr


  subroutine zqr(n,m, A, R)
    integer, intent(in):: n,m
    complex(8), intent(inout) :: A(n,m), R(min(n,m),m)
    integer, parameter :: nb = 256
    integer :: lwork 
    complex(8) :: tau(min(n,m))
    complex(8) :: work(nb*n)
    complex(8) :: ZERO
    parameter( ZERO=(0.0d0,0.0d0) )
    integer info
    integer rnew, k,j
    lwork = nb*n
    call zgeqrf(n, m, A, n, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'zqr: zgeqrf failed'
    end if
    rnew = min(n,m)
    R(:,:)=ZERO
    do j=1,m
       R(1:min(j,n),j)=A(1:min(j,n),j)
    end do
    call zungqr(n,rnew,rnew,A,n,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'zqr: zungqr failed'
    end if

  end subroutine zqr


  subroutine drow_add(m,n,k,A,B)
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
  end subroutine drow_add
  
  subroutine zrow_add(m,n,k,A,B)
    integer, intent(in) :: m, n, k
    complex(8), intent(inout) :: A(*)
    complex(8), intent(in) :: B(*)
    complex(8) swp(m)
    integer i

    do i=n,1,-1
       call zcopy(m, A(1+(i-1)*m), 1, swp, 1)
       call zcopy(m, swp, 1, A(1+(i-1)*(m+k)), 1)
       call zcopy(k, B(1+(i-1)*k), 1, A(m+1+(i-1)*(m+k)), 1)
    end do
  end subroutine zrow_add

  subroutine drow_cut(m,n,k,A)
    integer, intent(in) :: m, n, k
    real(8), intent(inout) :: A(*)
    integer i

    do i=2,n
        call dcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
    end do
  end subroutine drow_cut

  subroutine zrow_cut(m,n,k,A)
    integer, intent(in) :: m, n, k
    complex(8), intent(inout) :: A(*)
    integer i

    do i=2,n
       call zcopy(k, A(1+(i-1)*m), 1, A(1+(i-1)*k), 1)
    end do
  end subroutine zrow_cut

  subroutine d1_transp(n, m, A, B)
    integer, intent(in):: n,m
    real(8), intent(in):: A(:)
    real(8), intent(inout), optional, target :: B(:)
    if ( present(B) ) then
        call dtransp(n, m, A, B)
    else
        call dtransp(n, m, A)
    end if
  end subroutine d1_transp

  subroutine dtransp(n, m, A, B)
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
  end subroutine dtransp
  
  subroutine ztransp(n, m, A, B)
    integer, intent(in):: n,m
    complex(8), intent(in):: A(n,m)
    complex(8), intent(inout), optional, target :: B(m,n)
    complex(8), pointer :: C(:,:)
    integer i,j
    if ( present(B) ) then
       C => B
    else 
       allocate(C(m,n))
    end if
    do i=1,n
       call zcopy(m, A(i,1), n, C(1,i),1)
    end do
    if ( .not. present(B) ) then
       call zcopy(n*m, C, 1, A, 1)
       deallocate(C)
    end if
  end subroutine ztransp
  
  subroutine z1_transp(n, m, A, B)
    integer, intent(in):: n,m
    complex(8), intent(in):: A(:)
    complex(8), intent(inout), optional, target :: B(:)
    if ( present(B) ) then
        call ztransp(n, m, A, B)
    else
        call ztransp(n, m, A)
    end if

  end subroutine z1_transp

  subroutine dperm1324(n1,n2,n3,n4, A, B)
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

  end subroutine dperm1324
  
  subroutine zperm1324(n1,n2,n3,n4, A, B)
    integer, intent(in) :: n1,n2,n3,n4
    complex(8), intent (in) :: A(n1,n2,n3,n4)
    complex(8), intent (inout) :: B(n1,n3,n2,n4)

    integer i2,i3,i4

    do i4=1,n4
       do i3=1,n3
          do i2=1,n2
             call zcopy(n1, A(1,i2,i3,i4), 1, B(1,i3,i2,i4), 1)
          end do
       end do
    end do

  end subroutine zperm1324

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
    integer function my_chop3(n, s, eps)
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
  
  
  subroutine dperm321(n1, n2, n3, A, B)
    integer, intent(in):: n1,n2,n3
    real(8), intent(in):: A(n1,n2,n3)
    real(8), intent(inout), optional, target :: B(n3,n2,n1)
    double precision, pointer :: C(:,:,:)
    integer i1,i2,i3
    if ( present(B) ) then
       C => B
    else 
       allocate(C(n3,n2,n1))
    end if
    do i1 = 1,n1
       do i2 = 1,n2
          do i3 = 1,n3
             C(i3,i2,i1) = A(i1,i2,i3)
          end do 
       end do
    end do 
    if ( .not. present(B) ) then
       call dcopy(n1*n2*n3, C, 1, A, 1)
       deallocate(C)
    end if
  end subroutine dperm321
  
  subroutine zperm321(n1, n2, n3, A, B)
    integer, intent(in):: n1,n2,n3
    complex(8), intent(in):: A(n1,n2,n3)
    complex(8), intent(inout), optional, target :: B(n3,n2,n1)
    complex(8), pointer :: C(:,:,:)
    integer i1,i2,i3
    if ( present(B) ) then
       C => B
    else 
       allocate(C(n3,n2,n1))
    end if
    do i1 = 1,n1
       do i2 = 1,n2
          do i3 = 1,n3
             C(i3,i2,i1) = A(i1,i2,i3)
          end do 
       end do
    end do 
    if ( .not. present(B) ) then
       call zcopy(n1*n2*n3, C, 1, A, 1)
       deallocate(C)
    end if
  end subroutine zperm321

end module matrix_util


