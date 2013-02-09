module ttnodeop

  use iso_c_binding
  use tt_lib
  use ptype_lib

  implicit none

    ! work arrays for qrs, projections, etc.
!     real(8), allocatable :: dwork1(:), dwork2(:)
!     complex(8), allocatable :: zwork1(:), zwork2(:)

    ! current nodes. Just to allocate them on-the-fly and keep sizes
    ! all sizes not zero => tt_matrix, m==0 => tt_tensor
    type,public :: dttnode
      integer :: r1, n, m, r2
      real(8), pointer :: u(:)
    end type
    type,public :: zttnode
      integer :: r1, n, m, r2
      complex(8), pointer :: u(:)
    end type


contains

  ! 1324 - permutation. Necessary for tt-matrices, etc.
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

  ! 2d permutation
  subroutine dtransp(n, m, A, B)  bind(c)
    integer, intent(in):: n,m
    real(8), intent(in):: A(n,m)
    real(8), intent(inout):: B(m,n)
    integer i,j

    do i=1,n
       call dcopy(m, A(i,1), n, B(1,i),1)
    end do
  end subroutine dtransp

  subroutine ztransp(n, m, A, B)  bind(c)
    integer, intent(in):: n,m
    complex(8), intent(in):: A(n,m)
    complex(8), intent(inout):: B(m,n)
    integer i,j

    do i=1,n
       call zcopy(m, A(i,1), n, B(1,i),1)
    end do
  end subroutine ztransp


  !!!!
  ! QRs between two blocks
  !!!!!
  subroutine dqr(C1, C2)
  ! QR left-to-right, C1->C2, real
    type(dttnode) :: C1, C2
    integer :: lwork
    integer info
    integer rnew, k,j, m
    real(8),allocatable :: R(:,:), U(:,:), tau(:), work(:)

    ! allocate the auxiliary memory
    m = C1%r1*C1%n
    lwork = 256*m
    allocate(tau(min(m,C1%r2)), work(max(lwork, C2%r1*C2%n*C2%r2)))

    allocate(U(m, C1%r2))
    call dcopy(m*C1%r2, C1%u, 1, U, 1)

    ! orthogonalize the left node (C1)
    call dgeqrf(m, C1%r2, U, m, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'dqr: dgeqrf failed'
    end if
    rnew = min(m, C1%r2)
    allocate(R(rnew, m))
    do j=1,m
      do k=1,rnew
        R(k,j)=0.0d0
      end do
    end do
    do j=1,C1%r2
      do k=1,min(j,m)
        R(k,j) = U(k,j)
      end do
    end do
    call dorgqr(m,rnew,rnew,U,m,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'dqr: dorgqr failed'
    end if
    C1%r2=rnew

    ! Copy orth. factor to C1
    deallocate(C1%u)
    allocate(C1%u(m*rnew))
    call dcopy(m*rnew, U, 1, C1%u, 1)

    ! Move R to C2
    call dgemm('N', 'N', rnew, C2%n*C2%r2, C2%r1, 1d0, R, rnew, C2%u, C2%r1, 0d0, work, rnew)
    deallocate(C2%u)
    allocate(C2%u(rnew*C2%n*C2%r2))
    call dcopy(rnew*C2%n*C2%r2, work, 1, C2%u, 1)
    C2%r1 = rnew

    deallocate(R,U,tau,work)
  end subroutine dqr

  subroutine zqr(C1, C2)
  ! QR left-to-right, C1->C2, complex
    type(zttnode) :: C1, C2
    integer :: lwork
    integer info
    integer rnew, k,j, m
    complex(8),allocatable :: R(:,:), U(:,:), tau(:), work(:)

    ! check the auxiliary memory
    m = C1%r1*C1%n
    lwork = 256*m
    allocate(tau(min(m,C1%r2)), work(max(lwork, C2%r1*C2%n*C2%r2)))

    allocate(U(m, C1%r2))
    call zcopy(m*C1%r2, C1%u, 1, U, 1)

    ! orthogonalize the left node (C1)
    call zgeqrf(m, C1%r2, U, m, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'zqr: dgeqrf failed'
    end if
    rnew = min(m, C1%r2)
    allocate(R(rnew, m))
    do j=1,m
      do k=1,rnew
        R(k,j)=(0d0, 0d0)
      end do
    end do
    do j=1,C1%r2
      do k=1,min(j,m)
        R(k,j) = U(k,j)
      end do
    end do
    call zungqr(m,rnew,rnew,U,m,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'zqr: dorgqr failed'
    end if
    C1%r2=rnew

    ! Copy orth. factor to C1
    deallocate(C1%u)
    allocate(C1%u(m*rnew))
    call zcopy(m*rnew, U, 1, C1%u, 1)

    ! Move R to C2
    call zgemm('N', 'N', rnew, C2%n*C2%r2, C2%r1, (1.0d0,0.0d0), R, rnew, C2%u, C2%r1, (0.0d0,0.0d0), work, rnew)
    deallocate(C2%u)
    allocate(C2%u(rnew*C2%n*C2%r2))
    call zcopy(rnew*C2%n*C2%r2, work, 1, C2%u, 1)
    C2%r1 = rnew

    deallocate(R,U,tau,work)
  end subroutine zqr

  subroutine drq(C1, C2)
  ! QR right-to-left, C2->C1, real
    type(dttnode) :: C1, C2
    integer :: lwork
    integer info
    integer rnew, k,j, m
    real(8),allocatable :: R(:,:), U(:,:), tau(:), work(:)

    ! check the auxiliary memory
    m = C2%n*C2%r2
    lwork = 256*m
    allocate(tau(min(m,C2%r1)), work(max(lwork, C1%r1*C1%n*C1%r2)))

    allocate(U(m, C2%r1))
    call dtransp(C2%r1, m, C2%u, U)

    ! orthogonalize the left node (C1)
    call dgeqrf(m, C2%r1, U, m, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'dqr: dgeqrf failed'
    end if
    rnew = min(m, C2%r1)
    allocate(R(rnew, m))
    do j=1,m
      do k=1,rnew
        R(k,j)=0.0d0
      end do
    end do
    do j=1,C2%r1
      do k=1,min(j,m)
        R(k,j) = U(k,j)
      end do
    end do
    call dorgqr(m,rnew,rnew,U,m,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'dqr: dorgqr failed'
    end if
    C2%r1=rnew

    ! Copy orth. factor to C2
    deallocate(C2%u)
    allocate(C2%u(rnew*m))
    call dtransp(m, rnew, U, C2%u)

    ! Move R to C2
    call dgemm('N', 'T', C1%r1*C1%n, rnew, C1%r2, 1d0, C1%u, C1%r1*C1%n, R, rnew, 0d0, work, C1%r1*C1%n)
    deallocate(C1%u)
    allocate(C1%u(C1%r1*C1%n*rnew))
    call dcopy(C1%r1*C1%n*rnew, work, 1, C1%u, 1)
    C1%r2 = rnew

    deallocate(R,U,tau,work)
  end subroutine drq

  subroutine zrq(C1, C2)
  ! QR right-to-left, C2->C1, real
    type(zttnode) :: C1, C2
    integer :: lwork
    integer info
    integer rnew, k,j, m
    complex(8),allocatable :: R(:,:), U(:,:), tau(:), work(:)

    ! check the auxiliary memory
    m = C2%n*C2%r2
    lwork = 256*m
    allocate(tau(min(m,C2%r1)), work(max(lwork, C1%r1*C1%n*C1%r2)))

    allocate(U(m, C2%r1))
    call ztransp(C2%r1, m, C2%u, U)

    ! orthogonalize the left node (C1)
    call zgeqrf(m, C2%r1, U, m, tau, work,lwork,info)
    if (info.ne.0) then
       print *, 'dqr: dgeqrf failed'
    end if
    rnew = min(m, C2%r1)
    allocate(R(rnew, m))
    do j=1,m
      do k=1,rnew
        R(k,j)=(0.0d0,0.0d0)
      end do
    end do
    do j=1,C2%r1
      do k=1,min(j,m)
        R(k,j) = U(k,j)
      end do
    end do
    call zungqr(m,rnew,rnew,U,m,tau,work,lwork,info)
    if (info.ne.0) then
       print *, 'dqr: dorgqr failed'
    end if
    C2%r1=rnew

    ! Copy orth. factor to C2
    deallocate(C2%u)
    allocate(C2%u(rnew*m))
    call ztransp(m, rnew, U, C2%u)

    ! Move R to C2
    call zgemm('N', 'T', C1%r1*C1%n, rnew, C1%r2, (1.0d0,0.0d0), C1%u, C1%r1*C1%n, R, rnew, (0.0d0,0.0d0), work, C1%r1*C1%n)
    deallocate(C1%u)
    allocate(C1%u(C1%r1*C1%n*rnew))
    call zcopy(C1%r1*C1%n*rnew, work, 1, C1%u, 1)
    C1%r2 = rnew

    deallocate(R,U,tau,work)
  end subroutine zrq

end module ttnodeop