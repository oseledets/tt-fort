module python_conv_lib
  use tt_lib
  interface sdv_to_arrays
    module procedure dsdv_to_arrays, zsdv_to_arrays
  end interface sdv_to_arrays
  interface arrays_to_sdv
    module procedure darrays_to_sdv, zarrays_to_sdv
  end interface
contains
  !   Internal subroutine for the conversion of the sdv format to the array format
  !   cr is assumed to be not allocated
  subroutine dsdv_to_arrays(n,r,d,ps,cr,tt)
    implicit none
    integer, intent(in) :: d
    integer, intent(inout) :: n(d)
    integer, intent(out) :: r(d+1)
    integer, intent(out) :: ps(d+1)
    real(8), allocatable, intent(out) :: cr(:)
    type(dtt), intent(in) :: tt
    integer :: M,i
    r(1:d+1) = tt%r(tt%l-1:tt%m)
    n(1:d) = tt%n(1:d) 
    ps(1) = 1
    do i=1,d
       ps(i+1)=ps(i) + tt%n(i) * r(i) * r(i+1)
    end do
    M=ps(d+1)-1
    allocate(cr(M))
    do i=1,d
        call dcopy(r(i)*n(i)*r(i+1),tt%u(i)%p,1,cr(ps(i)),1)
    end do

  end subroutine dsdv_to_arrays

  subroutine zsdv_to_arrays(n,r,d,ps,cr,tt)
    implicit none
    integer, intent(in) :: d
    integer, intent(inout) :: n(d)
    integer, intent(out) :: r(d+1)
    integer, intent(out) :: ps(d+1)
    complex(8), allocatable, intent(out) :: cr(:)
    type(ztt), intent(in) :: tt
    integer :: M,i
    r=tt%r(tt%l-1:tt%m)
    n(1:d) = tt%n(1:d) 
    ps(1) = 1
    do i=1,d
       ps(i+1)=ps(i)+tt%n(i)*r(i)*r(i+1)
    end do
    M=ps(d+1)-1
    allocate(cr(M))
    do i=1,d
       call zcopy(r(i)*n(i)*r(i+1),tt%u(i)%p,1,cr(ps(i)),1)
    end do

  end subroutine zsdv_to_arrays


  subroutine darrays_to_sdv(n,r,d,ps,cr,tt)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: ps(d+1)
    real(8), intent(in) :: cr(:)
    type(dtt), intent(out) :: tt
    integer :: i

    tt%l = 1
    tt%m = d
    tt%r(0:d)= r(1:d+1)
    tt%n(1:d) = n(1:d)
    call alloc(tt)
    do i=1,d
       call dcopy(r(i)*n(i)*r(i+1),cr(ps(i)),1,tt%u(i)%p,1)
    end do
  end subroutine darrays_to_sdv
  
  subroutine zarrays_to_sdv(n,r,d,ps,cr,tt)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: ps(d+1)
    complex(8), intent(in) :: cr(:)
    type(ztt), intent(out) :: tt
    integer :: i

    tt%l = 1
    tt%m = d
    tt%r(0:d)= r(1:d+1)
    tt%n(1:d) = n(1:d)
    call alloc(tt)
    do i=1,d
       call zcopy(r(i)*n(i)*r(i+1),cr(ps(i)),1,tt%u(i)%p,1)
    end do
  end subroutine zarrays_to_sdv


end module python_conv_lib
