module tt_f90m
  use tt_lib
  use ttop_lib
  use time_lib
  use python_conv_lib

contains

  subroutine full_to_tt(a,n,d,eps,r,ps,core)
    implicit none
    integer, intent(in) :: d
    real(8), intent(in) :: eps
    integer, intent(in) :: n(:)
    integer, intent(out) :: r(d+1)
    integer, intent(out) :: ps(d+1)
    real(8), intent(in) :: a(:)
    real(8), intent(out),allocatable :: core(:)
    type(dtt) :: tt
    character(len=*),parameter :: subnam='full_to_tt'
    integer :: i,M,lwork,info
    !print *,'size(n)=',size(n),'size(a)=',size(a),'d=',d
    !Print the array
    !print *, a(1:32)
    call dtt_svd0(n,a,tt,eps)
    call sdv_to_arrays(n,r,d,ps,core,tt)
    !Bring it back to the pythonic line-format
    !r=tt%r(tt%l-1:tt%m)
    !Calculate core positions + memory for the cores
    !ps(1) = 1
    !do i=1,d
    !   ps(i+1)=ps(i)+tt%n(i)*r(i)*r(i+1)
    !end do
    !M=sum((tt%n(tt%l:tt%m)*r(1:d)*r(2:d+1)))
    !M=ps(d+1)-1
    !print *,'Memory for the cores:', M
    !print *,'Positions',ps
    !print *,'FSIZE:',sum((tt%n(tt%l:tt%m)*r(1:d)*r(2:d+1)))
    !allocate(core(M))
    !do i=1,d
    !   call dcopy(r(i)*n(i)*r(i+1),tt%u(i)%p,1,core(ps(i)),1)
    !end do
    call dealloc(tt)
  end subroutine full_to_tt

  !a should be preallocated, and filled by zeros
  subroutine tt_to_full(n,r,d,ps,cr,a)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: ps(d+1)
    real(8), intent(in) :: cr(:)
    real(8), intent(inout) :: a(:)
    type(dtt) :: tt
    !print *,'d=',d
    !print *,'modes:',n
    !print *,'ranks:',r
    !print *,'number of cells in the output:',size(a)
    call arrays_to_sdv(n,r,d,ps,cr,tt)
    !print *,'wmodes',tt%n(1:d)
    !print *,'wranks',tt%r(tt%l-1:tt%m)
    call dtt_full(tt,a)
    call dealloc(tt)
  end subroutine tt_to_full
  !  Internal subroutine for the conversion of the array format
  !  to the local fortran format


  subroutine tt_add(n,d,r1,r2,ps1,ps2,core1,core2,rres,psres,core3)
    implicit none
    integer, intent(in)  :: d
    integer, intent(in)  :: n(d)
    integer, intent(in)  :: r1(d+1)
    integer, intent(in)  :: r2(d+1)
    integer, intent(in)  :: ps1(d+1)
    integer, intent(in)  :: ps2(d+1)
    integer, intent(out) :: rres(d+1)
    integer, intent(out) :: psres(d+1)
    real(8), intent(out),allocatable :: core3(:)
    real(8),    intent(in) :: core1(:)
    real(8),    intent(in) :: core2(:)
    type(dtt) :: tt1, tt2, tt
    call arrays_to_sdv(n,r1,d,ps1,core1,tt1)
    call arrays_to_sdv(n,r2,d,ps2,core2,tt2)
    call dtt_axpy(1.d0,tt1,1.d0,tt2)
    !tt1 is the sum
    call sdv_to_arrays(n,rres,d,psres,core3,tt2)

    call dealloc(tt1)
    call dealloc(tt2)


  end subroutine tt_add


  subroutine tt_compr2(n,d,r,ps,cr,eps)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(inout) :: r(d+1)
    integer, intent(inout) :: ps(d+1)
    real(8), intent(inout),allocatable :: cr(:)
    real(8), intent(in) :: eps
    type(dtt) :: tt
    call arrays_to_sdv(n,r,d,ps,cr,tt)
    call dtt_svd(tt,eps)
    call sdv_to_arrays(n,r,d,ps,cr,tt)
  end subroutine tt_compr2

  subroutine tt_norm(n,d,r,ps,cr,nrm)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: ps(d+1)
    real(8), intent(in) :: cr(:)
    real(8), intent(out) :: nrm
    type(dtt) :: tt
    !real(8):: t1,t2
    integer :: i
    call arrays_to_sdv(n,r,d,ps,cr,tt)
    nrm=dtt_norm(tt)
    end subroutine tt_norm
    !Later on we will avoid allocation in + and hdm, where the result have a very specific
    !size, i.e., ranks core size can be precomputed
    subroutine tt_hdm(n,d,r1,r2,ps1,ps2,core1,core2,rres,psres,core3)
      integer, intent(in)  :: d
      integer, intent(in)  :: n(d)
      integer, intent(in)  :: r1(d+1)
      integer, intent(in)  :: r2(d+1)
      integer, intent(in)  :: ps1(d+1)
      integer, intent(in)  :: ps2(d+1)
      integer, intent(out) :: rres(d+1)
      integer, intent(out) :: psres(d+1)
      real(8), intent(out),allocatable :: core3(:)
      real(8),    intent(in) :: core1(:)
      real(8),    intent(in) :: core2(:)
      type(dtt) :: tt1, tt2, tt
      call arrays_to_sdv(n,r1,d,ps1,core1,tt1)
      call arrays_to_sdv(n,r2,d,ps2,core2,tt2)
      call dtt_ha1(tt1,tt2,tt)
      call dealloc(tt1)
      call dealloc(tt2)
      call sdv_to_arrays(n,rres,d,psres,core3,tt)
      call dealloc(tt)
    end subroutine tt_hdm

! Check, if we can call an external python function from Fortran



!!$  subroutine tt_test(n,d,r1,r2,ps1,ps2,core1,core2,rres,psres)
!!$    implicit none
!!$    integer, intent(in)  :: d
!!$    integer, intent(in)  :: n(d)
!!$    integer, intent(in)  :: r1(d+1)
!!$    integer, intent(in)  :: r2(d+1)
!!$    integer, intent(in)  :: ps1(d+1)
!!$    integer, intent(in)  :: ps2(d+1)
!!$    integer, intent(out) :: rres(d+1)
!!$    integer, intent(out) :: psres(d+1)
!!$    real(8),    intent(in) :: core1(:)
!!$    real(8),    intent(in) :: core2(:)
!!$    type(dtt) :: tt1, tt2, tt
!!$    call arrays_to_sdv(n,r1,d,ps1,core1,tt1)
!!$    call arrays_to_sdv(n,r2,d,ps2,core2,tt2)
!!$    call sdv_to_arrays(n,rres,d,psres,core,tt2)
!!$
!!$    call dealloc(tt1)
!!$    call dealloc(tt2)
!!$
!!$
!!$  end subroutine tt_test

end module tt_f90m
