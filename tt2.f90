module tt2
  real(8), allocatable :: core(:)
  type, public ::  pointd
    double precision, dimension(:), pointer :: p=>null()
  end type pointd

contains

  subroutine full_to_tt(a,asize,n,d,eps,r)
    use matrix_util
    implicit none
    integer, intent(in) :: d
    real(8), intent(in) :: eps
    integer, intent(in) :: n(d)
    integer, intent(in) :: asize
    integer, intent(out) :: r(d+1)
    real(8), intent(inout) :: a(asize)
    integer :: rmax = 150
    double precision, allocatable :: sv(:)
    type(pointd) :: cr(d+1) 
    integer :: i, rnew, M, k, nn
    integer :: lwork, info
    double precision, allocatable :: work(:)
    M = product(n(2:d))
    lwork = 0
    r(1) = 1
    r(d+1) = 1
    allocate(sv(rmax))
    print *,'d=',d
    print *,n
    do i = 1, d
        rnew = min(r(i)*n(i),M)
        if ( lwork < 256*max(M, r(i)*n(i)) .or. .not. (allocated(work))) then
           lwork = 256*max(M, r(i)*n(i))
           if ( allocated(work) ) then
              deallocate(work)
           end if
           allocate(work(lwork))
        end if
        if ( rnew > rmax ) then
           deallocate(sv)
           rmax = 2 * rnew 
           allocate(sv(rmax))
        end if
        allocate(cr(i)%p(r(i)*n(i)*rnew))
        print *,'i=',i,'M=',M,'rnew=',rnew
        call dgesvd('s','o',r(i)*n(i),M,a,r(i)*n(i),sv,cr(i)%p,r(i)*n(i), a, rnew, work, lwork, info)
        if ( info .ne. 0 ) then
           print *,'tt2 full_to_tt failed with info=',info
           exit
        end if
        r(i+1) = my_chop3(rnew,sv,eps/sqrt(1d0*(d-1)))
        do k = 1,r(i+1)
           call dscal(r(i)*n(i), sv(k), cr(i)%p((k-1)*r(i)*n(i)+1), 1)
        end do 
        M = (r(i+1)*M)/(n(i)*r(i))
    end do 
    print *,'HERE!' 
    !And parse the input
    nn = 0
    do i = 1,d
       nn = nn + r(i) * n(i) * r(i+1) 
    end do
    if ( allocated(core) ) then
        if ( size(core) < nn ) then
            deallocate(core)
            allocate(core(nn))
        end if
    else
        allocate(core(nn))
    end if
    nn = 1
    print *,'And here!'
    do i = 1,d
        call dcopy(r(i)*n(i)*r(i+1),cr(i)%p,1,core(nn),1)
        nn = nn + r(i)*n(i)*r(i+1)
        deallocate(cr(i)%p)
    end do
    deallocate(work)
    deallocate(sv)
    print *,'And(2) here!'
 end subroutine full_to_tt
  
  subroutine tt_dealloc()
    deallocate(core)
  end subroutine tt_dealloc
 

  subroutine tt_to_full(n,r,d,cr,crsize,a,asize)
    use matrix_util
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: asize
    integer, intent(in) :: crsize
    double precision, intent(in) :: cr(crsize)
    double precision, intent(out) :: a(asize)
    double precision :: tmp(asize)
    integer :: ps(d+1) 
    integer :: i, M
    call compute_ps(d,r,n,ps)
    !call dcopy(r(1)*n(1)*r(2), cr(ps(1)), 1, a, 1)
    call eye(r(1),a) !first core is M x r(i) x r(i) x n(i) x r(i+1)
    !M = r(1; M = r(i) x n(i); 
    M = r(1)
    do i = 1,d
      call dgemm('n','n',M,n(i)*r(i+1),r(i),1d0,a,M,cr(ps(i)),r(i),0d0,tmp,M)
      call dcopy(M*n(i)*r(i+1),tmp,1,a,1)
      M = M * n(i)
    end do 
  end subroutine tt_to_full


  subroutine tt_add(n,d,r1,r2,ps1,ps2,core1,core2,rres,psres)
    implicit none
    integer, intent(in)  :: d
    integer, intent(in)  :: n(d)
    integer, intent(in)  :: r1(d+1)
    integer, intent(in)  :: r2(d+1)
    integer, intent(in)  :: ps1(d+1)
    integer, intent(in)  :: ps2(d+1)
    integer, intent(out) :: rres(d+1)
    integer, intent(out) :: psres(d+1)
    real(8), intent(in) :: core1(:)
    real(8), intent(in) :: core2(:)


  end subroutine tt_add


  subroutine tt_compr2(n,d,r,ps,cr,eps)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(inout) :: r(d+1)
    integer, intent(inout) :: ps(d+1)
    real(8), intent(in) :: cr(:)
    real(8), intent(in) :: eps
  end subroutine tt_compr2

  subroutine tt_norm(n,d,r,ps,cr,nrm)
    implicit none
    integer, intent(in) :: d
    integer, intent(in) :: n(d)
    integer, intent(in) :: r(d+1)
    integer, intent(in) :: ps(d+1)
    real(8), intent(in) :: cr(:)
    real(8), intent(out) :: nrm
    integer :: i
  
  end subroutine tt_norm
    !Later on we will avoid allocation in + and hdm, where the result have a very specific
    !size, i.e., ranks core size can be precomputed
    subroutine tt_hdm(n,d,r1,r2,ps1,ps2,core1,core2,rres,psres)
      integer, intent(in)  :: d
      integer, intent(in)  :: n(d)
      integer, intent(in)  :: r1(d+1)
      integer, intent(in)  :: r2(d+1)
      integer, intent(in)  :: ps1(d+1)
      integer, intent(in)  :: ps2(d+1)
      integer, intent(out) :: rres(d+1)
      integer, intent(out) :: psres(d+1)
      real(8),    intent(in) :: core1(:)
      real(8),    intent(in) :: core2(:)
    end subroutine tt_hdm


end module tt2
