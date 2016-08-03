module bfun_diag_ksl
  !Here we need to have all parameters required for the matrix-by-vector product (to call bfun3)
  integer, private  ::  rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T
  real(8), pointer, private :: phi1T(:), phi2T(:),res1T(:), res2T(:), AT(:)
  complex(8), pointer, private :: zphi1T(:), zphi2T(:), zres1T(:), zres2T(:), zAT(:)
  integer,private ::  xsizeT, ysizeT
  
  type, public ::  dpointd
     real(8), dimension(:), pointer :: p=>null()
  end type dpointd

  type, public :: zpointd
     complex(8), dimension(:), pointer :: p=>null()
  end type zpointd

contains 

  subroutine dmatvec(x, y) !This would be very simple diagonal matvec by setting mT = nT = 1
    use ttals, only: dbfun3
    implicit none
    real(8), intent(in) :: x(*)
    real(8) :: y(*)
    call dbfun3(rx1T, 1, rx2T, ry1T, 1, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
  end subroutine dmatvec

  subroutine zmatvec(x, y)
    use ttals, only: zbfun3
    implicit none
    complex(8), intent(in) :: x(*)
    complex(8) :: y(*)
    call zbfun3(rx1T, 1, rx2T, ry1T, 1, ry2T, ra1T, ra2T, zphi1T, zAT, zphi2T, x, y)
 end subroutine zmatvec


  subroutine dmatvec_transp(x, y)
    use ttals, only: dbfun3_transp
    implicit none
    real(8), intent(in) :: x(*)
    real(8) :: y(*)
    call dmatvec(x,y)
    !call bfun3_transp(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
  end subroutine dmatvec_transp


  subroutine zmatvec_transp(x, y)
    use ttals, only: zbfun3_transp
    implicit none
    complex(8), intent(in) :: x(*)
    complex(8) :: y(*)
    call zbfun3_transp(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, zphi1T, zAT, zphi2T, x, y)
  end subroutine zmatvec_transp


  subroutine init_bfun_sizes(rx1, m, rx2, ry1, n, ry2, ra1, ra2, xsize, ysize)
    implicit none
    integer, intent(in) :: rx1, m, rx2, ry1, n, ry2, ra1, ra2, xsize, ysize
    rx1T = rx1
    rx2T = rx2
    mT = m
    nT = n
    ry1T = ry1
    ry2T = ry2
    ra1T = ra1
    ra2T = ra2
    xsizeT = xsize
    ysizeT = ysize
  end subroutine init_bfun_sizes


  subroutine dinit_bfun_main(phi1, A, phi2) !This is now an initialization of the diagonal block
    real(8), target :: phi1(:), phi2(:), A(:)
    phi1T => phi1
    phi2T => phi2
    AT => A
  end subroutine dinit_bfun_main


  subroutine zinit_bfun_main(phi1, A, phi2)
    complex(8), target :: phi1(:), phi2(:), A(:)
    zphi1T => phi1
    zphi2T => phi2
    zAT => A
  end subroutine zinit_bfun_main


end module bfun_diag_ksl

module sfun_diag_ksl

  integer, private :: rx1, rx2, ra, ry1, ry2
  real(8), private, pointer :: phi1(:), phi2(:)
  complex(8), private, pointer :: zphi1(:), zphi2(:)


contains

  subroutine dinit_sfun(rx1T, rx2T, raT, ry1T, ry2T, phi1T, phi2T)
    implicit none
    integer, intent(in) ::  rx1T, rx2T, raT, ry1T, ry2T 
    real(8), intent(in), target :: phi1T(:), phi2T(:)
    phi1 => phi1T
    phi2 => phi2T
    rx1 = rx1T
    rx2 = rx2T
    ra = raT
    ry1 = ry1T
    ry2 = ry2T
  end subroutine dinit_sfun

  subroutine zinit_sfun(rx1T, rx2T, raT, ry1T, ry2T, phi1T, phi2T)
    implicit none
    integer, intent(in) ::  rx1T, rx2T, raT, ry1T, ry2T 
    complex(8), intent(in), target :: phi1T(:), phi2T(:)
    zphi1 => phi1T
    zphi2 => phi2T
    rx1 = rx1T
    rx2 = rx2T
    ra = raT
    ry1 = ry1T
    ry2 = ry2T
  end subroutine zinit_sfun


  subroutine dsfun_matvec(Sx, Sy)
    use matrix_util
    real(8), intent(in) :: Sx(*)
    real(8), intent(out) :: Sy(*)
    real(8) :: res1(rx1,ra,ry2)
    !phi1(ry1,rx1,ra)*phi2(rx2,ra,ry2)*S(rx1,rx2); S(rx2,rx1)
    !S(rx1,rx2)*phi2(rx2,ra,ry2) = res(rx1,ra,ry2)*phi1(ry1,rx1,ra) 
    !call dtransp(rx2,rx1,Sx)
    call dgemm('n','n',rx1,ra*ry2,rx2,1d0,Sx,rx1,phi2,rx2,0d0, res1, rx1)
    call dgemm('n','n',ry1,ry2,rx1*ra,1d0,phi1,ry1,res1,rx1*ra,0d0,Sy,ry1)
    !call dtransp(rx1,rx2,Sx)
    !call dtransp(ry1,ry2,Sy)
  end subroutine dsfun_matvec

  subroutine dsfun_matvec_transp(Sy, Sx)
    use matrix_util,  only: dtransp
    real(8), intent(in) :: Sy(*)
    real(8), intent(out) :: Sx(*)
    real(8) res1(ry2,rx1,ra)
    !phi1(ry1,rx1,ra)*phi2(rx2,ra,ry2)*Sy(ry1,ry2) -> I think we will need at least one dtranspose :( 
    ! (ry1 * (rx1*ra)) -> ry2 x rx2 x ra
    call dgemm('t','n',ry2, rx1*ra, ry1, 1d0, Sy, ry1, phi1, ry1,0d0, res1, ry2)
    !res1 is ry2 x rx1 x ra, conv with phi2(rx2,ra,ry2) 
    call dtransp(ry2*rx1,ra, res1)
    !res1 is ra x ry2 x rx1
    call dgemm('t','t',rx1,rx2,ra*ry2,1d0, res1, ra*ry2, phi2, rx2, 0d0, Sx, rx1)
    !Seems to be OK.
  end subroutine dsfun_matvec_transp

  subroutine zsfun_matvec(Sx, Sy)
    use matrix_util
    complex(8), intent(in) :: Sx(*)
    complex(8), intent(out) :: Sy(*)
    complex(8) :: res1(rx1, ra, ry2)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    !phi1(ry1,rx1,ra)*phi2(rx2,ra,ry2)*S(rx1,rx2); S(rx2,rx1)
    !S(rx1,rx2)*phi2(rx2,ra,ry2) = res(rx1,ra,ry2)*phi1(ry1,rx1,ra) 
    call zgemm('n', 'n', rx1, ra*ry2, rx2, ONE, Sx, rx1, zphi2, rx2, ZERO, res1, rx1)
    call zgemm('n', 'n', ry1, ry2, rx1*ra, ONE, zphi1, ry1, res1, rx1*ra, ZERO, Sy, ry1)
    !call dtransp(rx1,rx2,Sx)
    !call dtransp(ry1,ry2,Sy)
  end subroutine zsfun_matvec

  subroutine zsfun_matvec_transp(Sy, Sx)
    use matrix_util,  only: ztransp
    complex(8), intent(in) :: Sy(*)
    complex(8), intent(out) :: Sx(*)
    complex(8) res1(ry2,rx1,ra)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    !phi1(ry1,rx1,ra)*phi2(rx2,ra,ry2)*Sy(ry1,ry2) -> I think we will need at least one dtranspose :( 
    ! (ry1 * (rx1*ra)) -> ry2 x rx2 x ra
    call zgemm('t','n',ry2, rx1*ra, ry1, ONE, Sy, ry1, zphi1, ry1,ZERO, res1, ry2)
    !res1 is ry2 x rx1 x ra, conv with phi2(rx2,ra,ry2) 
    call ztransp(ry2*rx1,ra, res1)
    !res1 is ra x ry2 x rx1
    call zgemm('t','t',rx1,rx2,ra*ry2,ONE, res1, ra*ry2, zphi2, rx2, ZERO, Sx, rx1)
    !Seems to be OK.
  end subroutine zsfun_matvec_transp

end module sfun_diag_ksl


module dyn_diag_tt
  implicit none
  real(8), allocatable :: dresult_core(:)
  complex(8), allocatable :: zresult_core(:)
contains
  subroutine deallocate_result
    if ( allocated(dresult_core) ) then
       deallocate(dresult_core)
    end if
    if ( allocated(zresult_core) ) then
       deallocate(zresult_core) 
    end if
  end subroutine deallocate_result
  
  
  subroutine dextract_slice(ra1, n, ra2, core, core_number, res)
    implicit none
    integer :: ra1, n, ra2, core_number
    real(8) :: core(ra1, n, ra2)
    real(8) :: res(ra1, ra2)
    res(:, :) = core(:, core_number, :)
  end subroutine 


  subroutine dset_slice(ra1, n, ra2, core, core_number, res)
    implicit none
    integer :: ra1, n, ra2, core_number
    real(8) :: core(ra1, n, ra2)
    real(8) :: res(ra1, ra2)
    core(:, core_number, :) = res(:, :)
  end subroutine 


  subroutine zextract_slice(ra1, n, ra2, core, core_number, res)
    implicit none
    integer :: ra1, n, ra2, core_number
    complex(8) :: core(ra1, n, ra2)
    complex(8) :: res(ra1, ra2)
    res(:, :) = core(:, core_number, :)
  end subroutine 


  subroutine zset_slice(ra1, n, ra2, core, core_number, res)
    implicit none
    integer :: ra1, n, ra2, core_number
    complex(8) :: core(ra1, n, ra2)
    complex(8) :: res(ra1, ra2)
    core(:, core_number, :) = res(:, :)
  end subroutine 


  subroutine dphi_diag_left(ry1, n, ry2, ra1, ra2, phi, A, x, res)
      use ttals
      implicit none
      integer :: ry1, n, ry2, ra1, ra2
      real(8) :: A(ra1, n, ra2), x(ry1, n, ry2), phi(*), res(ry2, ry2, ra2)
      real(8) :: tmp(ry2, ry2, ra2), Aslice(ra1, ra2), Xslice(ry1, ry2)
      integer :: k
      res(:, :, :) = 0d0    
      do k = 1, n
          Aslice = A(:, k, :)
          Xslice = x(:, k, :)
          call dphi_left(ry1, 1, ry2, ry1, 1, ry2, ra1,&
              ra2, phi, Aslice, Xslice, Xslice, &
              tmp)
          res(:, :, :) = res(:, :, :) + tmp(:, :, :)
      end do
  end subroutine dphi_diag_left
  
  subroutine zphi_diag_left(ry1, n, ry2, ra1, ra2, phi, A, x, res)
      use ttals
      implicit none
      integer :: ry1, n, ry2, ra1, ra2
      complex(8) :: A(ra1, n, ra2), x(ry1, n, ry2), phi(*), res(ry2, ry2, ra2)
      complex(8) :: tmp(ry2, ry2, ra2), Aslice(ra1, ra2), Xslice(ry1, ry2)
      integer :: k
      res(:, :, :) = (0d0, 0d0)    
      do k = 1, n
          Aslice = A(:, k, :)
          Xslice = x(:, k, :)
          call zphi_left(ry1, 1, ry2, ry1, 1, ry2, ra1,&
              ra2, phi, Aslice, Xslice, Xslice, &
              tmp)
          res(:, :, :) = res(:, :, :) + tmp(:, :, :)
      end do
  end subroutine zphi_diag_left
  
  subroutine dphi_diag_right(ry1, n, ry2, ra1, ra2, phi, A, x, res)
      use ttals
      implicit none
      integer :: ry1, n, ry2, ra1, ra2
      real(8) :: A(ra1, n, ra2), x(ry1, n, ry2), phi(*), res(ry1, ry1, ra1)
      real(8) :: tmp(ry1, ry1, ra1), Aslice(ra1, ra2), Xslice(ry1, ry2)
      integer :: k
      res(:, :, :) = 0d0    
      do k = 1, n
          Aslice = A(:, k, :)
          Xslice = x(:, k, :)
          call dphi_right(ry1, 1, ry2, ry1, 1, ry2, ra1,&
              ra2, phi, Aslice, Xslice, Xslice, &
              tmp)
          res(:, :, :) = res(:, :, :) + tmp(:, :, :)
      end do
  end subroutine dphi_diag_right
  
  subroutine zphi_diag_right(ry1, n, ry2, ra1, ra2, phi, A, x, res)
      use ttals
      implicit none
      integer :: ry1, n, ry2, ra1, ra2
      complex(8) :: A(ra1, n, ra2), x(ry1, n, ry2), phi(*), res(ry1, ry1, ra1)
      complex(8) :: tmp(ry1, ry1, ra1), Aslice(ra1, ra2), Xslice(ry1, ry2)
      integer :: k
      res(:, :, :) = (0d0, 0d0)    
      do k = 1, n
          Aslice = A(:, k, :)
          Xslice = x(:, k, :)
          call zphi_right(ry1, 1, ry2, ry1, 1, ry2, ra1,&
              ra2, phi, Aslice, Xslice, Xslice, &
              tmp)
          res(:, :, :) = res(:, :, :) + tmp(:, :, :)
      end do
  end subroutine zphi_diag_right

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! KLS-scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! What we have: we have a starting vector + a matrix (no vector X!) 
  subroutine dtt_diag_ksl(d, n, ra, crA, crY0, ry, tau, rmax, kickrank, nswp, verb, typ0, order0)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_diag_ksl !bfun declarations
    use sfun_diag_ksl !We have to do this
    use explib !Krylov exponential
    use rnd_lib
    implicit none
    integer,intent(in) :: d, n(d), ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb, typ0
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    real(8), intent(in) :: crA(*), crY0(*)
    real(8), intent(in) :: tau
    real(8) :: tau0
    type(dpointd) :: crnew(d+1)
    type(dpointd) :: phinew(d+1)
    real(8), allocatable, target :: curcr(:)
    real(8), allocatable, target :: slice(:), new_slice(:), matrix_slice(:)
    real(8), allocatable :: R(:)
    real(8) eps
    real(8), allocatable :: phitmp(:), Stmp(:)
    integer, allocatable :: pa(:)
    integer :: i, k, swp, dir, mm, nn, rnew, rmax2 
    integer :: max_matrix_core_size, max_phi_size, max_core_size, max_slice_size, max_R_size
    integer :: typ
    integer, optional, intent(in) :: order0
    integer :: order
    real(8) :: ermax, min_res
    real(8) anorm
    real(8) ZERO, ONE
    parameter( ZERO=0.d0, ONE=1.d0)
    typ = 2 ! 1 - KSL, 2 - KSL-symm scheme which is used by default
    if ( present(typ0) ) then
        typ = typ0 
    endif
    order = 8
    if ( present(order0) ) then
        order = order0
    endif
    min_res = 1d-1
    rmax2 = rmax 
    !Inner parameters
    eps = 1e-8 !For local solvers
    call disp('Solving a real-valued dynamical problem with tau='//tostring(tau))
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
    allocate(pa(d+1))
    call compute_ps(d, ra, n(1:d), pa)
    !Find memory for the temporary arrays
    max_phi_size = 0
    max_core_size = 0
    max_slice_size = 0
    max_matrix_core_size = 0
    max_R_size = 0
    do i = 1, d
        if (ry(i)*n(i)*ry(i+1) > max_core_size) then
            max_core_size = ry(i)*n(i)*ry(i+1)
        end if
        if (ry(i)*ry(i) > max_R_size) then
            max_R_size = ry(i)*ry(i)
        end if
        if (ra(i)*ra(i+1) > max_matrix_core_size) then
            max_matrix_core_size = ra(i)*ra(i+1)
        end if
        if (ry(i)*ry(i+1) > max_slice_size) then
            max_slice_size = ry(i)*ry(i+1)
        end if
        if (ry(i)*ry(i)*ra(i) > max_phi_size) then
            max_phi_size = ry(i)*ry(i)*ra(i)
        end if
    end do
    allocate(curcr(max_core_size))
    allocate(slice(max_slice_size), new_slice(max_slice_size), matrix_slice(max_matrix_core_size), R(max_R_size))
    allocate(Stmp(max_R_size), phitmp(max_phi_size))
    mm = 1
    do i=1,d
       allocate(crnew(i)%p(ry(i)*n(i)*ry(i+1)*2))
       call dcopy(ry(i)*n(i)*ry(i+1), crY0(mm), 1, crnew(i)%p, 1)
       mm = mm + ry(i)*n(i)*ry(i+1)
    end do
    allocate(phinew(1)%p(1))
    allocate(phinew(d+1)%p(1))
    phinew(1)%p(1) = ONE
    phinew(d+1)%p(1) = ONE
    !   QR, psi
    dir = 1
    i = 1
    do while (i < d)
       rnew = min(ry(i)*n(i), ry(i+1))
       call dqr(ry(i)*n(i), ry(i+1), crnew(i)%p, R)
       if ( i < d ) then
          call dgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), ONE, R, rnew, crnew(i+1)%p, ry(i+1), ZERO, curcr, rnew)
          call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p,1)
          ry(i+1) = rnew;
          !     Phir
          !phinew(i+1) is ry(i)*n(i)*ry(i+1)
          allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)))
          call dphi_diag_left(ry(i), n(i), ry(i+1), ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, phinew(i+1)%p)
       end if
       i = i + dir
    end do
    ermax = 0d0
    i = d
    dir = -1
    swp = 0
    call init_seed()
    ermax = 0d0
    if ( typ .eq. 1) then
        tau0 = tau
    else if ( typ .eq. 2  ) then
        tau0 = tau/2
    endif
    !The symmetrized scheme consists in: U1 - S1 - U2 - S2 - U3
    !Step in U3, S2, U2, S1, U1;
    !Backwards would be: U1 (again) S1 U2 S2 U3, i.e. - step in core; step in S 
    do while (swp < 1)
       !True iteration when started from the left:
       !move (US), move S, next core
       !and backwards move S, move (US), prev. core
       if ( dir < 0 ) then
           call init_bfun_sizes(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
               ra(i), ra(i+1), ry(i)*n(i)*ry(i+1), ry(i)*n(i)*ry(i+1))
           do k = 1, n(i) !Take block diagonal part and solve it separatedbly
               call dextract_slice(ra(i), n(i), ra(i+1), crA(pa(i):pa(i+1)-1), k, matrix_slice)
               call dinit_bfun_main(phinew(i)%p, matrix_slice, phinew(i+1)%p) !Additional setup for the zmatvec subroutine 
               !anorm = znormest(ry(i)*n(i)*ry(i+1),4, zmatvec, zmatvec_transp)
               anorm = 1d0 !This is the point we need to fix later with the norm estimate 
               call dextract_slice(ry(i), n(i), ry(i+1), crnew(i)%p, k, slice)
               call dexp_mv(ry(i)*ry(i+1), order, tau0, slice, new_slice, eps, anorm, dmatvec)
               call dset_slice(ry(i), n(i), ry(i+1), curcr, k, new_slice)
           end do
           if ( i > 1 ) then
               call dtransp(ry(i), n(i)*ry(i+1), curcr)
               rnew = min(n(i)*ry(i+1), ry(i))
               call dqr(n(i)*ry(i+1), ry(i), curcr, R) 
               call dtransp(n(i) * ry(i+1), rnew, curcr)
               call dcopy(rnew * n(i) * ry(i+1), curcr, 1, crnew(i)%p, 1)
               call dtransp(rnew, ry(i), R)
               !This should be also rewritten to the diagonal A
               call dphi_diag_right(rnew, n(i), ry(i+1), ra(i), ra(i+1), phinew(i+1)%p, crA(pa(i)), curcr, phitmp)
               !phitmp is now ry(i) x ra(i) x ry(i) 
               call dinit_sfun(ry(i), ry(i), ra(i), rnew, rnew, phinew(i)%p, phitmp)
               !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
               anorm = 1d0
               call dexp_mv(ry(i)*rnew, order, -tau0, R, Stmp, eps, anorm, dsfun_matvec)
               call dgemm('n', 'n', ry(i-1) * n(i-1), rnew, ry(i), ONE, crnew(i-1)%p,&
               ry(i-1)*n(i-1), Stmp, ry(i), ZERO, curcr, ry(i-1)*n(i-1))
               ry(i) = rnew
               call dcopy(ry(i-1)*n(i-1)*ry(i), curcr, 1, crnew(i-1)%p, 1)
               call dcopy(ry(i)*ra(i)*ry(i), phitmp, 1, phinew(i)%p, 1) !Update phi  
           else !i == 1 
               call dcopy(ry(i)*n(i)*ry(i+1),curcr,1,crnew(i)%p,1) 
           end if
       else   ! dir > 0
           call init_bfun_sizes(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
           ra(i+1), ry(i)*n(i)*ry(i+1), ry(i)*n(i)*ry(i+1))
           do k = 1, n(i)
               call dextract_slice(ra(i), n(i), ra(i+1), crA(pa(i):pa(i+1)-1), k, matrix_slice)
               call dinit_bfun_main(phinew(i)%p, matrix_slice, phinew(i+1)%p) !Additional setup for the zmatvec subroutine 
               call dextract_slice(ry(i), n(i), ry(i+1), crnew(i)%p, k, slice)
               call dexp_mv(ry(i)*ry(i+1), order, tau0, slice, new_slice, eps, anorm, dmatvec)
               call dset_slice(ry(i), n(i), ry(i+1), curcr, k, new_slice)
           end do
           if ( i < d ) then ! Update S
               call dqr(ry(i)*n(i), ry(i+1), curcr, R) 
               rnew = min(ry(i)*n(i), ry(i+1))
               call dcopy(ry(i)*n(i)*rnew, curcr, 1, crnew(i)%p, 1)
               call dphi_diag_left(ry(i), n(i), rnew, ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), curcr, phitmp)
               call dinit_sfun(rnew, rnew, ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
               !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
               anorm = 1d0
               call dexp_mv(rnew * ry(i+1), order, -tau0, R, Stmp, eps, anorm, dsfun_matvec)
               call dgemm('n', 'n', rnew, n(i+1) * ry(i+2), ry(i+1), ONE, &
               Stmp, rnew, crnew(i+1)%p, ry(i+1), &
               ZERO, curcr, rnew)
               call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p, 1) 
               ry(i+1) = rnew
               call dcopy(ry(i+1)*ra(i+1)*ry(i+1), phitmp, 1, phinew(i+1)%p, 1)
           else  ! Last core
               call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, crnew(i)%p, 1 ) 
           endif
       end if
       if ((dir>0) .and. (i==d )) then
           dir = -1
           i = d
           swp = swp + 1
           ermax = 0d0
       else if ((dir < 0) .and. (i == 1 )) then
           dir = 1
           i = 1
           if ( typ .eq. 1 ) then
               goto 105
           endif
       else
           i = i + dir
       end if
    end do
105 continue
    print *, 'Computation done'
    nn = sum(ry(2:d+1)*ry(1:d)*n(1:d))
    if ( allocated(dresult_core)) then
       if ( size(dresult_core) < nn ) then
          deallocate(dresult_core)
       end if
    end if
    if ( .not. allocated(dresult_core) ) then
       allocate(dresult_core(nn))
    end if
    nn = 1
    do i = 1, d
       call dcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, dresult_core(nn), 1)
       nn = nn + ry(i)*n(i)*ry(i+1)
    end do
    do i = 1, d
       if ( associated(crnew(i)%p)) then
          deallocate(crnew(i)%p)
       end if
    end do
    do i = 1, d+1
       if ( associated(phinew(i)%p)) then
          deallocate(phinew(i)%p)
       end if
    end do
    deallocate(R)
    deallocate(curcr)
    deallocate(slice)
    deallocate(new_slice)
    deallocate(matrix_slice)
    deallocate(pa)
  end subroutine dtt_diag_ksl

  subroutine ztt_diag_ksl(d, n, ra, crA, crY0, ry, tau, rmax, kickrank, nswp, verb, typ0, order0)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_diag_ksl !bfun declarations
    use sfun_diag_ksl !We have to do this
    use explib !Krylov exponential
    use rnd_lib
    implicit none
    integer,intent(in) :: d, n(d), ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb, typ0
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    complex(8), intent(in) :: crA(*), crY0(*)
    real(8), intent(in) :: tau
    real(8) :: tau0
    type(zpointd) :: crnew(d+1)
    type(zpointd) :: phinew(d+1)
    complex(8), allocatable, target :: curcr(:)
    complex(8), allocatable, target :: slice(:), new_slice(:), matrix_slice(:)
    complex(8), allocatable :: R(:)
    real(8) eps
    complex(8), allocatable :: phitmp(:), Stmp(:)
    integer, allocatable :: pa(:)
    integer :: i, k, swp, dir, mm, nn, rnew, rmax2
    integer :: max_matrix_core_size, max_phi_size, max_core_size, max_slice_size, max_R_size
    integer :: typ
    integer, optional, intent(in) :: order0
    integer :: order
    real(8) :: ermax, min_res
    real(8) anorm
    complex(8) ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    typ = 2 ! 1 - KSL, 2 - KSL-symm scheme which is used by default
    if ( present(typ0) ) then
        typ = typ0 
    endif
    order = 8
    if ( present(order0) ) then
        order = order0
    endif
    min_res = 1d-1
    rmax2 = rmax 
    !Inner parameters
    eps = 1e-8 !For local solvers
    call disp('Solving a complex-valued dynamical problem with tau='//tostring(tau))
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
    allocate(pa(d+1))
    call compute_ps(d, ra, n(1:d), pa)
    !Find memory for the temporary arrays
    max_phi_size = 0
    max_core_size = 0
    max_slice_size = 0
    max_matrix_core_size = 0
    max_R_size = 0
    do i = 1, d
        if (ry(i)*n(i)*ry(i+1) > max_core_size) then
            max_core_size = ry(i)*n(i)*ry(i+1)
        end if
        if (ry(i)*ry(i) > max_R_size) then
            max_R_size = ry(i)*ry(i)
        end if
        if (ra(i)*ra(i+1) > max_matrix_core_size) then
            max_matrix_core_size = ra(i)*ra(i+1)
        end if
        if (ry(i)*ry(i+1) > max_slice_size) then
            max_slice_size = ry(i)*ry(i+1)
        end if
        if (ry(i)*ry(i)*ra(i) > max_phi_size) then
            max_phi_size = ry(i)*ry(i)*ra(i)
        end if
    end do
    allocate(curcr(max_core_size))
    allocate(slice(max_slice_size), new_slice(max_slice_size), matrix_slice(max_matrix_core_size), R(max_R_size))
    allocate(Stmp(max_R_size), phitmp(max_phi_size))
    mm = 1
    do i = 1, d
       allocate(crnew(i)%p(ry(i)*n(i)*ry(i+1)))
       call zcopy(ry(i)*n(i)*ry(i+1), crY0(mm), 1, crnew(i)%p, 1)
       mm = mm + ry(i)*n(i)*ry(i+1)
    end do
    allocate(phinew(1)%p(1))
    allocate(phinew(d+1)%p(1))
    phinew(1)%p(1) = ONE
    phinew(d+1)%p(1) = ONE
    !   QR, psi
    dir = 1
    i = 1
    do while (i < d)
       rnew = min(ry(i)*n(i), ry(i+1))
       call zqr(ry(i)*n(i), ry(i+1), crnew(i)%p, R)
       if ( i < d ) then
          call zgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), ONE, R, rnew, crnew(i+1)%p, ry(i+1), ZERO, curcr, rnew)
          call zcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p, 1)
          ry(i+1) = rnew
          !     Phir
          !phinew(i+1) is ry(i)*n(i)*ry(i+1)
          allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)))
          call zphi_diag_left(ry(i), n(i), ry(i+1), ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, phinew(i+1)%p)
       end if
       i = i + dir
    end do
    ermax = 0d0
    i = d
    dir = -1
    swp = 0
    call init_seed()
    ermax = 0d0
    if ( typ .eq. 1) then
        tau0 = tau
    else if ( typ .eq. 2  ) then
        tau0 = tau/2
    endif
    !The symmetrized scheme consists in: U1 - S1 - U2 - S2 - U3
    !Step in U3, S2, U2, S1, U1;
    !Backwards would be: U1 (again) S1 U2 S2 U3, i.e. - step in core; step in S 
    do while (swp < 1)
       !True iteration when started from the left:
       !move (US), move S, next core
       !and backwards move S, move (US), prev. core
       if ( dir < 0 ) then
           call init_bfun_sizes(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
               ra(i), ra(i+1), ry(i)*n(i)*ry(i+1), ry(i)*n(i)*ry(i+1))
           do k = 1, n(i) !Take block diagonal part and solve it separatedbly
               call zextract_slice(ra(i), n(i), ra(i+1), crA(pa(i):pa(i+1)-1), k, matrix_slice)
               call zinit_bfun_main(phinew(i)%p, matrix_slice, phinew(i+1)%p) !Additional setup for the zmatvec subroutine 
               !anorm = znormest(ry(i)*n(i)*ry(i+1),4, zmatvec, zmatvec_transp)
               anorm = 1d0 !This is the point we need to fix later with the norm estimate 
               call zextract_slice(ry(i), n(i), ry(i+1), crnew(i)%p, k, slice)
               call zexp_mv(ry(i)*ry(i+1), min(order, ry(i)*ry(i+1)), tau0, slice, new_slice, eps, anorm, zmatvec)
               call zset_slice(ry(i), n(i), ry(i+1), curcr, k, new_slice)
           end do
           if ( i > 1 ) then
               call ztransp(ry(i), n(i)*ry(i+1), curcr)
               rnew = min(n(i)*ry(i+1), ry(i))
               call zqr(n(i)*ry(i+1), ry(i), curcr, R) 
               call ztransp(n(i) * ry(i+1), rnew, curcr)
               call zcopy(rnew*n(i)*ry(i+1), curcr, 1, crnew(i)%p, 1)
               call ztransp(rnew, ry(i), R)
               !This should be also rewritten to the diagonal A
               call zphi_diag_right(rnew, n(i), ry(i+1), ra(i), ra(i+1), phinew(i+1)%p, crA(pa(i)), curcr, phitmp)
               !phitmp is now ry(i) x ra(i) x ry(i) 
               call zinit_sfun(ry(i), ry(i), ra(i), rnew, rnew, phinew(i)%p, phitmp)
               !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
               anorm = 1d0
               call zexp_mv(ry(i)*rnew, order, -tau0, R, Stmp, eps, anorm, zsfun_matvec)
               call zgemm('n', 'n', ry(i-1) * n(i-1), rnew, ry(i), ONE, crnew(i-1)%p,&
               ry(i-1)*n(i-1), Stmp, ry(i), ZERO, curcr, ry(i-1)*n(i-1))
               ry(i) = rnew
               call zcopy(ry(i-1)*n(i-1)*ry(i), curcr, 1, crnew(i-1)%p, 1)
               call zcopy(ry(i)*ra(i)*ry(i), phitmp, 1, phinew(i)%p, 1) !Update phi  
           else !i == 1 
               call zcopy(ry(i)*n(i)*ry(i+1), curcr, 1, crnew(i)%p, 1) 
           end if
       else   ! dir > 0
           call init_bfun_sizes(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
           ra(i+1), ry(i)*n(i)*ry(i+1), ry(i)*n(i)*ry(i+1))
           do k = 1, n(i)
               call zextract_slice(ra(i), n(i), ra(i+1), crA(pa(i):pa(i+1)-1), k, matrix_slice)
               call zinit_bfun_main(phinew(i)%p, matrix_slice, phinew(i+1)%p) !Additional setup for the zmatvec subroutine 
               call zextract_slice(ry(i), n(i), ry(i+1), crnew(i)%p, k, slice)
               call zexp_mv(ry(i)*ry(i+1), order, tau0, slice, new_slice, eps, anorm, zmatvec)
               call zset_slice(ry(i), n(i), ry(i+1), curcr, k, new_slice)
           end do
           if ( i < d ) then ! Update S
               call zqr(ry(i)*n(i), ry(i+1), curcr, R) 
               rnew = min(ry(i)*n(i), ry(i+1))
               call zcopy(ry(i)*n(i)*rnew, curcr, 1, crnew(i)%p, 1)
               call zphi_diag_left(ry(i), n(i), rnew, ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), curcr, phitmp)
               call zinit_sfun(rnew, rnew, ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
               !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
               anorm = 1d0
               call zexp_mv(rnew*ry(i+1), order, -tau0, R, Stmp, eps, anorm, zsfun_matvec)
               call zgemm('n', 'n', rnew, n(i+1)*ry(i+2), ry(i+1), ONE, &
               Stmp, rnew, crnew(i+1)%p, ry(i+1), &
               ZERO, curcr, rnew)
               call zcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p, 1) 
               ry(i+1) = rnew
               call zcopy(ry(i+1) * ra(i+1)*ry(i+1), phitmp, 1, phinew(i+1)%p, 1)
           else  ! Last core
               call zcopy(ry(i)*n(i)*ry(i+1), curcr, 1, crnew(i)%p, 1) 
           endif
       end if
       if ((dir>0) .and. (i==d )) then
           dir = -1
           i = d
           swp = swp + 1
           ermax = 0d0
       else if ((dir < 0) .and. (i == 1 )) then
           dir = 1
           i = 1
           if ( typ .eq. 1 ) then
               goto 105
           endif
       else
           i = i + dir
       end if
    end do !While swp
105 continue
    print *, 'Computation done'
    nn = sum(ry(2:d+1)*ry(1:d)*n(1:d))
    if ( allocated(zresult_core)) then
       if ( size(zresult_core) < nn ) then
          deallocate(zresult_core)
       end if
    end if
    if ( .not. allocated(zresult_core) ) then
       allocate(zresult_core(nn))
    end if
    nn = 1
    do i = 1, d
       call zcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, zresult_core(nn), 1)
       nn = nn + ry(i)*n(i)*ry(i+1)
    end do
    do i = 1, d
       if ( associated(crnew(i)%p)) then
          deallocate(crnew(i)%p)
       end if
    end do
    do i = 1, d+1
       if ( associated(phinew(i)%p)) then
          deallocate(phinew(i)%p)
       end if
    end do
    deallocate(R)
    deallocate(curcr)
    deallocate(slice)
    deallocate(new_slice)
    deallocate(matrix_slice)
    deallocate(pa)
  end subroutine ztt_diag_ksl
end module dyn_diag_tt
