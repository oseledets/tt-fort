module bfun_ksl
  !Here we need to have all parameters required for the matrix-by-vector product (to call bfun3)
  integer, private  ::  rx1T,mT,rx2T,ry1T,nT,ry2T,ra1T,ra2T
  real(8), pointer, private :: phi1T(:), phi2T(:),res1T(:), res2T(:), AT(:)
  complex(8), pointer, private :: zphi1T(:), zphi2T(:), zres1T(:), zres2T(:), zAT(:)
  integer,private ::  xsizeT, ysizeT
  type, public ::  pointd
     real(8), dimension(:), pointer :: p=>null()
  end type pointd

  type, public :: zpointd
     complex(8), dimension(:), pointer :: p=>null()
  end type zpointd

contains 

  subroutine dmatvec(x,y)
    use ttals, only: dbfun3
    implicit none
    real(8), intent(in) :: x(*)
    real(8) :: y(*)
    call dbfun3(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
  end subroutine dmatvec

  subroutine zmatvec(x,y)
    use ttals, only: zbfun3
    implicit none
    complex(8), intent(in) :: x(*)
    complex(8) :: y(*)
    !call zcopy(rx1T*mT*rx2T, x, 1, y, 1)
    call zbfun3(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, zphi1T, zAT, zphi2T, x, y)
 end subroutine zmatvec


  subroutine dmatvec_transp(x,y)
    use ttals, only: dbfun3_transp
    implicit none
    real(8), intent(in) :: x(*)
    real(8) :: y(*)
    call dmatvec(x,y)
    !call bfun3_transp(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
  end subroutine dmatvec_transp

  subroutine zmatvec_transp(x,y)
    use ttals, only: zbfun3_transp
    implicit none
    complex(8), intent(in) :: x(*)
    complex(8) :: y(*)
    call zbfun3_transp(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, zphi1T, zAT, zphi2T, x, y)
  end subroutine zmatvec_transp


  subroutine init_bfun_sizes(rx1,m,rx2,ry1,n,ry2,ra1,ra2,xsize,ysize)
    implicit none
    integer, intent(in) :: rx1,m,rx2,ry1,n,ry2,ra1,ra2,xsize,ysize
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


  subroutine dinit_bfun_main(phi1,A,phi2)
    real(8), target :: phi1(:), phi2(:), A(:)
    phi1T => phi1
    phi2T => phi2
    AT => A
  end subroutine dinit_bfun_main

  subroutine zinit_bfun_main(phi1,A,phi2)
    complex(8), target :: phi1(:), phi2(:), A(:)
    zphi1T => phi1
    zphi2T => phi2
    zAT => A
  end subroutine zinit_bfun_main


end module bfun_ksl

module sfun_ksl

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


  subroutine dsfun_matvec(Sx,Sy)
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

  subroutine dsfun_matvec_transp(Sy,Sx)
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

  subroutine zsfun_matvec(Sx,Sy)
    use matrix_util
    complex(8), intent(in) :: Sx(*)
    complex(8), intent(out) :: Sy(*)
    complex(8) :: res1(rx1,ra,ry2)
    complex(8) ::  ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
    !phi1(ry1,rx1,ra)*phi2(rx2,ra,ry2)*S(rx1,rx2); S(rx2,rx1)
    !S(rx1,rx2)*phi2(rx2,ra,ry2) = res(rx1,ra,ry2)*phi1(ry1,rx1,ra) 
    !call dtransp(rx2,rx1,Sx)
    call zgemm('n','n',rx1,ra*ry2,rx2,ONE,Sx,rx1,zphi2,rx2,ZERO, res1, rx1)
    call zgemm('n','n',ry1,ry2,rx1*ra,ONE,zphi1,ry1,res1,rx1*ra,ZERO,Sy,ry1)
    !call dtransp(rx1,rx2,Sx)
    !call dtransp(ry1,ry2,Sy)
  end subroutine zsfun_matvec

  subroutine zsfun_matvec_transp(Sy,Sx)
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

end module sfun_ksl


module dyn_tt
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! KLS-scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! What we have: we have a starting vector + a matrix (no vector X!) 
  subroutine tt_ksl(d,n,m,ra,crA, crY0, ry, tau, rmax, kickrank, nswp, verb, typ0, order0)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_ksl !bfun declarations
    use sfun_ksl !We have to do this
    use explib !Krylov exponential
    implicit none
    integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    real(8), intent(in) :: crA(*), crY0(*)
    real(8), intent(in) :: tau
    type(pointd) :: crnew(d+1)
    type(pointd) :: phinew(d+1)
    real(8), allocatable, target :: curcr(:)
    real(8), allocatable :: R(:)
    real(8) eps
    real(8), allocatable :: phitmp(:), Stmp(:)
    integer, allocatable :: pa(:)
    integer :: i, swp, dir, mm, nn, rnew, rmax2 
    integer :: max_phi_size, max_core_size, max_R_size
    integer :: typ, order
    integer, optional :: typ0, order0
    real(8) :: ermax, min_res
    real(8) anorm
    typ = 2 ! 1 - KSL, 2 - KSL-symm
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
    call compute_ps(d,ra,n(1:d)*m(1:d),pa)
    !Find memory for the temporary arrays
    max_phi_size = 0
    max_core_size = 0
    max_R_size = 0
    do i = 1, d
        if (ry(i)*n(i)*ry(i+1) > max_core_size) then
            max_core_size = ry(i)*n(i)*ry(i+1)
        end if
        if (ry(i)*ry(i) > max_R_size) then
            max_R_size = ry(i)*ry(i)
        end if
        if (ry(i)*ry(i)*ra(i) > max_phi_size) then
            max_phi_size = ry(i)*ry(i)*ra(i)
        end if
    end do
    allocate(curcr(max_core_size))
    allocate(R(max_R_size))
    allocate(Stmp(max_R_size), phitmp(max_phi_size))
    mm = 1
    do i=1,d
       allocate(crnew(i)%p(ry(i)*n(i)*ry(i+1)*2))
       call dcopy(ry(i)*n(i)*ry(i+1), crY0(mm), 1, crnew(i)%p, 1)
       mm = mm + ry(i)*n(i)*ry(i+1)
    end do
    allocate(phinew(1)%p(1))
    allocate(phinew(d+1)%p(1))
    phinew(1)%p(1) = 1d0
    phinew(d+1)%p(1) = 1d0
    !   QR, psi
    dir = 1
    i = 1
    do while (i < d)
       rnew = min(ry(i)*n(i), ry(i+1))
       call dqr(ry(i)*n(i), ry(i+1), crnew(i) % p, R)
       if ( i < d ) then
          call dgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, R, rnew, crnew(i+1)%p, ry(i+1), 0d0, curcr, rnew)
          call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p,1)
          ry(i+1) = rnew;
          !     Phir
          !phinew(i+1) is ry(i)*n(i)*ry(i+1)
          allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
          call dphi_left(ry(i), m(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i),&
               ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, &
               phinew(i+1)%p)
       end if
       i = i+dir
    end do
    i = d
    ermax = 0d0
    i = d
    dir = -1
    swp = 1
    call init_seed()
    ermax = 0d0
    swp = 1
    do while (swp .eq. 1)
       !True iteration when started from the left:
       !move (US), move S, next core
       !and backwards move S, move (US), prev. core
       if ( dir < 0 ) then
          call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
          call dinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)
          anorm = normest(ry(i)*n(i)*ry(i+1),4, dmatvec, dmatvec_transp)
          call dexp_mv(ry(i)*n(i)*ry(i+1),order,tau/2,crnew(i)%p,curcr,eps,anorm,dmatvec)
          if ( i < d ) then
             !In this case, we have to put S(i+1) backwards in time (heh)
             call dqr(n(i)*ry(i), ry(i+1), curcr, R) 

             call dphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), & 
             ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), curcr, curcr, phitmp)
             call dinit_sfun(ry(i+1), ry(i+1), ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
             anorm = normest(ry(i+1)*ry(i+1), 4, dsfun_matvec, dsfun_matvec_transp)
             call dexp_mv(ry(i+1)*ry(i+1), order, -tau/2, R, Stmp, eps, anorm, dsfun_matvec)
             call dgemm('n','n',ry(i)*n(i),ry(i+1),ry(i+1),1d0,curcr, ry(i)*n(i), Stmp, ry(i+1), 0d0, crnew(i)%p, ry(i)*n(i))
             call dcopy(ry(i)*n(i)*ry(i+1),crnew(i)%p,1,curcr,1)
          end if

          if ( i > 1 ) then
             call dtransp(ry(i),n(i)*ry(i+1), curcr)
             call dqr(n(i)*ry(i+1), ry(i), curcr, R)
             call dtransp(n(i)*ry(i+1), ry(i), curcr, crnew(i)%p)
             call dgemm('n','t',ry(i-1)*n(i),ry(i), ry(i), 1d0, crnew(i-1)%p, ry(i-1)*n(i), R, ry(i), 0d0, curcr, ry(i-1)*n(i))
             call dcopy(ry(i-1)*n(i)*ry(i), curcr, 1, crnew(i-1)%p, 1)
             !And recompute phi
             if ( size(phinew(i)%p) < ry(i)*ry(i)*ra(i) ) then
                deallocate(phinew(i)%p)
                allocate(phinew(i)%p(ry(i)*ry(i)*ra(i)*2))
             end if
             call dphi_right(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
             ra(i+1), phinew(i+1)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phinew(i)%p)
          end if
       end if
       if ( dir > 0 ) then
          !We will rewrite to (an ordinary) LR-iteration that starts from the first core, but it seems to be not the second order

          !dqr ->
          !step S
          !update core
          !dqr ->
          !mul 
          if ( i < d ) then
             call dqr(ry(i)*n(i),ry(i+1),crnew(i)%p,R)
             !The size of the updated s would be ry(i+1) -> 
             call dphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
             ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phitmp)
             call dinit_sfun(ry(i+1), ry(i+1), ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
             anorm = normest(ry(i+1)*ry(i+1),4,dsfun_matvec,dsfun_matvec_transp)
             call dexp_mv(ry(i+1)*ry(i+1), order, -tau/2, R, Stmp, eps, anorm, dsfun_matvec)
             call dgemm('n','n',ry(i)*n(i),ry(i+1),ry(i+1),1d0,crnew(i)%p,ry(i)*n(i),Stmp,ry(i+1),0d0,curcr,ry(i)*n(i))
          else
             call dcopy(ry(i)*n(i)*ry(i+1),crnew(i)%p,1,curcr,1)
          end if
          call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
          call dinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)
          anorm = normest(ry(i)*n(i)*ry(i+1),4, dmatvec, dmatvec_transp)
          call dexp_mv(ry(i)*n(i)*ry(i+1), order, tau/2, curcr, crnew(i)%p, eps, anorm, dmatvec)
          if ( i < d ) then
             call dqr(ry(i)*n(i),ry(i+1),crnew(i)%p,R)
             call dgemm('n','n',ry(i+1),n(i+1)*ry(i+2),ry(i+1),1d0,R,ry(i+1),crnew(i+1)%p,ry(i+1),0d0,curcr,ry(i+1))
             call dcopy(ry(i+1)*n(i+1)*ry(i+2),curcr,1,crnew(i+1)%p,1)
             if ( size(phinew(i+1)%p) < ry(i+1)*ry(i+1)*ra(i+1) ) then
                deallocate(phinew(i+1)%p)
                allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
             end if
             call dphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
             ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phinew(i+1)%p)
          end if
       end if
       if ((dir>0) .and. (i==d )) then
          dir = -1
          i = d
          !call disp('swp: '//tostring(1d0*swp)//' er = '//tostring(ermax)//' rmax:'//tostring(1d0*maxval(ry(1:d))))
          swp = swp + 1
          ermax = 0d0
       else if ((dir < 0) .and. (i == 1 )) then
          dir = 1
          i = 1
          !goto 100
       else
          i = i + dir
       end if
    end do
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
    do i=1,d
       call dcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, dresult_core(nn), 1)
       nn = nn+ry(i)*n(i)*ry(i+1)
    end do
    do i = 1,d
       if ( associated(crnew(i)%p)) then
          deallocate(crnew(i)%p)
       end if
    end do
    do i = 1,d+1
       if ( associated(phinew(i)%p)) then
          deallocate(phinew(i)%p)
       end if
    end do
    deallocate(R)
    deallocate(curcr)
    deallocate(pa)


  end subroutine tt_ksl

  subroutine ztt_ksl(d,n,m,ra,crA, crY0, ry, tau, rmax, kickrank, nswp, verb, typ0, order0)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_ksl !bfun declarations
    use sfun_ksl !We have to do this
    use explib !Krylov exponential
    use rnd_lib
    implicit none
    integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb, typ0
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    complex(8), intent(in) :: crA(*), crY0(*)
    real(8), intent(in) :: tau
    real(8) :: tau0
    type(zpointd) :: crnew(d+1)
    type(zpointd) :: phinew(d+1)
    complex(8),allocatable, target :: curcr(:)
    complex(8),allocatable :: R(:)
    real(8) eps
    complex(8), allocatable :: phitmp(:), Stmp(:)
    integer,allocatable :: pa(:)
    integer :: i, swp, dir, mm, nn, rnew, rmax2 
    integer :: typ
    integer, optional, intent(in) :: order0
    integer :: max_phi_size, max_core_size, max_R_size
    integer :: order
    real(8) :: ermax, min_res
    real(8) anorm
    complex(8) ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

    typ = 2 ! 1 - KSL, 2 - KSL-symm
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
    call compute_ps(d,ra,n(1:d)*m(1:d),pa)
    max_phi_size = 0
    max_core_size = 0
    max_R_size = 0
    do i = 1, d
        if (ry(i)*n(i)*ry(i+1) > max_core_size) then
            max_core_size = ry(i)*n(i)*ry(i+1)
        end if
        if (ry(i)*ry(i) > max_R_size) then
            max_R_size = ry(i)*ry(i)
        end if
        if (ry(i)*ry(i)*ra(i) > max_phi_size) then
            max_phi_size = ry(i)*ry(i)*ra(i)
        end if
    end do
    allocate(curcr(max_core_size))
    allocate(R(max_R_size))
    allocate(Stmp(max_R_size), phitmp(max_phi_size))
    mm = 1
    do i=1,d
       allocate(crnew(i)%p(ry(i)*n(i)*ry(i+1)*2))
       call zcopy(ry(i)*n(i)*ry(i+1), crY0(mm), 1, crnew(i)%p, 1)
       mm = mm + ry(i)*n(i)*ry(i+1)
    end do
    allocate(phinew(1)%p(1))
    allocate(phinew(d+1)%p(1))
    phinew(1) % p(1) = ONE
    phinew(d+1)%p(1) = ONE
    !   QR, psi
    dir = 1
    i = 1
    do while (i < d)
       rnew = min(ry(i)*n(i), ry(i+1))
       call zqr(ry(i)*n(i), ry(i+1), crnew(i) % p, R)
       if ( i < d ) then
          call zgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), ONE, R, rnew, crnew(i+1)%p, ry(i+1), ZERO, curcr, rnew)
          call zcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p, 1)
          ry(i+1) = rnew
          !     Phir
          !phinew(i+1) is ry(i)*n(i)*ry(i+1)
          allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
          call zphi_left(ry(i), m(i), ry(i+1), ry(i), n(i), ry(i+1), &
                        ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, &
               phinew(i+1)%p)
       end if
       i = i + dir
    end do
    ermax = 0d0
    i = d
    dir = -1
    swp = 0
    call init_seed()
    ermax = 0d0
    if (typ .eq. 1) then
        tau0 = tau
    else if (typ .eq. 2) then
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
           call init_bfun_sizes(ry(i), n(i), ry(i + 1), ry(i), n(i), ry(i + 1), &
               ra(i), ra(i + 1), ry(i) * n(i) * ry(i + 1), ry(i) * n(i) * ry(i + 1))
           call zinit_bfun_main(phinew(i) % p, crA(pa(i):pa(i+1)-1), phinew(i+1) % p)
           !anorm = znormest(ry(i)*n(i)*ry(i+1),4, zmatvec, zmatvec_transp)
           anorm = 1d0
           call zexp_mv(ry(i) * n(i) * ry(i+1), order, tau0, &
                       crnew(i)%p, curcr, eps, anorm, zmatvec)
           if ( i > 1 ) then
               call ztransp(ry(i), n(i)*ry(i+1), curcr)
               rnew = min(n(i)*ry(i+1), ry(i))
               call zqr(n(i)*ry(i+1), ry(i), curcr, R) 
               call ztransp(n(i)*ry(i+1), rnew, curcr)
               call zcopy(rnew*n(i)*ry(i+1), curcr, 1, crnew(i)%p, 1)
               call ztransp(rnew, ry(i), R)
               call zphi_right(rnew, n(i), ry(i+1), rnew, n(i), ry(i+1), &
               ra(i), ra(i+1), phinew(i+1)%p, crA(pa(i)), curcr, curcr, phitmp)
               !phitmp is now ry(i) x ra(i) x ry(i) 
               call zinit_sfun(ry(i), ry(i), ra(i), rnew, rnew, phinew(i)%p, phitmp)
               !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
               anorm = 1d0
               call zexp_mv(ry(i)*rnew, order, -tau0, R, Stmp, eps, anorm, zsfun_matvec)
               call zgemm('n', 'n', ry(i-1) * n(i-1), rnew, ry(i), ONE, crnew(i-1)%p,&
               ry(i-1) * n(i-1), Stmp, ry(i), ZERO, curcr, ry(i-1) * n(i-1))
               ry(i) = rnew
               call zcopy(ry(i-1)*n(i-1)*ry(i), curcr, 1, crnew(i-1)%p, 1)
               call zcopy(ry(i)*ra(i)*ry(i),phitmp,1,phinew(i)%p,1) !Update phi 
           else !i == 1 
               call zcopy(ry(i)*n(i)*ry(i+1), curcr, 1, crnew(i)%p, 1) 
           end if
       else   ! dir > 0
           call init_bfun_sizes(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
           ra(i+1), ry(i)*n(i)*ry(i+1), ry(i)*n(i)*ry(i+1))
           call zinit_bfun_main(phinew(i)%p, crA(pa(i):pa(i+1)-1), phinew(i+1)%p)

           !anorm = znormest(ry(i) * n(i) * ry(i+1),4, zmatvec, zmatvec_transp)
           anorm = 1d0
           call zexp_mv(ry(i) * n(i) * ry(i+1), order, tau0, &
                       crnew(i) % p, curcr, eps, anorm, zmatvec)
           if ( i < d ) then ! Update S
               call zqr(ry(i) * n(i), ry(i+1), curcr, R) 
               rnew = min(ry(i) * n(i), ry(i+1))
               call zcopy(ry(i) * n(i) * rnew, curcr, 1, crnew(i)%p, 1)
               call zphi_left(ry(i), n(i), rnew, ry(i), n(i), rnew, &
               ra(i), ra(i+1), phinew(i) % p, crA(pa(i)), curcr, curcr, phitmp)
               call zinit_sfun(rnew, rnew, ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
               !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
               anorm = 1d0
               call zexp_mv(rnew * ry(i+1), order, -tau0, R, Stmp, eps, anorm, zsfun_matvec)
               call zgemm('n', 'n', rnew, n(i+1) * ry(i+2), ry(i+1), ONE, &
               Stmp, rnew, crnew(i+1) % p, ry(i+1), &
               ZERO, curcr, rnew)
               call zcopy(rnew * n(i+1) * ry(i+2), curcr, 1, crnew(i+1) % p, 1) 
               ry(i+1) = rnew
               call zcopy(ry(i+1) * ra(i+1) * ry(i+1), phitmp, 1, phinew(i+1) % p, 1)
           else  ! Last core
               call zcopy(ry(i) * n(i) * ry(i+1), curcr, 1, crnew(i) % p, 1 ) 
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
    do i = 1,d
       if ( associated(crnew(i)%p)) then
          deallocate(crnew(i)%p)
       end if
    end do
    do i = 1,d+1
       if ( associated(phinew(i)%p)) then
          deallocate(phinew(i)%p)
       end if
    end do
    deallocate(R)
    deallocate(curcr)
    deallocate(pa)
  end subroutine ztt_ksl
end module dyn_tt
