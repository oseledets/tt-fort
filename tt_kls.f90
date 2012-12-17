module bfun_kls
  !Here we need to have all parameters required for the matrix-by-vector product (to call bfun3)
  integer, private  ::  rx1T,mT,rx2T,ry1T,nT,ry2T,ra1T,ra2T
  double precision, pointer, private :: phi1T(:), phi2T(:),res1T(:), res2T(:), AT(:)
  complex(8), pointer, private :: zphi1T(:), zphi2T(:), zres1T(:), zres2T(:), zAT(:)
  integer,private ::  xsizeT, ysizeT
  type, public ::  pointd
     double precision, dimension(:), pointer :: p=>null()
  end type pointd

  type, public :: zpointd
     complex(8), dimension(:), pointer :: p=>null()
  end type zpointd

contains 

  subroutine dmatvec(x,y)
    use ttals, only: dbfun3
    implicit none
    double precision, intent(in) :: x(*)
    double precision :: y(*)
    call dbfun3(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
  end subroutine dmatvec

  subroutine zmatvec(x,y)
    use ttals, only: zbfun3
    implicit none
    complex(8), intent(in) :: x(*)
    complex(8) :: y(*)
    !call zcopy(rx1T*mT*ry1T, x, 1, y, 1)
    call zbfun3(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, zphi1T, zAT, zphi2T, x, y)
 end subroutine zmatvec


  subroutine dmatvec_transp(x,y)
    use ttals, only: dbfun3_transp
    implicit none
    double precision, intent(in) :: x(*)
    double precision :: y(*)
    call dmatvec(x,y)
    !call bfun3_transp(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
  end subroutine dmatvec_transp

  subroutine zmatvec_transp(x,y)
    use ttals, only: zbfun3_transp
    implicit none
    complex(8), intent(in) :: x(*)
    complex(8) :: y(*)
    call zmatvec(x,y)
    !call bfun3_transp(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x, y)
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
    double precision, target :: phi1(:), phi2(:), A(:)
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


end module bfun_kls

module sfun_kls

  integer, private :: rx1, rx2, ra, ry1, ry2
  double precision, private, pointer :: phi1(:), phi2(:)
  complex(8), private, pointer :: zphi1(:), zphi2(:)


contains

  subroutine dinit_sfun(rx1T, rx2T, raT, ry1T, ry2T, phi1T, phi2T)
    implicit none
    integer, intent(in) ::  rx1T, rx2T, raT, ry1T, ry2T 
    double precision, intent(in), target :: phi1T(:), phi2T(:)
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
    double precision, intent(in) :: Sx(*)
    double precision, intent(out) :: Sy(*)
    double precision :: res1(rx1,ra,ry2)
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
    double precision, intent(in) :: Sy(*)
    double precision, intent(out) :: Sx(*)
    double precision res1(ry2,rx1,ra)
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

end module sfun_kls


module dyn_tt
  implicit none
  real(8), allocatable :: result_core(:)
  complex(8), allocatable :: zresult_core(:)
contains
  subroutine deallocate_result
    if ( allocated(result_core) ) then
       deallocate(result_core)
    end if
    if ( allocated(zresult_core) ) then
       deallocate(zresult_core) 
    end if
  end subroutine deallocate_result



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! KLS-scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !! What we have: we have a starting vector + a matrix (no vector X!) 
  subroutine tt_kls(d,n,m,ra,crA, crY0, ry, tau, rmax, kickrank, nswp, verb)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_kls !bfun declarations
    use sfun_kls !We have to do this
    use explib !Krylov exponential
    implicit none
    integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    real(8), intent(in) :: crA(*), crY0(*)
    double precision, intent(in) :: tau
    type(pointd) :: crnew(d+1)
    type(pointd) :: phinew(d+1)
    real(8),allocatable, target :: curcr(:)
    real(8),allocatable, target :: work(:)
    real(8),allocatable :: R(:)
    real(8),allocatable :: full_core(:)
    double precision eps
    double precision :: sv(rmax)
    real(8), allocatable :: rnorms(:), W(:,:), X(:,:), Bmat(:,:), U(:,:), phitmp(:), Stmp(:)
    integer,allocatable :: pa(:)
    integer :: i,j, k, swp, dir, lwork, mm, nn, rnew, max_matvecs, rmax2 
    double precision :: err, ermax, res, res_old, min_res
    double precision anorm
    real(8) dznrm2

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

    call compute_ps(d,ra,n(1:d)*m(1:d),pa)



    lwork = rmax*maxval(n(1:d))*rmax

    nn = maxval(ra(1:d+1))*rmax*rmax
    allocate(curcr(lwork))
    nn = maxval(ra(1:d+1))*rmax*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
    allocate(work(nn))
    allocate(R(lwork))
    allocate(full_core(nn))
    allocate(phitmp(maxval(ra(1:d+1))*maxval(ry(1:d+1))**2))
    allocate(Stmp(maxval(ry(1:d+1))**2))

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
          call exp_mv(ry(i)*n(i)*ry(i+1),30,tau/2,crnew(i)%p,curcr,eps,anorm,dmatvec)
          if ( i < d ) then
             !In this case, we have to put S(i+1) backwards in time (heh)
             call dqr(n(i)*ry(i), ry(i+1), curcr, R) 

             call dphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), & 
             ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), curcr, curcr, phitmp)
             call dinit_sfun(ry(i+1), ry(i+1), ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
             anorm = normest(ry(i+1)*ry(i+1), 4, dsfun_matvec, dsfun_matvec_transp)
             call exp_mv(ry(i+1)*ry(i+1), 30, -tau/2, R, Stmp, eps, anorm, dsfun_matvec)
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
             call exp_mv(ry(i+1)*ry(i+1), 30, -tau/2, R, Stmp, eps, anorm, dsfun_matvec)
             call dgemm('n','n',ry(i)*n(i),ry(i+1),ry(i+1),1d0,crnew(i)%p,ry(i)*n(i),Stmp,ry(i+1),0d0,curcr,ry(i)*n(i))
          else
             call dcopy(ry(i)*n(i)*ry(i+1),crnew(i)%p,1,curcr,1)
          end if

          call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
          call dinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)

          anorm = normest(ry(i)*n(i)*ry(i+1),4, dmatvec, dmatvec_transp)
          call exp_mv(ry(i)*n(i)*ry(i+1),30,tau/2,curcr,crnew(i)%p,eps,anorm,dmatvec)

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

100 continue

    nn = sum(ry(2:d+1)*ry(1:d)*n(1:d))
    if ( allocated(result_core)) then
       if ( size(result_core) < nn ) then
          deallocate(result_core)
       end if
    end if
    if ( .not. allocated(result_core) ) then
       allocate(result_core(nn))
    end if

    nn = 1
    do i=1,d
       call dcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, result_core(nn), 1)
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
    deallocate(work)
    deallocate(curcr)
    deallocate(full_core)
    deallocate(pa)


  end subroutine tt_kls

  subroutine ztt_kls(d,n,m,ra,crA, crY0, ry, tau, rmax, kickrank, nswp, verb)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_kls !bfun declarations
    use sfun_kls !We have to do this
    use explib !Krylov exponential
    implicit none
    integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    complex(8), intent(in) :: crA(*), crY0(*)
    real(8), intent(in) :: tau
    type(zpointd) :: crnew(d+1)
    type(zpointd) :: phinew(d+1)
    complex(8),allocatable, target :: curcr(:)
    complex(8),allocatable, target :: work(:)
    complex(8),allocatable :: R(:)
    complex(8),allocatable :: full_core(:)
    double precision eps
    double precision :: sv(rmax)
    complex(8), allocatable :: rnorms(:), W(:,:), X(:,:), Bmat(:,:), U(:,:), phitmp(:), Stmp(:)
    integer,allocatable :: pa(:)
    integer :: i,j, k, swp, dir, lwork, mm, nn, rnew, max_matvecs, rmax2 
    double precision :: err, ermax, res, res_old, min_res
    double precision anorm
    real(8) dznrm2
    complex(8) ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

    !Debug print
    !print *,'d=',d,'n=',n,'m=',m,'ra=',ra,'ry=',ry
    !print *,cry0(1:2*d)  
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



    lwork = rmax*maxval(n(1:d))*rmax

    nn = maxval(ra(1:d+1))*rmax*rmax
    allocate(curcr(lwork))
    nn = maxval(ra(1:d+1))*rmax*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
    allocate(work(nn))
    allocate(R(lwork))
    allocate(full_core(nn))
    allocate(phitmp(maxval(ra(1:d+1))*maxval(ry(1:d+1))**2))
    allocate(Stmp(maxval(ry(1:d+1))**2))

    mm = 1
    do i=1,d
       allocate(crnew(i)%p(ry(i)*n(i)*ry(i+1)*2))
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
       call zqr(ry(i)*n(i), ry(i+1), crnew(i) % p, R)
       if ( i < d ) then
          call zgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), ONE, R, rnew, crnew(i+1)%p, ry(i+1), ZERO, curcr, rnew)
          call zcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p,1)
          ry(i+1) = rnew;
          !     Phir

          !phinew(i+1) is ry(i)*n(i)*ry(i+1)
          allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
          call zphi_left(ry(i), m(i), ry(i+1), ry(i), n(i), ry(i+1), &
                        ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, &
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
          call zinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)
          !anorm = znormest(ry(i)*n(i)*ry(i+1),4, zmatvec, zmatvec_transp)
          anorm = 1d0
          call zexp_mv(ry(i)*n(i)*ry(i+1),30,tau/2,crnew(i)%p,curcr,eps,anorm,zmatvec)
          if ( i < d ) then
             !In this case, we have to put S(i+1) backwards in time (heh)
             call zqr(n(i)*ry(i), ry(i+1), curcr, R) 

             call zphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
                           ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), curcr, curcr, phitmp)
             call zinit_sfun(ry(i+1), ry(i+1), ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
             !anorm = znormest(ry(i+1)*ry(i+1), 4, zsfun_matvec, zsfun_matvec_transp)
             anorm = 1d0
             call zexp_mv(ry(i+1)*ry(i+1), 30, -tau/2, R, Stmp, eps, anorm, zsfun_matvec)
             call zgemm('n','n',ry(i)*n(i),ry(i+1),ry(i+1),ONE,curcr, ry(i)*n(i), Stmp, ry(i+1), ZERO, crnew(i)%p, ry(i)*n(i))
             call zcopy(ry(i)*n(i)*ry(i+1),crnew(i)%p,1,curcr,1)
          end if

          if ( i > 1 ) then
             call ztransp(ry(i),n(i)*ry(i+1), curcr)
             call zqr(n(i)*ry(i+1), ry(i), curcr, R)
             call ztransp(n(i)*ry(i+1), ry(i), curcr, crnew(i)%p)
             call zgemm('n','t',ry(i-1)*n(i),ry(i), ry(i), ONE, crnew(i-1)%p, ry(i-1)*n(i), R, ry(i), ZERO, curcr, ry(i-1)*n(i))
             call zcopy(ry(i-1)*n(i)*ry(i), curcr, 1, crnew(i-1)%p, 1)

             !And recompute phi
             if ( size(phinew(i)%p) < ry(i)*ry(i)*ra(i) ) then
                deallocate(phinew(i)%p)
                allocate(phinew(i)%p(ry(i)*ry(i)*ra(i)*2))
             end if
             call zphi_right(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
                            ra(i), ra(i+1), phinew(i+1)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phinew(i)%p)
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
             call zqr(ry(i)*n(i),ry(i+1),crnew(i)%p,R)
             !The size of the updated s would be ry(i+1) -> 
             call zphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
             ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phitmp)
             call zinit_sfun(ry(i+1), ry(i+1), ra(i+1), ry(i+1), ry(i+1), phitmp, phinew(i+1)%p)
             !anorm = znormest(ry(i+1)*ry(i+1),4,zsfun_matvec,zsfun_matvec_transp)
             anorm = 1d0
             call zexp_mv(ry(i+1)*ry(i+1), 30, -tau/2, R, Stmp, eps, anorm, zsfun_matvec)
             call zgemm('n','n',ry(i)*n(i),ry(i+1),ry(i+1),ONE,crnew(i)%p,ry(i)*n(i),Stmp,ry(i+1),ZERO,curcr,ry(i)*n(i))
          else
             call zcopy(ry(i)*n(i)*ry(i+1),crnew(i)%p,1,curcr,1)
          end if

          call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
          call zinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)

          !anorm = znormest(ry(i)*n(i)*ry(i+1),4, zmatvec, zmatvec_transp)
          anorm = 1d0
          call zexp_mv(ry(i)*n(i)*ry(i+1),30,tau/2,curcr,crnew(i)%p,eps,anorm,zmatvec)

          if ( i < d ) then
             call zqr(ry(i)*n(i),ry(i+1),crnew(i)%p,R)
             call zgemm('n','n',ry(i+1),n(i+1)*ry(i+2),ry(i+1),ONE,R,ry(i+1),crnew(i+1)%p,ry(i+1),ZERO,curcr,ry(i+1))
             call zcopy(ry(i+1)*n(i+1)*ry(i+2),curcr,1,crnew(i+1)%p,1)
             if ( size(phinew(i+1)%p) < ry(i+1)*ry(i+1)*ra(i+1) ) then
                deallocate(phinew(i+1)%p)
                allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
             end if
             call zphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
             ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phinew(i+1)%p)
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
          !goto 105
       else
          i = i + dir
       end if
    end do

105 continue

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
    do i=1,d
       call zcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, zresult_core(nn), 1)
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
    deallocate(work)
    deallocate(curcr)
    deallocate(full_core)
    deallocate(pa)


  end subroutine ztt_kls
  !! What we have: we have a starting vector + a matrix (no vector X!) 
  
  subroutine tt_ksl(d,n,m,ra,crA, crY0, ry, tau, rmax, kickrank, nswp, verb)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_kls !bfun declarations
    use sfun_kls !We have to do this
    use explib !Krylov exponential
    implicit none
    integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    real(8), intent(in) :: crA(*), crY0(*)
    double precision, intent(in) :: tau
    type(pointd) :: crnew(d+1)
    type(pointd) :: phinew(d+1)
    real(8),allocatable, target :: curcr(:)
    real(8),allocatable, target :: work(:)
    real(8),allocatable :: R(:)
    real(8),allocatable :: full_core(:)
    double precision eps
    double precision :: sv(rmax)
    real(8), allocatable :: rnorms(:), W(:,:), X(:,:), Bmat(:,:), U(:,:), phitmp(:), Stmp(:)
    integer,allocatable :: pa(:)
    integer :: i,j, k, swp, dir, lwork, mm, nn, rnew, max_matvecs, rmax2 
    double precision :: err, ermax, res, res_old, min_res
    double precision anorm
    real(8) dznrm2

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

    call compute_ps(d,ra,n(1:d)*m(1:d),pa)



    lwork = rmax*maxval(n(1:d))*rmax

    nn = maxval(ra(1:d+1))*rmax*rmax
    allocate(curcr(lwork))
    nn = maxval(ra(1:d+1))*rmax*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
    allocate(work(nn))
    allocate(R(lwork))
    allocate(full_core(nn))
    allocate(phitmp(maxval(ra(1:d+1))*maxval(ry(1:d+1))**2))
    allocate(Stmp(maxval(ry(1:d+1))**2))

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
    dir = -1
    swp = 1

    call init_seed()
    ermax = 0d0
    swp = 1
    !   Algorithm is utterly funny: do a tau step in the core,     
    do i = d,1,-1
        call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
        call dinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)
        anorm = normest(ry(i)*n(i)*ry(i+1),4, dmatvec, dmatvec_transp)
        call exp_mv(ry(i)*n(i)*ry(i+1),30,tau,crnew(i)%p,curcr,eps,anorm,dmatvec)
        !Incorrect, we have to put other guy back in time
        if ( i > 1) then
            call dtransp(ry(i),n(i)*ry(i+1), curcr)
            call dqr(n(i)*ry(i+1), ry(i), curcr, R)
            call dtransp(n(i)*ry(i+1), ry(i), curcr, crnew(i)%p)
            call dphi_right(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
            ra(i+1), phinew(i+1)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phitmp)
             !We have to orthogonalize and put it back in time
              
            call dinit_sfun(ry(i), ry(i), ra(i), ry(i), ry(i), phinew(i)%p, phitmp)
            anorm = normest(ry(i)*ry(i), 4, dsfun_matvec, dsfun_matvec_transp)
            call dsfun_matvec(R,Stmp)
            call exp_mv(ry(i)*ry(i), 30, -1.0*tau, R, Stmp, eps, anorm, dsfun_matvec)
            call dgemm('n','t',ry(i-1)*n(i),ry(i), ry(i), 1d0, crnew(i-1)%p, ry(i-1)*n(i), Stmp, ry(i), 0d0, curcr, ry(i-1)*n(i))
            call dcopy(ry(i-1)*n(i)*ry(i), curcr, 1, crnew(i-1)%p, 1)
            call dcopy(ry(i) * ra(i) * ry(i), phitmp, 1, phinew(i) % p, 1)
        else
            call dcopy(ry(i) * n(i) * ry(i+1), curcr, 1, crnew(i) % p, 1)
        end if
    end do

100 continue

    nn = sum(ry(2:d+1)*ry(1:d)*n(1:d))
    if ( allocated(result_core)) then
       if ( size(result_core) < nn ) then
          deallocate(result_core)
       end if
    end if
    if ( .not. allocated(result_core) ) then
       allocate(result_core(nn))
    end if

    nn = 1
    do i=1,d
       call dcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, result_core(nn), 1)
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
    deallocate(work)
    deallocate(curcr)
    deallocate(full_core)
    deallocate(pa)


  end subroutine tt_ksl

  subroutine ztt_ksl(d,n,m,ra,crA, crY0, ry, tau, rmax, kickrank, nswp, verb)
    use dispmodule
    use matrix_util ! dqr
    use ttals   ! als stuff
    use estnorm  !1-norm estimation
    use bfun_kls !bfun declarations
    use sfun_kls !We have to do this
    use explib !Krylov exponential
    implicit none
    integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
    integer,intent(inout) :: ry(d+1)
    integer, intent(in), optional :: kickrank, nswp, verb
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    complex(8), intent(in) :: crA(*), crY0(*)
    real(8), intent(in) :: tau
    type(zpointd) :: crnew(d+1)
    type(zpointd) :: phinew(d+1)
    complex(8),allocatable, target :: curcr(:)
    complex(8),allocatable, target :: work(:)
    complex(8),allocatable :: R(:)
    complex(8),allocatable :: full_core(:)
    double precision eps
    double precision :: sv(rmax)
    complex(8), allocatable :: rnorms(:), W(:,:), X(:,:), Bmat(:,:), U(:,:), phitmp(:), Stmp(:)
    integer,allocatable :: pa(:)
    integer :: i,j, k, swp, dir, lwork, mm, nn, rnew, max_matvecs, rmax2 
    double precision :: err, ermax, res, res_old, min_res
    double precision anorm
    real(8) dznrm2
    complex(8) ZERO, ONE
    parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

    !Debug print
    !print *,'d=',d,'n=',n,'m=',m,'ra=',ra,'ry=',ry
    !print *,cry0(1:2*d)  
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



    lwork = rmax*maxval(n(1:d))*rmax

    nn = maxval(ra(1:d+1))*rmax*rmax
    allocate(curcr(lwork))
    nn = maxval(ra(1:d+1))*rmax*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
    allocate(work(nn))
    allocate(R(lwork))
    allocate(full_core(nn))
    allocate(phitmp(maxval(ra(1:d+1))*maxval(ry(1:d+1))**2))
    allocate(Stmp(maxval(ry(1:d+1))**2))

    mm = 1
    do i=1,d
       allocate(crnew(i)%p(ry(i)*n(i)*ry(i+1)*2))
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
       call zqr(ry(i)*n(i), ry(i+1), crnew(i) % p, R)
       if ( i < d ) then
          call zgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), ONE, R, rnew, crnew(i+1)%p, ry(i+1), ZERO, curcr, rnew)
          call zcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p,1)
          ry(i+1) = rnew;
          !     Phir

          !phinew(i+1) is ry(i)*n(i)*ry(i+1)
          allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
          call zphi_left(ry(i), m(i), ry(i+1), ry(i), n(i), ry(i+1), &
                        ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, &
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
    !   Algorithm is utterly funny: do a tau step in the core,     
    do i = d,1,-1
        call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
        call zinit_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)
        anorm = normest(ry(i)*n(i)*ry(i+1),4, zmatvec, zmatvec_transp)
        call zexp_mv(ry(i)*n(i)*ry(i+1),30,tau,crnew(i)%p,curcr,eps,anorm,dmatvec)
        !Incorrect, we have to put other guy back in time
        if ( i > 1) then
            call ztransp(ry(i),n(i)*ry(i+1), curcr)
            call zqr(n(i)*ry(i+1), ry(i), curcr, R)
            call ztransp(n(i)*ry(i+1), ry(i), curcr, crnew(i)%p)
            call zphi_right(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), &
            ra(i+1), phinew(i+1)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phitmp)
             !We have to orthogonalize and put it back in time
              
            call zinit_sfun(ry(i), ry(i), ra(i), ry(i), ry(i), phinew(i)%p, phitmp)
            anorm = normest(ry(i)*ry(i), 4, zsfun_matvec, zsfun_matvec_transp)
            call zsfun_matvec(R,Stmp)
            call zexp_mv(ry(i)*ry(i), 30, -1.0*tau, R, Stmp, eps, anorm, dsfun_matvec)
            call zgemm('n','t',ry(i-1)*n(i),ry(i), ry(i), ONE, crnew(i-1)%p, ry(i-1)*n(i), Stmp, ry(i), ZERO, curcr, ry(i-1)*n(i))
            call zcopy(ry(i-1)*n(i)*ry(i), curcr, 1, crnew(i-1)%p, 1)
            call zcopy(ry(i) * ra(i) * ry(i), phitmp, 1, phinew(i) % p, 1)
        else
            call zcopy(ry(i) * n(i) * ry(i+1), curcr, 1, crnew(i) % p, 1)
        end if
    end do

105 continue

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
    do i=1,d
       call zcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, zresult_core(nn), 1)
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
    deallocate(work)
    deallocate(curcr)
    deallocate(full_core)
    deallocate(pa)


  end subroutine ztt_ksl


end module dyn_tt
