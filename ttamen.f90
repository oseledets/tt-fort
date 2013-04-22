module ttamen_lib
!
!   This file contains AMEn solver for SPD problems in higher dimensions described in the following reference
!
!   Sergey V. Dolgov, Dmitry V. Savostyanov,
!   Alternating minimal energy methods for linear systems in higher dimensions. 
!   Part I: SPD systems, http://arxiv.org/abs/1301.6068,
!   Part II: Faster algorithm and application to nonsymmetric systems, http://arxiv.org/abs/1304.1222
!
!   Use {sergey.v.dolgov, dmitry.savostyanov}@gmail.com for feedback
!
!   Other experimental/development routines for approximation, eigenproblem, etc. are also included.
!
  use tt_lib
  use ttaux_lib
  use ptype_lib
  use ttlocsolve_lib
  use ttnodeop_lib
  use default_lib
!#ifdef MEXWRITE
!  use mexwrite_lib
!#endif

  implicit none
  double precision,parameter :: resid_damp=2d0
contains

! real solver
 subroutine dtt_amen_solve(A, y, tol, x, kickrank, nswp, local_prec, local_iters, local_restart, trunc_norm, max_full_size, verb)
  implicit none
  type(dtt),intent(in) :: A, y
  double precision,intent(in) :: tol
  type(dtt),intent(inout) :: x
  integer,intent(in), optional :: kickrank, local_iters, local_restart, nswp, trunc_norm, verb, max_full_size ! trunc_norm: 0 - fro, 1 - resid
  character,intent(in), optional :: local_prec ! 'n', 'l', 'c', 'r'
  character(len=*),parameter :: subnam='dtt_amen_solve'
  integer :: kickrank0, local_iters0, local_restart0, nswp0, trunc_norm0, verb0, max_full_size0
  character :: local_prec0
  type (dtt) :: z
  type(pointd) :: phixax(0:tt_size), phixy(0:tt_size), phizax(0:tt_size), phizy(0:tt_size)
  integer :: i, j, l, d, swp, rx1, rx2, ra1, ra2, info
  integer, pointer :: n(:)
  integer, allocatable :: ipiv(:)
  double precision,allocatable :: rhs(:), resid(:), B(:)
  double precision,pointer :: jacs(:), sol(:)
  double precision :: err, err_max, res_max, real_tol
  logical :: prc, jacsused
  double precision,external :: dnrm2
  logical,external :: lsame
  target :: A,y,x,z

  kickrank0 = default(4, kickrank)
  nswp0 = default(20, nswp)
  local_prec0 = default('n',local_prec)
  local_iters0 = default(2, local_iters)
  local_restart0 = default(40, local_restart)
  trunc_norm0 = default(1, trunc_norm)
  max_full_size0 = default(50, max_full_size)
  verb0 = default(1, verb)
  prc=.not.lsame(local_prec0,'n')

  l=y%l; d=y%m; n => y%n

  real_tol = tol/dsqrt(1.d0*(d-l))
  real_tol = real_tol/resid_damp

  ! prepare projections
  allocate(phixax(l-1)%p(1), phixy(l-1)%p(1), phixax(d)%p(1), phixy(d)%p(1))
  phixax(l-1)%p(1) = 1d0; phixax(d)%p(1) = 1d0
  phixy(l-1)%p(1) = 1d0; phixy(d)%p(1) = 1d0

  if (kickrank0.gt.0) then
  ! prepare z
    z%l = l; z%m = d; z%n(l:d)=n(l:d); z%r(l-1)=1; z%r(l:d-1)=kickrank0; z%r(d)=1
    call random(z)
    allocate(phizax(l-1)%p(1), phizy(l-1)%p(1), phizax(d)%p(1), phizy(d)%p(1))
    phizax(l-1)%p(1) = 1d0; phizax(d)%p(1) = 1d0
    phizy(l-1)%p(1) = 1d0; phizy(d)%p(1) = 1d0
  end if

  ! Main cycle
  do swp=1,nswp0
    ! orthogonalization
    do i=d,l+1,-1
      if ((kickrank0.gt.0).and.(swp.gt.1)) then
        ! update Z block
        call d2d_mv(x%r(i-1),n(i),x%r(i), z%r(i-1),n(i),z%r(i), A%r(i-1),A%r(i), phizax(i-1)%p, &
                    A%u(i)%p, phizax(i)%p, x%u(i)%p, z%u(i)%p)
        allocate(rhs(z%r(i-1)*n(i)*z%r(i)), resid(z%r(i-1)*n(i)*max(z%r(i),y%r(i))))
        call dgemm('n','n', z%r(i-1), n(i)*y%r(i), y%r(i-1), 1d0, phizy(i-1)%p, z%r(i-1), y%u(i)%p, y%r(i-1), 0d0, resid, z%r(i-1))
        call dgemm('n','n', z%r(i-1)*n(i), z%r(i), y%r(i), 1d0, resid, z%r(i-1)*n(i), phizy(i)%p, y%r(i), 0d0, rhs, z%r(i-1)*n(i))
        call daxpy(z%r(i-1)*n(i)*z%r(i), -1d0, rhs, 1, z%u(i)%p, 1)
        ! now z(i) = Z'(Ax-y)
        deallocate(rhs, resid)
      end if

      call dtt_qr(x, i, dir=-1)
      ! update X projections
      call dtt_YAX(phixax(i)%p, x, x, i, -1, phixax(i-1)%p, A=A)
      call dtt_YAX(phixy(i)%p, x, y, i, -1, phixy(i-1)%p)

      if (kickrank0>0) then
        call dtt_qr(z, i, dir=-1)
        ! update Z projections
        call dtt_YAX(phizax(i)%p, z, x, i, -1, phizax(i-1)%p, A=A)
        call dtt_YAX(phizy(i)%p, z, y, i, -1, phizy(i-1)%p)
      end if
    end do ! ort

    res_max = 0d0; err_max = 0d0
    ! optimization
    do i=l,d
      ! extract sizes
      rx1 = x%r(i-1); rx2 = x%r(i)
      ra1 = A%r(i-1); ra2 = A%r(i)

      ! prepare rhs
      ! first, memory allocations
      if (kickrank0.gt.0) then
        allocate(rhs(max(rx1,z%r(i-1))*n(i)*max(rx2,z%r(i))), resid(max(rx1,z%r(i-1))*n(i)*max(rx2,y%r(i),z%r(i))))
        allocate(sol(max(rx1,z%r(i-1))*n(i)*max(rx2,z%r(i))))
      else
        allocate(rhs(rx1*n(i)*rx2), resid(rx1*n(i)*max(rx2,y%r(i))))
        allocate(sol(rx1*n(i)*rx2))
      end if

      call dgemm('n','n', rx1, n(i)*y%r(i), y%r(i-1), 1d0, phixy(i-1)%p, rx1, y%u(i)%p, y%r(i-1), 0d0, resid, rx1)
      call dgemm('n','n', rx1*n(i), rx2, y%r(i), 1d0, resid, rx1*n(i), phixy(i)%p, y%r(i), 0d0, rhs, rx1*n(i))

      ! res_prev
      call d2d_mv(rx1, n(i), rx2, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, x%u(i)%p, resid)
      call daxpy(rx1*n(i)*rx2, -1d0, rhs, 1, resid, 1)
      err = dnrm2(rx1*n(i)*rx2, resid, 1)/dnrm2(rx1*n(i)*rx2, rhs, 1)
      res_max = dmax1(res_max,err)
      jacsused = .false.

      if (err.gt.real_tol) then
        if (rx1*n(i)*rx2<max_full_size0) then ! full solution
          allocate(B(rx1*n(i)*rx2*rx1*n(i)*rx2), ipiv(rx1*n(i)*rx2),stat=info)
          if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
          call d2d_fullmat(rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, B)
          call dcopy(rx1*n(i)*rx2, resid, 1, sol, 1)
          call dgesv(rx1*n(i)*rx2, 1, B, rx1*n(i)*rx2, ipiv, sol, rx1*n(i)*rx2, info)
          if(info.ne.0)then; write(*,*)subnam,': dgesv info: ',info; stop; end if
          if (verb0>1) write(*,"(A,I0,A)"), 'amen_solve: block=', i, ', direct local solver'
          deallocate(B, ipiv)
        else ! iter. solution
          ! allocate and compute jac prec, if necessary
          if (lsame(local_prec0,'l')) then
            allocate(jacs(rx1*rx1*n(i)*n(i)*max(ra2, rx2)))
            jacsused = .true.
          else if (lsame(local_prec0,'c'))then
            allocate(jacs(rx1*n(i)*n(i)*max(ra2, rx2)))
            jacsused = .true.
          else if (lsame(local_prec0,'r'))then
            allocate(jacs(rx1*n(i)*n(i)*rx2*max(ra2, rx2)))
            jacsused = .true.
          end if
          if (prc) call d2d_jac_gen(local_prec0, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, jacs)

          ! Run the solver
          call d2d_gmresr(phixax(i-1)%p, A%u(i)%p, phixax(i)%p, resid, rx1, n(i), rx2, ra1, ra2, &
                          local_restart0, real_tol/err, local_iters0, local_prec0, jacs, sol, max(verb0-1,0))

          if (prc) call d2d_jac_apply(local_prec0, rx1, n(i), rx2, jacs, sol, sol)
        end if
        ! Now, sol is the -correction, since the rhs was resid = Ax-y
        call daxpy(rx1*n(i)*rx2, -1d0, sol, 1, x%u(i)%p, 1)
        ! the updated core is ready
        err = dnrm2(rx1*n(i)*rx2, sol, 1)/dnrm2(rx1*n(i)*rx2, x%u(i)%p, 1)
        err_max = dmax1(err_max,err)
        ! compute res_new
        call d2d_mv(rx1, n(i), rx2, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, x%u(i)%p, resid)
        call daxpy(rx1*n(i)*rx2, -1d0, rhs, 1, resid, 1)
        err = dnrm2(rx1*n(i)*rx2, resid, 1)/dnrm2(rx1*n(i)*rx2, rhs, 1)
      end if


      if (i.lt.d) then
        if (kickrank0.gt.0) then
          if (trunc_norm0.eq.1) then ! truncation in the A-norm of error (i.e. residual)
            err = dmax1(err,real_tol*resid_damp)
            call dtt_trunc(x, i, err, 1, A=A, phi=phixax, rhs=rhs, xnew=sol)
          else  ! truncation in the frobenius norm of solution
            call dtt_trunc(x, i, real_tol*resid_damp, 1, xnew=sol)
          end if
          call dcopy(rx1*n(i)*rx2, sol, 1, rhs, 1)

          ! prepare the enrichments: both for X and for Z
          call d2d_mv(rx1, n(i), rx2, rx1, n(i), z%r(i), ra1, ra2, phixax(i-1)%p, A%u(i)%p, phizax(i)%p, rhs, sol)
          call d2d_mv(rx1,n(i),rx2, z%r(i-1),n(i),z%r(i), ra1,ra2, phizax(i-1)%p, A%u(i)%p, phizax(i)%p, rhs, z%u(i)%p)

          call dgemm('n','n', rx1, n(i)*y%r(i), y%r(i-1), 1d0, phixy(i-1)%p, rx1, y%u(i)%p, y%r(i-1), 0d0, resid, rx1)
          call dgemm('n','n', rx1*n(i), z%r(i), y%r(i), 1d0, resid, rx1*n(i), phizy(i)%p, y%r(i), 0d0, rhs, rx1*n(i))
          call daxpy(rx1*n(i)*z%r(i), -1d0, rhs, 1, sol, 1)
          ! now sol is the enrichment to x
          ! update the block for z
          call dgemm('n','n', z%r(i-1), n(i)*y%r(i), y%r(i-1), 1d0, phizy(i-1)%p, z%r(i-1), &
                              y%u(i)%p, y%r(i-1), 0d0, resid, z%r(i-1))
          call dgemm('n','n', z%r(i-1)*n(i), z%r(i), y%r(i), 1d0, resid, z%r(i-1)*n(i), phizy(i)%p, y%r(i), 0d0, rhs, z%r(i-1)*n(i))
          call daxpy(z%r(i-1)*n(i)*z%r(i), -1d0, rhs, 1, z%u(i)%p, 1)
          ! now z(i) = Z'(Ax-y)

          ! enrichment
          call dtt_enrich(x, i, z%r(i), left=sol)
          ! orthogonalization
          call dtt_qr(x, i, dir=+1)
          call dtt_qr(z, i, dir=+1)
          ! update projections
          call dtt_YAX(phizax(i-1)%p, z, x, i, 1, phizax(i)%p, A=A)
          call dtt_YAX(phizy(i-1)%p, z, y, i, 1, phizy(i)%p)
        else
          call dtt_qr(x, i, dir=+1)
        end if

        ! update projections
        call dtt_YAX(phixax(i-1)%p, x, x, i, 1, phixax(i)%p, A=A)
        call dtt_YAX(phixy(i-1)%p, x, y, i, 1, phixy(i)%p)

      end if ! i<d

      ! deallocate all work arrays
      if (prc)then
        if (jacsused)deallocate(jacs)
      end if
      deallocate(rhs, resid, sol)
    end do ! i

    if (verb0>0) then
!#ifdef MEXWRITE
!      write(mxline,"(A,I0,A,ES10.3,A,ES10.3,A,I0,A)"), 'amen_solve: swp=', swp, ', max_dx=', err_max, ', max_res=', res_max, ', max_rank=', maxval(x%r(l-1:d)), char(0)
!      call mexputstr()
!#else
      write(*,"(A,I0,A,ES10.3,A,ES10.3,A,I0)"), 'amen_solve: swp=', swp, ', max_dx=', err_max, &
                                                ', max_res=', res_max, ', max_rank=', maxval(x%r(l-1:d))
!#endif
    end if

    if (trunc_norm0.eq.1) then
      if (res_max.lt.tol) exit
    else
      if (err_max.lt.tol) exit
    end if

  end do ! swp

  do i=l-1,d
    deallocate(phixax(i)%p, phixy(i)%p)
    if (kickrank0>0)deallocate(phizax(i)%p, phizy(i)%p)
  end do
  if (kickrank0>0)call dealloc(z)
 end subroutine

! complex solver
 subroutine ztt_amen_solve(A, y, tol, x, kickrank, nswp, local_prec, local_iters, local_restart, trunc_norm, max_full_size, verb)
  implicit none
  type(ztt),intent(in) :: A, y
  double precision,intent(in) :: tol
  type(ztt),intent(inout) :: x
  integer,intent(in), optional :: kickrank, local_iters, local_restart, nswp, trunc_norm, verb, max_full_size ! trunc_norm: 0 - fro, 1 - resid
  character,intent(in), optional :: local_prec ! 'n', 'l', 'c', 'r'
  character(len=*),parameter :: subnam='ztt_amen_solve'
  double complex,parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
  integer :: kickrank0, local_iters0, local_restart0, nswp0, trunc_norm0, verb0, max_full_size0
  character :: local_prec0
  type(ztt) :: z
  type(pointz) :: phixax(0:tt_size), phixy(0:tt_size), phizax(0:tt_size), phizy(0:tt_size)
  integer :: i, j, l, d, swp, rx1, rx2, ra1, ra2, info
  integer, pointer :: n(:)
  double complex,allocatable :: rhs(:), resid(:), B(:)
  integer, allocatable :: ipiv(:)
  double complex,pointer :: jacs(:), sol(:)
  double precision :: err, err_max, res_max, real_tol
  logical :: prc, jacsused
  double precision,external :: dznrm2
  logical,external :: lsame
  target :: A,y,x,z

  kickrank0 = default(4, kickrank)
  nswp0 = default(20, nswp)
  local_prec0 = default('n',local_prec)
  local_iters0 = default(2, local_iters)
  local_restart0 = default(40, local_restart)
  max_full_size0 = default(50, max_full_size)
  trunc_norm0 = default(1, trunc_norm)
  verb0 = default(1, verb)
  prc=.not.lsame(local_prec0,'n')

  l=y%l; d=y%m; n => y%n

  real_tol = tol/dsqrt(1.d0*(d-l))
  real_tol = real_tol/resid_damp

  ! prepare projections
  allocate(phixax(l-1)%p(1), phixy(l-1)%p(1), phixax(d)%p(1), phixy(d)%p(1))
  phixax(l-1)%p(1)=one; phixax(d)%p(1)=one
  phixy(l-1)%p(1)=one; phixy(d)%p(1)=one

  if (kickrank0.gt.0) then
  ! prepare z
    z%l = l; z%m = d; z%n(l:d)=n(l:d); z%r(l-1)=1; z%r(l:d-1)=kickrank0; z%r(d)=1
    call random(z)
    allocate(phizax(l-1)%p(1), phizy(l-1)%p(1), phizax(d)%p(1), phizy(d)%p(1))
    phizax(l-1)%p(1)=one; phizax(d)%p(1)=one
    phizy(l-1)%p(1)=one; phizy(d)%p(1)=one
  end if

  ! Main cycle
  do swp=1,nswp0
    ! orthogonalization
    do i=d,l+1,-1
      if ((kickrank0.gt.0).and.(swp.gt.1)) then
        ! update Z block
        call z2d_mv(x%r(i-1),n(i),x%r(i), z%r(i-1),n(i),z%r(i), A%r(i-1),A%r(i), phizax(i-1)%p, &
                    A%u(i)%p, phizax(i)%p, x%u(i)%p, z%u(i)%p)
        allocate(rhs(z%r(i-1)*n(i)*z%r(i)), resid(z%r(i-1)*n(i)*max(z%r(i),y%r(i))))
        call zgemm('n','n', z%r(i-1),n(i)*y%r(i),y%r(i-1), one, phizy(i-1)%p,z%r(i-1), y%u(i)%p,y%r(i-1), zero, resid,z%r(i-1))
        call zgemm('n','n', z%r(i-1)*n(i),z%r(i),y%r(i), one, resid,z%r(i-1)*n(i), phizy(i)%p,y%r(i), zero, rhs,z%r(i-1)*n(i))
        call zaxpy(z%r(i-1)*n(i)*z%r(i), -one, rhs,1, z%u(i)%p,1)
        ! now z(i) = Z'(Ax-y)
        deallocate(rhs, resid)
      end if

      call ztt_qr(x, i, dir=-1)
      ! update X projections
      call ztt_YAX(phixax(i)%p, x, x, i, -1, phixax(i-1)%p, A=A)
      call ztt_YAX(phixy(i)%p, x, y, i, -1, phixy(i-1)%p)

      if (kickrank0>0) then
        call ztt_qr(z, i, dir=-1)
        ! update Z projections
        call ztt_YAX(phizax(i)%p, z, x, i, -1, phizax(i-1)%p, A=A)
        call ztt_YAX(phizy(i)%p, z, y, i, -1, phizy(i-1)%p)
      end if
    end do ! ort

    res_max = 0d0; err_max = 0d0
    ! optimization
    do i=l,d
      ! extract sizes
      rx1 = x%r(i-1); rx2 = x%r(i)
      ra1 = A%r(i-1); ra2 = A%r(i)

      ! prepare rhs
      ! first, memory allocations
      if (kickrank0.gt.0) then
        allocate(rhs(max(rx1,z%r(i-1))*n(i)*max(rx2,z%r(i))), resid(max(rx1,z%r(i-1))*n(i)*max(rx2,y%r(i),z%r(i))))
        allocate(sol(max(rx1,z%r(i-1))*n(i)*max(rx2,z%r(i))))
      else
        allocate(rhs(rx1*n(i)*rx2), resid(rx1*n(i)*max(rx2,y%r(i))))
        allocate(sol(rx1*n(i)*rx2))
      end if

      call zgemm('n','n', rx1, n(i)*y%r(i), y%r(i-1), (1d0,0d0), phixy(i-1)%p, rx1, y%u(i)%p, y%r(i-1), (0d0,0d0), resid, rx1)
      call zgemm('n','n', rx1*n(i), rx2, y%r(i), (1d0,0d0), resid, rx1*n(i), phixy(i)%p, y%r(i), (0d0,0d0), rhs, rx1*n(i))

      ! res_prev
      call z2d_mv(rx1, n(i), rx2, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, x%u(i)%p, resid)
      call zaxpy(rx1*n(i)*rx2, (-1d0,0d0), rhs, 1, resid, 1)
      err = dznrm2(rx1*n(i)*rx2, resid, 1)/dznrm2(rx1*n(i)*rx2, rhs, 1)
      res_max = dmax1(res_max,err)
      jacsused = .false.

      if (err.gt.real_tol) then
        if (rx1*n(i)*rx2<max_full_size0) then ! full solution
          allocate(B(rx1*n(i)*rx2*rx1*n(i)*rx2), ipiv(rx1*n(i)*rx2),stat=info)
          if(info.ne.0)then;write(*,*)subnam,': cannot allocate';stop;endif
          call z2d_fullmat(rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, B)
          call zcopy(rx1*n(i)*rx2, resid, 1, sol, 1)
          call zgesv(rx1*n(i)*rx2, 1, B, rx1*n(i)*rx2, ipiv, sol, rx1*n(i)*rx2, info)
          if(info.ne.0)then; write(*,*)subnam,': dgesv info: ',info; stop; end if
          if (verb0>1) write(*,"(A,I0,A)"), 'amen_solve: block=', i, ', direct local solver'
          deallocate(B, ipiv)
        else ! iterative solution
          ! allocate and compute jac prec, if necessary
          if (lsame(local_prec0,'l')) then
            allocate(jacs(rx1*rx1*n(i)*n(i)*max(ra2, rx2)))
            jacsused = .true.
          else if (lsame(local_prec0,'c'))then
            allocate(jacs(rx1*n(i)*n(i)*max(ra2, rx2)))
            jacsused = .true.
          else if (lsame(local_prec0,'r'))then
            allocate(jacs(rx1*n(i)*n(i)*rx2*max(ra2, rx2)))
            jacsused = .true.
          end if
          if (prc) call z2d_jac_gen(local_prec0, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, jacs)

          ! Run the solver
          call z2d_gmresr(phixax(i-1)%p, A%u(i)%p, phixax(i)%p, resid, rx1, n(i), rx2, ra1, ra2, &
                          local_restart0, real_tol/err, local_iters0, local_prec0, jacs, sol, max(verb0-1,0))

          if (prc) call z2d_jac_apply(local_prec0, rx1, n(i), rx2, jacs, sol, sol)
        end if ! max_full_size
        ! Now, sol is the -correction, since the rhs was resid = Ax-y
        call zaxpy(rx1*n(i)*rx2, (-1d0,0d0), sol, 1, x%u(i)%p, 1)
        ! the updated core is ready
        err = dznrm2(rx1*n(i)*rx2, sol, 1)/dznrm2(rx1*n(i)*rx2, x%u(i)%p, 1)
        err_max = dmax1(err_max,err)
        ! compute res_new
        call z2d_mv(rx1, n(i), rx2, rx1, n(i), rx2, ra1, ra2, phixax(i-1)%p, A%u(i)%p, phixax(i)%p, x%u(i)%p, resid)
        call zaxpy(rx1*n(i)*rx2, (-1d0,0d0), rhs, 1, resid, 1)
        err = dznrm2(rx1*n(i)*rx2, resid, 1)/dznrm2(rx1*n(i)*rx2, rhs, 1)
      end if

      if (i.lt.d) then
        if (kickrank0.gt.0) then
          if (trunc_norm0.eq.1) then ! truncation in the A-norm of error (i.e. residual)
            err = dmax1(err,real_tol*resid_damp)
            call ztt_trunc(x, i, err, 1, A=A, phi=phixax, rhs=rhs, xnew=sol)
          else  ! truncation in the frobenius norm of solution
            call ztt_trunc(x, i, real_tol*resid_damp, 1, xnew=sol)
          end if
          call zcopy(rx1*n(i)*rx2, sol, 1, rhs, 1)

          ! prepare the enrichments: both for X and for Z
          call z2d_mv(rx1, n(i), rx2, rx1, n(i), z%r(i), ra1, ra2, phixax(i-1)%p, A%u(i)%p, phizax(i)%p, rhs, sol)
          call z2d_mv(rx1,n(i),rx2, z%r(i-1),n(i),z%r(i), ra1,ra2, phizax(i-1)%p, A%u(i)%p, phizax(i)%p, rhs, z%u(i)%p)

          call zgemm('n','n', rx1, n(i)*y%r(i), y%r(i-1), (1d0,0d0), phixy(i-1)%p, rx1, y%u(i)%p, y%r(i-1), (0d0,0d0), resid, rx1)
          call zgemm('n','n', rx1*n(i), z%r(i), y%r(i), (1d0,0d0), resid, rx1*n(i), phizy(i)%p, y%r(i), (0d0,0d0), rhs, rx1*n(i))
          call zaxpy(rx1*n(i)*z%r(i), (-1d0,0d0), rhs, 1, sol, 1)
          ! now sol is the enrichment to x
          ! update the block for z
          call zgemm('n','n', z%r(i-1), n(i)*y%r(i), y%r(i-1), (1d0,0d0), phizy(i-1)%p, &
                              z%r(i-1), y%u(i)%p, y%r(i-1), (0d0,0d0), resid, z%r(i-1))
          call zgemm('n','n', z%r(i-1)*n(i), z%r(i), y%r(i), (1d0,0d0), resid, &
                              z%r(i-1)*n(i), phizy(i)%p, y%r(i), (0d0,0d0), rhs, z%r(i-1)*n(i))
          call zaxpy(z%r(i-1)*n(i)*z%r(i), (-1d0,0d0), rhs, 1, z%u(i)%p, 1)
          ! now z(i) = Z'(Ax-y)

          ! enrichment
          call ztt_enrich(x, i, z%r(i), left=sol)
          ! orthogonalization
          call ztt_qr(x, i, dir=+1)
          call ztt_qr(z, i, dir=+1)
          ! update projections
          call ztt_YAX(phizax(i-1)%p, z, x, i, 1, phizax(i)%p, A=A)
          call ztt_YAX(phizy(i-1)%p, z, y, i, 1, phizy(i)%p)
        else
          call ztt_qr(x, i, dir=+1)
        end if

        ! update projections
        call ztt_YAX(phixax(i-1)%p, x, x, i, 1, phixax(i)%p, A=A)
        call ztt_YAX(phixy(i-1)%p, x, y, i, 1, phixy(i)%p)

      end if ! i<d

      ! deallocate all work arrays
      if (prc)then
        if (jacsused)deallocate(jacs)
      end if
      deallocate(rhs, resid, sol)
    end do ! i

    if (verb0>0) then
!#ifdef MEXWRITE
!      write(mxline,"(A,I0,A,ES10.3,A,ES10.3,A,I0,A)"), 'amen_solve: swp=', swp, ', max_dx=', err_max, ', max_res=', res_max, ', max_rank=', maxval(x%r(l-1:d)), char(0)
!      call mexputstr()
!#else
      write(*,"(A,I0,A,ES10.3,A,ES10.3,A,I0)"), 'amen_solve: swp=', swp, ', max_dx=', err_max, &
                                                ', max_res=', res_max, ', max_rank=', maxval(x%r(l-1:d))
!#endif
    end if

    if (trunc_norm0.eq.1) then
      if (res_max.lt.tol) exit
    else
      if (err_max.lt.tol) exit
    end if

  end do ! swp

  do i=l-1,d
    deallocate(phixax(i)%p, phixy(i)%p)
    if (kickrank0>0)deallocate(phizax(i)%p, phizy(i)%p)
  end do
  if (kickrank0>0)call dealloc(z)
 end subroutine
end module
