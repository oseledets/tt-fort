module mv_bfun3
use ttals
use matrix_util
use dispmodule
!Here we need to have all parameters required for the matrix-by-vector product (to call bfun3)
 integer, private  ::  rx1T,mT,rx2T,ry1T,nT,ry2T,ra1T,ra2T
 real(8), pointer :: phi1T(:), phi2T(:),res1T(:), res2T(:), AT(:)
 integer,private ::  xsizeT, ysizeT
 type, public ::  pointd
    real(8), dimension(:), pointer :: p=>null()
 end type pointd
 integer, parameter :: primme_kind = 8

contains
 subroutine primme_matvec(x,y,k,primme)
     use ttals, only: dbfun3
     implicit none
     integer(kind=primme_kind) :: primme
     integer, intent(in) :: k
     real(8), intent(in) :: x(*)
     real(8) :: y(*)
     integer :: i
     do i = 1,k
        call dbfun3(rx1T, mT, rx2T, ry1T, nT, ry2T, ra1T, ra2T, phi1T, AT, phi2T, x((i-1)*xsizeT+1), y((i-1)*ysizeT+1))
     end do

 end subroutine primme_matvec

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

 subroutine init_bfun_main(phi1,A,phi2)
   real(8), target :: phi1(:), phi2(:), A(:)
   phi1T => phi1
   phi2T => phi2
   AT => A
 end subroutine init_bfun_main
end module

module tt_block_eig
use ttals
use matrix_util
use trans_lib
implicit none
real(8), allocatable :: result_core(:)

contains
subroutine deallocate_result
  if ( allocated(result_core) ) then
    deallocate(result_core)
  end if
end subroutine deallocate_result

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Primme matvec stuff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AUX routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function compute_residue(mv, M, B, X) result(res)
  use dispmodule
  integer, intent(in) :: M,B
  real(8) X(M,B)
  real(8)  :: Y(M,B), PHI(B,B)
  real(8) res, dnrm2
  external mv

  call mv(X, Y, B, 0)
  call dgemm('t','n', B, B, M, 1d0, X, M, Y, M, 0d0, PHI, B)
  call dgemm('n','n', M, B, B, -1d0, X, M, PHI, B, 1d0, Y, M)
  res = dnrm2(M*B, Y, 1)
end function compute_residue

!!! Generate full local matrix
  subroutine d2d_fullmat(rx1, m, rx2, ra1, ra2, phi1, A, phi2, B)
    integer, intent(in) :: rx1, m, rx2, ra1, ra2
    real(8), intent(in) :: phi1(*), A(*), phi2(*)
    real(8), intent(inout) :: B(*)
    real(8), allocatable :: res1(:), res2(:), phi2t(:)

    allocate(res1(rx1*m*rx1*m*max(ra2, rx2*rx2)), res2(rx1*m*rx1*m*max(ra2, rx2*rx2)), phi2t(ra2*rx2*rx2))

    ! phi1: ry1,rx1,ra1
    call dgemm('N', 'N', rx1*rx1, m*m*ra2, ra1, 1d0, phi1, rx1*rx1, A, ra1, 0d0, res1, rx1*rx1)
    ! res1: ry1,rx1,n,m,ra2
    call dperm1324(rx1, rx1, m, m*ra2, res1, res2)
    ! phi2: rx2,ra2,ry2
    call trans2d(rx2, ra2*rx2, phi2, phi2t)
    call dgemm('N', 'N', rx1*m*rx1*m, rx2*rx2, ra2, 1d0, res2, rx1*m*rx1*m, phi2t, ra2, 0d0, res1, rx1*m*rx1*m);
    call dperm1324(rx1*m, rx1*m, rx2, rx2, res1, B)
    deallocate(res1, res2, phi2t)
  end subroutine



!! What we have: we have a starting vector + a matrix (no vector X!)
subroutine tt_eigb(d,n,m,ra,crA, crY0, ry, eps, rmax, lambda, B, kickrank, nswp, max_full_size, verb)
use matrix_util
use ttals
use mv_bfun3
integer,intent(in) :: d,n(d),m(d),ra(d+1), rmax
integer,intent(inout) :: ry(d+1)
integer, intent(in), optional :: kickrank, nswp, verb, max_full_size
! verb: 0: off; 1: matlab; 2: console
integer :: kickrank0, nswp0, verb0, max_full_size0, max_basis_size
integer :: B
real(8), intent(in) :: crA(*),eps, crY0(*)
real(8) eps2
type(pointd) :: crnew(d+1)
type(pointd) :: phinew(d+1)
real(8),allocatable, target :: curcr(:), locmat(:), work(:)
real(8) erloc, resid_damp
real(8) :: sv(B*rmax)
real(8), intent(out) :: lambda(B)
real(8), allocatable :: rnorms(:), U(:,:)
real(8) fv, fvold !Functional
integer(kind=primme_kind) ::  primme, ierr, num_matvecs
include 'primme_f77.h'
integer, allocatable :: pa(:)
integer :: i,k, swp, dir, lwork, mm, nn, rnew, max_matvecs, rmax2, sz
integer info
integer total_mv
real(8) :: ermax, res,  min_res
real(8) dnrm2
INTEGER :: i_, n_, clock_
INTEGER, DIMENSION(:), ALLOCATABLE :: seed_

CALL RANDOM_SEED(size = n_)
ALLOCATE(seed_(n_))
CALL SYSTEM_CLOCK(COUNT=clock_)
seed_ = clock_ + 37 * (/ (i_ - 1, i_ = 1, n_) /)
CALL RANDOM_SEED(PUT = seed_)
DEALLOCATE(seed_)
fvold = 0d0
total_mv = 0
min_res = 1d-1
rmax2 = rmax
!Inner parameters
max_matvecs = 500
max_basis_size = 100
resid_damp = 10d0

call disp('Solving a block eigenvalue problem')
if ( B .ne. ry(d+1) ) then
   print *,'Inconsistent number of eigenvalues sought and initial guess last rank'
   return
end if
 allocate(rnorms(B)) !For the residual norms
eps2 = eps/sqrt(real(d,8))
call disp('Looking for '//tostring(B*1d0)//' eigenvalues with accuracy '//tostring(eps))

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
  max_full_size0 = 2000
  if (present(max_full_size)) then
    max_full_size0 = max_full_size
  end if
  allocate(pa(d+1))
  call compute_ps(d,ra,n(1:d)*m(1:d),pa)
! Work arrays for MPS blocks
  sz = rmax*maxval(n(1:d))*rmax*max(B+1, max_basis_size)
  lwork = max(rmax*max(rmax,max(B,maxval(n(1:d)))), max_full_size0)*256
  allocate(curcr(sz), work(lwork), stat=info)
  if (info .ne. 0) then
    print *, 'cannot allocate'
    stop
  end if
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
  do while (i <= d)
    rnew = min(ry(i)*n(i), ry(i+1))
    call dqr(ry(i)*n(i), ry(i+1), crnew(i) % p, work)
    if ( i < d ) then
       call dgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, work, rnew, crnew(i+1)%p, ry(i+1), 0d0, curcr, rnew)
       call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, crnew(i+1)%p,1)
       ry(i+1) = rnew;
       !     Phir
       allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
       call dphi_left(ry(i), m(i), ry(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, &
            phinew(i+1)%p)
    end if
    i = i+dir
  end do
  i = d
  !goto 100
  ry(d+1) = 1 !Actually, this is the true stuff.

  ermax = 0d0
  i = d
  dir = -1
  swp = 0
  swp = swp + 1

  ermax = 0d0
  lambda(1:B) = 0d0
  do while (swp .le. nswp0)
      call init_bfun_sizes(ry(i),n(i),ry(i+1),ry(i),n(i),ry(i+1),ra(i),ra(i+1),ry(i)*n(i)*ry(i+1),ry(i)*n(i)*ry(i+1))
      call init_bfun_main(phinew(i)%p,crA(pa(i):pa(i+1)-1),phinew(i+1)%p)

      if (ry(i)*n(i)*ry(i+1)<max_full_size0) then
          allocate(locmat(ry(i)*n(i)*ry(i+1)*ry(i)*n(i)*ry(i+1)))
          call d2d_fullmat(ry(i), n(i), ry(i+1), ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), phinew(i+1)%p, locmat)
          ! Full EIG
          call dsyev('V', 'U', ry(i)*n(i)*ry(i+1), locmat, ry(i)*n(i)*ry(i+1), curcr, work, lwork, info)
          if ( info .ne. 0 ) then
              print *,'tt_eigb: dsyev failed', info
              return
          end if

          call dcopy(B, curcr, 1, lambda, 1)
          call dcopy(ry(i)*n(i)*ry(i+1)*B, locmat, 1, curcr,1)

          res = 0d0
          num_matvecs = 0
          deallocate(locmat)
      else
          !Initialization of the primme stuff;
          call primme_initialize_f77(primme)
          call primme_set_member_f77(primme, PRIMMEF77_numEvals, B)
          !        call primme_set_method_f77(primme,PRIMMEF77_LOBPCG_OrthoBasis, ierr)  !Select LOBPCG
          call primme_set_method_f77(primme, PRIMMEF77_DYNAMIC, ierr)

          call primme_set_member_f77(primme, PRIMMEF77_n, ry(i)*n(i)*ry(i+1))

          call primme_set_member_f77(primme, PRIMMEF77_matrixMatvec, primme_matvec)
          call primme_set_member_f77(primme, PRIMMEF77_printLevel, 1)
          !if ( dir < 0 ) then
          !   call primme_set_member_f77(primme, PRIMMEF77_maxMatvecs, 500)
          !end if
          !Uncomment if dynamic
          !call primme_set_member_f77(primme, PRIMMEF77_eps, min(min_res, res_old/(1000)))
          !call primme_display_params_f77(primme)
          !call primme_set_member_f77(primme, PRIMMEF77_initSize, B)
          !call primme_set_member_f77(primme, PRIMMEF77_maxBlockSize, ry(i)*n(i)*ry(i+1))
          call primme_set_member_f77(primme, PRIMMEF77_maxMatvecs, max_matvecs)
          !Calculate initial residue
          call primme_set_member_f77(primme, PRIMMEF77_minRestartSize,min(B+1,ry(i)*n(i)*ry(i+1)-1))
          call primme_set_member_f77(primme, PRIMMEF77_maxBasisSize, max_basis_size)
          call primme_set_member_f77(primme, PRIMMEF77_eps, eps2/resid_damp)
          call primme_set_member_f77(primme, PRIMMEF77_initSize, B)
          ! if ( (swp .eq. 4) .and. (i .eq. 35) ) then
          !call primme_set_member_f77(primme, PRIMMEF77_printLevel,4)
          ! end if
          !Solver part
          !return
          call dcopy(ry(i)*n(i)*ry(i+1)*B,crnew(i)%p,1,curcr,1)
          call dprimme_f77(lambda, curcr, rnorms, primme, ierr)
          call primmetop_get_member_f77(primme, PRIMMEF77_stats_numMatvecs, num_matvecs)
          call primme_Free_f77(primme)
          !Compute the local residue
          if ( ierr < 0 .and. (ierr .ne. -3) ) then
              print *,'tt_eigb, Primme failed with ierr=',ierr
              return
          end if

          res = dnrm2(B,rnorms,1)

      end if
      fv = sum(lambda(1:B))
      erloc = (fvold - fv)/abs(fv)
      ermax = max(ermax, erloc)
      if (verb0 > 0) then
      call disp('swp: '//tostring(1d0*swp)//' i: ['//tostring(1d0*i)//'/'//&
          tostring(1d0*d)//'] loc_size: '//tostring(1d0*ry(i)*n(i)*ry(i+1))//&
          ' matvecs: '//tostring(1d0*num_matvecs)//' res: ' //tostring(res)//&
          ' ermax: '//tostring(ermax)//' dfv:'//tostring(fvold-fv))
      end if
      total_mv = total_mv + num_matvecs
      fvold = fv

      if ( (dir < 0) .and.  (i > 1) ) then
          !the curcr is ry(i)*n(i)*ry(i+1)*B,
          call dtransp(ry(i)*n(i)*ry(i+1),B,curcr)
          !curcr [B x ry(i)] x [n(i) x ry(i+1)]
          rnew = min(B*ry(i),n(i)*ry(i+1))
          allocate(U(B*ry(i),rnew))

          call dgesvd('S','O',B*ry(i),n(i)*ry(i+1),curcr,B*ry(i),sv,U,B*ry(i), curcr, rnew , work, lwork, info)
          if ( info .ne. 0 ) then
              print *,'tt_eigb: dgesvd failed'
              return
          end if

          do k = 1,rnew
          call dscal(B*ry(i),sv(k), U(1,k), 1)
          end do
          nn = my_chop3(rnew,sv,eps2)
          nn = min(nn, rmax)
          call drow_cut(B*ry(i), n(i)*ry(i+1), nn, curcr)
          rnew = nn
          call dcopy(rnew*n(i)*ry(i+1),curcr,1,crnew(i)%p,1)


          !call dcopy(rnew*n(i)*ry(i+1),curcr,1,cr(1,i),1)
          !And the transformation matrix is stored in U with is B*ry(i) * rnew
          call dtransp(B,ry(i)*rnew, U)

          !U is ry(i)*rnew*B
          call dgemm('n','n',ry(i-1)*n(i-1),rnew*B,ry(i),1d0,crnew(i-1)%p,ry(i-1)*n(i-1),U,ry(i),0d0,curcr,ry(i-1)*n(i-1))
          deallocate(U)
          !curcr ry(i-1)*n(i-1)*rnew*B
          if ( ry(i-1)*n(i-1)*rnew < B ) then
              print *,'tt_eigb: Something stupid happened, have to think!'
              return
          end if

          if ( size(crnew(i-1)%p) < ry(i-1)*n(i-1)*rnew*B ) then
              deallocate(crnew(i-1)%p)
              allocate(crnew(i-1)%p(ry(i-1)*n(i-1)*rnew*B*2))
          end if
          !call dqr(ry(i-1)*n(i-1)*rnew,B,curcr,R,work,lwork,tau) !Here we keep orthogonality
          call dcopy(ry(i-1)*n(i-1)*rnew*B,curcr,1,crnew(i-1)%p,1)
          ry(i) = rnew
          !And recompute phi
          if ( size(phinew(i)%p) < ry(i)*ry(i)*ra(i) ) then
              deallocate(phinew(i)%p)
              allocate(phinew(i)%p(ry(i)*ry(i)*ra(i)*2))
          end if
          call dphi_right(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
              ra(i), ra(i+1), phinew(i+1)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phinew(i)%p)

      end if
      if ( (dir > 0) .and. (i < d) ) then
          !curcr is ry(i)*n(i)*ry(i+1)*B
          rnew = min(ry(i)*n(i),ry(i+1)*B)
          allocate(U(rnew, ry(i+1)*B))
          call dgesvd('O', 'S',ry(i)*n(i),ry(i+1)*B,curcr, ry(i)*n(i), sv, curcr, ry(i)*n(i), U, rnew, work, lwork, info)
          nn = my_chop3(rnew,sv,eps2)
          nn = min(nn, rmax)
          call dcopy(ry(i)*n(i)*nn,curcr,1,crnew(i)%p,1)
          !Cut rows
          call drow_cut(rnew, ry(i+1)*B,  nn, U)
          do k = 1,nn
          call dscal(ry(i+1)*B,sv(k), U(k,1), nn)
          end do

          rnew = nn
          !U is (rnew x ry(i+1) x B) x (ry(i+1) * n(i+1) * ry(i+2)) -> rnew x n(i+1) x ry(i+2) x B
          ! (rnew x ry(i+1)) * B  -> B * rnew * ry(i+1)
          call dtransp(rnew*ry(i+1),B,U)

          call dgemm('n','n',B*rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, U, B*rnew, crnew(i+1)%p, ry(i+1), 0d0, curcr, B*rnew)

          deallocate(U)
          call dtransp(B, rnew*n(i+1)*ry(i+2), curcr)
          ry(i+1) = rnew
          if ( size(crnew(i+1)%p) < ry(i+1)*n(i+1)*ry(i+2)*B) then
              deallocate(crnew(i+1)%p)
              allocate(crnew(i+1)%p(ry(i+1)*n(i+1)*ry(i+2)*B*2))
          end if
          call dcopy(ry(i+1)*n(i+1)*ry(i+2)*B, curcr, 1, crnew(i+1)%p, 1)
          !call dqr(ry(i+1)*n(i+1)*ry(i+2),B,crnew(i+1)%p,R,work,lwork,tau) !Here we keep orthogonality

          if ( size(phinew(i+1)%p) < ry(i+1)*ry(i+1)*ra(i+1) ) then
              deallocate(phinew(i+1)%p)
              allocate(phinew(i+1)%p(ry(i+1)*ry(i+1)*ra(i+1)*2))
          end if
          call dphi_left(ry(i), n(i), ry(i+1), ry(i), n(i), ry(i+1), &
              ra(i), ra(i+1), phinew(i)%p, crA(pa(i)), crnew(i)%p, crnew(i)%p, phinew(i+1)%p)
      end if


      if ((dir>0) .and. (i==d - 1)) then
          dir = -1
          i = d
          call disp('swp: '//tostring(1d0*swp)//' er = '//tostring(ermax)//' rmax:'//tostring(1d0*maxval(ry(1:d))))
          swp = swp + 1
          if (ermax<eps) then
              exit
          end if
          ermax = 0d0
      else if ((dir < 0) .and. (i == 2 )) then
          dir = 1
          i = 1
      else
          i = i + dir
      end if
 end do
  ry(d+1) = B
  nn = sum(ry(2:d+1)*ry(1:d)*n(1:d))
  if ( .not. allocated(result_core) ) then
     allocate(result_core(nn))
  end if

  nn = 1
  do i=1, d
    call dcopy(ry(i)*n(i)*ry(i+1), crnew(i)%p, 1, result_core(nn), 1)
    nn = nn+ry(i)*n(i)*ry(i+1)
  end do
  do i = 1, d+1
     if ( associated(crnew(i)%p)) then
        deallocate(crnew(i)%p)
     end if
     if ( associated(phinew(i)%p)) then
        deallocate(phinew(i)%p)
     end if
  end do
  call disp('Total number of matvecs: '//tostring(total_mv*1d0))
   deallocate(work)
   deallocate(curcr)
   deallocate(pa)
end subroutine



end module
