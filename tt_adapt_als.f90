module tt_adapt_als

  use iso_c_binding

  implicit none
  real(8), allocatable :: result_core(:)
!   character*120 :: matlab_ist_dumme_kuh

  type :: dpoint
     !   double precision,dimension(:),pointer :: p=>null()
     real(8), allocatable :: p(:)
  end type dpoint

contains


subroutine deallocate_result
  if ( allocated(result_core) ) then
    deallocate(result_core)
  end if
end subroutine deallocate_result


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! AMEn MatVec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine tt_mvk4(d,n,m,rx,ra,crA, crX, crY0, ry, eps, rmax, kickrank, nswp, verb)
    use tt_linalg
    use dispmodule

    integer,intent(in) :: d,n(*),m(*),rx(*),ra(*), rmax
    integer,intent(inout) :: ry(*)
    integer, intent(in), optional :: kickrank, nswp, verb
    ! verb: 0: off; 1: matlab; 2: console
    integer :: kickrank0, nswp0, verb0
    real(8), intent(in) :: crA(*), crX(*), eps, crY0(*)
    real(8) eps2
    ! integer cr_size
    ! real(8),allocatable ::  curcr(:), work(:), tau(:), R(:)
    ! real(8),allocatable :: crl(:,d), crr(:,d), phil(:,d), phir(:,d)
    real(8),allocatable :: cr(:,:)
    ! real(8) :: crr(1000)
    real(8),allocatable :: curcr(:)
    real(8),allocatable :: work(:)
    real(8),allocatable :: tau(:)
    real(8),allocatable :: R(:)
    ! real(8) :: phil(1000)
    real(8),allocatable :: phi(:,:)
    real(8),allocatable :: full_core(:)

    integer,allocatable :: pa(:), px(:)
    integer :: i,j,k, swp, dir, lwork, mm, nn, rnew
    integer info
    real(8) :: err, err_max

    real(8) dnrm2
    real(8) sqrt

    INTEGER :: i_, n_, clock_
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed_

    CALL RANDOM_SEED(size = n_)
    ALLOCATE(seed_(n_))
    CALL SYSTEM_CLOCK(COUNT=clock_)
    seed_ = clock_ + 37 * (/ (i_ - 1, i_ = 1, n_) /)
    CALL RANDOM_SEED(PUT = seed_)
    DEALLOCATE(seed_)


    eps2 = eps/sqrt(real(d,8))
    ! print *, eps2

    ! real(4) etime
    ! real(4), dimension(2) :: tarray
    ! real(4) :: t0,t1
    !
    !   print *, rmax
    !   print *, d
    !   print *, rmax*maxval(n(1:d))*rmax
    !   print *, rx
    !   print *, ra
    !   print *, ry

    kickrank0 = 5;
    if (present(kickrank)) then
       kickrank0 = kickrank
    end if
    !   if (.not.present(lwork_core)) then
    !     lwork_phi = sum(ry(2:d+1)*ry(1:d)*n(1:d))*2
    !   end if
    !   if (.not.present(lwork_phi)) then
    !     lwork_phi = sum(rx*ry*ra)*2
    !   end if
    nswp0 = 20
    if (present(nswp)) then
       nswp0 = nswp
    end if
    verb0 = 1
    if (present(verb)) then
       verb0 = verb
    end if

    allocate(pa(d+1), px(d+1))
    !   print *, nswp0

    call compute_ps(d,ra,n(1:d)*m(1:d),pa)
    call compute_ps(d,rx,m,px)
    !   px(1)=1;
    !   do i=1,d
    !     px(i+1) = px(i) + rx(i)*m(i)*rx(i+1)
    !   end do


    !   print *, pa
    !   print *, px

    !   working array: position is at i+1
    !   previous array: position is at i
    !   phi: at i
    !   pyr = sum(ry(2:d+1)*ry(1:d)*n(1:d))+1
    ! !   print *, pyr
    !   pyl = pyr-ry(d)*n(d)*ry(d+1);
    !   ppr = sum(rx*ry*ra) + 1
    !   ppl = ppr;
    lwork = rmax*maxval(n(1:d))*rmax

    allocate(cr(lwork, d))
    nn = maxval(rx(1:d+1))*maxval(ra(1:d+1))*rmax
    allocate(phi(nn, d+1))
    !   print *, lwork
    !   allocate(crr(pyr-1))
    !   allocate(phir(ppr-1))
    !   allocate(crl(pyr-1))
    !   allocate(phil(ppl-1))
    allocate(curcr(lwork))
    !   allocate(work(lwork))
    nn = maxval(ra(1:d+1))*maxval(rx(1:d+1))*max(maxval(n(1:d)),maxval(m(1:d)))*rmax
    allocate(work(nn))
    allocate(R(lwork))
    allocate(tau(lwork))
    allocate(full_core(nn))

    !   cr(:,:)=0d0
    mm = 1
    do i=1,d
       call dcopy(ry(i)*n(i)*ry(i+1), crY0(mm), 1, cr(1,i), 1)
       mm = mm + ry(i)*n(i)*ry(i+1)
    end do

    phi(1,1)=1d0
    phi(1,d+1)=1d0

    !   t0 =  etime(tarray)
    !   QR, psi
    dir = 1
    i = 1
    do while (i<d)
       !     current core ppr-r1*n*r2:ppr-1. Transpose it and QR
       !     mm=n(i)*ry(i+1); nn=ry(i)
       rnew = min(ry(i)*n(i), ry(i+1))
       !     call dcopy(nn*mm, crr(pyr-ry(i)*n(i)*ry(i+1)), 1, curcr, 1)
       !     print *, curcr(1:nn*mm)
       !     call dtransp(nn, mm, cr(1,i), curcr)
       call dqr(ry(i)*n(i), ry(i+1), cr(1,i), R, work, lwork, tau)
       !     print *, curcr(1:nn*mm)

       !     call dtransp(mm, rnew, curcr, cr(1,i))

       call dgemm('N', 'N', rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, R, rnew, cr(1,i+1), ry(i+1), 0d0, curcr, rnew)
       !     call dgemm('N', 'T', ry(i-1)*n(i-1), rnew, ry(i), 1d0, cr(1,i-1), ry(i-1)*n(i-1), R, rnew, 0d0, curcr, ry(i-1)*n(i-1))
       call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, cr(1,i+1),1)
       ry(i+1) = rnew;

       !     Phir
       !     call dphi_right(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), phi(1,i+1), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i))
       call dphi_left(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), &
       phi(1,i), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i+1), work, full_core)

       !     write(*,"(A,I0)"), 'als_fort:  qr:', i

       i = i+dir
    end do

    !   t1 = etime(tarray)
    !   print *, 'qr time: ', t1-t0

    !   print *, phi(1:2,1)
    !   print *, phi(1:2,2)
    !   print *, phi(1:2,3)
    !   print *, phi(1:2,4)

    if (verb0==2) then
       print *, ''
    end if
    err_max = 0d0
    i = d;
    dir = -1;
    swp = 1
    !   t0 = etime(tarray)
    do while (swp .le. nswp0)

       !     write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  dbfun3 started:[', swp, ',', i, ',', dir, ']'
       call dbfun3(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), ra(i+1), &
       phi(1,i), crA(pa(i)), phi(1,i+1), crX(px(i)), curcr, work, full_core)
       !     write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  dbfun3 done:[', swp, ',', i, ',', dir, ']'

       !     print *, 'cr: ', curcr(1:ry(i)*n(i)*ry(i+1)), 'ry1: ', ry(i), 'ry2: ', ry(i+1)

       ! curcr is of size: r1,n,r2. we need r1n,(r2+kick)

       !     print *, curcr(1:ry(i)*n(i)*ry(i+1))
       !    err = dnrm2(ry(i)*n(i)*ry(i+1), curcr(1:ry(i)*n(i)*ry(i+1)) - cr(1:ry(i)*n(i)*ry(i+1),i), 1) / dnrm2(ry(i)*n(i)*ry(i+1), curcr, 1)
       call dcopy(ry(i)*n(i)*ry(i+1), cr(1,i), 1, tau, 1)
       call daxpy(ry(i)*n(i)*ry(i+1), -1d0, curcr, 1, tau, 1)
       err = dnrm2(ry(i)*n(i)*ry(i+1), tau, 1) / dnrm2(ry(i)*n(i)*ry(i+1), curcr, 1)
       if (verb0>=2) then ! verb on, carriage_return, no matlab. Otherwise - a lot of shit comes...
          call disp('als_fort: iteration:['// tostring(swp*1d0) // ',' // &
          tostring(i*1d0) // ']  err:' // tostring(err) // '  ry(i):' // &
          tostring(ry(i)*1d0) // '  err_max:' // tostring(err_max))
          !      write(*,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3,20X$)"), 'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max
          !!       write(*,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3)")  'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max
       end if
       !     write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3)"),  'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max


       err_max = max(err_max,err)


       if ((dir>0) .and. (i<d)) then
          rnew = min(ry(i+1), ry(i)*n(i))
          if (kickrank0>0) then
             call dgesvd('O', 'S', ry(i)*n(i), ry(i+1), curcr, ry(i)*n(i), tau, 0, 1, R, rnew, work, lwork, info)
             ! curcr is U and has the size r1*n,rnew; R is of size rnew,r2
             nn = my_chop3(rnew, tau, eps2)
             ! 	nn = rnew
             call drow_cut(rnew, ry(i+1), nn, R)
             do j=1,nn
                call dscal(ry(i+1), tau(j), R(j), nn)
             end do

             ! 	write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  svd:[', swp, ',', i, ',', dir, ']'

             rnew = nn
             ! 	print *, 'nn: ', nn

             ! 	print *, curcr(1:ry(i)*n(i)*nn)
             ! 	print *, '^curcr, vR'
             ! 	print *, R(1:nn*ry(i+1))

             ! 	call random_number(curcr(ry(i)*n(i)*nn+1:ry(i)*n(i)*(nn+kickrank0)))
             ! 	call random_number(work(1:ry(i)*n(i)*kickrank0))
             ! 	call dcopy(ry(i)*n(i)*kickrank0, work, 1, curcr(ry(i)*n(i)*nn+1), 1)
             ! 	tau(1:kickrank0*ry(i+1))=0d0
             ! 	print *, 'R: ', R(1:nn*ry(i+1))
             ! 	call drow_add(nn,ry(i+1),kickrank0,R,tau)
             ! 	print *, 'R2: ', R(1:(nn+kickrank0)*ry(i+1))
             ! 	print *, nn, kickrank0, ry(i+1)
             ! 	print *, R(1:(nn+kickrank0)*ry(i+1))
             ! 	print *, R(1:nn*ry(i+1))
             ! 	print *, ''
             ! 	call dcopy((nn+kickrank0)*ry(i+1), work, 1, R, 1)

             ! 	call dgemm('N', 'N', ry(i)*n(i), ry(i+1), nn+kickrank0, 1d0, curcr, ry(i)*n(i), R, nn+kickrank0, 0d0, work, ry(i)*n(i))
             ! 	print *, 'cr2: ', work(1:ry(i)*n(i)*ry(i+1))


             ! 	call dqr(ry(i)*n(i), nn+kickrank0, curcr, cr(1,i), work, lwork, tau)
             ! 	rnew = min(ry(i)*n(i), nn+kickrank0)
             ! 	print *, '^R  vR2'
             ! 	print *, cr(1:rnew*(nn+kickrank0),i)
             ! 	print *, '^R2  v curcr'
             ! 	print *, curcr(1:ry(i)*n(i)*rnew)
             ! 	print *, 'rnew ', rnew
             ! 	call dgemm('N', 'N', rnew, ry(i+1), nn+kickrank0, 1d0, cr(1,i), rnew, R, nn+kickrank0, 0d0, work, rnew)
             ! 	call dcopy(rnew*ry(i+1), work, 1, R, 1)

             ! 	call dgemm('N', 'N', ry(i)*n(i), ry(i+1), rnew, 1d0, curcr, ry(i)*n(i), R, rnew, 0d0, work, ry(i)*n(i))
             ! ! 	call dgemm('T', 'N', rnew, rnew, ry(i)*n(i), 1d0, curcr, ry(i)*n(i), curcr, ry(i)*n(i), 0d0, work, rnew)
             ! 	print *, 'cr2: ', work(1:ry(i)*n(i)*ry(i+1))
             ! 	print *, R(1:rnew*ry(i+1))
             ! 	print *, ''

             ! 	return
             !         call dtransp(rnew, ry(i+1), R, work)
             !         do j=1,nn
             !           call dscal(ry(i+1), tau(j), work(1+(j-1)*ry(i+1)), 1)
             !         end do
             !         call dtransp(ry(i+1), nn, work, R)
          else
             call dqr(ry(i)*n(i), ry(i+1), curcr, R, work, lwork, tau)
          end if

          call dcopy(ry(i)*n(i)*rnew, curcr, 1, cr(1,i), 1)
          call dgemm('N','N', rnew, n(i+1)*ry(i+2), ry(i+1), 1d0, R, rnew, cr(1,i+1), ry(i+1), 0d0, curcr, rnew)
          call dcopy(rnew*n(i+1)*ry(i+2), curcr, 1, cr(1,i+1), 1)

          ry(i+1) = rnew
          call dphi_left(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), &
          ra(i+1), phi(1,i), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i+1), work, full_core)

          !       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  phi:[', swp, ',', i, ',', dir, ']'
       end if

       !     print *, phi(1:2,1)
       !     print *, phi(1:2,2)
       !     print *, phi(1:2,3)
       !     print *, phi(1:2,4)

       if ((dir<0) .and. (i>1)) then
          rnew = min(n(i)*ry(i+1), ry(i))
          if (kickrank0>0) then
             call dgesvd('S', 'O', ry(i), n(i)*ry(i+1), curcr, ry(i), tau, R, ry(i), 0, 1, work, lwork, info)
             ! 	print *, 'VT: ', curcr(1:ry(i)*n(i)*ry(i+1))
             ! 	print *, 'R: ', R(1:ry(i)*rnew)
             ! 	print *, tau(1:rnew)
             ! curcr is V and has the size rnew,n*r2; R is U and of size r1,rnew
             nn = my_chop3(rnew, tau, eps2)
             call dtransp(ry(i), n(i)*ry(i+1), curcr, work)
             call dcopy(n(i)*ry(i+1)*nn, work, 1, curcr, 1) ! now, n1,ry2,rnew
             call dtransp(ry(i), nn, R, work)
             call dcopy(nn*ry(i), work,1, R,1) ! rnew,r1
             do j=1,nn
                call dscal(ry(i), tau(j), R(j), nn)
             end do
             !         call drow_cut(ry(i), n(i)*ry(i+1), nn, curcr) ! its size is correct

             ! 	write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  svd:[', swp, ',', i, ',', dir, ']'
             !         print *, 'calling full_core...'
             !         if (i==d) then
             !           print *, crA(pa(i):pa(i+1)-1)
             !           print *, crX(px(i):px(i+1)-1)
             !         end if
             ! kick
             call dbfun3_right(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), &
             ra(i+1), crA(pa(i)), phi(1,i+1), crX(px(i)), full_core, full_core, work)
             !         if (i==d) then
             !           print *, full_core(1:rx(i)*ra(i)*n(i)*ry(i+1))
             !         end if
             ! 	! full_core is of size rx1*ra1, n*ry2
             call duchol_fort('N', rx(i)*ra(i), n(i)*ry(i+1), full_core,  min(kickrank0, n(i)*ry(i+1)), work, tau, mm);
             ! 	if (i==d) then
             !           print *, work(1:n(i)*ry(i+1)*kickrank0)
             !         end if
             ! 	print *, mm
             call dcopy(n(i)*ry(i+1)*mm, work, 1, curcr(n(i)*ry(i+1)*nn+1), 1)
             ! 	! work=du: n(i)*ry(i+1),kickrank

             !         mm = kickrank0
             ! 	call random_number(curcr(n(i)*ry(i+1)*nn+1:n(i)*ry(i+1)*(nn+kickrank0)))
             do j=1,mm*ry(i)
                tau(j)=0d0
             end do
             call drow_add(nn, ry(i), mm, R, tau)

             !         print *, 'VT2: ', curcr(1:nn*n(i)*ry(i+1))
             !         call dtransp(ry(i), nn, R, work) ! R has to be transposed
             ! 	call dcopy(nn*ry(i), work, 1, R, 1)

             call dqr(n(i)*ry(i+1), nn+mm, curcr, cr(1,i), work, lwork, tau)
             rnew = min(n(i)*ry(i+1), nn+mm)
             call dgemm('N', 'N', rnew, ry(i), nn+mm, 1d0, cr(1,i), rnew, R, nn+mm, 0d0, work, rnew)
             call dcopy(rnew*ry(i), work, 1, R, 1)

             ! 	write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  kick:[', swp, ',', i, ',', dir, ']'

             !         rnew = nn
             !         print *, rnew
          else
             call dtransp(ry(i), n(i)*ry(i+1), curcr, cr(1,i))
             call dqr(n(i)*ry(i+1), ry(i), cr(1,i), R, work, lwork, tau)
             call dcopy(n(i)*ry(i+1)*ry(i), cr(1,i), 1, curcr, 1)
             !         call dtransp(n(i)*ry(i+1), rnew, cr(1,i), curcr)
          end if

          call dtransp(n(i)*ry(i+1), rnew, curcr, cr(1,i))
          !       call dcopy(rnew*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
          !       call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
          call dgemm('N','T', ry(i-1)*n(i-1), rnew, ry(i), 1d0, cr(1,i-1), ry(i-1)*n(i-1), R, rnew, 0d0, curcr, ry(i-1)*n(i-1))
          call dcopy(ry(i-1)*n(i-1)*rnew, curcr, 1, cr(1,i-1), 1)

          ry(i) = rnew
          call dphi_right(rx(i), m(i), rx(i+1), ry(i), n(i), ry(i+1), ra(i), &
          ra(i+1), phi(1,i+1), crA(pa(i)), crX(px(i)), cr(1,i), phi(1,i), work, full_core)
          !       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  phi:[', swp, ',', i, ',', dir, ']'
       end if


       if ((dir>0) .and. (i==d)) then
          call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
          if (verb0>=1) then
            call disp('als_fort: iteration:'// tostring(swp*1d0) // '(' // &
            tostring(dir*1d0) // ')  max(ry):' // tostring(maxval(ry(1:d+1))*1d0) &
            // ' err_max:' // tostring(err_max));
             !!    write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,I0,A,ES10.3)") 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry(1:d+1)), '  err_max:', err_max
             ! 	write(*,"(A,I0,A,I0,A,I0,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry), '  err_max:', err_max
             ! 	call mexPrintf(matlab_ist_dumme_kuh//achar(10))
          end if
          dir = -1
          i = i+1
          swp = swp+1
          if (err_max<eps) then
             exit
          end if
          err_max = 0d0
          !       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  sweep_reverse:[', swp, ',', i, ',', dir, ']'
       end if

       if ((dir<0) .and. (i==1)) then
          call dcopy(ry(i)*n(i)*ry(i+1), curcr, 1, cr(1,i), 1)
          if (verb0>=1) then
            call disp('als_fort: iteration:'// tostring(swp*1d0) // '(' // tostring(dir*1d0) &
            // ')  max(ry):' // tostring(maxval(ry(1:d+1))*1d0) // ' err_max:' // tostring(err_max));
             !!    write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,I0,A,ES10.3)") 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry(1:d+1)), '  err_max:', err_max
             ! 	write(*,"(A,I0,A,I0,A,I0,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(ry), '  err_max:', err_max
             ! 	call mexPrintf(matlab_ist_dumme_kuh//achar(10))
          end if
          dir = 1
          i = i-1
          !       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  sweep_reverse:[', swp, ',', i, ',', dir, ']'
       end if

       i=i+dir

    end do

    if (verb0==2) then
       print *, ''
    end if

    !   t1 = etime(tarray)
    !   print *, 'sweep time: ', t1-t0


    nn = sum(ry(2:d+1)*ry(1:d)*n(1:d));
    allocate(result_core(nn))

    nn = 1;
    do i=1,d
       call dcopy(ry(i)*n(i)*ry(i+1), cr(1,i), 1, result_core(nn), 1)
       nn = nn+ry(i)*n(i)*ry(i+1)
    end do

    !   write(*,"(A)"), 'als_fort:  result_core dcopy'
    !   print *, ry
    !   call dcopy(sz, crr(pyr), 1, result_core, 1)

    deallocate(cr)
    !    deallocate(crl)
    deallocate(R)
    deallocate(work)
    deallocate(curcr)
    deallocate(phi)
    !    deallocate(phil)
    deallocate(tau)
    deallocate(full_core)
    deallocate(pa,px)

    !   print *, pyl

  end subroutine tt_mvk4


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! AMEn lin. solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tt_amen_solve
! Solution of linear systems in the TT-format
! Input:
! integer d, real(8) arrays n(d), m(d)
! Matrix A: 
! Multilevel matrix with mode sizes n(1)n(2)...n(d) x m(1) m(2) ... m(d)

 subroutine tt_amen_solve(d,n,m,ry,ra,crA, crY, crX0, rx, eps, kickrank, nswp, verb, prec, nrestart, niters)
    use tt_linalg
    use dispmodule

    integer,intent(in) :: d,n(*),m(*),ry(*),ra(*)
    integer,intent(inout) :: rx(*)
    integer, intent(in), optional :: kickrank, nswp, verb, nrestart, niters
    ! verb: 0: off; 1: matlab; 2: console
    character, intent(in), optional :: prec ! 'n', 'l', 'r', 'c'
    integer :: kickrank0, nswp0, verb0, nrestart0, niters0
    character :: prec0
    real(8), intent(in) ::  crY(*), eps, crX0(*)
    real(8), target, intent(in) :: crA(*)
    real(8) eps2
    ! real(8),allocatable :: cr(:,:)
    type(dpoint), allocatable :: cr(:)

    real(8),allocatable,target :: curcr(:)
    ! real(8) curcr2(65536)
    real(8),allocatable :: rhs(:)
    real(8),allocatable :: jacs(:)
    real(8),allocatable :: work(:)
    real(8),allocatable :: tau(:)
    real(8),allocatable :: R(:)
    real(8),allocatable :: res1(:), res2(:)
    type(dpoint), allocatable :: phiA(:), phiy(:)
    real(8),allocatable :: kick_block(:)

    integer,allocatable :: pa(:), py(:)
    integer :: i,j,k, swp, dir, lwork, mm, nn, rnew
    integer info
    real(8) :: err, err_max, res_max, res_new, norm_rhs

    real(8) dnrm2
    real(8) sqrt

    kickrank0 = 4;
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
    prec0 = 'n'
    if (present(prec)) then
       prec0 = prec
    end if
    nrestart0 = 40
    if (present(nrestart)) then
       nrestart0 = nrestart
    end if
    niters0 = 2
    if (present(niters)) then
       niters0 = niters
    end if

    !   print *, n(1:d), ra(1:d+1), ry(1:d+1), rx(1:d+1), nrestart0, niters0, kickrank0, prec0

    eps2 = eps/sqrt(real(d,8))

    allocate(pa(d+1), py(d+1), cr(d), phiA(d+1), phiy(d+1))

    call compute_ps(d,ra,n(1:d)*m(1:d),pa)
    call compute_ps(d,ry,n,py)

    lwork = maxval(rx(1:d+1))*maxval(m(1:d))*maxval(rx(1:d+1))*5
    !   lwork = 256*maxval(m(1:d))*256*5
    allocate(curcr(lwork))
    if (.not.(prec0=='n')) then
       allocate(jacs(lwork))
    end if
    if (prec0=='n') then
       allocate(jacs(1))
    end if
    allocate(rhs(lwork))
    allocate(work(lwork))
    allocate(R(lwork))
    allocate(tau(lwork))
    allocate(kick_block(lwork))


    do i=1,d+1
       allocate(phiA(i)%p(rx(i)*ra(i)*rx(i)), phiy(i)%p(rx(i)*ry(i)))
       !     allocate(phiA(i)%p(256*ra(i)*256), phiy(i)%p(256*ry(i)))
    end do

    mm = 1
    do i=1,d
       allocate(cr(i)%p(rx(i)*m(i)*rx(i+1)))
       !     allocate(cr(i)%p(256*m(i)*256))
       call dcopy(rx(i)*m(i)*rx(i+1), crX0(mm), 1, cr(i)%p, 1)
       mm = mm + rx(i)*m(i)*rx(i+1)
    end do
    phiy(1)%p(1)=1d0
    phiy(d+1)%p(1)=1d0
    phiA(1)%p(1)=1d0
    phiA(d+1)%p(1)=1d0

    allocate(res1(rx(2)*rx(2)*m(1)*ra(2)), res2(rx(2)*rx(2)*m(1)*ra(2)))
    !   allocate(res1(256*256*m(1)*maxval(ra(1:d))), res2(256*256*m(1)*maxval(ra(1:d))))
    ! Initial QR
    dir = 1
    i = 1
    do while (i<d)
       rnew = min(rx(i)*n(i), rx(i+1))
       call dqr(rx(i)*n(i), rx(i+1), cr(i)%p, R, work, lwork, tau)
       call dgemm('N', 'N', rnew, n(i+1)*rx(i+2), rx(i+1), 1d0, R, rnew, cr(i+1)%p, rx(i+1), 0d0, curcr, rnew)
       call dcopy(rnew*n(i+1)*rx(i+2), curcr, 1, cr(i+1)%p,1)
       rx(i+1) = rnew;

       ! size update and phiA
       nn = max(rx(i),rx(i+1))*max(rx(i),rx(i+1))*m(i)*max(ra(i),ra(i+1))
       if (size(res1)<nn) then
          !       print *, size(res1), nn
          deallocate(res1, res2)
          allocate(res1(nn), res2(nn))
       end if
       call dphi_left(rx(i), m(i), rx(i+1), rx(i), m(i), rx(i+1), ra(i), &
       ra(i+1), phiA(i)%p, crA(pa(i)), cr(i)%p, cr(i)%p, phiA(i+1)%p, res1, res2)
       ! size update and phiy
       nn = rx(i)*m(i)*ry(i+1)
       if (size(work)<nn) then
          deallocate(work)
          allocate(work(nn))
       end if
       call dphi2_left(ry(i), ry(i+1), rx(i), m(i), rx(i+1), phiy(i)%p, crY(py(i)), cr(i)%p, phiy(i+1)%p, work)
       i = i+dir
    end do


    ! SWEEPS
    err_max = 0d0
    res_max = 0d0
    i = d;
    dir = -1;
    swp = 1
    do while (swp .le. nswp0)

       ! Generate the RHS
       if (dir<0) then
          if (size(work)<ry(i)*m(i)*rx(i+1)) then
             deallocate(work)
             allocate(work(ry(i)*m(i)*rx(i+1)))
          end if
          call dgemm('N','N', ry(i)*n(i), rx(i+1), ry(i+1), 1d0, crY(py(i)), &
          ry(i)*n(i), phiy(i+1)%p, ry(i+1), 0d0, work, ry(i)*n(i))
          !     print *, 'rhs1'
          ! save for kick
          nn = (rx(i)*ra(i)+ry(i))*n(i)*rx(i+1)
          if (size(kick_block)<nn) then
             deallocate(kick_block)
             allocate(kick_block(nn))
          end if
          call dcopy(ry(i)*m(i)*rx(i+1), work, 1, kick_block, 1)
          !       print *, 'rhs_kick1.5'
          if (size(rhs)<rx(i)*m(i)*rx(i+1)) then
             deallocate(rhs)
             allocate(rhs(rx(i)*m(i)*rx(i+1)))
          end if
          call dgemm('N', 'N', rx(i), n(i)*rx(i+1), ry(i), 1d0, phiy(i)%p, rx(i), work, ry(i), 0d0, rhs, rx(i))
       else
          if (size(work)<rx(i)*m(i)*ry(i+1)) then
             deallocate(work)
             allocate(work(rx(i)*m(i)*ry(i+1)))
          end if
          call dgemm('N', 'N', rx(i), n(i)*ry(i+1), ry(i), 1d0, phiy(i)%p, rx(i), crY(py(i)), ry(i), 0d0, work, rx(i))
          nn = rx(i)*n(i)*(rx(i+1)*ra(i+1)+ry(i+1))
          if (size(kick_block)<nn) then
             deallocate(kick_block)
             allocate(kick_block(nn))
          end if
          call dcopy(rx(i)*m(i)*ry(i+1), work, 1, kick_block, 1)
          if (size(rhs)<rx(i)*m(i)*rx(i+1)) then
             deallocate(rhs)
             allocate(rhs(rx(i)*m(i)*rx(i+1)))
          end if
          call dgemm('N','N', rx(i)*n(i), rx(i+1), ry(i+1), 1d0, work, rx(i)*n(i), phiy(i+1)%p, ry(i+1), 0d0, rhs, rx(i)*n(i))
       end if
       !     print *, 'rhs2'

       ! allocations, initial guess, update, correction, etc
       nn = max(rx(i)*m(i)*ra(i+1)*rx(i+1), rx(i)*ra(i)*n(i)*rx(i+1))
       if (size(res1)<nn) then
          deallocate(res1, res2)
          allocate(res1(nn), res2(nn))
       end if
       nn = rx(i)*m(i)*rx(i+1)
       if (dir<0) then
          nn = nn+(m(i)*rx(i+1)*kickrank0)
       else
          nn = nn+(rx(i)*m(i)*kickrank0)
       end if
       if (size(curcr)<nn) then
          deallocate(curcr)
          allocate(curcr(nn))
       end if
       if (lwork<nn*5) then
          lwork = nn*5
       end if
       if (size(work)<lwork) then ! for svd
          deallocate(work)
          allocate(work(lwork))
       end if
       if (size(tau)<lwork) then ! for svd
          deallocate(tau)
          allocate(tau(lwork))
       end if
       if (size(R)<lwork) then ! for svd
          deallocate(R)
          allocate(R(lwork))
       end if
       nn = rx(i)*m(i)*rx(i+1)

       ! R0
!        call dbfun3(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), ra(i+1), phiA(i)%p, crA(pa(i)), phiA(i+1)%p, cr(i)%p, work, res1, res2)
       call dbfun32(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), ra(i+1), &
       phiA(i)%p, crA(pa(i)), phiA(i+1)%p, cr(i)%p, work)
       !     print *, 'R0, curcr size:', size(curcr)
       call daxpy(nn, -1d0, rhs, 1, work, 1)
       norm_rhs = dnrm2(nn, rhs, 1)
       err = dnrm2(nn, work, 1)/norm_rhs
       !     print *, 'oldres', dnrm2(nn, work, 1), dnrm2(nn, rhs, 1)
       res_max = max(res_max,err)

       tau(1) = eps2/(err*2d0)

       if (tau(1)<1d0) then
          ! Generate PREC, if necessary
          if (prec0=='c') then
             mm = rx(i)*m(i)*m(i)*max(ra(i+1), rx(i+1))
             if (size(jacs)<mm) then
                deallocate(jacs)
                allocate(jacs(mm))
             end if
             call djac_gen(prec0, rx(i), m(i), rx(i+1), ra(i), ra(i+1), phiA(i)%p, crA(pa(i)), phiA(i+1)%p, jacs)
          end if
          if (prec0=='l') then
             mm = rx(i)*rx(i)*m(i)*m(i)*max(ra(i+1), rx(i+1))
             if (size(jacs)<mm) then
                deallocate(jacs)
                allocate(jacs(mm))
             end if
             call djac_gen(prec0, rx(i), m(i), rx(i+1), ra(i), ra(i+1), phiA(i)%p, crA(pa(i)), phiA(i+1)%p, jacs)
          end if
          if (prec0=='r') then
             mm = max(rx(i)*m(i)*m(i)*rx(i+1)*rx(i+1), rx(i)*m(i)*m(i)*ra(i+1))
             if (size(jacs)<mm) then
                deallocate(jacs)
                allocate(jacs(mm))
             end if
             call djac_gen(prec0, rx(i), m(i), rx(i+1), ra(i), ra(i+1), phiA(i)%p, crA(pa(i)), phiA(i+1)%p, jacs)
          end if

          call dgmresr_hh_fort(phiA(i)%p(1), crA(pa(i)), phiA(i+1)%p(1), &
          work(1), rx(i), m(i), rx(i+1), ra(i), ra(i+1), min(nrestart0,nn), &
          tau(1), niters0, prec0, jacs, curcr, max(verb0-1,0))

          ! prec and correction
          if (.not.(prec0=='n')) then
             call djac_apply(prec0, rx(i), n(i), rx(i+1), jacs, curcr, work, tau);
          else
             call dcopy(nn, curcr, 1, work, 1)
          end if
          call daxpy(nn, -1d0, work, 1, cr(i)%p, 1)

!           call dbfun3(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), ra(i+1), phiA(i)%p, crA(pa(i)), phiA(i+1)%p, cr(i)%p, tau, res1, res2)
          call dbfun32(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), &
          ra(i+1), phiA(i)%p, crA(pa(i)), phiA(i+1)%p, cr(i)%p, tau)
          call daxpy(nn, -1d0, rhs, 1, tau, 1)
          res_new = dnrm2(nn, tau, 1)/norm_rhs

          ! dx
          err = dnrm2(nn, work, 1)/dnrm2(nn, cr(i)%p, 1)
          err_max = max(err_max,err)
       else
          res_new = err
          err = 0d0
       end if
       !     print *, 'newres', dnrm2(nn, R, 1), dnrm2(nn, rhs(1:nn)-cr(i)%p(1:nn), 1)

       if (verb0>=2) then ! verb on, carriage_return, no matlab. Otherwise - a lot of shit comes...
         call disp('als_fort: iteration:['// tostring(swp*1d0) // ',' // tostring(i*1d0) &
         // ']  err:' // tostring(err) // '  ry(i):' // tostring(ry(i)*1d0) // '  err_max:' // &
         tostring(err_max) // ' res_max:' // tostring(res_max))
!                 write(*,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3,20X,A$)"), 'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', ry(i), '  err_max:', err_max, 13
          !!      write(*,"(A,I0,A,I0,A,ES10.3,A,I0,A,ES10.3,A,ES10.3)")  'als_fort:  iteration:[', swp, ',', i, ']  err:', err, '  ry(i):', rx(i), '  err_max:', err_max, '  res_max:', res_max
       end if

       ! Rank truncation and kick
       call dcopy(nn, cr(i)%p, 1, curcr, 1)

       if ((dir>0) .and. (i<d)) then
          rnew = min(rx(i+1), rx(i)*n(i))
          if ((kickrank0>0).or.(kickrank0==-1)) then
             call dgesvd('O', 'S', rx(i)*n(i), rx(i+1), curcr, rx(i)*n(i), tau, 0, 1, R, rnew, work, lwork, info)
             ! curcr is U and has the size r1*n,rnew; R is of size rnew,r2
             !         nn = my_chop3(rnew, tau, eps2)
             do j=1,rnew
                call dscal(rx(i+1), tau(j), R(j), rnew)
             end do
             ! residual truncation
             do nn=rnew-1,1,-1
                ! 	  call dcopy(rx(i+1), R(nn), rnew, tau, 1)
                ! 	  call dscal(rx(i+1), 0d0, R(nn), rnew)
                call dgemm('N', 'N', rx(i)*n(i), rx(i+1), nn, 1d0, curcr, rx(i)*n(i), R, rnew, 0d0, work, rx(i)*n(i))
!                 call dbfun3(rx(i),n(i),rx(i+1),rx(i),n(i),rx(i+1),ra(i),ra(i+1),phiA(i)%p, crA(pa(i)), phiA(i+1)%p, work, work, res1, res2)
                call dbfun32(rx(i),n(i),rx(i+1),rx(i),n(i),rx(i+1),ra(i),ra(i+1),phiA(i)%p, crA(pa(i)), phiA(i+1)%p, work, work)
                call daxpy(rx(i)*n(i)*rx(i+1), -1d0, rhs, 1, work, 1)
                err = dnrm2(rx(i)*n(i)*rx(i+1), work, 1)/norm_rhs
                if ((err>res_new*2d0).and.(err>eps2)) then
                   exit
                end if
             end do
             ! 	call dcopy(rx(i+1), tau, 1, R(nn), rnew)
             nn=nn+1
             call drow_cut(rnew, rx(i+1), nn, R)

             rnew = nn
             ! kick
             if (kickrank0>0) then
                mm = max(ry(i)*ra(i)*rx(i), max(ra(i)*m(i),n(i)*ra(i+1))*ry(i)*rx(i+1), ra(i)*n(i)*m(i)*ra(i+1))
                if (size(res1)<mm) then
                   deallocate(res1, res2)
                   allocate(res1(mm), res2(mm))
                end if
                if (size(tau)<rx(i)*n(i)*ra(i+1)*rx(i+1)) then
                   deallocate(tau)
                   allocate(tau(rx(i)*n(i)*ra(i+1)*rx(i+1)))
                end if
                !     call dgemm('N', 'N', rx(i)*n(i), rx(i+1), nn, 1d0, curcr, rx(i)*n(i), R, rnew, 0d0, work, rx(i)*n(i))
                call dbfun3_left(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), &
                ra(i+1), phiA(i)%p, crA(pa(i)), cr(i)%p, tau, res1, res2)
                call dcopy(rx(i)*n(i)*rx(i+1)*ra(i+1), tau, 1, kick_block(rx(i)*n(i)*ry(i+1)+1), 1)
                ! kick_block is of size rx1*ra1+ry1, n*rx2
                call duchol_fort('T', rx(i+1)*ra(i+1)+ry(i+1), rx(i)*n(i), kick_block,  &
                min(kickrank0, rx(i)*n(i)), work, tau, mm)

                call dcopy(rx(i)*n(i)*mm, work, 1, curcr(rx(i)*n(i)*nn+1), 1)

                do j=1,mm*rx(i+1)
                   tau(j)=0d0
                end do
                call drow_add(nn, rx(i+1), mm, R, tau)

                if (size(cr(i)%p)<rx(i)*n(i)*(nn+mm)) then
                   deallocate(cr(i)%p)
                   allocate(cr(i)%p(rx(i)*n(i)*(nn+mm)))
                end if

                rnew = min(rx(i)*n(i), nn+mm)

                call dqr(rx(i)*n(i), nn+mm, curcr, cr(i)%p, work, lwork, tau)
                ! 	print *, 'QR'
                call dgemm('N', 'N', rnew, rx(i+1), nn+mm, 1d0, cr(i)%p, rnew, R, nn+mm, 0d0, work, rnew)
                ! 	print *, 'R*R'
                call dcopy(rnew*rx(i+1), work, 1, R, 1)
             end if
          else
             call dqr(rx(i)*n(i), rx(i+1), curcr, R, work, lwork, tau)
          end if

          call dcopy(rx(i)*n(i)*rnew, curcr, 1, cr(i)%p, 1)
          if (size(curcr)<rnew*n(i+1)*rx(i+2)) then
             deallocate(curcr)
             allocate(curcr(rnew*n(i+1)*rx(i+2)))
          end if
          call dgemm('N','N', rnew, n(i+1)*rx(i+2), rx(i+1), 1d0, R, rnew, cr(i+1)%p, rx(i+1), 0d0, curcr, rnew)
          if (size(cr(i+1)%p)<rnew*n(i+1)*rx(i+2)) then
             deallocate(cr(i+1)%p)
             allocate(cr(i+1)%p(rnew*n(i+1)*rx(i+2)))
          end if
          call dcopy(rnew*n(i+1)*rx(i+2), curcr, 1, cr(i+1)%p, 1)

          rx(i+1) = rnew

          ! New allocations and phis
          nn = max(rx(i)*n(i)*ra(i)*rx(i+1), rx(i)*ra(i+1)*m(i)*rx(i+1), rx(i+1)*rx(i+1)*ra(i+1))
          if (size(res1)<nn) then
             deallocate(res1, res2)
             allocate(res1(nn), res2(nn))
          end if
          nn = rx(i+1)*ra(i+1)*rx(i+1)
          if (size(phiA(i+1)%p)<nn) then
             deallocate(phiA(i+1)%p)
             allocate(phiA(i+1)%p(nn))
          end if
          call dphi_left(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), &
          ra(i+1), phiA(i)%p, crA(pa(i)), cr(i)%p, cr(i)%p, phiA(i+1)%p, res1, res2)

          nn = ry(i)*n(i)*rx(i+1)
          if (size(work)<nn) then
             deallocate(work)
             allocate(work(nn))
          end if
          nn = ry(i+1)*rx(i+1)
          if (size(phiy(i+1)%p)<nn) then
             deallocate(phiy(i+1)%p)
             allocate(phiy(i+1)%p(nn))
          end if
          call dphi2_left(ry(i), ry(i+1), rx(i), n(i), rx(i+1), phiy(i)%p, crY(py(i)), cr(i)%p, phiy(i+1)%p, work)

       end if


       if ((dir<0) .and. (i>1)) then
          rnew = min(n(i)*rx(i+1), rx(i))
          if ((kickrank0>0).or.(kickrank0==-1)) then
             ! 	mm=-1
             ! 	call dgesvd('S', 'O', rx(i), n(i)*rx(i+1), curcr, rx(i), tau, R, rx(i), 0, 1, work, mm, info)
             ! 	mm = work(1)
             ! 	if (lwork<
             ! 	print *, lwork, mm
             call dgesvd('S', 'O', rx(i), n(i)*rx(i+1), curcr, rx(i), tau, R, rx(i), 0, 1, work, lwork, info)
             ! 	print *, 'svd', lwork
             ! curcr is V and has the size rnew,n*r2, but lda=r1; R is U and of size r1,rnew
             do j=1,rnew
                call dscal(rx(i), tau(j), R((j-1)*rx(i)+1), 1)
             end do
             ! residual truncation
             do nn=rnew-1,1,-1
                ! 	  call dcopy(rx(i+1), R(nn), rnew, tau, 1)
                ! 	  call dscal(rx(i+1), 0d0, R(nn), rnew)
                call dgemm('N', 'N', rx(i), n(i)*rx(i+1), nn, 1d0, R, rx(i), curcr, rx(i), 0d0, work, rx(i))
!                 call dbfun3(rx(i),n(i),rx(i+1),rx(i),n(i),rx(i+1),ra(i),ra(i+1),phiA(i)%p, crA(pa(i)), phiA(i+1)%p, work, work, res1, res2)
                call dbfun32(rx(i),n(i),rx(i+1),rx(i),n(i),rx(i+1),ra(i),&
                            ra(i+1),phiA(i)%p, crA(pa(i)), phiA(i+1)%p, work, work)
                call daxpy(rx(i)*n(i)*rx(i+1), -1d0, rhs, 1, work, 1)
                err = dnrm2(rx(i)*n(i)*rx(i+1), work, 1)/norm_rhs
                if ((err>res_new*2d0).and.(err>eps2)) then
                   exit
                end if
             end do
             ! 	call dcopy(rx(i+1), tau, 1, R(nn), rnew)
             nn=nn+1
             !         nn = my_chop3(rnew, tau, eps2)
             call dtransp(rx(i), n(i)*rx(i+1), curcr, work)
             call dcopy(n(i)*rx(i+1)*nn, work, 1, curcr, 1) ! now, n1,ry2,rnew
             call dtransp(rx(i), nn, R, work)
             call dcopy(nn*rx(i), work,1, R,1) ! rnew,r1
             !         do j=1,nn
             !           call dscal(rx(i), tau(j), R(j), nn)
             !         end do
             !         print *, 'curcr, R'
             rnew = nn;

             ! kick
             if (kickrank0>0) then
                if (size(res1)<max(rx(i)*m(i)*ra(i+1)*rx(i+1), rx(i)*ra(i)*n(i)*rx(i+1))) then
                   deallocate(res1, res2)
                   allocate(res1(max(rx(i)*m(i)*ra(i+1)*rx(i+1), rx(i)*ra(i)*n(i)*rx(i+1))), &
                   res2(max(rx(i)*m(i)*ra(i+1)*rx(i+1), rx(i)*ra(i)*n(i)*rx(i+1))))
                end if
                if (size(tau)<rx(i)*ra(i)*n(i)*rx(i+1)) then
                   deallocate(tau)
                   allocate(tau(rx(i)*ra(i)*n(i)*rx(i+1)))
                end if
                ! 	  call dgemm('N', 'N', rx(i), n(i)*rx(i+1), nn, 1d0, R, rx(i), curcr, rx(i), 0d0, work, rx(i))
                call dbfun3_right(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), &
                ra(i+1), crA(pa(i)), phiA(i+1)%p, cr(i)%p, tau, res1, res2)
                call drow_add(ry(i),n(i)*rx(i+1),rx(i)*ra(i),kick_block,tau)
                call duchol_fort('N', rx(i)*ra(i)+ry(i), n(i)*rx(i+1), kick_block,  &
                min(kickrank0, n(i)*rx(i+1)), work, tau, mm);
                call dcopy(n(i)*rx(i+1)*mm, work, 1, curcr(n(i)*rx(i+1)*nn+1), 1)

                do j=1,mm*rx(i)
                   tau(j)=0d0
                end do
                call drow_add(nn, rx(i), mm, R, tau)

                if (size(cr(i)%p)<(nn+mm)*n(i)*rx(i+1)) then
                   deallocate(cr(i)%p)
                   allocate(cr(i)%p((nn+mm)*n(i)*rx(i+1)))
                end if

                rnew = min(n(i)*rx(i+1), nn+mm)

                call dqr(n(i)*rx(i+1), nn+mm, curcr, cr(i)%p, work, lwork, tau)
                ! 	print *, 'QR'
                call dgemm('N', 'N', rnew, rx(i), nn+mm, 1d0, cr(i)%p, rnew, R, nn+mm, 0d0, work, rnew)
                ! 	print *, 'R*R'
                call dcopy(rnew*rx(i), work, 1, R, 1)
             end if
          else
             call dtransp(rx(i), n(i)*rx(i+1), curcr, cr(i)%p)
             call dqr(n(i)*rx(i+1), rx(i), cr(i)%p, R, work, lwork, tau)
             call dcopy(n(i)*rx(i+1)*rx(i), cr(i)%p, 1, curcr, 1)
          end if

          !       call dgemm('T', 'N', rnew, rnew, n(i)*rx(i+1), 1d0, curcr, n(i)*rx(i+1), curcr, n(i)*rx(i+1), 0d0, work, rnew)
          !       print *, work(1:rnew*rnew)

          call dtransp(n(i)*rx(i+1), rnew, curcr, cr(i)%p)
          !       print *, 'cr(i)'
          if (size(curcr)<rx(i-1)*n(i-1)*rnew) then
             deallocate(curcr)
             allocate(curcr(rx(i-1)*n(i-1)*rnew))
          end if
          call dgemm('N','T', rx(i-1)*n(i-1), rnew, rx(i), 1d0, cr(i-1)%p, rx(i-1)*n(i-1), R, rnew, 0d0, curcr, rx(i-1)*n(i-1))
          !       print *, 'cr(i-1)1'
          if (size(cr(i-1)%p)<rx(i-1)*n(i-1)*rnew) then
             deallocate(cr(i-1)%p)
             allocate(cr(i-1)%p(rx(i-1)*n(i-1)*rnew))
          end if
          call dcopy(rx(i-1)*n(i-1)*rnew, curcr, 1, cr(i-1)%p, 1)
          !       print *, 'cr(i-1)2'

          rx(i) = rnew

          ! Phis
          nn = max(rx(i)*n(i)*ra(i+1)*rx(i+1), rx(i)*ra(i)*m(i)*rx(i+1), rx(i)*rx(i)*ra(i))
          if (size(res1)<nn) then
             deallocate(res1, res2)
             allocate(res1(nn), res2(nn))
          end if
          nn = rx(i)*ra(i)*rx(i)
          if (size(phiA(i)%p)<nn) then
             deallocate(phiA(i)%p)
             allocate(phiA(i)%p(nn))
          end if
          call dphi_right(rx(i), m(i), rx(i+1), rx(i), n(i), rx(i+1), ra(i), ra(i+1), phiA(i+1)%p, &
          crA(pa(i)), cr(i)%p, cr(i)%p, phiA(i)%p, res1, res2)
          !       print *, 'phi_r'

          nn = rx(i)*n(i)*ry(i+1)
          if (size(work)<nn) then
             deallocate(work)
             allocate(work(nn))
          end if
          nn = ry(i)*rx(i)
          if (size(phiy(i)%p)<nn) then
             deallocate(phiy(i)%p)
             allocate(phiy(i)%p(nn))
          end if
          call dphi2_right(ry(i), ry(i+1), rx(i), n(i), rx(i+1), phiy(i+1)%p, crY(py(i)), cr(i)%p, phiy(i)%p, work)
          !       print *, 'phi2_r'
       end if


       if ((dir>0) .and. (i==d)) then
          !       call dcopy(rx(i)*n(i)*rx(i+1), curcr, 1, cr(i)%p, 1) ! done at gmres
          if (verb0>=1) then
            call disp('als_fort: iteration:'// tostring(swp*1d0) // '(' // &
            tostring(dir*1d0) // ')  max(ry):' // tostring(maxval(rx(1:d+1))*1d0) &
            // ' err_max:' // tostring(err_max) // ' res_max:' // tostring(res_max))
             !         write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(rx(1:d+1)), '  err_max:', err_max, '  res_max:', res_max
             !!        write(*,"(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3)") 'als_fort:  iteration:', swp, '(', dir, ')  max(ry):', maxval(rx(1:d+1)), '  err_max:', err_max, '  res_max:', res_max

          end if
          dir = -1
          i = i+1
          if (kickrank0<0) then
             exit
          end if
          if (res_max<eps) then
             kickrank0 = -1
          end if
          swp = swp+1
          err_max = 0d0
          res_max = 0d0
          !       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  sweep_reverse:[', swp, ',', i, ',', dir, ']'
       end if

       if ((dir<0) .and. (i==1)) then
          !       call dcopy(rx(i)*n(i)*rx(i+1), curcr, 1, cr(i)%p, 1) ! done at gmres
          if (verb0>=1) then
            call disp('als_fort: iteration:'// tostring(swp*1d0) // '(' // &
            tostring(dir*1d0) // ')  max(ry):' // tostring(maxval(rx(1:d+1))*1d0) &
            // ' err_max:' // tostring(err_max) // ' res_max:' // tostring(res_max))
             !!        write(*,"(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3)") 'als_fort:  iteration:', swp, '(', dir, ')  max(rx):', maxval(rx(1:d+1)), '  err_max:', err_max, '  res_max:', res_max
             !         write(matlab_ist_dumme_kuh,"(A,I0,A,I0,A,I0,A,ES10.3,A,ES10.3$)"), 'als_fort:  iteration:', swp, '(', dir, ')  max(rx):', maxval(rx(1:d+1)), '  err_max:', err_max, '  res_max:', res_max

          end if
          if (kickrank0<0) then
             exit
          end if
          if (res_max<eps) then
             kickrank0 = -1
          end if
          dir = 1
          i = i-1
          err_max = 0d0
          res_max = 0d0
          !       write(*,"(A,I0,A,I0,A,I0,A)"), 'als_fort:  sweep_reverse:[', swp, ',', i, ',', dir, ']'
       end if

       !     print *, 'new core is ready:'
       !     print *, cr(i)%p(1:rx(i)*n(i)*rx(i+1))

       i = i+dir
    end do   ! sweep

    ! copy result_core
    nn = sum(rx(2:d+1)*rx(1:d)*n(1:d));
    allocate(result_core(nn))
    nn = 1;
    do i=1,d
       call dcopy(rx(i)*n(i)*rx(i+1), cr(i)%p, 1, result_core(nn), 1)
       nn = nn+rx(i)*n(i)*rx(i+1)
    end do

    deallocate(work, R, tau, curcr, rhs, kick_block, res1, res2)
    !   if (.not.(prec0=='n')) then
    deallocate(jacs)
    !   end if
    deallocate(pa,py)

    do i=1,d+1
       if (i<=d) then
          deallocate(cr(i)%p)
       end if
       deallocate(phiA(i)%p, phiy(i)%p)
    end do
    deallocate(cr, phiA, phiy)

  end subroutine tt_amen_solve

end module tt_adapt_als
