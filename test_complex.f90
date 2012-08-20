      program main
      use dyn_tt
      integer :: d, nn, i 
      integer, allocatable :: ry(:), ra(:), n(:)
      complex(8), allocatable :: cra(:), cry(:)
      real(8) :: tmpx, tmpy
      INTEGER :: i_, n_, clock_
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed_

        CALL RANDOM_SEED(size = n_)
        ALLOCATE(seed_(n_))
        CALL SYSTEM_CLOCK(COUNT=clock_)
        seed_ = clock_ + 37 * (/ (i_ - 1, i_ = 1, n_) /)
        CALL RANDOM_SEED(PUT = seed_)
        DEALLOCATE(seed_)
      d = 5
      allocate(n(d))
      allocate(ra(d+1))
      allocate(ry(d+1))
      ra(1) = 1
      ra(d+1) = 1
      ry(1) = 1
      ry(d+1) = 1
      do i = 1,d
        n(i) = 2
      end do 
      do i = 2,d
        ry(i) = 2
        ra(i) = 2
      end do
      nn = 0
      do i = 1,d
         nn = nn + ry(i)*n(i)*ry(i+1)
      end do 
      allocate(cry(nn))
      do i = 1, nn
         call random_number(tmpx)
         call random_number(tmpy)
         cry(i) = tmpx + (0d0,1d0)*tmpy
      end do 
      !do i = 1,nn
      !   cry(i) = 1d0
      !end do 
      nn = 0
      do i = 1,d
         nn = nn + ra(i)*n(i)*n(i)*ra(i+1)
      end do
      allocate(cra(nn))
      do i = 1, nn
         call random_number(tmpx)
         call random_number(tmpy)
         cra(i) = tmpx + (0d0,1d0)*tmpy
      end do 
      !do i = 1,d
      !  cra(i) = 1d0
      !end do 
      call  ztt_kls(d,n,n,ra,cra,cry,ry,1d-2,150,5,20, 0)
      deallocate(n)
      deallocate(ra)
      deallocate(ry)
      end program main
