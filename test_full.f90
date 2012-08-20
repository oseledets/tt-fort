      program main
      use dyn_tt
      integer :: d, nn, i 
      integer, allocatable :: ry(:), ra(:), n(:)
      double precision, allocatable :: cra(:), cry(:)
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
      call random_number(cry(1:nn))
      !do i = 1,nn
      !   cry(i) = 1d0
      !end do 
      nn = 0
      do i = 1,d
         nn = nn + ra(i)*n(i)*n(i)*ra(i+1)
      end do
      allocate(cra(nn))
      call random_number(cra(1:nn))
      !do i = 1,d
      !  cra(i) = 1d0
      !end do 
      call  tt_kls(d,n,n,ra,cra,cry,ry,1d-2,150,5,20, 0)
      deallocate(n)
      deallocate(ra)
      deallocate(ry)
      end program main
