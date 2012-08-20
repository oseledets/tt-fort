      subroutine matvec(x,y)
      complex(8) :: x(4), y(4)
      !call zcopy(4,x,1,y,1)
      y(:) = x(:)
      end subroutine matvec
      
      
      program main
      use explib
      integer i, n
      complex(8) :: v(4), w(4)
      complex(8) :: zdotc, hij, ONE, ZERO
      real(8) :: tol
      external matvec
      ONE = (1d0,0d0)
      ZERO = (0d0, 0d0)
      do i = 1,4
        v(i) = (0d0, 1d0)
        w(i) = (0d0, 1d0)
      end do 
      n = 4
      !hij = zdotc(n,v,1,w,1)
      !call zgemm('c','n',1,1,n,ONE,v,n,w,n,ZERO,hij,1)
      !print *, n 
      !call temp()
      tol = 1e-8
      call zexp_mv(4,4,1d-2,v,w,tol,1d0,matvec)
      end program main
