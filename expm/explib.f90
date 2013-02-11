module explib
 contains 
   subroutine exp_mv(n,m,tau,v,w,tol,anorm,matvec)
     integer, intent(in) :: n,m
     double precision, intent(in) :: tau, tol, anorm
     double precision, intent(in) :: v(n)
     double precision, intent(out) :: w(n)
     integer :: m1
     integer :: ideg = 6
     integer :: lwsp 
     integer :: liwsp 
     double precision, allocatable :: wsp(:)
     integer, allocatable :: iwsp(:)
     integer :: itrace = 0
     integer iflag
     external matvec
     lwsp = (n*(m+1)+n+(m+2)**2+4*(m+2)**2+ideg+1)
     lwsp = max(lwsp,200)
     liwsp = max((m+2),200)
     allocate(wsp(lwsp))
     allocate(iwsp(liwsp))
     m1 = min(m,n) !Fix sizes
     if ( m1 .eq. 0 ) then !The matrix is 1 x 1 :)
        call matvec(1d0,w)
        w(1) = v(1) * exp(w(1) * tau)
     else
         call dgexpv(n,m1,tau,v,w,tol,anorm,wsp,lwsp,iwsp,liwsp,matvec,itrace,iflag)
        if ( iflag .ne. 0 ) then
           print *,'exp_mv failed with iflag=',iflag
        end if
     end if
     deallocate(wsp)
     deallocate(iwsp)
   end subroutine exp_mv
   

   subroutine zexp_mv(n,m,tau,v,w,tol,anorm,matvec)
     implicit none
     integer, intent(in) :: n,m
     double precision, intent(in) :: tau, tol, anorm
     complex(8), intent(in) :: v(n)
     complex(8), intent(out) :: w(n)
     integer :: m1
     integer :: ideg = 6
     integer :: lwsp 
     integer :: liwsp 
     integer, allocatable :: iwsp(:)
     integer :: itrace = 0
     integer iflag
     complex(8), allocatable :: wsp(:)
     external matvec
     lwsp = (n*(m+1)+n+(m+2)**2+4*(m+2)**2+ideg+1)
     lwsp = max(lwsp,200)
     liwsp = max((m+2),200)
     allocate(wsp(lwsp))
     allocate(iwsp(liwsp))
     m1 = min(m,n) !Fix sizes
     if ( m1 .eq. 1 ) then !The matrix is 1 x 1 :)
        call matvec((1d0,0d0),w)
        w(1) = v(1) * exp(w(1) * tau)
     else
        call zgexpv(n,m1,tau,v,w,tol,anorm,wsp,lwsp,iwsp,liwsp,matvec,itrace,iflag)
        if ( iflag .ne. 0 ) then
           print *,'exp_mv failed with iflag=',iflag
        end if
     end if
     deallocate(wsp)
     deallocate(iwsp)
   end subroutine zexp_mv

end module explib
