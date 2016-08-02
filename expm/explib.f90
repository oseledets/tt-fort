module explib
 contains 
   subroutine dexp_mv(n, m, tau, v, w, tol, anorm, matvec, verb)
     integer, intent(in) :: n,m
     double precision, intent(in) :: tau, tol, anorm
     double precision, intent(in) :: v(n)
     double precision, intent(out) :: w(n)
     integer, intent(in), optional :: verb
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
     if ( present(verb) ) then 
        write(*, '(a, I5)') 'matvecs: ', iwsp(1) 
     end if
     deallocate(wsp)
     deallocate(iwsp)
 end subroutine dexp_mv
   

   subroutine zexp_mv(n, m, tau, v, w, tol, anorm, matvec, verb)
     ! Thin layer around EXPOKIT package
     ! Input arguments:
     !
     ! integer n - the problem size
     !
     ! integer m - maximum size of the Krylov basis
     !
     ! real(8) tau - the time step
     !
     ! complexl(8) v(n), input - the starting vector
     !
     ! complex(8) w(n), output - the result
     !
     ! real(8) tol ---  solution tolerance
     !
     ! real(8) anorm --- norm estimate 
     !
     ! external matvec --- matrix-by-vector procedure, call matvec(x,y)
     !
     ! integer verb --- optional, verbosity level


     implicit none
     integer, intent(in) :: n, m
     double precision, intent(in) :: tau, tol, anorm
     integer, intent(in), optional :: verb
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
        call matvec((1d0,0d0), w)
        w(1) = v(1) * exp(w(1) * tau)
     else
        call zgexpv(n, m1, tau, v, w, tol, anorm, wsp, &
                   lwsp, iwsp, liwsp, matvec, itrace, iflag)
        if ( iflag .ne. 0 ) then
           print *,'exp_mv failed with iflag=',iflag
        end if
     end if
     if ( present(verb) ) then 
        write(*, '(a, I5)') 'matvecs: ', iwsp(1) 
        print *,'time covered:', wsp(8)
        print *,'error:', wsp(6) 
     end if
     deallocate(wsp)
     deallocate(iwsp)
   end subroutine zexp_mv

end module explib
