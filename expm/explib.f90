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
     !TMP SECTION
     double precision, allocatable :: tmp(:,:),vtmp(:)
     integer :: k
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
     complex(8), allocatable :: wsp(:)
     !TMP SECTION
     integer :: k
     integer, allocatable :: iwsp(:)
     integer :: itrace = 0
     integer iflag
     !Test area
     integer i,j, iv, j1v, mh, ih
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
     !   call zcopy( n, v,1, w,1 )
     !   mh = m1 + 2
     !   iv = 1
     !   j1v = 5
     !   ih = iv + n*(m+1) + n
     !   do i = 1,n
     !      wsp(iv + i-1) = w(i)
     !   enddo
     !   do i = 1,mh*mh
     !      wsp(ih+i-1) = ZERO
     !   enddo
        !call matvec( wsp(j1v-n), wsp(j1v) )
        !call zcopy(n,v,1,w,1)
        !call matvec(v,w)
        !print *,'old core:', v(1:2)
        call zgexpv(n,m1,tau,v,w,tol,anorm,wsp,lwsp,iwsp,liwsp,matvec,itrace,iflag)
        !print *,'new core:', w(1:2)
        !print *,'true sol:', exp(tau)*v(1:2)
        !read(*,*)
        if ( iflag .ne. 0 ) then
           print *,'exp_mv failed with iflag=',iflag
        end if
     end if
     deallocate(wsp)
     deallocate(iwsp)
   end subroutine zexp_mv

end module explib
