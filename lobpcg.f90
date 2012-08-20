module blopex
! See http://arxiv.org/pdf/0705.2626.pdf
  use dispmodule
contains
 subroutine get_sbm(n,m,ind,x,y)
 integer, intent(in) :: n,m
 integer, intent(in) :: ind(*)
 double precision, intent(in) :: x(n,*)
 double precision, intent(inout) :: y(n,*)
 integer :: i,j
 do i = 1,n
    do j = 1,m
       y(i,j) = x(i,ind(j))
    end do
 end do
 end subroutine get_sbm
 subroutine dzero(n,x)
   integer, intent(in) :: n
   double precision, intent(inout) :: x(*)
   integer :: i
   do i = 1,n
      x(i) = 0d0
   end do
 end subroutine dzero

 subroutine qr(n,m, A, R)
   integer, intent(in):: n,m
   double precision, intent(inout) :: A(n,m)
   double precision, intent(out) :: R(m,min(n,m))
   double precision, allocatable :: work(:)
   double precision :: Q1(n,m)
   double precision ::  work1(1), tau(m)
   
   integer lwork
   integer info
   integer rnew, k,j
   call dcopy(n*m,A,1,Q1,1)
   call dgeqrf(n, m, Q1, n, tau, work1,-1,info)
   lwork = work1(1)
   allocate(work(lwork))
   call dgeqrf(n, m, Q1, n, tau, work,lwork,info)
   deallocate(work)
   if (info.ne.0) then
      print *, 'qr: dgeqrf failed'
   end if
   rnew = min(n,m)
   R(:,:)=0d0
   do j=1,m
      R(1:min(j,n),j)=Q1(1:min(j,n),j)
   end do
   
   call dorgqr(n,rnew,rnew,Q1,n,tau,work1,-1,info)
   lwork = work1(1)
   allocate(work(lwork))
   call dorgqr(n,rnew,rnew,Q1,n,tau,work,lwork,info)
   deallocate(work)
   if (info.ne.0) then
      print *, 'qr: dorgqr failed'
   end if
   call dcopy(n*rnew,Q1,1,A,1)
end subroutine

! Some special lobpcg-based stuff that we will need to use

 subroutine compute_gramA(n,m,mcur,lambda,X,AW,AP,W,P,gramA,gsize)
   implicit none 
   integer, intent(in) :: n,m,mcur,gsize
   double precision, intent(in) :: lambda(3*m)
   double precision, intent(in) :: X(n,m), AW(n,mcur), AP(n,mcur), W(n,mcur), P(n,mcur)
   double precision, intent(inout) :: gramA(:,:)
   integer :: i
   double precision tmp_mxp(m,mcur)
   double precision tmp_pxp(mcur,mcur)
   call dzero(gsize*gsize,gramA)
   do i = 1, m
      gramA(i,i) = lambda(i)
   end do
   call dgemm('t','n',m,mcur,n,1d0,X,n,AW,n,0d0,tmp_mxp,m)
   gramA(1:m,m+1:m+mcur) = tmp_mxp(1:m,1:mcur)
   call dgemm('t','n',m,mcur,n,1d0,X,n,AP,n,0d0,tmp_mxp,m)
   gramA(1:m,m+mcur+1:m+2*mcur) = tmp_mxp(1:m,1:mcur)
   call dgemm('t','n',mcur,mcur,n,1d0,W,n,AW,n,0d0,tmp_pxp,mcur)
   gramA(m+1:m+mcur,m+1:m+mcur) = tmp_pxp(1:mcur,1:mcur)
   call dgemm('t','n',mcur,mcur,n,1d0,W,n,AP,n,0d0,tmp_pxp,mcur)
   gramA(m+1:m+mcur,m+mcur+1:m+2*mcur) = tmp_pxp(1:mcur,1:mcur)
   call dgemm('t','n',mcur,mcur,n,1d0,P,n,AP,n,0d0,tmp_pxp,mcur)
   gramA(m+mcur+1:m+2*mcur,m+mcur+1:m+2*mcur) = tmp_pxp(1:mcur,1:mcur)

   
 end subroutine compute_gramA

 subroutine compute_gramB(n,m,mcur,X,W,P,gramB,Gsize)
 implicit none
 integer, intent(in) :: n,m,mcur, Gsize
 double precision, intent(in) :: X(n,m), W(n,mcur), P(n,mcur)
 double precision, intent(inout) :: gramB(Gsize,Gsize)
 double precision :: tmp1(m,mcur), tmp2(mcur,mcur)
 integer i
 call dzero(Gsize*Gsize,gramB)
 do i = 1,Gsize
    gramB(i,i) = 1d0
 end do
 call dgemm('t','n',m,mcur,n,1d0,X,n,W,n,0d0,tmp1,m)
 gramB(1:m,m+1:m+mcur) = tmp1(1:m,1:mcur)

 call dgemm('t','n',m,mcur,n,1d0,X,n,P,n,0d0,tmp1,m)
 gramB(1:m,m+mcur+1:m+2*mcur) = tmp1(1:m,1:mcur)

call dgemm('t','n',mcur,mcur,n,1d0,W,n,P,n,0d0,tmp2,m)
 gramB(m+1:m+mcur,m+mcur+1:m+2*mcur) = tmp2(1:mcur,1:mcur)
 
 end subroutine compute_gramB


 subroutine lobpcg2(n,m,avec,bvec,tvec,X,maxiter,tol,lambda)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: maxiter
    integer k,max_restarts
    double precision, intent(in) :: tol
    double precision, intent(inout) :: X(n,m)
    double precision, intent(out) :: lambda(3*m)
    integer info 
    external avec, bvec, tvec
    max_restarts = 5
    info = 1
    do k = 1,max_restarts
       call lobpcg(n,m,avec,bvec,tvec,X,maxiter,tol,lambda,info)
       if ( info .eq. 0 ) then
          return 
       end if
    end do
    end subroutine lobpcg2

 subroutine lobpcg(n,m, avec,bvec,tvec,X,maxiter,tol,lambda,info)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: maxiter
    double precision, intent(in) :: tol
    double precision, intent(inout) :: X(n,m)
    double precision, intent(out) :: lambda(3*m)
    integer, intent(out) :: info
    integer :: lwork,i,mcur,p0,k, ind(m)
    logical restart_flag
    double precision, allocatable :: work(:)
    double precision dnrm2, tmp4(3*m,3*m)
    !Memory
    double precision :: W(n,m), TW(n,m), P(n,m), Q(n,m), AX(n,m), AW(n,m), TMP0(n,m), TMP2(n,m)
    !These are for the selected vectors
    double precision :: P_loc(n,m), AP_loc(n,m), vec(n)
    double precision :: AP(n,m), BX(n,m), BW(n,m), BP(n,m), TMP1(m)
    double precision :: work1(1)
    double precision, allocatable :: tmp_pxp(:,:), tmp_mxp(:,:),R(:,:), sol(:,:), tmp_pxm(:,:)
    double precision, allocatable :: gramA(:,:), gramB(:,:)
    external avec, bvec, tvec
    call dzero(n*m,P)
    !B-orthonormalize B
    allocate(tmp_pxp(m,m))
    call qr(n,m,X,tmp_pxp)
    !do i = 1,n
    !   print *,X(i,1), X(i,2)
    !end do 
    !return

    !call bvec(n,m,X,BX)
    !call dgemm('t','n',m,m,n,1d0,X,n,BX,n,0d0,tmp_pxp,m)
    !call dpotrf('l',m,tmp_pxp,m,info)
    !if ( info .ne. 0 ) then
    !   print *,'lobpcg: dpotrf failed with info=',info
    !   return
    !end if
    !do i = 1, n
    !   
    !   call dcopy(m,X(i,1),n,TMP1,1)
    !   call dtrtrs('l','n','n',m,1,tmp_pxp,m,TMP1,m,info)
    !   call dcopy(m,TMP1,1,X(i,1),n)
    !   if ( info .ne. 0 ) then
    !      print *,'lobpcg: dtrtrs failed with info=',info
    !      return
    !   end if
    !   call dcopy(m,BX(i,1),n,TMP1,1)
    !   call dtrtrs('l','n','n',m,1,tmp_pxp,m,TMP1,m,info)
    !   call dcopy(m,TMP1,1,BX(i,1),n)
    !   if ( info .ne. 0 ) then
    !      print *,'lobpcg: dtrtrs failed with info=',info
    !      return
    !   end if
    !end do 

    ! X = X * R^{-1}
    mcur = m
    deallocate(tmp_pxp)

       allocate(tmp_mxp(m,mcur))
       allocate(tmp_pxp(mcur,mcur))
       allocate(tmp_pxm(mcur,m))

    restart_flag = .false.
    call avec(n,m,X,AX)
    
    call dgemm('t','n',m,m,n,1d0,X,n,AX,n,0.d0,tmp_pxp,m)
    
    !Compute the initial Ritz vector
    !Query workspace
    call dsyev('v','l',m,tmp_pxp,m,lambda,work1,-1,info)
    lwork = work1(1)
    allocate(work(lwork))
    call dsyev('v','l',m,tmp_pxp,m,lambda,work,lwork,info)
    deallocate(work)
    call dgemm('n','n',n,m,m,1d0,X,n,tmp_pxp,m,0d0,TMP0,n)

    call dcopy(n*m,TMP0,1,X,1)

    call dgemm('n','n',n,m,m,1d0,AX,n,tmp_pxp,m,0d0,TMP0,n)
    call dcopy(n*m,TMP0,1,AX,1)

    if ( info .ne. 0 ) then
       print *,'Lobpcg, dysev failed with info=',info
       return
    end if

    do k = 1,m
       ind(k) = k
    end do

    do k = 1, maxiter
       print *,21
       !write(*,'(AI3)') 'Iteration:',k
       print *,22,'mcur=',mcur, 'ind=',ind
       !call get_sbm(n,mcur,ind,X,TMP0)
       !call avec(n,mcur,TMP0,AX)
       p0 = 0
       do i = 1, mcur
          vec(1:n) = AX(1:n,ind(i)) - X(1:n,ind(i)) * lambda(i)
          print *,dnrm2(n,vec,1)
          if ( dnrm2(n,vec,1) > tol ) then
            ind(p0+1) = ind(i) !Hope it will not kill ind
            p0 = p0 + 1
            W(:,p0) = vec(:)
            if ( k > 1 ) then
               P_loc(:,p0) = P(:,ind(i)) 
               AP_loc(:,p0) = AP(:,ind(i))
            end if
          end if
       end do

       mcur = p0 !New number of active constaints
       if ( mcur == 0 ) then
          write(*,'(A,1X,I3,1X)') 'Done in', k,'iterations' 
          print *,'Eigenvalues:',lambda(1:m)
          return
       end if


       call tvec(n,mcur,W,TW)
       call dcopy(n*mcur,TW,1,W,1)
       !Orthogonalize W
       call qr(n,mcur,W,tmp_pxp)
       
       

       call avec(n,mcur,W,AW)
       if ( k > 1 ) then
          !B-ort of P_loc
          call qr(n,mcur,P_loc,tmp_pxp)
          !Temporary fix
          call avec(n,mcur,P_loc,AP_loc) 
       end if

       !Rayleigh-Ritz procedure
       if ( k > 1  ) then
          allocate(gramA(m+2*mcur,m+2*mcur))
          allocate(gramB(m+2*mcur,m+2*mcur))
          print *,'bA'
          call compute_gramA(n,m,mcur,lambda,X,AW,AP_loc,W,P_loc,gramA,m+2*mcur)
          call disp(gramA)
          print *,'aa'
          print *,'bB'
          call compute_gramB(n,m,mcur,X,W,P_loc,gramB,m+2*mcur)
          call dsyev('n','u',m+2*mcur, gramB, m+2*mcur, vec, work1,-1,info)
          lwork = work1(1)
          allocate(work(lwork))
          call dsyev('n','u',m+2*mcur, gramB, m+2*mcur, vec, work,lwork,info)
          deallocate(work)
          call compute_gramB(n,m,mcur,X,W,P_loc,gramB,m+2*mcur)

          print *,'Gram B:'
          call disp(gramB)
          print *,'Eigenvalues:'
          call disp(vec(1:m+2*mcur))
          !call dcopy((m+2*mcur)*(m+2*mcur),gramB,1,tmp4,1)
          
          !call dpstrf('u',m+2*mcur,
          !call dpotrf('u',m+2*mcur,gramB,m+2*mcur,info)
          
          print *,'Matrix A:'
          call disp(gramA)
          !What it was before:
          call dsygv(1,'v','u',m+2*mcur,gramA,m+2*mcur,gramB,m+2*mcur,lambda, work1, -1, info)
          lwork = work1(1)
          allocate(work(lwork))
          call dsygv(1,'v','u',m+2*mcur,gramA,m+2*mcur,gramB,m+2*mcur,lambda, work, lwork, info)
          deallocate(work)

          print *,lambda(1:m)
          if ( info .ne. 0 ) then
             !Just restart
             !restart_flag = .true.
             print *,'restarting'
             deallocate(gramA)
             deallocate(gramB)
             deallocate(tmp_pxp)
             deallocate(tmp_mxp)
             deallocate(tmp_pxm)
             return
             !print *,'Generalized eigenproblem in lobpcg failed with info = ', info
             !print *,'mcur = ',mcur,'m=',m
             !print *,'X:'
             !call disp(X(1:n,1:2))
             !print *,'W:'
             !call disp(W(1:n,1:mcur))
             !print *,'P:'
             !call disp(P_loc(1:n,1:mcur))
             !deallocate(tmp_pxp)
             !deallocate(tmp_mxp)
             !deallocate(tmp_pxm)
             !call compute_gramB(n,m,mcur,X,W,P_loc,lambda,m+2*mcur)
             !print *,'HERE!'
             !lambda(:) = 0d0
             !do i = 1,9
             !   lambda(i) = tmp3(i)
             !end do
             !call dcopy(9,tmp3,1,lambda,1)
             !call dgemm('t','n',mcur,mcur,n,1d0,P_loc,n,P_loc,n,0d0,tmp_pxp,mcur)
             !print *,tmp_pxp(1:mcur,1:mcur)
             !return
        end if
          
       else
          allocate(gramA(m+mcur,m+mcur))
          allocate(gramB(m+mcur,m+mcur))
          gramA(:,:) = 0d0

          do i = 1, m
             gramA(i,i) = lambda(i)
          end do

          call dgemm('t','n',m,mcur,n,1d0,X,n,AW,n,0d0,tmp_mxp,m)

          gramA(1:m,m+1:m+mcur) = tmp_mxp(1:m,1:mcur)
          call dgemm('t','n',mcur,mcur,n,1d0,W,n,AW,n,0d0,tmp_pxp,mcur)
          gramA(m+1:m+mcur,m+1:m+mcur) = tmp_pxp(1:mcur,1:mcur)
          gramB(:,:) = 0d0
          do i = 1,m+mcur
             gramB(i,i) = 1d0
          end do
          call dgemm('t','n',m,mcur,n,1d0,X,n,W,n,0d0,tmp_mxp,m)
          gramB(1:m,m+1:m+mcur) = tmp_mxp(1:m,1:mcur)
          call dsygv(1,'v','u',m+mcur,gramA,m+mcur,gramB,m+mcur,lambda, work1, -1, info)
          lwork = work1(1)
          allocate(work(lwork))
          call dsygv(1,'v','u',m+mcur,gramA,m+mcur,gramB,m+mcur,lambda, work, lwork, info)
          deallocate(work)
          if ( info .ne. 0 ) then
             print *,'Generalized eigenproblem in lobpcg failed!'
          end if
       end if
       print *,7
       if ( k > 1 ) then
          call dgemm('n','n',n,m,mcur,1d0,W,n,gramA(m+1:m+mcur,1:m),mcur,0d0,TMP0,n)
          call dgemm('n','n',n,m,mcur,1d0,P_loc,n,gramA(m+mcur+1:m+2*mcur,1:m),mcur,1d0,TMP0,n)
          P(:,:) = TMP0(:,:)
          call dgemm('n','n',n,m,mcur,1d0,AW,n,gramA(m+1:m+mcur,1:m),mcur,0d0,TMP0,n)
          tmp_pxm(1:mcur,1:m) = gramA(m+mcur+1:m+2*mcur,1:m)
          call dgemm('n','n',n,m,mcur,1d0,AP_loc,n,tmp_pxm,mcur,1d0,TMP0,n)

          AP(:,:) = TMP0(:,:)


       else 
          call dgemm('n','n',n,m,mcur,1d0,W,n,gramA(m+1:m+mcur,1:m),mcur,0d0,P,n)
          call dgemm('n','n',n,m,mcur,1d0,AW,n,gramA(m+1:m+mcur,1:m),mcur,0d0,AP,n)
       endif
       call dgemm('n','n',n,m,m,1d0,X,n,gramA(1:m,1:m),m,0d0,TMP0,n)
       X(:,:) = P(:,:) + TMP0(:,:)
       print *,17
       call dgemm('n','n',n,m,m,1d0,AX,n,gramA(1:m,1:m),m,0d0,TMP0,n)
       print *,18
       AX(:,:) = AP(:,:)! + TMP0(:,:)
       deallocate(gramA)
       deallocate(gramB)

    end do
       deallocate(tmp_pxm)
       deallocate(tmp_pxp)
       deallocate(tmp_mxp)


  end subroutine lobpcg
end module blopex
