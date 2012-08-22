module dmrg_lib
 use tt_lib
 use ttaux_lib
 use ttop_lib
 use rnd_lib
 double precision,parameter,private :: tpi=6.28318530717958647692528676655900577d0
contains 
 subroutine dtt_dmrg(arg, tol, rhs, maxiter, coresize, demo)
  ![tt] fit rhs using dmrg
  implicit none
  type(dtt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  double precision,intent(in) :: rhs(*)
  integer,intent(in),optional :: maxiter,coresize,demo
  character(len=*),parameter :: subnam='dtt_dmrg'
  integer :: l,m,dd,it,itt,maxit,g,nn,nu,nv,ng,ld,rd,mn,lwork,info,i,j,k,ii,jj,kk,bit(64),p
  integer :: idum(tt_size)
  integer,pointer :: r(:),n(:)
  double precision,allocatable :: work(:),u(:,:),v(:,:),sv(:),a(:),d(:)
  double precision :: eps,nrm,err,hx,hy
  double precision,external :: dnrm2
  type(dtt) :: tmp
  character(len=64) :: str
  
  l=arg%l; m=arg%m; if(l.gt.m)return; dd=m-l+1
  if(arg%r(l-1).ne.1)then;write(*,*)subnam,': left border rank should be 1';stop;endif
  if(arg%r(m).ne.1)then;write(*,*)subnam,': right border rank should be 1';stop;endif
  !if(.not.all(arg%n(l:m)==2))then;write(*,*)subnam,': not bitt on input';stop;endif
  r=>arg%r; n=>arg%n
  eps=tol/dsqrt(dble(dd))
  maxit=16;if(present(maxiter))maxit=maxiter
  g=1;if(present(coresize))g=max(coresize-1,1)
  if(m-l.le.g)then
   nn=product(n(l:m))
   call svd(n(l:m),rhs,arg,tol)
   return
  end if
  
  call dtt_rev(arg)
  call ort(arg)
  call dtt_rev(arg)
       
       !nn=product(n(l:m)) 
       !allocate(work(nn))
       !call full(arg,work)
       !nrm=maxval(work)
       !call dscal(nn,1.d0/nrm,work,1)
       !hx=1.d0/nn
       !write(str,'(a,i2.2,a)') '_sin_m',m,'_init'
       !open(10,file=str)
       !do i=1,nn
       ! write(10,'(2(f15.12,1x))')i*hx-hx/2,work(i)
       !end do
       !close(10)
       !STOP

  do it=1,maxit
   write(*,'(a,i2,a)') ' DMRG  --(',it,')-- '
   
   k=l-1
   do itt=-(m-g-l),m-g-l
    if(itt.le.0)then;k=k+1;else;k=k-1;endif
    write(*,'(a,3(i3,a))') 'itt: ',itt,' supercore [',k,':',k+g,']'
    nu=max(1,product(n(l:k-1))); ng=product(n(k:k+g)); nv=max(1,product(n(k+g+1:m)))
    allocate(d(r(k-1)*ng*r(k+g)),u(nu,r(k-1)),v(r(k+g),nv),a(r(k-1)*ng*nv),stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot allocate ';stop;endif
    
    ! form supercore
    if(k.gt.l)then; call dtt_full(arg,u,part=(/l,k-1/)); else; u=1.d0; endif
    if(k+g.lt.m)then; call dtt_full(arg,v,part=(/k+g+1,m/)); else; v=1.d0; endif  
    call dgemm('t','n',r(k-1),ng*nv,nu,1.d0,u,nu,rhs,nu,0.d0,a,r(k-1))
    call dgemm('n','t',r(k-1)*ng,r(k+g),nv,1.d0,a,r(k-1)*ng,v,r(k+g),0.d0,d,r(k-1)*ng)
    deallocate(u,v,a,stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot deallocate';stop;endif

    ! decompose the supercore and recover train format
    if(g.eq.1)then
     ld=r(k-1)*n(k); rd=product(n(k+1:k+g))*r(k+g); mn=min(ld,rd); lwork=128*mn*maxval(r(l:m))
     allocate(u(ld,mn),v(mn,rd),sv(mn),work(lwork), stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate ';stop;endif
     call dgesvd('s','s',ld,rd,d,ld,sv,u,ld,v,mn,work,lwork,info)
     if(info.ne.0)then; write(*,*)subnam,': dgesvd info: ',info; stop; end if
     r(k)=chop(sv,tol=eps)
     deallocate(arg%u(k)%p,arg%u(k+1)%p)
     allocate(arg%u(k)%p(r(k-1),n(k),r(k)),arg%u(k+1)%p(r(k),n(k+1),r(k+1)))
     call dcopy(r(k-1)*n(k)*r(k),u,1,arg%u(k)%p,1)
     call d2submat(r(k),n(k+1)*r(k+1),v,mn,arg%u(k+1)%p)
     if(itt.lt.0)then
      write(*,*) 'ort->'
      forall(i=1:r(k))arg%u(k+1)%p(i,:,:)=sv(i)*arg%u(k+1)%p(i,:,:)
     else 
      write(*,*) '<-ort'
      forall(i=1:r(k))arg%u(k)%p(:,:,i)=sv(i)*arg%u(k)%p(:,:,i)
     endif
     deallocate(u,v,sv,work)
    else
     idum(1)=r(k-1)*n(k);idum(2:g)=n(k+1:k+g-1);idum(g+1)=n(k+g)*r(k+g); tmp%l=1
     call svd(idum(1:g+1),d,tmp,tol=eps)
     if(itt.lt.0)call ort(tmp)
     r(k:k+g-1)=tmp%r(1:g)
     do j=k,k+g
      deallocate(arg%u(j)%p)
      allocate(arg%u(j)%p(r(j-1),n(j),r(j)))
      call dcopy(size(tmp%u(j-k+1)%p),tmp%u(j-k+1)%p,1,arg%u(j)%p,1)
     enddo
     call dealloc(tmp)
    endif 

    deallocate(d,stat=info)
    if(info.ne.0)then;write(*,*)subnam,': cannot deallocate';stop;endif
    call say(arg)
!----
    nn=product(n(l:m)) 
    allocate(work(nn))
    call full(arg,work)
    
    if(present(demo))then
     select case(demo)
      case(1)
       hx=1.d0/nn
       write(str,'(3(a,i2.2))') '_sin_m',m,'_',it,'.',itt+m-g-l
       open(10,file=str)
       do i=1,nn
        write(10,'(2(f15.12,1x))')i*hx-hx/2,work(i)
       end do
       close(10)

      case(2)
       p=2**(m/2); hx=1.d0/p; hy=1.d0/p; allocate(u(p,p))
       do i=1,nn
        do j=0,31
         bit(j)=ibits(i-1,j,1)
        enddo
        ii=1;jj=1
        do kk=0,m/2-1
         ii=ii+bit(2*kk)*(2**kk)
         jj=jj+bit(2*kk+1)*(2**kk)
        enddo
        u(ii,jj)=work(i)
       enddo
       
       write(str,'(3(a,i2.2))') '_simplex_k2_m',m/2,'_',it,'.',itt+m-g-l
       open(10,file=str)
       do ii=0,p
        do jj=0,p
         if(ii.gt.0.and.jj.gt.0)then
          write(10,'(3(f15.12,1x))')ii*hx,jj*hy,u(ii,jj)
         else 
          write(10,'(3(f15.12,1x))')ii*hx,jj*hy,0.d0
         endif 
        end do
        write(10,*)
       end do 
       close(10)
       deallocate(u)
     end select
    end if 
    
    nrm=dnrm2(nn,work,1)
    call daxpy(nn,-1.d0,rhs,1,work,1)
    err=dnrm2(nn,work,1)
    write(*,*) err/nrm
    deallocate(work)
!----
   end do   
   
  end do 
 end subroutine
 

 subroutine dtt_dmrgf(arg, tol, maxiter, coresize, kick, par)
  ![tt] fit rhs using dmrg and maxvol
  use maxvol_lib
  use sort_lib
  use ort_lib
  implicit none
  type(dtt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  integer,intent(in),optional :: maxiter,coresize,kick
  double precision,intent(in),optional :: par(:)
  character(len=*),parameter :: subnam='dtt_dmrgf'
  integer,pointer :: r(:),n(:)
  type(pointi2) :: vip(0:tt_size,2)
  type(pointd2) :: mat(0:tt_size,2),inv(0:tt_size,2)
  type(ttind) :: ind
  double precision,allocatable :: b(:),c(:),work(:),u(:,:),v(:,:),sv(:)
  integer :: ijk,i,j,k,p,pp,q,qq,l,m,mm,dd,d1,nn,n1,t,s,rl,rr,g,mn,ld,rd,lp,rp,dir,lft,rgt,it,maxit,lwork,info
  double precision :: hx,eps,er2,err,nrm,mean,nval
  character*60 :: frm, sdir
  double precision,external :: dnrm2
  logical :: done
  
  l=arg%l; m=arg%m; if(l.gt.m)return
  dd=m-l+1; hx=1.d0/dsqrt(2.d0**dd); d1=dd/3; hx=1.d0/(2.d0**d1); n1=2**d1
  if(arg%r(l-1).ne.1)then;write(*,*)subnam,': left border rank should be 1';stop;endif
  if(arg%r(m).ne.1)then;write(*,*)subnam,': right border rank should be 1';stop;endif
  !if(.not.all(arg%n(l:m)==2))then;write(*,*)subnam,': not bitt on input';stop;endif
  r=>arg%r; n=>arg%n
  maxit=16;if(present(maxiter))maxit=maxiter
  eps=dmax1(tol/(10*dd),1.d-13)
  g=2;if(present(coresize))g=max(coresize,g)
  q=1;if(present(kick))q=max(kick,0)
  write(frm,'(a,i2.2,a)')'(',dd,'i1,f20.14)'
  if(m-l.le.g)then;write(*,*)'tensor smaller than coresize: ',m-l,coresize;return;end if
  
  ! maxvol ->
  lft=1; rgt=2
  allocate(vip(0,lft)%p(r(l-1),1)); forall(i=1:r(l-1))vip(0,lft)%p(i,1)=i
  allocate(mat(0,lft)%p(r(l-1),r(l-1))); call eye(mat(0,lft)%p)
  allocate(inv(0,lft)%p(r(l-1),r(l-1))); call eye(inv(0,lft)%p)
  ! maxvol <-
  call reverse(arg)
  call ort(arg); call maxvol(arg,vip(:,rgt),mat(:,rgt))
  call reverse(arg)
  do i=0,dd-1
   allocate(inv(i,rgt)%p(r(m-i),r(m-i)))
   call matinv(mat(i,rgt)%p,inv(i,rgt)%p)
  enddo
  
  ! iterate
  done=.false.; it=0; nval=0.d0
  do while(.not.done)
   it=it+1; er2=0.d0
   do dir=1,2
    lft=dir;rgt=3-dir
    if(dir.eq.1)then;sdir='>>';else if(dir.eq.2)then;sdir='<<';else;sdir='??';endif
    do p=l,m-g+1
     r=>arg%r; n=>arg%n
     !call say(arg)
     pp=p-l+1; lp=pp-1; rp=dd-lp-g
     rl=r(p-1); rr=r(p+g-1)
     mm=product(n(p:p+g-1))

      !write(*,*)'(a)'
     ! compute long indices and form rhs
     allocate(b(rl*mm*rr),c(rl*mm*rr),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate rhs';stop;endif

!$OMP PARALLEL DO SHARED(b) PRIVATE(i,j,k,ind,t,s)
     do ijk=1,rl*mm*rr
      k=(ijk-1)/(mm*rl)+1
      j=(ijk-1)/rl-(k-1)*mm+1
      i=ijk-(k-1)*mm*rl-(j-1)*rl
      ind%p=0; ind%n(1:dd)=arg%n(l:m); ind%m=dd
      t=i; do s=lp,1,-1; ind%p(s)=vip(s,lft)%p(t,2); t=vip(s,lft)%p(t,1); enddo
      call mindex(j,n(p:p+g-1),ind%p(pp:pp+g-1))
      t=k; do s=rp,1,-1; ind%p(dd-s+1)=vip(s,rgt)%p(t,2); t=vip(s,rgt)%p(t,1); enddo
      if(dir.eq.2)call reverse(ind)
      b(ijk)=dfun(ind,par)
      !write(*,frm)ind%p(1:dd),b(i+rl*(j-1)+rl*mm*(k-1))
     end do
!$OMP END PARALLEL DO
     nval=nval+rl*mm*rr 
      !write(*,*)'(b)'
     
     ! compute supercore
     ld=rl*mm; rd=rr; call dgemm('n','t',ld,rd,rr,1.d0,b,ld,inv(rp,rgt)%p,rr,0.d0,c,ld)
     ld=rl; rd=mm*rr; call dgemm('n','n',ld,rd,rl,1.d0,inv(lp,lft)%p,rl,c,ld,0.d0,b,ld)
      
      !write(*,*)'(c)'
     
     ! check step size
     ld=rl*n(p); rd=n(p+1)*rr; call dcopy(ld*rd,b,1,c,1); nrm=dnrm2(ld*rd,c,1)
     call dgemm('n','n',ld,rd,r(p),-1.d0,arg%u(p)%p,ld,arg%u(p+1)%p,r(p),1.d0,c,ld)
     err=dnrm2(ld*rd,c,1); er2=er2+err*err
      
      !write(*,*)'(d)'
     
     ! decompose the supercore and recover train format
     if(g.eq.1)then
      stop
     else if(g.eq.2)then
      ld=rl*n(p); rd=n(p+1)*rr; mn=min(ld,rd); lwork=128*mn*max(ld,rd)
      allocate(u(ld,mn+q),v(mn,rd),sv(mn),work(lwork), stat=info)
      if(info.ne.0)then;write(*,*)subnam,': cannot allocate ';stop;endif
      call dgesvd('s','s',ld,rd,b,ld,sv,u,ld,v,mn,work,lwork,info)
      if(info.ne.0)then; write(*,*)subnam,': dgesvd info: ',info; stop; end if
      r(p)=chop(sv,tol=eps)
      qq=min(q,ld-r(p))
      if(qq.gt.0)then
       call random(u(:,r(p)+1:r(p)+qq))
       call orto(u(:,1:r(p)),u(:,r(p)+1:r(p)+qq))
      end if 
      deallocate(arg%u(p)%p,arg%u(p+1)%p)
      allocate(arg%u(p)%p(r(p-1),n(p),r(p)+qq),arg%u(p+1)%p(r(p)+qq,n(p+1),r(p+1)))
      call dcopy(r(p-1)*n(p)*(r(p)+qq),u,1,arg%u(p)%p,1)
      forall(i=1:r(p),        j=1:n(p+1),k=1:r(p+1)) arg%u(p+1)%p(i,j,k)=sv(i)*v(i,j+(k-1)*n(p+1))
      forall(i=r(p)+1:r(p)+qq,j=1:n(p+1),k=1:r(p+1)) arg%u(p+1)%p(i,j,k)=0.d0
      r(p)=r(p)+qq
      deallocate(u,v,sv,work)
     else
      stop
     end if
     deallocate(b,c)
     
      !write(*,*)'(e)'

     mean=sumall(arg)/numel(arg)
     !write(*,'(i2,a2,i3,a,f7.3,2(a,e14.7),a,e10.3)') it,sdir,p,' rank ',rank(arg),' mean ',mean,' nrm ',nrm,' err ',err/nrm
     
     ! compute new maxvol positions
     allocate(u(r(p-1)*n(p),r(p)),v(r(p-1)*n(p),r(p)),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate u,v';stop;endif
      
      !write(*,*)'(f)'

     call dgemm('n','n',r(p-1),n(p)*r(p),r(p-1),1.d0,mat(pp-1,lft)%p,r(p-1),arg%u(p)%p,r(p-1), 0.d0,u,r(p-1))
     if(associated(vip(pp,lft)%p))deallocate(vip(pp,lft)%p); allocate(vip(pp,lft)%p(r(p),2)); vip(pp,lft)%p=0
     if(associated(mat(pp,lft)%p))deallocate(mat(pp,lft)%p); allocate(mat(pp,lft)%p(r(p),r(p)))
     if(associated(inv(pp,lft)%p))deallocate(inv(pp,lft)%p); allocate(inv(pp,lft)%p(r(p),r(p)))
      
     call ort0(u,v); call maxvol(v,vip(pp,lft)%p(:,1)); call sort(vip(pp,lft)%p(:,1))
     forall(i=1:r(p),j=1:r(p))mat(pp,lft)%p(i,j)=u(vip(pp,lft)%p(i,1),j)
     call matinv(mat(pp,lft)%p,inv(pp,lft)%p)
     
     vip(pp,lft)%p(:,2)=(vip(pp,lft)%p(:,1)-1)/r(p-1)
     vip(pp,lft)%p(:,1)=vip(pp,lft)%p(:,1)-vip(pp,lft)%p(:,2)*r(p-1)
     vip(pp,lft)%p(:,2)=vip(pp,lft)%p(:,2)+1
     deallocate(u,v)
    end do
    call reverse(arg)
   end do
   mean=sumall(arg)/numel(arg); err=dsqrt(er2); nrm=norm(arg); done=err.le.tol*nrm .or. it.ge.maxit
   write(*,'(i3,a,f7.3,2(a,e14.7),2(a,e10.3))')it,' rank ',rank(arg),' mean ',mean,' nrm ',nrm/dsqrt(numel(arg)),' nval ',nval,' err ',err/nrm
  end do
  call svd(arg,tol)
 end subroutine 
 subroutine ztt_dmrgf(arg, tol, maxiter, coresize, kick, par)
  ![tt] fit rhs using dmrg and maxvol
  use maxvol_lib
  use sort_lib
  use ort_lib
  implicit none
  type(ztt),intent(inout),target :: arg
  double precision,intent(in) :: tol
  integer,intent(in),optional :: maxiter,coresize,kick
  double complex,intent(in),optional :: par(:)
  character(len=*),parameter :: subnam='ztt_dmrgf'
  double complex,parameter :: zero=(0.d0,0.d0),one=(1.d0,0.d0)
  integer,pointer :: r(:),n(:)
  type(pointi2) :: vip(0:tt_size,2)
  type(pointz2) :: mat(0:tt_size,2),inv(0:tt_size,2)
  type(ttind) :: ind
  double complex,allocatable :: b(:),c(:),work(:),u(:,:),v(:,:)
  double precision,allocatable :: sv(:),rwork(:)
  integer :: ijk,i,j,k,p,pp,q,qq,l,m,mm,dd,d1,nn,n1,t,s,rl,rr,g,mn,ld,rd,lp,rp,dir,lft,rgt,it,maxit,lwork,info
  double complex :: mean
  double precision :: hx,eps,er2,err,nrm
  character*60 :: frm, sdir
  double precision,external :: dznrm2
  logical :: done
  
  l=arg%l; m=arg%m; if(l.gt.m)return
  dd=m-l+1
  if(arg%r(l-1).ne.1)then;write(*,*)subnam,': left border rank should be 1';stop;endif
  if(arg%r(m).ne.1)then;write(*,*)subnam,': right border rank should be 1';stop;endif
  !if(.not.all(arg%n(l:m)==2))then;write(*,*)subnam,': not bitt on input';stop;endif
  r=>arg%r; n=>arg%n
  maxit=16;if(present(maxiter))maxit=maxiter
  eps=dmax1(tol/(10*dd),1.d-13)
  g=2;if(present(coresize))g=max(coresize,g)
  q=1;if(present(kick))q=max(kick,0)
  write(frm,'(a,i2.2,a)')'(',dd,'i1,f20.14)'
  if(m-l.le.g)then;write(*,*)'tensor smaller than coresize: ',m-l,coresize;return;end if
   
  ! maxvol ->
  lft=1; rgt=2
  allocate(vip(0,lft)%p(r(l-1),1)); forall(i=1:r(l-1))vip(0,lft)%p(i,1)=i
  allocate(mat(0,lft)%p(r(l-1),r(l-1))); call eye(mat(0,lft)%p)
  allocate(inv(0,lft)%p(r(l-1),r(l-1))); call eye(inv(0,lft)%p)
  ! maxvol <-
  call reverse(arg)
  call ort(arg); call maxvol(arg,vip(:,rgt),mat(:,rgt))
  call reverse(arg)
  do i=0,dd-1
   allocate(inv(i,rgt)%p(r(m-i),r(m-i)))
   call matinv(mat(i,rgt)%p,inv(i,rgt)%p)
  enddo
  
  ! iterate
  done=.false.; it=0
  do while(.not.done)
   it=it+1; er2=0.d0
   do dir=1,2
    lft=dir;rgt=3-dir
    if(dir.eq.1)then;sdir='>>';else if(dir.eq.2)then;sdir='<<';else;sdir='??';endif
    do p=l,m-g+1
     r=>arg%r; n=>arg%n
     !call say(arg)
     pp=p-l+1; lp=pp-1; rp=dd-lp-g
     rl=r(p-1); rr=r(p+g-1)
     mm=product(n(p:p+g-1))
     
     ! compute long indices and form rhs
     allocate(b(rl*mm*rr),c(rl*mm*rr),stat=info)
     if(info.ne.0)then;write(*,*)subnam,': cannot allocate rhs';stop;endif

!$OMP PARALLEL DO SHARED(b) PRIVATE(i,j,k,ind,t,s)
     do ijk=1,rl*mm*rr
      k=(ijk-1)/(mm*rl)+1
      j=(ijk-1)/rl-(k-1)*mm+1
      i=ijk-(k-1)*mm*rl-(j-1)*rl
      ind%p=0; ind%n(1:dd)=arg%n(l:m); ind%m=dd
      t=i; do s=lp,1,-1; ind%p(s)=vip(s,lft)%p(t,2); t=vip(s,lft)%p(t,1); enddo
      call mindex(j,n(p:p+g-1),ind%p(pp:pp+g-1))
      t=k; do s=rp,1,-1; ind%p(dd-s+1)=vip(s,rgt)%p(t,2); t=vip(s,rgt)%p(t,1); enddo
      if(dir.eq.2)call reverse(ind)
      b(ijk)=zfun(ind,par)
      !write(*,frm)ind%p(1:dd),b(i+rl*(j-1)+rl*mm*(k-1))
     end do
!$OMP END PARALLEL DO
    
     ! compute supercore
     ld=rl*mm; rd=rr; call zgemm('n','t',ld,rd,rr,one,b,ld,inv(rp,rgt)%p,rr,zero,c,ld)
     ld=rl; rd=mm*rr; call zgemm('n','n',ld,rd,rl,one,inv(lp,lft)%p,rl,c,ld,zero,b,ld)
     
     ! check step size
     ld=rl*n(p); rd=n(p+1)*rr; call zcopy(ld*rd,b,1,c,1); nrm=dznrm2(ld*rd,c,1)
     call zgemm('n','n',ld,rd,r(p),-one,arg%u(p)%p,ld,arg%u(p+1)%p,r(p),one,c,ld)
     err=dznrm2(ld*rd,c,1); er2=er2+err*err
     
     ! decompose the supercore and recover train format
     if(g.eq.1)then
      stop
     else if(g.eq.2)then
      ld=rl*n(p); rd=n(p+1)*rr; mn=min(ld,rd); lwork=128*mn*max(rl,rr)
      allocate(u(ld,mn+q),v(mn,rd),sv(mn),work(lwork),rwork(8*mn), stat=info)
      if(info.ne.0)then;write(*,*)subnam,': cannot allocate ';stop;endif
      call zgesvd('s','s',ld,rd,b,ld,sv,u,ld,v,mn,work,lwork,rwork,info)
      if(info.ne.0)then; write(*,*)subnam,': zgesvd info: ',info; stop; end if
      r(p)=chop(sv,tol=eps)
      qq=min(q,ld-r(p))
      if(qq.gt.0)then
       call random(u(:,r(p)+1:r(p)+qq))
       call orto(u(:,1:r(p)),u(:,r(p)+1:r(p)+qq))
      end if 
      deallocate(arg%u(p)%p,arg%u(p+1)%p)
      allocate(arg%u(p)%p(r(p-1),n(p),r(p)+qq),arg%u(p+1)%p(r(p)+qq,n(p+1),r(p+1)))
      call zcopy(r(p-1)*n(p)*(r(p)+qq),u,1,arg%u(p)%p,1)
      forall(i=1:r(p),        j=1:n(p+1),k=1:r(p+1)) arg%u(p+1)%p(i,j,k)=sv(i)*v(i,j+(k-1)*n(p+1))
      forall(i=r(p)+1:r(p)+qq,j=1:n(p+1),k=1:r(p+1)) arg%u(p+1)%p(i,j,k)=zero
      r(p)=r(p)+qq
      deallocate(u,v,sv,work,rwork)
     else
      stop
     end if
     deallocate(b,c)

     mean=sumall(arg)/numel(arg)
     !write(*,'(i2,a2,i3,a,f7.3,2(a,e14.7),a,e10.3)') it,sdir,p,' rank ',rank(arg),' mean ',cdabs(mean),' nrm ',nrm,' err ',err/nrm
    
     ! compute new maxvol positions
     allocate(u(r(p-1)*n(p),r(p)),v(r(p-1)*n(p),r(p)))
     call zgemm('n','n',r(p-1),n(p)*r(p),r(p-1),one,mat(pp-1,lft)%p,r(p-1),arg%u(p)%p,r(p-1), zero,u,r(p-1))
     if(associated(vip(pp,lft)%p))deallocate(vip(pp,lft)%p); allocate(vip(pp,lft)%p(r(p),2)); vip(pp,lft)%p=0
     if(associated(mat(pp,lft)%p))deallocate(mat(pp,lft)%p); allocate(mat(pp,lft)%p(r(p),r(p)))
     if(associated(inv(pp,lft)%p))deallocate(inv(pp,lft)%p); allocate(inv(pp,lft)%p(r(p),r(p)))
     call ort0(u,v); call maxvol(v,vip(pp,lft)%p(:,1)); call sort(vip(pp,lft)%p(:,1))
     forall(i=1:r(p),j=1:r(p))mat(pp,lft)%p(i,j)=u(vip(pp,lft)%p(i,1),j)
     call matinv(mat(pp,lft)%p,inv(pp,lft)%p)
    
     vip(pp,lft)%p(:,2)=(vip(pp,lft)%p(:,1)-1)/r(p-1)
     vip(pp,lft)%p(:,1)=vip(pp,lft)%p(:,1)-vip(pp,lft)%p(:,2)*r(p-1)
     vip(pp,lft)%p(:,2)=vip(pp,lft)%p(:,2)+1
     deallocate(u,v)
    end do
    call reverse(arg)
   end do
   mean=sumall(arg)/numel(arg); err=dsqrt(er2); nrm=norm(arg); done=err.le.tol*nrm .or. it.ge.maxit
   write(*,'(i3,a,f7.3,a,e14.7,a,e10.3)')it,' rank ',rank(arg),' nrm ',nrm/dsqrt(numel(arg)),' err ',err/nrm
  end do
  call svd(arg,tol)
 end subroutine 
 


end module
