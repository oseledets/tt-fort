program test_eigb_f
 use tt_lib
 use ttaux_lib
 use ttop_lib
 use ttio_lib
!  use ttm_lib
 use time_lib
 use tt_block_eig
 implicit none

 double precision,parameter :: tpi=6.28318530717958647692528676655900577d0

 type(dtt) :: H,x
 integer :: i, r1,j,r2, d, rmax, B, sz, info, Bmax, max_full_size
 integer, allocatable :: n(:), ra(:), rx(:)
 double precision :: tol,t1,t2
 double precision, allocatable :: crX(:), crA(:), theta(:), theta_ex(:), theta_lam(:)
 character(len=40) :: fnam

  INTEGER :: i_, n_, clock_
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed_

  CALL RANDOM_SEED(size = n_)
  ALLOCATE(seed_(n_))
  CALL SYSTEM_CLOCK(COUNT=clock_)
  seed_ = clock_ + 37 * (/ (i_ - 1, i_ = 1, n_) /)
  CALL RANDOM_SEED(PUT = seed_)
  DEALLOCATE(seed_)


  open(10,file='test_eigb.inp')
  read (10,*) B
  read (10,*) Bmax
  read (10,*) rmax
  read (10,*) max_full_size
  read (10,*) tol
  read (10,*) fnam
!  B = 2
!  rmax = 200
!  tol = 1.0d-5
  close(10)

  call OMP_SET_NUM_THREADS(2)

!  read(*,*) j, np
 call dtt_read(H,fnam,info)
 if(info.ne.0)then;write(*,*) 'read fails: ',info;stop;endif
 call say(H)

!  call dtt_read(x,'Heisen_x0.sdv',info)
!  if(info.ne.0)then;write(*,*) 'read fails: ',info;stop;endif
!  call say(x)

 d = H%m - H%l + 1
 print *, 'd=', d

 sz = sum(H%r(0:d-1)*H%n(1:d)*H%r(1:d))
 allocate(crA(sz))
 sz = 0
 do i=1,d
  forall(r1=1:H%r(i-1),j=1:H%n(i),r2=1:H%r(i)) crA(sz+r1+(j-1)*H%r(i-1)+(r2-1)*H%r(i-1)*H%n(i))=H%u(i)%p(r1,j,r2)
  sz=sz+H%r(i-1)*H%n(i)*H%r(i)
 end do

!  sz = sum(x%r(0:d-1)*x%n(1:d)*x%r(1:d))
!  allocate(crX(sz))
!  sz = 0
!  do i=1,d
!   forall(r1=1:x%r(i-1),j=1:x%n(i),r2=1:x%r(i)) crX(sz+r1+(j-1)*x%r(i-1)+(r2-1)*x%r(i-1)*x%n(i))=x%u(i)%p(r1,j,r2)
!   sz=sz+x%r(i-1)*x%n(i)*x%r(i)
!  end do

 allocate (n(d), ra(d+1), rx(d+1))
 ra(1) = 1
 ra(d+1) = 1
 do i=1,d
   n(i) = nint(sqrt(real(H%n(i))))
   ra(i+1) = H%r(i)
   rx(i+1) = B
 end do
 rx(1) = 1

 sz = n(1)*B + sum(n(2:d))*B*B
 allocate(crX(sz))

 do i=1,sz
   call random_number(crX(i))
 end do

 allocate(theta(B))

 t1=timef()
 call tt_eigb(d,n,n,ra,crA, crX, rx, tol, rmax, theta, B, max_full_size=max_full_size,nswp=4)
 t2=timef()-t1


 allocate(theta_lam(4),theta_ex(1+d*(d+3)/2+d*(d-1)*(d-2)/6))
 do i=1,4
   theta_lam(i) = 4*(sin(tpi*dble(i)/(4*(n(1)+1))))**2
 end do
 theta_ex(1)=d*theta_lam(1)
 forall(i=2:d+1) theta_ex(i)=(d-1)*theta_lam(1) + theta_lam(2)
 forall(i=2+d:1+d*(d+1)/2) theta_ex(i)=(d-2)*theta_lam(1) + 2*theta_lam(2)
 forall(i=2+d*(d+1)/2:1+d*(d+3)/2)theta_ex(i)=(d-1)*theta_lam(1) + 1*theta_lam(3)
 forall(i=2+d*(d+3)/2:1+d*(d+3)/2+d*(d-1)*(d-2)/6)theta_ex(i)=(d-3)*theta_lam(1) + 3*theta_lam(2)


 write(*,*) 'Energy levels (compared to Laplace)'
 do i=1,B; write(*,'(i4,1x,e20.12,1x,e12.3)') i-1,theta(i),theta_ex(i)-theta(i); enddo
 print *, 'Wall Time = ', t2
 print *, 'Max Rank = ', maxval(rx(1:d+1))

 open(10,file='test_eigb.out')
 write(10,'(I0)') n(1)
 write(10,'(I0)') d
 write(10,'(I0)') B
 write(10,'(ES8.1)') tol
 write(10,'(ES12.5)') t2
 write(10,'(I0)') maxval(rx(1:d+1))
 do i=1,B
  write(10,*) theta(i)
 end do
 do i=B+1,Bmax
   write(10,'(I0)') 0
 end do
 close(10)

 call deallocate_result

!  m=128; n=np*m
!  write(fnam,'(a,i4.4)')'ubiq_out.',j
!  open(10,file=fnam)
!  do i=1,m
!   ind=i +(j-1)*m
!   freq=200.d0+2800.d0 * (ind-1)/(n-1)
!   LL=H2
!   call dtt_axpy(2.d0*tpi*freq,HH,1.d0,LL)
!   call dtt_axpy(((tpi*freq)**2) + (damp**2),id,1.d0,LL)
!
!   if(i.eq.1) xx=rho0
!   t1=timef()
!   call dtt_amen_solve(LL,rho0,1.d-3,xx, kickrank=5,nswp=100,local_prec='l',local_iters=2,local_restart=40,max_full_size=10000,trunc_norm=1,verb=1)
!   t2=timef()
!   fid=damp*dot_product(coil,xx)
!   write(*,'(a,i5,a,i5,a,e12.5,a,f6.1,a,e10.3,a,f9.2)')'[',ind,'/',n,'] freq ',freq,' rank ',rank(xx),' fid',fid,' time ',t2-t1
!
!   write(10,'(e15.7,3x,e15.7,3x,f9.2)')freq,fid,t2-t1
!  end do
!  close(10)

 deallocate (crA, crX, ra, rx, n, theta, theta_ex)
end