! ifort -O2 test_ksl_cme.f90 tt_ksl.o expokit.o explib.o normest.o dlacn1.o zlacn1.o dlarpc.o dlapst.o mytt.a -o test_ksl_cme -mkl

program test_eigb_f
 use tt_lib
 use ttaux_lib
 use ttop_lib
 use ttio_lib
 use time_lib
 use dyn_tt
 implicit none

 double precision,parameter :: tpi=6.28318530717958647692528676655900577d0

 type(ztt) :: A,u
 integer :: i, r1,j,r2, d, sz, info, Nt
 integer, allocatable :: n(:), ra(:), ru(:)
 double precision :: tol,t1,t2, tau, T
 double complex, allocatable :: crU(:), crA(:)
 character(len=40) :: Anam, unam, unam_res

  INTEGER :: i_, n_, clock_
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed_

  CALL RANDOM_SEED(size = n_)
  ALLOCATE(seed_(n_))
  CALL SYSTEM_CLOCK(COUNT=clock_)
  seed_ = clock_ + 37 * (/ (i_ - 1, i_ = 1, n_) /)
  CALL RANDOM_SEED(PUT = seed_)
  DEALLOCATE(seed_)


  open(10,file='test_ksl_cme.inp')
  read (10,*) Anam
  read (10,*) unam
  read (10,*) unam_res
  read (10,*) T
  read (10,*) Nt
  close(10)

  tau = T/Nt

  call OMP_SET_NUM_THREADS(1)

! Read dtt data
  call ztt_read(A,Anam,info)
  if(info.ne.0)then;write(*,*) 'read fails: ',info;stop;endif
  call say(A)

  call ztt_read(u,unam,info)
  if(info.ne.0)then;write(*,*) 'read fails: ',info;stop;endif
  call say(u)


  d = u%m - u%l + 1
  print *, 'd=', d

! Convert into the linear storage
  sz = sum(A%r(0:d-1)*A%n(1:d)*A%r(1:d))
  allocate(crA(sz))
  sz = 0
  do i=1,d
    forall(r1=1:A%r(i-1),j=1:A%n(i),r2=1:A%r(i)) crA(sz+r1+(j-1)*A%r(i-1)+(r2-1)*A%r(i-1)*A%n(i))=A%u(i)%p(r1,j,r2)
    sz=sz+A%r(i-1)*A%n(i)*A%r(i)
  end do

  sz = sum(u%r(0:d-1)*u%n(1:d)*u%r(1:d))
  allocate(crU(sz))
  sz = 0
  do i=1,d
    forall(r1=1:u%r(i-1),j=1:u%n(i),r2=1:u%r(i)) crU(sz+r1+(j-1)*u%r(i-1)+(r2-1)*u%r(i-1)*u%n(i))=u%u(i)%p(r1,j,r2)
    sz=sz+u%r(i-1)*u%n(i)*u%r(i)
  end do

 allocate (n(d), ra(d+1), ru(d+1))
 ra(1) = 1
 ra(d+1) = 1
 do i=1,d
   n(i) = nint(sqrt(real(A%n(i))))
   ra(i+1) = A%r(i)
   ru(i+1) = u%r(i)
 end do
 ru(1) = 1


! Run the solver and measure the time
 t1=timef()
 do i=1,Nt
  call ztt_ksl(d,n,n,ra,crA, crU, ru, tau, 500)
!   call tt_ksl(d,n,n,ra,crA, crU, ru, tau, 500)
  call zcopy(sum(ru(1:d)*n(1:d)*ru(2:d+1)), zresult_core, 1, crU, 1)
  print *, 'step ', i, 'done'
 end do
 t2=timef()-t1


 print *, ru(1:d+1)
 print *, 'Wall Time = ', t2

  open(10,file='test_ksl_cme.out')
  do i=1,sum(ru(1:d)*n(1:d)*ru(2:d+1))
  write (10,*) zresult_core(i)
  end do
  close(10)

 ! Cast result_core back into dtt and save
  call dealloc(u)
  u%l = 1
  u%m = d
  u%r(0)=1
  do i=1,d
    u%n(i)=n(i)
    u%r(i)=ru(i+1)
  end do
  call alloc(u)
  sz = 0
  do i=1,d
    forall(r1=1:u%r(i-1),j=1:u%n(i),r2=1:u%r(i)) u%u(i)%p(r1,j,r2) = zresult_core(sz+r1+(j-1)*u%r(i-1)+(r2-1)*u%r(i-1)*u%n(i))
    sz=sz+u%r(i-1)*u%n(i)*u%r(i)
  end do

! deallocate result_core
 call deallocate_result

 call ztt_write(u,unam_res,info)
 if(info.ne.0)then;write(*,*) 'write fails: ',info;stop;endif

 call dealloc(u)

 deallocate (crA, crU, ra, ru, n)
end