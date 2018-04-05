 program main
 use dyn_tt
 integer :: d, rmax, nswp, verb, kickrank, Asize, Ysize, ntrials
 real(8) :: tau
 integer, allocatable :: n(:), m(:), ra(:), ry(:)
 complex(8), allocatable ::  crA(:), crY(:)
 double precision, allocatable ::  drA(:), drY(:)
 integer :: i
    open(unit=10,status='old',file='test_ksl.dat',form='unformatted',access='stream')
    !open(unit=10,status='old',file='test_eye_ksl.dat',form='unformatted',access='stream')
    read(10) d
    print *,'d=',d
    allocate(n(d))
    allocate(m(d))
    read(10) n(1:d)
    read(10) m(1:d)
    allocate(ra(d+1))
    allocate(ry(d+1))
    read(10) ra(1:d+1)
    read(10) ry(1:d+1)

    read(10) Asize
    allocate(crA(Asize))
    read(10) crA(1:Asize)

    allocate(drA(Asize))
    drA(:) = 0.0
    drA(:) = -1.0
    
    read(10) Ysize
    allocate(crY(Ysize))
    read(10) crY(1:Ysize)
    read(10) tau,rmax,kickrank,nswp,verb
    close(10)

    allocate(drY(Ysize))
    drY(:) = 1.0
    
    !Test if we read all correctly
    print *,'n=',n(1:d)
    print *,'m=',m(1:d)
    print *,'ra=',ra(1:d+1)
    print *,'ry=',ry(1:d+1)
    print *,'tau=',tau, 'rmax=',rmax,'kickrank=',kickrank,'nswp=',nswp,'verb=',verb

    ntrials = 10
    do while ( ntrials > 0 )
       call ztt_ksl(d,n,m,ra,crA, crY, ry, tau, rmax, kickrank, nswp, verb)
       ntrials = ntrials - 1
    end do
    ! Test real code as well
    ntrials = 10
    do while ( ntrials > 0 )
       call tt_ksl(d,n,m,ra,drA, drY, ry, tau, rmax, kickrank, nswp, verb)
       ntrials = ntrials - 1
    end do
    
    !open(unit=10,status='replace',file='test_ksl.dat',form='unformatted')
    !write(10) d,n,m,ra,ry,pa(d+1)-1,crA(1:(pa(d+1)-1)),mm-1,crY0(1:mm-1),tau,rmax,kickrank,nswp,verb 
    !close(10)
    !return
 end program
