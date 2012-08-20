module quad_lib
 implicit none
contains
 subroutine quad_rinv1(n,q)
  implicit none
  ! 1/t = \sum w_p exp(-a_p t^2), t >= 1
  ! w_p = 2 h c_p / sqrt(pi) (1+exp(-s_p)); a_p = log^2(1+exp(s_p))
  ! c_p = cosh t_p, s_p = sinh t_p, t_p = h*(p-n/2-1), p=1,n
  ! [Hackbusch, Khoromskij, 2006 part 1, (5.3)]
  integer,intent(inout) :: n
  double precision :: q(2,n)
  character(len=*),parameter :: subnam='quad_rinv1'
  double precision,parameter :: tpi=6.28318530717958647692528676655900577d0
  integer :: nq,i,m
  double precision :: s,c,es,loghuge,wq,aq,ht,t
  
  t=0.d0
  loghuge=dlog(huge(t))
  nq=(n-3)/2;m=1
  ht=dlog(tpi*nq)/nq
  q=0.d0
  do i=-nq,nq
   t=dble(i)*ht; s=dsinh(t); c=dcosh(t)
   if(s.lt.-loghuge)then
    !q(1,1)=q(1,1)+2*c*ht/dsqrt(tpi/2)
    !q(1,1)=0.d0
   ! write(*,*) i, ' small'
   else if(s.gt.loghuge)then
   ! write(*,*) i, ' big'
   else
    m=m+1
    q(1,m)=2*c*ht/(dsqrt(tpi/2)*(1.d0+dexp(-s)))
    q(2,m)=(dlog(1.d0+dexp(s)))**2
   end if 
   !write(*,'(2e25.15)') q(:,m)
  end do
  n=m
 end subroutine

 subroutine testquad_rinv(nq,q,a,b,n,err)
  implicit none
  integer,intent(in) :: nq
  double precision,intent(in) :: q(2,nq)
  double precision,intent(in) :: a,b
  integer,intent(in) :: n
  double precision,intent(out),optional :: err
  character(len=*),parameter :: subnam='testquad_rinv'
  double precision :: t,x,y,lt,lx,ly,val,err1,errm
  integer :: i,j
  if(a.le.0 .or. b.le.0)then; write(*,*)subnam,': illegal interval: ',a,b;stop;endif
  if(b.lt.1)then;write(*,*)subnam,': too few step numbers: ',n;stop;endif
  x=dmin1(a,b);y=dmax1(a,b);lx=dlog(x);ly=dlog(y); errm=0.d0
  open(10,file='_testquad.rinv',access='sequential')
  do i=0,n-1
   lt=lx+dble(i)*(ly-lx)/(n-1)
   t=dexp(lt)
   val=0.d0
   do j=1,nq
    val=val+q(1,j)*dexp(-q(2,j)*t*t)
   end do
   err1=t*dabs(1.d0/t-val)
   errm=dmax1(err1,errm)
   write(10,'(e12.5,3x,2f25.15,3x,e10.2)')t,1.d0/t,val,err1
  end do
  close(10)
  if(present(err))err=errm
 end subroutine
end module
