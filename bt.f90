module bt_lib
 implicit none
! y''(t) + b D_*^{3/2} y(t) + c y(t) = f(t); y(0)=p;y'(0)=q
contains

 pure double precision function dbt(b,c,p,q,t,t1,t2,f1,f2) result (ft)
 !pure double precision function dbt(b,c,p,q,t,t1,t2,t3) result (ft)
  implicit none
  double precision,intent(in) :: b,c,p,q,t, t1,t2,f1,f2!,t3
  integer :: m,n,d,i,j,k,info
  double precision :: a,ht,aq0,tk
  double precision,allocatable :: y(:,:),f(:),y0(:),aq(:),bq(:),rhs(:)
  double precision,allocatable :: mat(:,:),hmat(:,:),mmat(:,:)
 
  m=4;d=10;n=2**d; a=.5d0; ht=t/n; !b=1.d0/2;c=1.d0/2;p=0.d0;q=0.d0
  allocate(y(0:m-1,1:n),y0(0:m-1),rhs(0:m-1),f(0:n),aq(0:n),bq(0:n),mat(m,m),hmat(m,m),mmat(m,m),stat=info)
  !if(info.ne.0)then;write(*,*)'cannot allocate';stop;endif
  
  y0(0)=p;y0(1)=0.d0;y0(2)=q;y0(3)=0.d0
  do k=0,n;tk=k*ht+ht/2;f(k)=0.d0;if(tk.le.1.d0)f(k)=8.d0;if(dabs(tk-t1).le..25d0)f(k)=f1;if(dabs(tk-t2).le..25d0)f(k)=f2;enddo
  !forall(k=0:n)f(k)=dbt_force(k*ht+ht/2)
  forall(k=0:n)aq(k)=((ht**a)/(gamma(a)*a*(a+1)))*(dble(k+2)**(a+1)-2*dble(k+1)**(a+1)+dble(k)**(a+1))
  forall(k=0:n)bq(k)=((ht**a)/(gamma(a)*a))*(dble(k+1)**a-dble(k)**a)

  aq0=ht**a/(gamma(a)*a*(a+1))
  mat=0.d0; mat(1,1)=1.d0; mat(2,2)=1.d0; mat(3,3)=1.d0; mat(4,4)=1.d0
  hmat=0.d0; hmat(4,4)=-b*aq0;hmat(4,1)=-c*aq0; hmat(1,2)=aq0; hmat(2,3)=aq0; hmat(3,4)=aq0
  mmat=hmat
  do i=1,16
   mat=mat+mmat
   mmat=matmul(mmat,hmat)
  end do

  do k=0,n-1
   aq0=(dble(k)**(a+1)-(k-a)*(dble(k+1)**a))*(ht**a)/(gamma(a)*a*(a+1))
   rhs(0)=y0(0)+aq0*y0(1)
   rhs(1)=y0(1)+aq0*y0(2)
   rhs(2)=y0(2)+aq0*y0(3)
   rhs(3)=y0(3)+aq0*(-c*y0(0)-b*y0(3))

   do j=1,k
    rhs(0)=rhs(0)+aq(k-j)*y(1,j)
    rhs(1)=rhs(1)+aq(k-j)*y(2,j)
    rhs(2)=rhs(2)+aq(k-j)*y(3,j)
    rhs(3)=rhs(3)+aq(k-j)*(-c*y(0,j)-b*y(3,j))
   end do
   do j=0,k
    rhs(3)=rhs(3)+bq(k-j)*f(j)
   end do
   y(:,k+1)=matmul(mat,rhs)
  end do
  ft=y(0,n)
  deallocate(y,y0,rhs,f,aq,bq,mat,mmat,hmat)
  return
 end function

 pure double precision function dbt_force(t) result (f)
  implicit none
  double precision,intent(in) :: t
  !f=dcos(t)/(t+1)
  !f=0.d0; if(t.le.1.d0)f=8.d0 ! podlubny
  f=t+1
 end function
 
 pure double precision function dbt_ttoep2(x) result (f)
  implicit none
  ! f(x)= (4/3) * [(x-1)^{3/2} - 2x^{3/2} + (x+1)^{3/2}]
  double precision,intent(in) :: x
  integer :: i,n
  double precision :: df,err,top,bot
  f=0.d0
  if(x.le.0.d0)return
  if(x.lt.4.d0)then
   top=(x-1)**1.5d0 + (x+1)**1.5d0
   bot=2*(x**1.5d0)
   f=4*(top-bot)/3
   err=dlog(1.e15*(top-bot)/top)/dlog(10.d0)
  else 
   err=0.d0;f=1.d0;df=f
   n=int(8.d0*dlog(10.d0)/dlog(x))+2
   do i=1,n
    df=df*(1.d0-2.5d0/(2*i+1))*(1.d0-2.5d0/(2*i+2))/(x*x)
    f=f+df
   end do
   f=f/dsqrt(x)
  end if
  !if(present(snr))snr=err
  return
 end function
end module
