!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The program celle subdivide the space of catalogue in cell containing each one neve events !!!!!!!!!
!!!! 		   It looks for the largest earthquake j not yet assigned to a cell             !!!!!!!!!
!!!!	   Then looks for all the events around the j_th event within a circle of radius d      !!!!!!!!!
!!!!        If the number of events is larger than neve, d is reduced by a quantity delta       !!!!!!!!!
!!!!      If the number of events is smaller than neve, d is increased by a quantity delta      !!!!!!!!!
!!!!  delta is adjusted on the base of the average interdistance between the evnts in the cell  !!!!!!!!!
!!!!                     Each cell is charcterized by its index ke			        !!!!!!!!!
!!!!      Each event in the catalogue is marked with the cell index through the vector iflag    !!!!!!!!!
!!!!  The procedure continues up to when all the events in the catalogue are assignd to a cell  !!!!!!!!!
!!!!            For reasons of speedness the procedure have a tolerance tol                     !!!!!!!!!
!!!!  The exit of the program is a new catalogue of events with their coordinates 		!!!!!!!!!
!!!!                       the cell index and the magnitude					!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



parameter(k=3000000)
dimension x(k),y(k),q(k),iflag(k)
double precision d,dr,x,y,p,dr1,dr2,dd,df,xr,yr,xx,yy
character*30 a(k),nome

p=3.1415926535879/180.
p=p*1.d+0

open(80,file='california.dat')
open(90,file='celle_cali.dat')

!!!! sets the starting parameters and the number of events per cell  !!!!

neve=2000
delta=.1
ntol=50
tol=.01

!!!!  reads the data from the opened file 
!!!!  sets iflag(i) and find the maximum event in the catalogue 

qmax=-10.
do i=1,npt !3000000
  read(80,'(a25,1x,i7,2x,f8.5,1x,f10.5,1x,f7.3,2x,f4.2)',end=99)ia,m,ig,io,mi,s,x(i),y(i),z,q(i)   
  iflag(i)=0
  if(q(i)>qmax)then
    imax=i
    qmax=q(i)
  end if
end do
99 npt=i-1

iff=0
ke=0
ntot=0
do ic=1,100000000
  ieve=imax
  if(iff>0)print*,iff,ne
  ke=ke+1
  iflag(ieve)=ke
  xr=x(ieve)*p
  yr=y(ieve)*p
  d=10.d+0   !!!! d is the initial radius !!!!!
  
!!!! looks for the events in a cell of radius d around the imax event
  
  do iv=1,200
    ne=0
    ds=0.
    do j=1,npt
      if(iflag(j).ne.0)cycle
      xx=x(j)*p
      yy=y(j)*p
      df=abs(yy-yr)
      dr=dacos(dsin(xx)*dsin(xr)+dcos(xx)*dcos(xr)*dcos(df))
      dr=dr*6370.d+0
      if(dr<d)then
        ds=ds+dr 
        ne=ne+1
      end if
    end do
    if(ne==0)then
      d=2*d
      cycle
    end if
    dm=ds/float(ne)
    s=dm*delta
    if(ne==1)s=5.
    if(neve-ntol<=ne.and.ne<=neve+ntol)then   !!!!  exits from the loop when ne=neve+-ntol
      dd=d
      d=10.d+0
      delta=.1
      exit
    end if
    if(ne<neve)d=d+s 
    if(ne>neve)d=d-s  
  end do
  
!!!! This is for speedness, if the search of the optimal radius is too much long the program moves to the next largest event
  
  if(iff>10)then   
    qmax=0.
    do j=1,npt
      if(iflag(j).ne.0)cycle
      if(q(j)>qmax)then
        qmax=q(j)
        imax=j
      end if
    end do
    iff=0
    delta=.1
    if(ll==7)print*,imax
    cycle
  end if

  if(ne<10.and.iff>5)then
    iff=0
    exit
  end if
  if(ne<neve-50.or.ne>neve+50)then 
    iff=iff+1
    delta=delta*.5
    cycle
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! assigns all the events to the cell !!!!!

  l=0
  do j=1,npt
    if(iflag(j).ne.0)cycle
    xx=x(j)*p
    yy=y(j)*p
    df=abs(yy-yr)
    dr=dacos(dsin(xx)*dsin(xr)+dcos(xx)*dcos(xr)*dcos(df))
    dr=dr*6370.D+0
    if(dr<=dd)then
      iflag(j)=ke
      l=l+1
    end if
  end do

  ntot=ntot+ne


  if(float(npt-ntot)/float(npt)<tol)exit !!!! if all the events are assigned to a cell with a tolerance tol ends the program

!!!!  looks for the nex largest event

  qmax=-10.
  do j=1,npt
    if(iflag(j).ne.0)cycle
    if(q(j)>qmax)then
      qmax=q(j)
      imax=j
    end if
  end do
  iflag(imax)=ke
end do

!!!!!  writes the data to the opened file

do i=1,npt
  if(iflag(i)==99999)iflag(i)=0
  write(90,'(2(e12.6,1x),i4,1x,e10.4)')x(i),y(i),iflag(i),q(i)
end do

end
