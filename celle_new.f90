parameter(k=3000000)
dimension x(k),y(k),q(k),iflag(k),ip(k),ngr(200)
double precision d,dr,x,y,p,dr1,dr2,dd,df,xr,yr,xx,yy

p=3.1415926535879/180.
p=p*1.d+0

open(80,file='nz.dat')
open(90,file='nz_bmap_500_new.dat')

neve=500

kk=0
do i=1,3000000
  read(80,*,end=99)ia,m,ig,io,mi,s,y(i),x(i),pp,q(i)   !,id   t,x(i),y(i),z,q(i)   !
  iflag(i)=0
end do
99 npt=i-1

print*,npt


ke=1
ntot=0
do ic=1,100000000
  ieve=int(rand()*float(npt))
  if(iflag(ieve).ne.0)cycle
  kk=kk+1
  if(kk==6626)exit
  xr=x(ieve)*p
  yr=y(ieve)*p
  d=10.d+0
  do iv=1,500
    do j=1,k
      ip(j)=0
    end do
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
        if(ne>k)print*,ne,'*****'
        ip(ne)=j
      end if
    end do
    if(ne==0)exit
    dm=ds/float(ne)
    s=dm*.1
    if(ne==1)s=5.
    if(neve-50<=ne.and.ne<=neve+50)then
      dd=d
      d=10.d+0
      exit
    end if
    if(ne<neve)d=d+s !1.d-2
    if(ne>neve)d=d-s  !1.d-2
  end do
  if(ne<neve-50.or.ne>neve+50)cycle
  do j=1,200
    ngr(j)=0
  end do
  qmax=-20.
  do j=1,ne
    if(q(ip(j))>qmax)qmax=q(ip(j))
    ll=int(q(ip(j))*10.)
    if(ll>200)print*,ll,'ll'
    ngr(ll)=ngr(ll)+1
  end do
  nmax=0
  do j=1,200
    if(ngr(j)>nmax)then
      nmax=ngr(j)
      jmax=j
    end if
  end do
  qc=jmax*.1
  if(qmax-qc<2.)cycle
  print*,ic,iv,ne,ke,kk,ntot
  do j=1,ne
    iflag(ip(j))=ke
  end do
  ke=ke+1
  ntot=ntot+ne

  if(float(npt-ntot)/float(npt)<.01)exit
end do
print*,npt-ntot,float(npt-ntot)/float(npt)

do i=1,npt
  write(90,'(2(e10.4,1x),i4,1x,f3.1)')x(i),y(i),iflag(i),q(i)
end do

end
