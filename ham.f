       implicit none
       integer NDAT,N1,N2,nn1,nn2,i,mm,k,l,ii,j,lk,jk,ll,jj
       parameter (N1=33,N2=33,NDAT=N1*N2)
       double precision umu,D0,alpha,xe,dx,xmax,xmin,dk,T,term,pot1,KE,s
     &,hbar,x(N1),zpi,grid_length,y(N2),ymax,ymin,dy,dk1,xx(NDAT),
     &xy(NDAT),grid_length2,KE1,T1,KE2,s1,pot2,V1(NDAT)

       double precision H0(N1,N2,N1,N2),V(NDAT)

       umu=916.675d0
       D0=0.1744d0
       alpha=1.02764d0
       xe=1.40201d0
       zpi=4.0d0*atan(1.0d0) 
       nn1=(N1-1)/2.0d0
       nn2=(N2-1)/2.0d0
       hbar=1.0d0
       xmax=2.0d0
       xmin=0.0d0
       dx=(xmax-xmin)/(N1)
       grid_length=N1*dx 
       dk=(2.0d0*zpi)/grid_length
	ymax=2.0d0
       ymin=0.0d0
       dy=(ymax-ymin)/(N2)
       grid_length2=N2*dy
       dk1=(2.0d0*zpi)/grid_length2       

       do k=1,N1
       x(k)=xmin+k*dx
       y(k)=ymin+k*dy
       end do

       k=1
       i=1
23     j=1
33     xx(k)=x(i)
       xy(k)=y(j)
    
       if (j .lt.N2)then
       k=k+1
       j=j+1
       go to 33
       end if

       if (i .lt. N1)then
       k=k+1
       i=i+1
       go to 23
       end if
       
       do mm=0,N1-1
       do ii=0,N1-1
       k=1+mm
       l=1+ii
       do ll=0,N2-1
       do jj=0,N2-1
       lk=1+ll
       jk=1+jj

       V(k)=D0*((1-dexp(-alpha*(x(k)-xe)))**2) 
       V1(lk)=D0*((1-dexp(-alpha*(y(lk)-xe)))**2) 

       s=0.0d0
       do i=1,nn1
       T=(2.0d0/umu)*(((hbar*zpi*i)/grid_length)**2)
       term=dcos((i*2.0d0*zpi*(mm-ii))/N1)*T
       s=s+term
       end do
       KE1=(2.0d0*s)/N1

       s1=0.0d0
       do i=1,nn2
       T1=(2.0d0/umu)*(((hbar*zpi*i)/grid_length2)**2)
       term=dcos((i*2.0d0*zpi*(ll-jj))/N2)*T1
       s1=s1+term
       end do
       KE2=(2.0d0*s1)/N2

       if (k .eq. l)then
       pot1=V(k)
       else
       pot1=0.0d0
       end if

       if (lk .eq. jk)then
       pot2=V1(lk)
       else
       pot2=0.0d0
       end if

       H0(k,l,lk,jk)=KE1+KE2+pot1+pot2
       write(231,*) H0(k,l,lk,jk)
       end do
       end do
       end do
       end do 

       return
       end
             
