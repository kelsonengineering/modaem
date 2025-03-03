c***********************************************************************
      subroutine contur (zp,idmz,nx,ny,rgrx1,rgry1,rdx,rdy,
     .                   tour,ison,iOutLU)
c     contur is adapted from subroutine cntour (c) 1976
c     university of minnesota
c
c  zp        2-d real array  containing z-values to be contoured
c  idmz      integer first dimension of array zp  idimz>nx
c  nx        integer number of grid points in x direction
c  ny        integer number of grid points in y direction
c  rgrx1,rgry1 coordinates of the lower-left corner
c  rdx       horizontal grid spacing.
c  rdy       vertical grid spacing.
c  tour      real value of contour level to be plotted
c  ison      1-d integer array used for temporary storage of
c            flag bits. must be dimensined to a minimum
c            size of 1+(nx*ny-1)/31

      integer*4 idmz,nx,ny,ison,ipen,ix,iy,irtn,jrtn,iOutLU
      logical lbtest
      dimension ison(*),zp(idmz,ny)

c  clear switches
      j=1+(nx*ny-1)/31
      do 80 i=1,j
 80     ison(i)=0

c  Search left edge.
 90   i=1
 91   ii=ny-i
      if(zp(1,ii).le.tour) goto 160
      j=ii+1
      if(zp(1,j) .gt.tour) goto 160
      call toura(1,ii,1,j,rgrx1,rgry1,rdx,rdy,tour,idmz,nx,ny,zp,ison,
     &           iOutLU)
 160  i=i+1
      if(i.lt.ny) goto 91

c  Search bottom edge.
      i=2
 96   if (zp(i,1).le.tour) goto 100
      j=i-1
      if (zp(j,1).gt.tour) goto 100
      call toura(i,1,j,1,rgrx1,rgry1,rdx,rdy,tour,idmz,nx,ny,zp,ison,
     &           iOutLU)
 100  i=i+1
      if(i.le.nx) goto 96

c  Search right edge.
      i=2
 116  if (zp(nx,i).le.tour) goto 120
      j=i-1
      if(zp(nx,j).gt.tour) goto 120
      call toura(nx,i,nx,j,rgrx1,rgry1,rdx,rdy,tour,idmz,nx,ny,zp,ison,
     &           iOutLU)
 120  i=i+1
      if(i.le.ny) goto 116

c  Search top edge.
      i=1
 136  ii=nx-i
      if (zp(ii,ny).le.tour) goto 140
      j=ii+1
      if (zp(j,ny).gt.tour) goto 140
      call toura(ii,ny,j,ny,rgrx1,rgry1,rdx,rdy,tour,idmz,nx,ny,zp,ison,
     &           iOutLU)
 140  i=i+1
      if(i.lt.nx) goto 136

c  Search the rest of the array.
      i=2
 166  j=2
 176  jj=j-1
      if(zp(j,i).le.tour) goto 180
      if(zp(jj,i).gt.tour) goto 180
      ibit=(i-1)*nx+(j-1)
      iword=ibit/31
      ibit=ibit-iword*31
      iword=iword+1
      if(lbtest(ison(iword),ibit)) goto 180
      call toura(j,i,jj,i,rgrx1,rgry1,rdx,rdy,tour,idmz,nx,ny,zp,ison,
     &           iOutLU)
 180  j=j+1
      if(j.lt.nx) goto 176
      i=i+1
      if(i.lt.ny) goto 166
      return
      end


c***********************************************************************
c
      subroutine toura(irpp,jrpp,ispp,jspp,rgrx1,rgry1,rdx,rdy,
     &                 tour,idmz,nx,ny,zp,ison,iOutLU)
      integer*4 irp,jrp,isp,jsp,iOutLU
      real*4    rgrx1,rgry1,rdx,rdy
      integer*4 idmz,nx,ny
      real*4    zp(idmz,ny)
      integer*4 ison(*)

      integer*4 locate,it(9),jt(9)
      logical dtest(9),chkk,lbset
      data it/0,1,1,0,9,0,-1,-1,0/,jt/-1,0,0,-1,9,1,0,0,1/
      data dtest/.false.,.true.,.false.,.true.,.false.,.true.,.false.,
     &           .true.,.false./

      irp=irpp
      jrp=jrpp
      isp=ispp
      jsp=jspp
      ipen=0
      write (iOutLU,999)tour
 999  format("LEVEL ",e14.6)
      goto 1020

 1000 irp=in
      jrp=jn
c  Switch points.
 1020 xr=float(irp-1)*rdx+rgrx1
      yr=float(jrp-1)*rdy+rgry1
      hr= zp(irp,jrp)
      xs=float(isp-1)*rdx+rgrx1
      ys=float(jsp-1)*rdy+rgry1
      hs=zp(isp,jsp)
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)

c  Find the next point to check by looking in table.
 1030 locate=3*(jrp-jsp)+irp-isp+5
      in=isp+it(locate)
      jn=jsp+jt(locate)
c  Test for an edge.
      if (in.gt.nx.or.in.lt.1.or.jn.gt.ny.or.jn.lt.1) return

c  It may be a diagonal.
      if (locate.eq.6) then
        ibit=(jrp-1)*nx+(irp-1)
        iword=ibit/31
        ibit=ibit-iword*31
        iword=iword+1
        chkk=lbset(ison(iword),ibit)
        if (chkk) return
      endif

      if (dtest(locate)) goto 1060
c  Determine that it is a contour or switch points.
      zpp=zp(in,jn)
      if (zpp.gt.tour) goto 1000
      isp=in
      jsp=jn
      goto 1020

c  Diagonals get special treatment.
c  Calculate height and location of the midpoint.
 1060 vx=0.5*rdx*float(irp+in-2)
      vy=0.5*rdy*float(jrp+jn-2)
      locate=3*(jrp-jn)+irp-in+5
      inn=in+it(locate)
      jnn=jn+jt(locate)
      htm=(zp(irp,jrp)+zp(isp,jsp)+zp(in,jn)+zp(inn,jnn))/4.
      if(htm.gt.tour) goto 1120
c  midpoint less than contour
      xr=float(irp-1)*rdx+rgrx1
      yr=float(jrp-1)*rdy+rgry1
      hr=zp(irp,jrp)
      xs=vx+rgrx1
      ys=vy+rgry1
      hs=htm
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)

 1070 if (zp(inn,jnn).gt.tour) goto 1080
c  turn off sharp right
      isp=inn
      jsp=jnn
      goto 1020
 1080 xr=float(inn-1)*rdx+rgrx1
      yr=float(jnn-1)*rdy+rgry1
      hr=zp(inn,jnn)
      xs=vx+rgrx1
      ys=vy+rgry1
      hs=htm
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)

 1090 if (zp(in,jn).gt.tour) goto 1100
c  continue straight through
 1140 irp=inn
      jrp=jnn
      isp=in
      jsp=jn
      goto 1020
c  wide left turn
 1100 xr=float(in-1)*rdx+rgrx1
      yr=float(jn-1)*rdy+rgry1
      hr=zp(in,jn)
      xs=vx+rgrx1
      ys=vy+rgry1
      hs=htm
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)
      goto 1000

c.......................................................................
c  midpoint greater than contour
 1120 xr=vx+rgrx1
      yr=vy+rgry1
      hr=htm
      xs=float(isp-1)*rdx+rgrx1
      ys=float(jsp-1)*rdy+rgry1
      hs=zp(isp,jsp)
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)

c  it may be a sharp left turn
 1125 if (zp(in,jn).gt.tour) goto 1000
      xr=vx+rgrx1
      yr=vy+rgry1
      hr=htm
      xs=float(in-1)*rdx+rgrx1
      ys=float(jn-1)*rdy+rgry1
      hs=zp(in,jn)
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)

c  wide right turn
 1130 if(zp(inn,jnn).gt.tour) goto 1140
      xr=vx+rgrx1
      yr=vy+rgry1
      hr=htm
      xs=float(inn-1)*rdx+rgrx1
      ys=float(jnn-1)*rdy+rgry1
      hs=zp(inn,jnn)
      call tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)
 1135 isp=inn
      jsp=jnn
      goto 1020
      end


c***********************************************************************
c  Plot line to next point.
      subroutine tourbi(xr,yr,hr,xs,ys,hs,tour,ipen,iOutLU)
      integer*4 ipen,ix,iy,iOutLU
      frac=(hr-tour)/(hr-hs)
      xri=xr-frac*(xr-xs)
      yri=yr-frac*(yr-ys)

      write(iOutLU,100) xri,yri
 100  format("      ",2e14.6)  

      ipen=1
      return
      end


c***********************************************************************
      logical function lbtest (iword,ibit)
      data ihigh /1073741824/
      ifac=ihigh
      idum=iword
      do 1000 i=2,32
        j=32-i
        kbyte=idum/ifac
        if(kbyte.eq.1) idum=idum-ifac
        if(j.eq.ibit) goto 2000
        ifac=ifac/2
 1000 continue
 2000 lbtest=kbyte.eq.1
      return
      end


c***********************************************************************
      logical function lbset(iword,ibit)
      logical lbtest
      lbset=lbtest(iword,ibit)
      if(lbset) return
      iword=iword+2**ibit
      return
      end

