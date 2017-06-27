      Program hal3gw

* Two dimensional hexagonal plane!! 2 bands!!
* Compute orbital magnetization IC contrubution as well as Chern number!!

* rmu calcalted as in the middle of the gap !!
* Read Greens function from 116 or 117 !!
* nphi same as in hgw4.f !!


      implicit real*8 (a-h,o-z)
      parameter (maxp=1500,kdim=2,maxk=10000) !maxk is new (Just an exaggeratred size for the Green's function for each phi)

      double precision eA,eB,t1,t2,phi
      double precision rtemp1,rtemp2,rtemp3 


      double precision pi,twopi 
      double precision uxm,uym,uzm,numk,numkbound
      double precision kx,ky,dkx,dky
      double precision ux,uy,uz,dux,duy,duz, kmaxGK,kmaxKKP

      double precision nos(maxp),dos(maxp)
      double precision tnos(maxp),tdos(maxp),x1(2),f1(2)
      double precision emin,emax,ebot,etop,e(2,8),ee(4),bb(4)
      double precision freq(maxp),rdos(maxp),cdos(maxp),
     .                 v2(maxp),v2d(maxp),rgamma(maxp)

      double precision rsp(maxp),csp(maxp),sumr(maxp),sumc(maxp)
      double precision h1(2,2,2), o1(2,2,2), z1(2,2,2), w1(2,11), 
     .                 e1(2)

      double precision e1xp(2)  !kx+d
      double precision e1xm(2)  !kx-d
      double precision e1yp(2)  !ky+d
      double precision e1ym(2)  !ky-d

      double precision Gnnp(2,2,maxk)  ! Gnnp is new 

      double precision z1xp(2,2,2)  !kx+d
      double precision z1xm(2,2,2)  !kx-d
      double precision z1yp(2,2,2)  !ky+d
      double precision z1ym(2,2,2)  !ky-d
      double precision phiC(maxp), tx(2), ty(2), ekn(2) 

      complex*16 H0, Hx, Hy, Hz
      complex*16 expphi, ctemp1, ctemp2


      complex*16 ukn(2,2), xi, csum, csum1, csum2, csum3, csumx, csumy
      complex*16 csum5, csum6, csum7
      complex*16 ukxp(2)
      complex*16 ukxm(2)
      complex*16 ukyp(2)
      complex*16 ukym(2)
      complex*16 derxuk(2)
      complex*16 deryuk(2)

      integer nkabc(3), p(4,6), iw1(2)
      integer itest 



      data p/1,2,4,5,4,5,7,8,2,4,5,7,1,3,4,5,4,5,6,8,3,4,5,6/

      pi    = 4.d0*atan(1.d0)
      twopi = 8.d0*atan(1.d0)
      xi    = dcmplx(0.d0,1.d0)

* Open files.
      open(unit=1,file='haldane3.d',status='OLD')
      open(unit=2,file='Morb3-gw.d',status='UNKNOWN')
      open(unit=3,file='Chern3-gw.d',status='UNKNOWN')
      open(unit=116,file='fort.116',status='UNKNOWN') !116 is new (NOT USED)
      open(unit=117,file='fort.117',status='UNKNOWN') !117 is new

      read(1,*)ipp,itotal
      read(1,*)nkabc
      read(1,*)npts,emin,emax
      read(1,*)t1,t2
      read(1,*)phi
      read(1,*)dx,dy
      read(1,*)rmu
      read(1,*)nphi
      close(1,status='keep')
      if(npts.gt.maxp)then
      print*,' Dimension too small, increase maxp!!'
      stop
      endif

* Control parameters
      write(6,*)'kdim:',kdim
      write(6,*)'nkabc(i)',nkabc(1),nkabc(2),nkabc(3)
      write(6,'(a20,f10.3)')'t1:',t1
      write(6,'(a20,f10.3)')'t2:',t2
      write(6,'(a20,f10.3)')'phi (units of pi):',phi
      write(6,'(a20,f10.3)')'dx:',dx
      write(6,'(a20,f10.3)')'dy:',dy
      write(6,'(a20,f10.3)')'Chemical potential:',rmu
      write(6,'(a20,i3)')'nphi:',nphi

      phi = phi*pi
      write(6,'(a20,f10.3)')'phi (in radians):',phi

* Basis!
      iphase = 0
      if(iphase.eq.1)then
      tx(1) = 0.d0
      tx(2) = 0.d0
      ty(1) = 0.d0
      ty(2) = 1.d0
      write(6,*)'With phasefactors!!'
      else
      tx(1) = 0.d0
      tx(2) = 0.d0
      ty(1) = 0.d0
      ty(2) = 0.d0
      write(6,*)'No phasefactors!!'
      endif



* Calculate M.
* If boundary k-points is used !!
      volwgt = ( (8.d0*pi*pi)/(3.d0*dsqrt(3.d0)) )*
     .         ( 1.d0/((nkabc(1)+1)*(nkabc(2)+1)) )
* If no boundary k-points is used !!
*      volwgt = ( (8.d0*pi*pi)/(3.d0*dsqrt(3.d0)) )*
*     .         ( 1.d0/(nkabc(1)*nkabc(2)) )

      numkbound = (nkabc(1)+1)*(nkabc(2)+1)
      write(6,*)'Number k (with boundaries):',numkbound
      numk = nkabc(1)*nkabc(2)
      write(6,*)'Number k (no boundaries):',numk

      uxm=(2.d0/3.d0)*2.d0*pi
      uym=uxm
      dux=dble(uxm/nkabc(1))
      duy=dble(uym/nkabc(2))

* Choose yy = Delta/t2 !!
      yy = 6.00d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 5.25d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 4.00d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 3.67d0
*      yy = 3.5d0
!      yy = 3.00d0     !Ceresoli PRB74, 024408
*      yy = 1.00d0
*      yy = 0.00d0

      beta = 1.d0/0.05d0  !Inverse temperature; see Ceresoli PRB74, 024408

      eA = -t2*yy
      eB = -eA

      ! new: removed writing if normal or chern
      write(6,'(a20,f10.3)')'eA:',eA
      write(6,'(a20,f10.3)')'eB:',eB
      write(6,'(a20,f10.3)')'Delta/t2:',dabs(eA)/t2
      write(6,'(a20,f10.3)')'beta:',beta

!!!!! Beginning of new part

* Test k-mesh !!
* Without boundary we have (2*nkabc(1)+1) fewer points compared with the boundary case!!
      write(6,*)'k-mesh!'
      ux=-dux

      icount = 0

      
      do ix = 1, nkabc(1) +1  !Also boundary
*      do ix = 1, nkabc(1)      !No boundary 2->3 and 3->4 

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) +1    !Also boundary 
*      do iy = 1, nkabc(2)        !No boundary 2->3 and 3->4

      uy=uy+duy


* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy

      icount = icount + 1
      write(6,'(2i4,2f10.6)')ix,iy,kx/twopi, ky/twopi


* Find index for Dirac point k=K.
      tolk = 0.000001d0
      if(dabs(kx-twopi*(2.d0/3.d0)*(1.d0/dsqrt(3.d0))).lt.tolk)then
      if(dabs(ky).lt.tolk)then
      write(6,*)'Index K:',ix,iy
      endif
      endif

      enddo
      enddo

      write(6,*)'icount:',icount
*      stop 

!!!!! End of new part


      dphi = pi/dble(nphi-1)
      do i = 1,nphi
      phiC(i) = dphi*dble(i-1)
      enddo


* Only one rewind !!
*      rewind(116)
      rewind(117)      ! This is new
      tol = 0.0001d0   ! This is new


      do ip = 1, nphi 
*      do ip = 11, 11 

      phi = phiC(ip) 

*      read (116,'(i4,5f12.6)') ip0, phi0
      read (117,'(i4,5f12.6)') ip0, phi0           ! This is new
      write (115,'(i4,5f12.6)') ip0, phi0          ! This is new
      if(dabs(phi0-phi).ge.tol) stop 'Error phi'   ! This is new


* Find rmu for each phi !!
* Gap is always at K-point!
      kx = twopi*(2.d0/3.d0)*(1.d0/dsqrt(3.d0))
      ky = 0.d0

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =   h1(1,2,1)
      h1(2,1,2)  =  -h1(1,2,2)

* For kx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      rmu = (e1(2) + e1(1))/2.d0
      gap = (e1(2) - e1(1))/2.d0


      write(69,'(5f10.6)')phi/pi, rmu, e1(1),e1(2),gap 


! This is where the interesting code starts

* For Chern number and IC contribution!
      csum1 = dcmplx(0.d0,0.d0)
      csum2 = dcmplx(0.d0,0.d0)
      csum3 = dcmplx(0.d0,0.d0)
      csum7 = dcmplx(0.d0,0.d0)

      icount = 0
      ux=-dux

      do ix = 1, nkabc(1) +1   !With boundary
*      do ix = 1, nkabc(1)       !No boundary

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) +1   !With boundary
*      do iy = 1, nkabc(2)       !No boundary 

      uy=uy+duy

      icount = icount + 1

* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy

!!!!! Beginning of new part

*      read (116,'(i4,6f12.6)') iq, tempkx,tempky,
*     .                         Gnnp(1,1,iq),Gnnp(1,2,iq),
*     .                         Gnnp(2,1,iq),Gnnp(2,2,iq)
      read (117,'(i4,6f12.6)') iq, tempkx,tempky,
     .                         Gnnp(1,1,iq),Gnnp(1,2,iq),
     .                         Gnnp(2,1,iq),Gnnp(2,2,iq)
      write (115,'(i4,6f12.6)') iq, tempkx,tempky,
     .                         Gnnp(1,1,iq),Gnnp(1,2,iq),
     .                         Gnnp(2,1,iq),Gnnp(2,2,iq)

!      Gnnp(1,1,icount) = 1.d0 This should probably be removed (so I
!      comment it)

*      if(icount.ne.iq)   stop
*      if(kx.ne.tempkx)   stop
*      if(ky.ne.tempky)   stop

!!!!! End of new part

* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)

* For kx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

* Fermifactor.
      ff    = 1.d0/(1.d0+dexp( beta*(e1(1)-rmu) ))


      do ib      = 1,kdim  !Band
      ekn(ib) = e1(ib) - rmu
      do i1      = 1,kdim
      ukn(i1,ib) = dcmplx(z1(i1,ib,1),z1(i1,ib,2))
      enddo
      enddo

! new: Comment on the above: For unshifted stuff we already include all bands.
!      We want to do that also for shifted things (at the same time as we use 
!      the proper Green function). I guess we should have two entries on
!      the ukxp etc.

* For kx+dx,ky!!
      kx  =  ux*dsqrt(3.d0)/2.d0 + dx
      ky  = -ux/2.d0 + uy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx+dx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1xp,e1xp)


* Determine dual state! Band 1 only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,1,1),-z1  (i1,1,2))*
     .               dcmplx(z1xp(i1,1,1), z1xp(i1,1,2))*
     .               cdexp(xi*dx*( tx(i1)))

      enddo
      do i1    = 1,kdim
      ukxp(i1) = dcmplx(z1xp(i1,1,1),z1xp(i1,1,2))/csum
      enddo

!!!!! Beginning of new (commented) code
  
* Determine dual state! Band 2 occupied.
* Resta!
*      csum  = dcmplx(0.d0,0.d0)
*      do i1 = 1,kdim
*      csum  = csum + dcmplx(z1  (i1,2,1),-z1  (i1,2,2))*
*     .               dcmplx(z1xp(i1,2,1), z1xp(i1,2,2))*
*     .               cdexp(xi*dx*( tx(i1)))
*
*      enddo
*
*      do i1    = 1,kdim
*      ukxp(i1) = dcmplx(z1xp(i1,2,1),z1xp(i1,2,2))/csum
*      enddo

!!!!! End of new (commented) code

* For kx-dx,ky!!
      kx  =  ux*dsqrt(3.d0)/2.d0 - dx
      ky  = -ux/2.d0 + uy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx-dx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1xm,e1xm)


* Determine dual state! Band 1 only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,1,1),-z1  (i1,1,2))*
     .               dcmplx(z1xm(i1,1,1), z1xm(i1,1,2))*
     .               cdexp(xi*dx*(-tx(i1)))
      enddo


      do i1    = 1,kdim
      ukxm(i1) = dcmplx(z1xm(i1,1,1),z1xm(i1,1,2))/csum
      enddo

!!!!! Beginning of new (commented) code

* Determine dual state! Band 2 occupied.
* Resta!
*      csum  = dcmplx(0.d0,0.d0)
*      do i1 = 1,kdim
*      csum  = csum + dcmplx(z1  (i1,2,1),-z1  (i1,2,2))*
*     .               dcmplx(z1xm(i1,2,1), z1xm(i1,2,2))*
*     .               cdexp(xi*dx*(-tx(i1)))
*      enddo
*
*      do i1    = 1,kdim
*      ukxm(i1) = dcmplx(z1xm(i1,2,1),z1xm(i1,2,2))/csum
*      enddo

!!!!! End of new (commented) code


* For kx,ky+dy!!
      kx  =  ux*dsqrt(3.d0)/2.d0 
      ky  = -ux/2.d0 + uy + dy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx,ky+dy!!
      call diagno(kdim,h1,o1,w1,iw1,z1yp,e1yp)


* Determine dual state! Band 1 only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,1,1),-z1  (i1,1,2))*
     .               dcmplx(z1yp(i1,1,1), z1yp(i1,1,2))*
     .               cdexp(xi*dy*( ty(i1)))
      enddo


      do i1    = 1,kdim
      ukyp(i1) = dcmplx(z1yp(i1,1,1),z1yp(i1,1,2))/csum
      enddo

!!!!! Beginning of new (commented) code

* Determine dual state! Band 2 ounoccupied.
* Resta!
*      csum  = dcmplx(0.d0,0.d0)
*      do i1 = 1,kdim
*      csum  = csum + dcmplx(z1  (i1,2,1),-z1  (i1,2,2))*
*     .               dcmplx(z1yp(i1,2,1), z1yp(i1,2,2))*
*     .               cdexp(xi*dy*( ty(i1)))
*      enddo
*
*      do i1    = 1,kdim
*      ukyp(i1) = dcmplx(z1yp(i1,2,1),z1yp(i1,2,2))/csum
*      enddo

!!!!! End of new (commented) code


* For kx,ky-dy!!
      kx  =  ux*dsqrt(3.d0)/2.d0
      ky  = -ux/2.d0 + uy - dy
* Fix diagonal elements.
      do i = 1,kdim
      do j = 1,kdim
      h1(i,j,1)  = 0.d0
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0
      enddo

* Diagonal elements H11 and H22 first!
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      h1(1,1,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eA

      rtemp1     = dcos(dsqrt(3.d0)*kx + phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 + phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 + phi)

      h1(2,2,1)  =  -2.d0*t2*(rtemp1 + rtemp2 + rtemp3) + eB

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      
      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* Hopping H12.
      rtemp1     = dcos(ky)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dcos( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,1)  =  -t1*(rtemp1 + rtemp2 + rtemp3)

      rtemp1     = dsin(ky)
      rtemp2     = dsin(-dsqrt(3.d0)*kx/2.d0 - ky/2.d0)
      rtemp3     = dsin( dsqrt(3.d0)*kx/2.d0 - ky/2.d0)

      h1(1,2,2)  =   t1*(rtemp1 + rtemp2 + rtemp3)

* H21.
      h1(2,1,1)  =    h1(1,2,1)
      h1(2,1,2)  =  - h1(1,2,2)


* For kx,ky-dy!!
      call diagno(kdim,h1,o1,w1,iw1,z1ym,e1ym)


* Determine dual state! Band 1 only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,1,1),-z1  (i1,1,2))*
     .               dcmplx(z1ym(i1,1,1), z1ym(i1,1,2))*
     .               cdexp(xi*dy*(-ty(i1)))
      enddo


      do i1    = 1,kdim
      ukym(i1) = dcmplx(z1ym(i1,1,1),z1ym(i1,1,2))/csum
      enddo

!!!!! Beginning of new (commented) code

* Determine dual state! Band 2 unoccupied.                
* Resta! 
*      csum  = dcmplx(0.d0,0.d0)
*      do i1 = 1,kdim
*      csum  = csum + dcmplx(z1  (i1,2,1),-z1  (i1,2,2))*     
*     .               dcmplx(z1ym(i1,2,1), z1ym(i1,2,2))*
*     .               cdexp(xi*dy*(-ty(i1)))
*      enddo
*
*      do i1    = 1,kdim
*      ukym(i1) = dcmplx(z1ym(i1,2,1),z1ym(i1,2,2))/csum
*      enddo

!!!!! End of new (commented) code


* Determine the derivatives!
      do i2      = 1,kdim
      derxuk(i2) = ( ukxp(i2) - ukxm(i2) )/(2.d0*dx) 
      deryuk(i2) = ( ukyp(i2) - ukym(i2) )/(2.d0*dy) 
      enddo

!!!!! Beginning of new csum1

      derxEn     = (e1xp(1)-e1xm(1))/(2.d0*dx)
      deryEn     = (e1yp(1)-e1ym(1))/(2.d0*dy)

      ctemp1 =          ( dconjg(derxuk(1))*deryuk(1)    !IC contribution
     .               +    dconjg(derxuk(2))*deryuk(2)    ! 
     .               -    dconjg(deryuk(1))*derxuk(1)    ! old part
     .               -    dconjg(deryuk(2))*derxuk(2) )* !
     .                    (e1(1)-rmu)                    ! 
     .
     .               +  ( dconjg(deryuk(1))*ukn(1,1)     ! 
     .               +    dconjg(deryuk(2))*ukn(2,1) )*  !
     .                    derxEn                         ! new part
     .               -  ( dconjg(derxuk(1))*ukn(1,1)     !
     .               +    dconjg(derxuk(2))*ukn(2,1) )*  !
     .                    deryEn                         !

      csum1 = csum1 + Gnnp(1,1,icount)*ctemp1            ! new (I guess this has to be modified 
                                                         ! since we only use Gnp(1,1,icount) here 

!!!!! End of new csum1 (old code commented below)
!!!!!                  (csum2 same as before)

!      csum1 = csum1  +  ( dconjg(derxuk(1))*deryuk(1)    !IC contribution
!     .               +    dconjg(derxuk(2))*deryuk(2) 
!     .               -    dconjg(deryuk(1))*derxuk(1) 
!     .               -    dconjg(deryuk(2))*derxuk(2) )*
!     .                    (e1(1)-rmu) 



      csum2 = csum2  +  dconjg(derxuk(1))*deryuk(1)   !Chern number
     .               +  dconjg(derxuk(2))*deryuk(2)
     .               -  dconjg(deryuk(1))*derxuk(1)
     .               -  dconjg(deryuk(2))*derxuk(2)


!!!!! Beginning of new code

      rtemp1 = dreal(     dconjg(deryuk(1))*ukn(1,1)
     .               +    dconjg(deryuk(2))*ukn(2,1)  
     .               -    dconjg(derxuk(1))*ukn(1,1)
     .               -    dconjg(derxuk(2))*ukn(2,1) ) 
      rtemp2 = dimag(     dconjg(deryuk(1))*ukn(1,1)
     .               +    dconjg(deryuk(2))*ukn(2,1)  
     .               -    dconjg(derxuk(1))*ukn(1,1)
     .               -    dconjg(derxuk(2))*ukn(2,1) ) 


      tol1 = 0.000000001d0
      if(rtemp1.gt.tol1) stop 'Error !!'
      if(rtemp2.gt.tol1) stop 'Error !!'


*      write(6,'(a30,4f10.6)')'<der u(n=1) | u(n=1) >:',kx,ky,
*     .                        rtemp1, rtemp2

!!!!! End of new code

      enddo
      enddo     !End k-sums


*      write(6,*)'icount and numk:',icount,numk


      fact = 1.d0/(twopi)**2 ! new: Before, we divided by 2


      if(ip.eq.1) write(6,'(a30,3f10.2)')'IC contribution to M:'
      if(ip.eq.1) write(2,'(a30,3f10.2)')'IC contribution to M:'

      write(6,'(a30,3f10.4)')'M:',phi/pi,
     .-dimag(csum1*volwgt*fact) ! new: This change cancels the change in fact

      write(2,'(3f10.4)') phi/pi, -dimag(csum1*volwgt*fact) ! new: This
change cancels the change in fact 

      if(ip.eq.1) write(3,'(a30,3f10.2)')'Chern number:'
      write(3,'(5f10.3)')phi/pi,csum2*volwgt*xi/twopi, 
     . 3.d0*dsqrt(3.d0)*dsin(phi),
     .-3.d0*dsqrt(3.d0)*dsin(phi)


      enddo     !End phi-loop


      write(6,*)'Done!'

      stop
      end

      subroutine diagno(ndim,h,o,wk,iwk,z,eb)

C- Diagonalize secular equation with overlap
C----------------------------------------------------------------------
Ci Inputs
Ci    ndim: dimension of problem
Ci    h,o:  hamiltonian, overlap matrices
Ci    wk, work array length at least 11*ndim; iwk, work array
Co Outputs
Co    z:    eigenvectors; eb, eigenvalues
Cr Remarks
Cr    h,o,z are dimensioned (ndim,ndim,2) (real followed by imag. parts)
Cr    h,o are OVERWRITTEN in this routine
C----------------------------------------------------------------------
C Passed parameters
      integer ndim, iwk(ndim)
      double precision h(ndim,ndim,2),o(ndim,ndim,2),z(ndim,ndim,2),
     .         eb(ndim),wk(ndim,11)

C Local variables:
      integer ierr

C --- Make unit eigenvector matrix ---
      call zinit(z,ndim**2)
      call dvcpy(1.d0,0,z,ndim+1,ndim)

C --- Find the eigenvalues and vectors of O^-1/2  H  O^-1/2 ---
      call bchd(ndim,o,o(1,1,2),wk,ierr)
      if (ierr .ne. 0) stop 'DIAGNO: error in bchd'
      call bred(ndim,h,h(1,1,2),o,o(1,1,2),z,z(1,1,2),wk)
      call btridi(ndim,h,h(1,1,2),wk(1,2),wk(1,3),wk(1,4),wk(1,5))
      call imtqlv(ndim,wk(1,2),wk(1,3),wk(1,4),eb,iwk,ierr,wk(1,7))
      if (ierr .ne. 0) stop 'DIAGNO: error in imtqlv'


      call binvit(ndim,wk(1,2),wk(1,3),wk(1,4),ndim,eb,iwk,
     .           z,ierr,wk(1,7),wk(1,8),wk(1,9),wk(1,10),wk(1,11))
ckk (comment write statement)
      if (ierr .ne. 0) write(6,*)'Not converged eigenvector:',iabs(ierr)
      if (ierr .ne. 0) stop 'DIAGNO: error in binvit'

      call btribk(ndim,h,h(1,1,2),wk(1,5),ndim,z,z(1,1,2))
C --- Get the eigenvectors of H - E O ---
      call bbak(ndim,o,o(1,1,2),wk,ndim,z,z(1,1,2))

      end


      double precision function alagr2 (x,xi,fi)

c 92.03.02
c 92.04.10 from alagr3
c two-point interpolation
c given a function fi at two points xi, the routine interpolates
c the function at x
c f(x) = [ (x-x2)/(x1-x2) ] f1
c      + [ (x-x1)/(x2-x1) ] f2

c x  = the point at which the function is to be interpolated
c xi(2) = points where the function is given
c fi(2) = the function at xi

      implicit real*8 (a-h,o-z)
      dimension xi(2),fi(2)

      xx1        = x-xi(1)
      xx2        = x-xi(2)
      x12        = xi(1)-xi(2)
      alagr2     = (xx2*fi(1) - xx1*fi(2))/x12

      return
      end                       
