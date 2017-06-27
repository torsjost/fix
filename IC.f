      Program haldane3

* Two dimensional hexagonal plane!! 2 bands!!
* Compute orbital magnetization IC contrubution as well as Chern number!!
* See PRL95, 137205 (2005). Compare with Fig. 2.


      implicit real*8 (a-h,o-z)
      parameter (maxp=1500,kdim=2)

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
      open(unit=1,file='input.d',status='OLD')
      open(unit=2,file='Morb3.d',status='UNKNOWN')
      open(unit=3,file='Chern3.d',status='UNKNOWN')

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

      uxm=(2.d0/3.d0)*2.d0*pi ! This is the dimension of the 4 sided version of BZ
      uym=uxm
      dux=dble(uxm/nkabc(1))
      duy=dble(uym/nkabc(2))

* Choose yy = Delta/t2 !!
      yy = 6.00d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 5.25d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 4.00d0    !Thonhauser PRL 2005 E0=2 and t2=1/3
*      yy = 3.67d0    ! <- This is at the transition for phi=pi/4
*      yy = 3.5d0
*      yy = 3.00d0     !Ceresoli PRB74, 024408
*      yy = 1.00d0
*      yy = 0.00d0

      beta = 1.d0/0.05d0  !Inverse temperature; see Ceresoli PRB74, 024408

      eA = -t2*yy
      eB = -eA

      ! Below, 3.67=3*sqrt(3)*sin(phi=pi/4)... Just a choice
      if(yy.gt.3.67d0) write(6,'(a20,f10.3)')'Normal insulator !!'
      if(yy.le.3.67d0) write(6,'(a20,f10.3)')'Chern insulator !!'
      write(6,'(a20,f10.3)')'eA:',eA
      write(6,'(a20,f10.3)')'eB:',eB
      write(6,'(a20,f10.3)')'Delta/t2:',dabs(eA)/t2
      write(6,'(a20,f10.3)')'beta:',beta



      dphi = pi/dble(nphi-1)
      do i = 1,nphi
      phiC(i) = dphi*dble(i-1)
      enddo

      do ip = 1, nphi 

      phi = phiC(ip) 


* Find rmu for each phi !!
* Gap is always at K-point!
      kx = twopi*(2.d0/3.d0)*(1.d0/dsqrt(3.d0))
      ky = 0.d0

* Fix diagonal elements.
      do i = 1,kdim ! kdim=2
      do j = 1,kdim
      ! h1 is the Hamiltonian
      h1(i,j,1)  = 0.d0 ! I assume 1 and 2 here are real and imag
      h1(i,j,2)  = 0.d0
      o1(i,j,1)  = 0.d0
      o1(i,j,2)  = 0.d0
      enddo
      enddo

      do i = 1,kdim
      o1(i,i,1)  = 1.d0 ! Put the diagonal "atom" components of o1 to 1
      enddo

* Diagonal elements H11 and H22 first!
      ! Recall (Kai Sun)
      rtemp1     = dcos(dsqrt(3.d0)*kx - phi)
      rtemp2     = dcos(-dsqrt(3.d0)*kx/2.d0 + 3.d0*ky/2.d0 - phi)
      rtemp3     = dcos(-dsqrt(3.d0)*kx/2.d0 - 3.d0*ky/2.d0 - phi)

      ! Atom 1 = A, Atom 2 = B
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

! The Hamiltonian matrix is now completed

* For kx,ky!!
      call diagno(kdim,h1,o1,w1,iw1,z1,e1)

      ! We put chem. pot. as the average of energies at the K point
      ! (since it gives the right value), and the gap as the difference
      ! in energy at K. 
      rmu = (e1(2) + e1(1))/2.d0
      gap = (e1(2) - e1(1))/2.d0


      write(69,'(5f10.6)')phi/pi, rmu, e1(1),e1(2),gap 



! Note that we are still at a certain ip (phi value)
* For Chern number and IC contribution!
      csum1 = dcmplx(0.d0,0.d0)
      csum2 = dcmplx(0.d0,0.d0)
      csum3 = dcmplx(0.d0,0.d0)
      csum7 = dcmplx(0.d0,0.d0)

      icount = 0
      ux=-dux

      ! With nkabc(1)=25 this is 1->26 (we thus enter neighboring zones
      ! in k space
      do ix = 1, nkabc(1) +1   !With boundary
*      do ix = 1, nkabc(1)       !No boundary

      ux=ux+dux
      uy=-duy

      do iy = 1, nkabc(2) +1   !With boundary
*      do iy = 1, nkabc(2)       !No boundary 

      uy=uy+duy

      ! icount counts the k index
      icount = icount + 1

* For fix (ux,uy) solve for (kx,ky).
      kx  =  ux*dsqrt(3.d0)/2.d0
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

! Same as before but for this specific k (not K, which we used when
! calculating the chemical potential and the gap)

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

! Note that the chemical potential enters ff (we never make use of the
! gap)

* Fermifactor.
      ff    = 1.d0/(1.d0+dexp( beta*(e1(1)-rmu) ))



      do ib      = 1,kdim  !Band
      ekn(ib)    = e1(ib) - rmu ! For this particular k point and band
                                ! ib this is e-mu
      do i1      = 1,kdim
      ukn(i1,ib) = dcmplx(z1(i1,ib,1),z1(i1,ib,2)) ! |u_{kn}> is the
                                                   ! eigenvector and is complex valued
                                                   ! i1 is the atom
                                                   ! index
      enddo
      enddo


! Note, we are still in the do ix, do iy 
* For kx+dx,ky!!
      ! dx is for the shifted mesh needed for derivatives
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

! This is important (Dual states, if we use covariant derivative)
! Later: I will change the method (nice to first do here, and then with GW)

* Determine dual state! Band 1 only occupied.
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim ! sum over atoms but middle index (band) = 1 (occ)
      csum  = csum + dcmplx(z1  (i1,1,1),-z1  (i1,1,2))*
     .               dcmplx(z1xp(i1,1,1), z1xp(i1,1,2))*
     .               cdexp(xi*dx*( tx(i1)))

      enddo

*      csum = csum/cdabs(csum)

      do i1    = 1,kdim
      ukxp(i1) = dcmplx(z1xp(i1,1,1),z1xp(i1,1,2))/csum
      enddo 

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

*      csum = csum/cdabs(csum)

      do i1    = 1,kdim
      ukxm(i1) = dcmplx(z1xm(i1,1,1),z1xm(i1,1,2))/csum
      enddo


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

*      csum = csum/cdabs(csum)

      do i1    = 1,kdim
      ukyp(i1) = dcmplx(z1yp(i1,1,1),z1yp(i1,1,2))/csum
      enddo


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


* Determine dual state! Band 1 only occupied. ! Modify this later on
      csum  = dcmplx(0.d0,0.d0)
      do i1 = 1,kdim
      csum  = csum + dcmplx(z1  (i1,1,1),-z1  (i1,1,2))*
     .               dcmplx(z1ym(i1,1,1), z1ym(i1,1,2))*
     .               cdexp(xi*dy*(-ty(i1)))
      enddo

*      csum = csum/cdabs(csum)

      do i1    = 1,kdim
      ukym(i1) = dcmplx(z1ym(i1,1,1),z1ym(i1,1,2))/csum
      enddo


* Determine the derivatives!
      do i2      = 1,kdim
      derxuk(i2) = ( ukxp(i2) - ukxm(i2) )/(2.d0*dx) 
      deryuk(i2) = ( ukyp(i2) - ukym(i2) )/(2.d0*dy) 
      enddo

! Here we use the chemical potential again (csum1=0 before the first k
! loop)
      csum1 = csum1  +  ( dconjg(derxuk(1))*deryuk(1)    !IC contribution
     .               +    dconjg(derxuk(2))*deryuk(2) 
     .               -    dconjg(deryuk(1))*derxuk(1) 
     .               -    dconjg(deryuk(2))*derxuk(2) )*
     .                    (e1(1)-rmu) 


      csum2 = csum2  +  dconjg(derxuk(1))*deryuk(1)   !Chern number
     .               +  dconjg(derxuk(2))*deryuk(2)
     .               -  dconjg(deryuk(1))*derxuk(1)
     .               -  dconjg(deryuk(2))*derxuk(2)



      enddo
      enddo     !End k-sums





      fact = 1.d0/(twopi)**2
      fact = fact/2.d0


      if(ip.eq.1) write(6,'(a30,3f10.2)')'IC contribution to M:'
      if(ip.eq.1) write(2,'(a30,3f10.2)')'IC contribution to M:'

      write(6,'(a30,3f10.4)')'M:',phi/pi,
     .dimag(csum1*volwgt*fact)*(-2.d0)

*      write(6,'(a30,3f10.4)')'M:',phi/pi,
*     .dimag(csum1*volwgt*fact)*(-2.d0)  
*      write(6,'(a30,3f10.4)')'M:',phi/pi, 
*     .dimag(csum1*volwgt*fact)*(-2.d0)*2.d0  !FA extra factor 2 

*      write(2,'(3f10.4)') phi/pi, dimag(csum7*volwgt*fact)*( 2.d0)*2.d0    !Need extra factor 2 for Resta
      write(2,'(3f10.4)') phi/pi, dimag(csum1*volwgt*fact)*(-2.d0) 
*      write(2,'(3f10.4)') phi/pi, dimag(csum1*volwgt*fact)*(-2.d0)*2.d0  !FA extra factor 2 

      if(ip.eq.1) write(3,'(a30,3f10.2)')'Chern number:'
      write(3,'(5f10.3)')phi/pi,csum2*volwgt*xi/twopi, 
     .3.d0*dsqrt(3.d0)*dsin(phi),
     .-3.d0*dsqrt(3.d0)*dsin(phi)


      enddo     !End phi-loop


      write(6,*)'Done!'

      stop
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

* eispack1.f below.
C#define NEWVERS
C#ifdef NEWVERS
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


      subroutine bchd(n,br,bi,dl,ierr)
C- Cholesky decomposition of hermetian matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   n:  order of the matrix system.  If b has already been Cholesky
Ci       decomposed, n should be prefixed with a minus sign.
Ci   b:  nonorthogonality matrix (only full upper triangle is used)
Ci       If n is negative, both b and dl must be altered
Co Outputs
Co   b:  strict lower triangle of Cholesky factor l of b
Co       The strict upper triangle is the transposed of the lower one
Co   dl: diagonal elements of l.
Co   ierr:nonzero if decomposition finds matrix not pos. def.
Cr Remarks
Cr   Adapted from hchd which is adapted from eispack reduc for
Cr   hermetian matrices. bchd is written for square matrices.
Cr   At the beginning the strict lower triangle is set to the
Cr   transposed of the strict upper one. Then it is possible to
Cr   change the indices in particular for the innermost loop so
Cr   that the first array index equals the loop variable.
C ----------------------------------------------------------------
C Passed parameters
      integer ierr,n
      double precision br(n,n),bi(n,n),dl(n)
C Local parameters
      integer i,im1,j,k
      double precision xr,xi,y

      do 40 j = 1, n
        do 40 i = 1, j
          br(j,i) = br(i,j)
          bi(j,i) = bi(i,j)
   40 continue

c     .......... form l in the arrays b and dl ..........
      do 80 i = 1, n
         im1 = i - 1
         do 80 j = i, n
            xr =  br(j,i)
            xi = -bi(j,i)
            do 20 k = 1, im1
c             x = x - b*(i,k) * b(j,k)
              xr = xr - (br(k,i)*br(k,j) + bi(k,i)*bi(k,j))
              xi = xi - (br(k,i)*bi(k,j) - bi(k,i)*br(k,j))
   20       continue
            if (j .ne. i) go to 60
            ierr = i
            if (xr .le. 0.0d0) return
            y = dsqrt(xr)
            dl(i) = y
            go to 80
   60       br(j,i) = xr / y
            bi(j,i) = xi / y
            br(i,j) = br(j,i)
            bi(i,j) = bi(j,i)
   80 continue
      ierr = 0
      end
      subroutine bred(n,ar,ai,br,bi,atr,ati,dl)
C- Reduction of nonorthogonal hermetian matrix to orthogonal form
C ----------------------------------------------------------------
Ci Inputs
Ci   n:  order of the matrix system.
Ci   a:  hamiltonian matrix(only full upper triangle is used)
Ci   b:  Cholesky-decomposed nonorthogonality matrix
Ci       (strict upper triangle; see hchd)
Ci   at: work array
Ci   dl: diagonal elements of l (see hchd).
Co Outputs
Co   a:  reduced hermetian matrix:  a <-  ldag^-1 a l^-1
Co       The strict upper triangle
Cr Remarks
Cr   Adapted from hred which is adapted from eispack reduc for
Cr   hermetian matrices. bred is written for square matrices.
Cr   At the beginning the strict upper triangle of at is set to the
Cr   transposed of the strict lower one of a. Then it is possible to
Cr   change the indices in particular for the innermost loop so
Cr   that the first array index equals the loop variable.
C ----------------------------------------------------------------
C Passed parameters
      integer n
      double precision ar(n,n),ai(n,n),br(n,n),bi(n,n),dl(n),
     .                 atr(n,n),ati(n,n)
C Local parameters
      integer i,im1,j,j1,k
      double precision xr,xi,y
c
      do 40 j = 1, n
        do 40 i = 1, j
          atr(i,j) = ar(i,j)
          ati(i,j) = -ai(i,j)
   40 continue
c
c     .......... form the transpose of the upper triangle of inv(l)*a
c                in the lower triangle of the array a ..........
      do 200 i = 1, n
         im1 = i - 1
         y = dl(i)
c
C  (l^-1 * a)_ij  =  (a_ij - sum_k<i  l_ik * (l^-1)_kj) / l_ii
         do 200 j = i, n
            xr = ar(i,j)
            xi = ai(i,j)
            do 160 k = 1, im1
              xr = xr - (br(k,i)*atr(k,j) - bi(k,i)*ati(k,j))
              xi = xi - (br(k,i)*ati(k,j) + bi(k,i)*atr(k,j))
  160       continue
            ar(j,i) = xr / y
            ai(j,i) = xi / y
            atr(i,j) = ar(j,i)
            ati(i,j) = ai(j,i)
  200 continue
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 j = 1, n
         j1 = j - 1
c
         do 300 i = j, n
            xr =  ar(i,j)
            xi = -ai(i,j)
            im1 = i - 1
            do 220 k = j, im1
c             x = x - a(k,j) * b(i,k)
              xr = xr - (atr(k,j)*br(k,i) - ati(k,j)*bi(k,i))
              xi = xi - (atr(k,j)*bi(k,i) + ati(k,j)*br(k,i))
  220       continue
            do 260 k = 1, j1
c             x = x - a*(j,k) * b(i,k)
              xr = xr - (ar(k,j)*br(k,i) + ai(k,j)*bi(k,i))
              xi = xi - (ar(k,j)*bi(k,i) - ai(k,j)*br(k,i))
  260       continue
            ar(j,i) = xr / dl(i)
            ai(j,i) = xi / dl(i)
            atr(i,j) = ar(j,i)
            ati(i,j) = ai(j,i)
  300 continue
      end
      subroutine bbak(n,br,bi,dl,m,zr,zi)
C- Back-transforms eigenvectors to nonorthogonal representation
C ----------------------------------------------------------------
Ci Inputs
Ci   n:  order of the matrix system.
Ci   b:  Cholesky-decomposed nonorthogonality matrix, as decomposed
Ci       by reduch (in its strict lower triangle).
Ci   dl: contains further information about the transformation.
Ci   m:  number of eigenvectors to be back transformed.
Ci   z:  eigenvectors to be back transformed (first m columns)
Co Outputs
Co   z transformed eigenvectors
Cr Remarks
Cr   Adapted from rebakh, eispack
C ----------------------------------------------------------------
C Passed parameters
      integer m,n
      double precision br(n,n),bi(n,n),dl(n),zr(n,m),zi(n,m)
C Local parameters
      integer i,j,k
      double precision xr,xi

      do 100 j = 1, m
         do 100 i = n, 1, -1
            xr = zr(i,j)
            xi = zi(i,j)
            do 60 k = i+1, n
              xr = xr - (br(k,i)*zr(k,j) + bi(k,i)*zi(k,j))
              xi = xi - (br(k,i)*zi(k,j) - bi(k,i)*zr(k,j))
   60       continue
            zr(i,j) = xr/dl(i)
            zi(i,j) = xi/dl(i)
  100 continue
      end
C#elseC
C      subroutine diagno(ndim,h,o,wk,iwk,z,eb)
C
CC- Diagonalize secular equation with overlap
CC----------------------------------------------------------------------
CCi Inputs
CCi    ndim: dimension of problem
CCi    h,o:  hamiltonian, overlap matrices
CCi    wk, work array length at least 5*ndim
CCo Outputs
CCo    z:    eigenvectors; eb, eigenvalues
CCr Remarks
CCr    h,o,z are dimensioned (ndim,ndim,2) (real followed by imag. parts
CCr    h,o are OVERWRITTEN in this routine
CC----------------------------------------------------------------------
CC Passed parameters
C      integer ndim, iwk(ndim)
C      double precision h(ndim,ndim,2),o(ndim,ndim,2),z(ndim,ndim,2),
C     .         eb(ndim),wk(ndim,5)
C
CC Local variables:
C      integer ierr
C
CC --- Make unit eigenvector matrix ---
C      call zinit(z,ndim**2)
C      call dvcpy(1.d0,0,z,ndim+1,ndim)
C
CC --- Find the eigenvalues and vectors of O^-1/2  H  O^-1/2 ---
C      call hchd(ndim,ndim,o,o(1,1,2),wk,ierr)
C      if (ierr .ne. 0) stop 'DIAGNO: error in hchd'
C      call hred(ndim,ndim,h,h(1,1,2),o,o(1,1,2),wk)
C      call htridi(ndim,ndim,h,h(1,1,2),eb,wk(1,2),wk(1,3),wk(1,4))
C      call imtql2(ndim,ndim,eb,wk(1,2),z,ierr)
C      if (ierr .ne. 0) stop 'DIAGNO: error in imtql2'
C      call htribk(ndim,ndim,h,h(1,1,2),wk(1,4),ndim,z,z(1,1,2))
C
CC --- Get the eigenvectors of H - E O ---
C      call hbak(ndim,ndim,o,o(1,1,2),wk,ndim,z,z(1,1,2))
C      end
C      subroutine hchd(nm,n,br,bi,dl,ierr)
CC- Cholesky decomposition of hermetian matrix
CC ----------------------------------------------------------------
CCi Inputs
CCi   nm: true row dimension of arrays b and z, as declared by caller
CCi   n:  order of the matrix system.  If b has already been Cholesky
CCi       decomposed, n should be prefixed with a minus sign.
CCi   b:  nonorthogonality matrix (only full upper triangle is used)
CCi       If n is negative, both b and dl must be altered
CCo Outputs
CCo   b:  strict lower triangle of Cholesky factor l of b
CCo       The strict upper triangle unaltered.
CCo   dl: diagonal elements of l.
CCo   ierr:nonzero if decomposition finds matrix not pos. def.
CCr Remarks
CCr   Adapted from eispack reduc for hermetian matrices
CC ----------------------------------------------------------------
CC Passed parameters
C      integer ierr,n,nm
C      double precision br(nm,n),bi(nm,n),dl(n)
CC Local parameters
C      integer i,im1,j,k
C      double precision xr,xi,y
C
Cc     .......... form l in the arrays b and dl ..........
C      do 80 i = 1, n
C         im1 = i - 1
C         do 80 j = i, n
C            xr =  br(i,j)
C            xi = -bi(i,j)
C            do 20 k = 1, im1
Cc             x = x - b*(i,k) * b(j,k)
C              xr = xr - (br(i,k)*br(j,k) + bi(i,k)*bi(j,k))
C              xi = xi - (br(i,k)*bi(j,k) - bi(i,k)*br(j,k))
C   20       continue
C            if (j .ne. i) go to 60
C            ierr = i
C            if (xr .le. 0.0d0) return
C            y = dsqrt(xr)
C            dl(i) = y
C            go to 80
C   60       br(j,i) = xr / y
C            bi(j,i) = xi / y
C   80 continue
C      ierr = 0
C      end
C      subroutine hred(nm,n,ar,ai,br,bi,dl)
CC- Reduction of nonorthogonal hermetian matrix to orthogonal form
CC ----------------------------------------------------------------
CCi Inputs
CCi   nm: true row dimension of arrays b and z, as declared by caller
CCi   n:  order of the matrix system.
CCi   a:  hamiltonian matrix (only full upper triangle is used)
CCi   b:  Cholesky-decomposed nonorthogonality matrix
CCi       (strict lower triangle; see hchd)
CCi   dl: diagonal elements of l (see hchd).
CCo Outputs
CCo   a:  reduced hermetian matrix:  a <-  ldag^-1 a l^-1
CCo       The strict upper triangle unaltered.
CCr Remarks
CCr   Adapted from eispack reduc for hermetian matrices
CC ----------------------------------------------------------------
CC Passed parameters
C      integer n,nm
C      double precision ar(nm,n),ai(nm,n),br(nm,n),bi(nm,n),dl(n)
CC Local parameters
C      integer i,im1,j,j1,k
C      double precision xr,xi,y
C
Cc     .......... form the transpose of the upper triangle of inv(l)*a
Cc                in the lower triangle of the array a ..........
C      do 200 i = 1, n
C         im1 = i - 1
C         y = dl(i)
Cc
CC  (l^-1 * a)_ij  =  (a_ij - sum_k<i  l_ik * (l^-1)_kj) / l_ii
C         do 200 j = i, n
C            xr = ar(i,j)
C            xi = ai(i,j)
C            do 160 k = 1, im1
C              xr = xr - (br(i,k)*ar(j,k) - bi(i,k)*ai(j,k))
C              xi = xi - (br(i,k)*ai(j,k) + bi(i,k)*ar(j,k))
C  160       continue
C            ar(j,i) = xr / y
C            ai(j,i) = xi / y
C  200 continue
Cc     .......... pre-multiply by inv(l) and overwrite ..........
C      do 300 j = 1, n
C         j1 = j - 1
Cc
C         do 300 i = j, n
C            xr =  ar(i,j)
C            xi = -ai(i,j)
C            im1 = i - 1
C            do 220 k = j, im1
Cc             x = x - a(k,j) * b(i,k)
C              xr = xr - (ar(k,j)*br(i,k) - ai(k,j)*bi(i,k))
C              xi = xi - (ar(k,j)*bi(i,k) + ai(k,j)*br(i,k))
C  220       continue
C            do 260 k = 1, j1
Cc             x = x - a*(j,k) * b(i,k)
C              xr = xr - (ar(j,k)*br(i,k) + ai(j,k)*bi(i,k))
C              xi = xi - (ar(j,k)*bi(i,k) - ai(j,k)*br(i,k))
C  260       continue
C            ar(i,j) = xr / dl(i)
C            ai(i,j) = xi / dl(i)
C  300 continue
C      end
C      subroutine hbak(nm,n,br,bi,dl,m,zr,zi)
CC- Back-transforms eigenvectors to nonorthogonal representation
CC ----------------------------------------------------------------
CCi Inputs
CCi   nm: true row dimension of arrays b and z
CCi   n:  order of the matrix system.
CCi   b:  Cholesky-decomposed nonorthogonality matrix, as decomposed
CCi       by reduch (in its strict lower triangle).
CCi   dl: contains further information about the transformation.
CCi   m:  number of eigenvectors to be back transformed.
CCi   z:  eigenvectors to be back transformed (first m columns)
CCo Outputs
CCo   z transformed eigenvectors
CCr Remarks
CCr   Adapted from rebakh, eispack
CC ----------------------------------------------------------------
CC Passed parameters
C      integer m,n,nm
C      double precision br(nm,n),bi(nm,n),dl(n),zr(nm,m),zi(nm,m)
CC Local parameters
C      integer i,j,k
C      double precision xr,xi
C
C      do 100 j = 1, m
C         do 100 i = n, 1, -1
C            xr = zr(i,j)
C            xi = zi(i,j)
C            do 60 k = i+1, n
C              xr = xr - (br(k,i)*zr(k,j) + bi(k,i)*zi(k,j))
C              xi = xi - (br(k,i)*zi(k,j) - bi(k,i)*zr(k,j))
C   60       continue
C            zr(i,j) = xr/dl(i)
C            zi(i,j) = xi/dl(i)
C  100 continue
C      end
C#endif
      subroutine hcinv(nm,n,br,bi,dl,lopt)
C- Obtain inverse from Cholesky decomposed matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   nm: true row dimension of arrays b and z, as declared by caller
Ci   n:  order of the matrix system.
Ci   b:  Cholesky-decomposed nonorthogonality matrix
Ci       (strict lower triangle; see hchd)
Ci   dl: diagonal elements of l (see hchd).
Ci   lopt: false, only upper half of b\-1 is generated, and lower
Ci         half continues to hold the c.d. of b
Ci         true,  lower half of b also filled in
Co Outputs
Co   b:  Inverse of Cholesky matrix
Cr Remarks
Cr   Adapted from eispack reduc for hermetian matrices
C ----------------------------------------------------------------
C Passed parameters
      integer n,nm
      logical lopt
      double precision br(nm,n),bi(nm,n),dl(n)
C Local parameters
      integer i,im1,j,jm1,k
      double precision xr,xi,y
C#ifdefC BLAS
C      double precision ddot
C#endif

C --- form inv(l) and store transpose in full upper triangle of b ---
C  (l\-1)_ij  =  (del_ij - sum_k<i  l_ik * (l\-1)_kj) / l_ii
      do 200 j = 1, n
        do 200 i = j, n
          y = dl(i)
C#ifdefC BLAS
C          if (j .eq. i) then
C            xr = 1d0
C            xi = 0
C          else
C            xr = - ddot(i-j,br(i,j),nm,br(j,j),nm)
C     .           + ddot(i-j,bi(i,j),nm,bi(j,j),nm)
C            xi = - ddot(i-j,br(i,j),nm,bi(j,j),nm)
C     .           - ddot(i-j,bi(i,j),nm,br(j,j),nm)
C          endif
C#else
          xr = 0
          xi = 0
          if (j .eq. i) xr = 1d0
          im1 = i-1
          do 160 k = j, im1
            xr = xr - (br(i,k)*br(j,k) - bi(i,k)*bi(j,k))
            xi = xi - (br(i,k)*bi(j,k) + bi(i,k)*br(j,k))
  160     continue
C#endif
          br(j,i) = xr / y
          bi(j,i) = xi / y
  200 continue
c
c     .......... pre-multiply by inv(l) and overwrite ..........
      do 300 i = 1, n
         do 300 j = i, n
C           x = x - l\dag_ik * l_kj; here i <= j <= k
C#ifdefC BLAS
C           br(i,j) = ddot(n+1-j,br(i,j),nm,br(j,j),nm) +
C     .               ddot(n+1-j,bi(i,j),nm,bi(j,j),nm)
C           bi(i,j) = ddot(n+1-j,br(i,j),nm,bi(j,j),nm) -
C     .               ddot(n+1-j,bi(i,j),nm,br(j,j),nm)
C#else
            xr = 0
            xi = 0
            do 260 k = j, n
              xr = xr + (br(i,k)*br(j,k) + bi(i,k)*bi(j,k))
              xi = xi + (br(i,k)*bi(j,k) - bi(i,k)*br(j,k))
  260       continue
            br(i,j) = xr
            bi(i,j) = xi
C#endif
  300 continue

      if (.not. lopt) return
      do  400  i = 1, n
       do  400  j = i+1, n
         br(j,i) =  br(i,j)
         bi(j,i) = -bi(i,j)
  400 continue
      end

      subroutine yhifa(ar,ai,lda,n,kpvt,info)
C- Factor an hermetian matrix with symmetric pivoting
      integer lda,n,kpvt(1),info
      double precision ar(lda,1),ai(lda,1)
c
c     to solve  a*x = b , follow yhifa by yhisl.
c     to compute  inverse(a)*c , follow yhifa by yhisl.
c     to compute  determinant(a) , follow yhifa by yhidi.
c     to compute  inertia(a) , follow yhifa by yhidi.
c     to compute  inverse(a) , follow yhifa by yhidi.
c
c     on entry
c
c        ar,ai   the hermitian matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*ctrans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , ctrans(u) is the
c                conjugate transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that yhisl or yhidi may
c                     divide by zero if called.
c
c     Adapted from linpack's version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     internal variables
c
      double precision mulk(2),mulkm1(2),denom(2),ak(2),akm1(2),
     .  bk(2),bkm1(2),t(2)

      double precision absakk,alpha,colmax,rowmax,dcabs1
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,iyamax
      logical swap

c
C --- Initialize ---
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
      info = 0

C --- Main loop on k, which goes from n to 1. ---
      k = n
   10 continue

C ---    Leave the loop if k=0 or k=1. ---
        if (k .eq. 0) go to 200
        if (k .gt. 1) go to 20
          kpvt(1) = 1
          if (dcabs1(ar,ai) .eq. 0.0d0) info = 1
          go to 200
   20   continue
c
c     This section of code determines the kind of
c     elimination to be performed.  when it is completed,
c     kstep will be set to the size of the pivot block, and
c     swap will be set to .true. if an interchange is
c     required.

        km1 = k - 1
        absakk = dcabs1(ar(k,k),ai(k,k))

C ---  Determine the largest off-diagonal element in column k ---
        imax = iyamax(k-1,ar(1,k),ai(1,k),1)
        colmax = dcabs1(ar(imax,k),ai(imax,k))
        if (absakk .lt. alpha*colmax) go to 30
          kstep = 1
          swap = .false.
          go to 90
   30   continue

C ---  Determine the largest off-diagonal element in row imax. ---
        rowmax = 0.0d0
        imaxp1 = imax + 1
        do 40 j = imaxp1, k
          rowmax = dmax1(rowmax,dcabs1(ar(imax,j),ai(imax,j)))
   40   continue
        if (imax .eq. 1) go to 50
          jmax = iyamax(imax-1,ar(1,imax),ai(1,imax),1)
          rowmax=dmax1(rowmax,dcabs1(ar(jmax,imax),ai(jmax,imax)))
   50   continue
        if (dcabs1(ar(imax,imax),ai(imax,imax)) .lt. alpha*rowmax)
     .    go to 60
          kstep = 1
          swap = .true.
          go to 80
   60   continue
        if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
          kstep = 1
          swap = .false.
          go to 80
   70   continue
        kstep = 2
        swap = imax .ne. km1
   80   continue
   90   continue
        if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100

C ---   Column k is zero.  set info and iterate the loop. ---
        kpvt(k) = k
        info = k
        go to 190
  100   continue
        if (kstep .eq. 2) go to 140

C ---      1 x 1 pivot block. ---
          if (.not.swap) go to 120

C ---       Perform an interchange. ---
            call dswap(imax,ar(1,imax),1,ar(1,k),1)
            call dswap(imax,ai(1,imax),1,ai(1,k),1)
            do 110 jj = imax, k
              j = k + imax - jj
              t(1) = ar(j,k)
              t(2) = -ai(j,k)
              ar(j,k) =  ar(imax,j)
              ai(j,k) = -ai(imax,j)
              ar(imax,j) = t(1)
              ai(imax,j) = t(2)
  110       continue
  120       continue

C ---       Perform the elimination. ---
            do 130 jj = 1, km1
              j = k - jj
C     complex divide ... mulk = -a(j,k)/a(k,k) ...
              call cdiv(-ar(j,k),-ai(j,k),
     .          ar(k,k),ai(k,k),mulk(1),mulk(2))
c     call zaxpy(j,mulk*,ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(2),ai(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ai(1,k),1,ai(1,j),1)
              call daxpy(j,-mulk(2),ar(1,k),1,ai(1,j),1)
              ai(j,j) = 0
              ar(j,k) = mulk(1)
              ai(j,k) = mulk(2)
  130       continue

C --- Set the pivot array. ---
            kpvt(k) = k
            if (swap) kpvt(k) = imax
            go to 190
  140     continue

C --- 2 x 2 pivot block. ---
          if (.not.swap) go to 160
c
C ---      Perform an interchange. ---
            call dswap(imax,ar(1,imax),1,ar(1,k-1),1)
            call dswap(imax,ai(1,imax),1,ai(1,k-1),1)
            do 150 jj = imax, km1
              j = km1 + imax - jj
              t(1) = ar(j,k-1)
              t(2) = -ai(j,k-1)
              ar(j,k-1) =  ar(imax,j)
              ai(j,k-1) = -ai(imax,j)
              ar(imax,j) = t(1)
              ai(imax,j) = t(2)
  150       continue
            t(1) = ar(k-1,k)
            t(2) = ai(k-1,k)
            ar(k-1,k) = ar(imax,k)
            ai(k-1,k) = ai(imax,k)
            ar(imax,k) = t(1)
            ai(imax,k) = t(2)
  160     continue

C ---       Perform the elimination. ---
          km2 = k - 2
          if (km2 .eq. 0) go to 180
c ... ak = a(k,k)/a(k-1,k),  akm1 = a(k-1,k-1)/dconjg(a(k-1,k)) ...
            call cdiv(ar(k,k),ai(k,k),ar(k-1,k),ai(k-1,k),ak(1),ak(2))
            call cdiv(ar(k-1,k-1),ai(k-1,k-1),ar(k-1,k),-ai(k-1,k),
     .        akm1(1),akm1(2))
c ... denom = 1.0d0 - ak*akm1 ...
            denom(1) = 1 - ak(1)*akm1(1) + ak(2)*akm1(2)
            denom(2) =   - ak(1)*akm1(2) - ak(2)*akm1(1)
            do 170 jj = 1, km2
              j = km1 - jj
c ... bk = a(j,k)/a(k-1,k),  bkm1 = a(j,k-1)/dconjg(a(k-1,k))
              call cdiv(ar(j,k),ai(j,k),
     .          ar(k-1,k),ai(k-1,k),bk(1),bk(2))
              call cdiv(ar(j,k-1),ai(j,k-1),
     .          ar(k-1,k),-ai(k-1,k),bkm1(1),bkm1(2))
c ... mulk = (akm1*bk - bkm1)/denom
              call cpy(akm1(1),akm1(2),bk(1),bk(2),mulk(1),mulk(2))
              mulk(1) = mulk(1) - bkm1(1)
              mulk(2) = mulk(2) - bkm1(2)
              call cdiv(mulk(1),mulk(2),denom(1),denom(2),
     .          mulk(1),mulk(2))
c ... mulkm1 = (ak*bkm1 - bk)/denom
              call cpy(ak(1),ak(2),bkm1(1),bkm1(2),mulkm1(1),mulkm1(2))
              mulkm1(1) = mulkm1(1) - bk(1)
              mulkm1(2) = mulkm1(2) - bk(2)
              call cdiv(mulkm1(1),mulkm1(2),denom(1),denom(2),
     .          mulkm1(1),mulkm1(2))
c     call zaxpy(j,mulk*,ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ar(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(2),ai(1,k),1,ar(1,j),1)
              call daxpy(j,mulk(1),ai(1,k),1,ai(1,j),1)
              call daxpy(j,-mulk(2),ar(1,k),1,ai(1,j),1)
c     call zaxpy(j,mulkm1*,ar(1,k-1),1,ar(1,j),1)
              call daxpy(j,mulkm1(1),ar(1,k-1),1,ar(1,j),1)
              call daxpy(j,mulkm1(2),ai(1,k-1),1,ar(1,j),1)
              call daxpy(j,mulkm1(1),ai(1,k-1),1,ai(1,j),1)
              call daxpy(j,-mulkm1(2),ar(1,k-1),1,ai(1,j),1)
              ar(j,k) = mulk(1)
              ai(j,k) = mulk(2)
              ar(j,k-1) = mulkm1(1)
              ai(j,k-1) = mulkm1(2)
              ai(j,j) = 0
  170       continue
  180     continue

C ---    Set the pivot array. ---
          kpvt(k) = 1 - k
          if (swap) kpvt(k) = -imax
          kpvt(k-1) = kpvt(k)
  190   continue
        k = k - kstep
        go to 10
  200 continue
      end
      subroutine yhidi(ar,ai,lda,n,kpvt,det,inert,work,job)
C- Compute determinant, inertia and inverse using factors from yhifa
      integer lda,n,job
      double precision ar(lda,lda),ai(lda,lda),work(2,1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     on entry
c
c        ar,ai   output from yhifa
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from yhifa.
c
c        work    complex*16(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    double precision(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. dabs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  yhico  has set rcond .eq. 0.0
c        or  yhifa  has set  info .ne. 0 .
c
c     Adapted from linpack's version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     internal variables.
c
      double precision ten,d,t,ak,akp1,akkp1(2),temp,cdabs2,ddot
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
        if (noert) go to 10
          inert(1) = 0
          inert(2) = 0
          inert(3) = 0
   10   continue
        if (nodet) go to 20
          det(1) = 1.0d0
          det(2) = 0.0d0
          ten = 10.0d0
   20   continue
        t = 0.0d0
        do 130 k = 1, n
          d = ar(k,k)

C ---       Check if 1 by 1 ---
          if (kpvt(k) .gt. 0) go to 50

c     2 by 2 block
c     use det (d  s)  =  (d/t * c - t) * t  ,  t = cdabs(s)
c                        (s  c)
c     to avoid underflow/overflow troubles.
c     take two passes through scaling.  use  t  for flag.
c
            if (t .ne. 0.0d0) go to 30
              t = cdabs2(ar(k,k+1),ai(k,k+1))
              d = (d/t)*ar(k+1,k+1) - t
              go to 40
   30       continue
            d = t
            t = 0.0d0
   40       continue
   50     continue
c
          if (noert) go to 60
            if (d .gt. 0.0d0) inert(1) = inert(1) + 1
            if (d .lt. 0.0d0) inert(2) = inert(2) + 1
            if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60     continue
c
          if (nodet) go to 120
            det(1) = d*det(1)
            if (det(1) .eq. 0.0d0) go to 110
   70       if (dabs(det(1)) .ge. 1.0d0) go to 80
              det(1) = ten*det(1)
              det(2) = det(2) - 1.0d0
              go to 70
   80       continue
   90       if (dabs(det(1)) .lt. ten) go to 100
              det(1) = det(1)/ten
              det(2) = det(2) + 1.0d0
              go to 90
  100       continue
  110       continue
  120       continue
  130    continue
  140 continue
c
C --- Compute inverse(a) ---
      if (noinv) go to 270
      k = 1
  150 if (k .gt. n) go to 260
        km1 = k - 1
        if (kpvt(k) .lt. 0) go to 180
c
C --- 1 by 1 ---
        ar(k,k) = 1/ar(k,k)
        ai(k,k) = 0
        if (km1 .lt. 1) go to 170
          call dcopy(km1,ar(1,k),1,work,2)
          call dcopy(km1,ai(1,k),1,work(2,1),2)
          do 160 j = 1, km1
            ar(j,k) = ddot(j,ar(1,j),1,work,2)
     .        +ddot(j,ai(1,j),1,work(2,1),2)
            ai(j,k) = ddot(j,ar(1,j),1,work(2,1),2)
     .        -ddot(j,ai(1,j),1,work,2)
c     call zaxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,-work(2,j),ai(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ai(1,j),1,ai(1,k),1)
            call daxpy(j-1,work(2,j),ar(1,j),1,ai(1,k),1)
  160     continue
          ar(k,k) = ar(k,k)
     .      +ddot(km1,work,2,ar(1,k),1)
     .      +ddot(km1,work(2,1),2,ai(1,k),1)
          ai(k,k) = 0
  170   continue
        kstep = 1
        go to 220
  180   continue

C  ---2 by 2 ---
        t = cdabs2(ar(k,k+1),ai(k,k+1))
        ak = ar(k,k)/t
        akp1 = ar(k+1,k+1)/t
c ... akkp1 = a(k,k+1)/t
        akkp1(1) = ar(k,k+1)/t
        akkp1(2) = ai(k,k+1)/t
        d = t*(ak*akp1 - 1.0d0)
        ar(k,k) = akp1/d
        ai(k,k) = 0
        ar(k+1,k+1) = ak/d
        ai(k+1,k+1) = 0
        ar(k,k+1) = -akkp1(1)/d
        ai(k,k+1) = -akkp1(2)/d
        if (km1 .lt. 1) go to 210
          call dcopy(km1,ar(1,k+1),1,work,2)
          call dcopy(km1,ai(1,k+1),1,work(2,1),2)
          do 190 j = 1, km1
            ar(j,k+1) =
     .        ddot(j,ar(1,j),1,work,2) +
     .        ddot(j,ai(1,j),1,work(2,1),2)
            ai(j,k+1) =
     .        ddot(j,ar(1,j),1,work(2,1),2) -
     .        ddot(j,ai(1,j),1,work,2)
c     call zaxpy(j-1,work(1,j),ar(1,j),1,ar(1,k+1),1)
            call daxpy(j-1,work(1,j),ar(1,j),1,ar(1,k+1),1)
            call daxpy(j-1,-work(2,j),ai(1,j),1,ar(1,k+1),1)
            call daxpy(j-1,work(1,j),ai(1,j),1,ai(1,k+1),1)
            call daxpy(j-1,work(2,j),ar(1,j),1,ai(1,k+1),1)
  190     continue
          ar(k+1,k+1) = ar(k+1,k+1)
     .      +ddot(km1,work,2,ar(1,k+1),1)
     .      +ddot(km1,work(2,1),2,ai(1,k+1),1)
          ai(k+1,k+1) = 0
          ar(k,k+1) = ar(k,k+1)
     .      +ddot(km1,ar(1,k),1,ar(1,k+1),1)
     .      +ddot(km1,ai(1,k),1,ai(1,k+1),1)
          ai(k,k+1) = ai(k,k+1)
     .      +ddot(km1,ar(1,k),1,ai(1,k+1),1)
     .      -ddot(km1,ai(1,k),1,ar(1,k+1),1)
          call dcopy(km1,ar(1,k),1,work,2)
          call dcopy(km1,ai(1,k),1,work(2,1),2)
          do 200 j = 1, km1
            ar(j,k) =
     .        ddot(j,ar(1,j),1,work,2)
     .        +ddot(j,ai(1,j),1,work(2,1),2)
            ai(j,k) =
     .        ddot(j,ar(1,j),1,work(2,1),2)
     .        -ddot(j,ai(1,j),1,work,2)
c     call zaxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ar(1,j),1,ar(1,k),1)
            call daxpy(j-1,-work(2,j),ai(1,j),1,ar(1,k),1)
            call daxpy(j-1,work(1,j),ai(1,j),1,ai(1,k),1)
            call daxpy(j-1,work(2,j),ar(1,j),1,ai(1,k),1)
  200     continue
          ar(k,k) = ar(k,k)
     .      +ddot(km1,work,2,ar(1,k),1)
     .      +ddot(km1,work(2,1),2,ai(1,k),1)
          ai(k,k) = 0
  210   continue
        kstep = 2
  220   continue

C --- Swap ---
        ks = iabs(kpvt(k))
        if (ks .eq. k) go to 250
          call dswap(ks,ar(1,ks),1,ar(1,k),1)
          call dswap(ks,ai(1,ks),1,ai(1,k),1)
          do 230 jb = ks, k
            j = k + ks - jb
            temp = ar(j,k)
            ar(j,k) = ar(ks,j)
            ar(ks,j) = temp
            temp = -ai(j,k)
            ai(j,k) = -ai(ks,j)
            ai(ks,j) = temp
  230     continue
          if (kstep .eq. 1) go to 240
            temp = ar(ks,k+1)
            ar(ks,k+1) = ar(k,k+1)
            ar(k,k+1) = temp
            temp = ai(ks,k+1)
            ai(ks,k+1) = ai(k,k+1)
            ai(k,k+1) = temp
  240     continue
  250     continue
          k = k + kstep
          go to 150
  260   continue
  270 continue
      end
      integer function iyamax(n,zxr,zxi,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c

      double precision zxr(1),zxi(1)
      double precision smax
      double precision cdabs2
      double precision dcabs1 
*      real*8 zxr(1),zxi(1)
*      real*8 smax
*      real*8 cdabs2
c
      iyamax = 1
      if(n.le.1)return
      if(incx.eq.1)go to 20
c
c        code for increments not equal to 1
c
      ix = 1
      smax = cdabs2(zxr,zxi)
      ix = ix + incx
      do 10 i = 2,n
         if(cdabs2(zxr(ix),zxi(ix)).le.smax) go to 5
         iyamax = i
         smax = cdabs2(zxr(ix),zxi(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increments equal to 1
c
   20 smax = dcabs1(zxr,zxi)

      do 30 i = 2,n
         if(cdabs2(zxr(i),zxi(i)).le.smax) go to 30
         iyamax = i
         smax = cdabs2(zxr(i),zxi(i))
   30 continue
      return
      end

      double precision function dcabs1(zr,zi)
      double precision zr,zi
      dcabs1 = dabs(zr) + dabs(zi)
      end
      double precision function cdabs2(x1,x2)
      double precision x1,x2
      cdabs2 = dsqrt(x1**2 + x2**2)
      end

      subroutine cpy(tr,ti,dr,di,t1,t2)
C- complex multiply (t1,t2) = (tr,ti) * (dr,di)
      double precision tr,ti,dr,di,t1,t2
      double precision tmp
      tmp = tr*dr - ti*di
      t2  = tr*di + ti*dr
      t1 = tmp
      end
      subroutine cdiv(tr,ti,dr,di,t1,t2)
C- complex divide (t1,t2) = (tr,ti) / (dr,di)
Cr Remarks
Cr   Adapted from eispack.
Cr   It is permissible for (t1,t2) to occupy the same address space
Cr   as either (tr,ti) or (dr,di)
      double precision tr,ti,dr,di,t1,t2
      double precision rr,d,tmp

      if (dabs(di) .gt. dabs(dr)) then
        rr = dr / di
        d = dr * rr + di
        tmp = (tr * rr + ti) / d
        t2 = (ti * rr - tr) / d
        t1 = tmp
      else
        rr = di / dr
        d = dr + di * rr
        tmp = (tr + ti * rr) / d
        t2 = (ti - tr * rr) / d
        t1 = tmp
      endif
      end

C#ifndef SLATEC
C#ifdef NEWVERS
      subroutine btridi(n,ar,ai,d,e,e2,tau)
c
      integer i,j,k,l,n,ii,jp1
      double precision ar(n,n),ai(n,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the upper triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full upper
c          triangles.  their strict lower triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version is adapted from the eispack htridi. In contrast
c     only the upper triangle of a is used for input and output
c
c     -----------------------------------------------------------------
c
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c
      do 100 i = 1, n
  100 d(i) = ar(i,i)
c     .......... for i=n step -1 until 1 do -- ..........
      do 300 ii = 1, n
         i = n + 1 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(ar(k,i)) + dabs(ai(k,i))
c
         if (scale .ne. 0.0d0) go to 140
         tau(1,l) = 1.0d0
         tau(2,l) = 0.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do 150 k = 1, l
            ar(k,i) = ar(k,i) / scale
            ai(k,i) = ai(k,i) / scale
            h = h + ar(k,i) * ar(k,i) + ai(k,i) * ai(k,i)
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(l,i),ai(l,i))
c     .......... form next diagonal element of matrix t ......... .
         if (f .eq. 0.0d0) go to 160
         tau(1,l) = (ai(l,i) * tau(2,i) - ar(l,i) * tau(1,i)) / f
         si = (ar(l,i) * tau(2,i) + ai(l,i) * tau(1,i)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(l,i) = g * ar(l,i)
         ai(l,i) = g * ai(l,i)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(1,l) = -tau(1,i)
         si = tau(2,i)
         ar(l,i) = g
  170    f = 0.0d0
c
         do 240 j = 1, l
            g = 0.0d0
            gi = 0.0d0
c     .......... form element of a*u ..........
            do 180 k = 1, j
               g = g + ar(k,j) * ar(k,i) + ai(k,j) * ai(k,i)
               gi = gi - ar(k,j) * ai(k,i) + ai(k,j) * ar(k,i)
  180       continue
c
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + ar(j,k) * ar(k,i) - ai(j,k) * ai(k,i)
               gi = gi - ar(j,k) * ai(k,i) - ai(j,k) * ar(k,i)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(2,j) = gi / h
            f = f + e(j) * ar(j,i) - tau(2,j) * ai(j,i)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do 260 j = 1, l
            f = ar(j,i)
            g = e(j) - hh * f
            e(j) = g
            fi = -ai(j,i)
            gi = tau(2,j) - hh * fi
            tau(2,j) = -gi
c
            do 260 k = 1, j
               ar(k,j) = ar(k,j) - f * e(k) - g * ar(k,i)
     x                           + fi * tau(2,k) + gi * ai(k,i)
               ai(k,j) = ai(k,j) - f * tau(2,k) - g * ai(k,i)
     x                           - fi * e(k) - gi * ar(k,i)
  260    continue
c
  270    do 280 k = 1, l
            ar(k,i) = scale * ar(k,i)
            ai(k,i) = scale * ai(k,i)
  280    continue
c
         tau(2,l) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)
C
      INTEGER I,J,K,L,M,N,II,MML,TAG,IERR
      DOUBLE PRECISION D(N),E(N),E2(N),W(N),RV1(N)
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2,PYTHAG
      INTEGER IND(N)
C
C     THIS ROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION OF
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS ROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM
C     THEIR CORRESPONDING SUBMATRIX INDICES.
C
C     ON INPUT
C
C        N IS THE ORDER OF THE MATRIX.
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C     ON OUTPUT
C
C        D AND E ARE UNALTERED.
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO.
C
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES.
C
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC..
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      K = 0
      TAG = 0
C
      DO 100 I = 1, N
         W(I) = D(I)
         IF (I .NE. 1) RV1(I-1) = E(I)
  100 CONTINUE
C
      E2(1) = 0.0D0
      RV1(N) = 0.0D0
C
      DO 290 L = 1, N
         J = 0
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = DABS(W(M)) + DABS(W(M+1))
            TST2 = TST1 + DABS(RV1(M))
            IF (TST2 .EQ. TST1) GO TO 120
C     .......... GUARD AGAINST UNDERFLOWED ELEMENT OF E2 ..........
            IF (E2(M+1) .EQ. 0.0D0) GO TO 125
  110    CONTINUE
C
  120    IF (M .LE. K) GO TO 130
         IF (M .NE. N) E2(M+1) = 0.0D0
  125    K = M
         TAG = TAG + 1
  130    P = W(L)
         IF (M .EQ. L) GO TO 215
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         G = (W(L+1) - P) / (2.0D0 * RV1(L))
         R = PYTHAG(G,1.0D0)
         G = W(M) - P + RV1(L) / (G + DSIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            I = M - II
            F = S * RV1(I)
            B = C * RV1(I)
            R = PYTHAG(F,G)
            RV1(I+1) = R
            IF (R .EQ. 0.0D0) GO TO 210
            S = F / R
            C = G / R
            G = W(I+1) - P
            R = (W(I) - G) * S + 2.0D0 * C * B
            P = S * R
            W(I+1) = G + P
            G = C * R - B
  200    CONTINUE
C
         W(L) = W(L) - P
         RV1(L) = G
         RV1(M) = 0.0D0
         GO TO 105
C     .......... RECOVER FROM UNDERFLOW ..........
  210    W(I+1) = W(I+1) - P
         RV1(M) = 0.0D0
         GO TO 105
C     .......... ORDER EIGENVALUES ..........
  215    IF (L .EQ. 1) GO TO 250
C     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
         DO 230 II = 2, L
            I = L + 2 - II
            IF (P .GE. W(I-1)) GO TO 270
            W(I) = W(I-1)
            IND(I) = IND(I-1)
  230    CONTINUE
C
  250    I = 1
  270    W(I) = P
         IND(I) = TAG
  290 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
      SUBROUTINE BINVIT(N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
C*
C*
C*    AUTHORS -
C*       THIS ROUTINE IS A MODIFICATION OF EISPACK ROUTINE TINVIT
C*       WHICH IS A TRANSLATION OF THE INVERSE ITERATION TECHNIQUE
C*       IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C*       HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C*       THE EISPACK VERSION WAS MODIFIED BY C. MOLER AND D. SPANGLER
C*       (NRCC - 1 APR 1979).  THIS VERSION IS BY S. T. ELBERT
C*       (AMES LABORATORY - USDOE), WHO FURTHER STREAMLINED THE CODE.
C*
C*    PURPOSE -
C*       THIS ROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C*       SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C*       USING INVERSE ITERATION.
C*
C*    METHOD -
C*       THE CALCULATIONS PROCEED AS FOLLOWS.  FIRST, THE  E2  ARRAY
C*       IS INSPECTED FOR THE PRESENCE OF A ZERO ELEMENT DEFINING A
C*       SUBMATRIX.  THE EIGENVALUES BELONGING TO THIS SUBMATRIX ARE
C*       IDENTIFIED BY THEIR COMMON SUBMATRIX INDEX IN  IND.
C*
C*       THE EIGENVECTORS OF THE SUBMATRIX ARE THEN COMPUTED BY
C*       INVERSE ITERATION.  FIRST, THE  LU  DECOMPOSITION OF THE SUB-
C*       MATRIX WITH AN APPROXIMATE EIGENVALUE SUBTRACTED FROM ITS
C*       DIAGONAL ELEMENTS IS ACHIEVED BY GAUSSIAN ELIMINATION USING
C*       PARTIAL PIVOTING.  THE MULTIPLIERS DEFINING THE LOWER
C*       TRIANGULAR MATRIX  L  ARE STORED IN THE TEMPORARY ARRAY  RV4
C*       AND THE UPPER TRIANGULAR MATRIX  U  IS STORED IN THE THREE
C*       TEMPORARY ARRAYS  RV1,  RV2, AND  RV3.  SAVING THESE
C*       QUANTITIES IN RV1,  RV2,  RV3, AND  RV4 AVOIDS REPEATING
C*       THE  LU  DECOMPOSITION IF FURTHER ITERATIONS ARE REQUIRED.
C*       AN APPROXIMATE VECTOR, STORED IN  RV6, IS COMPUTED STARTING
C*       FROM AN INITIAL VECTOR, AND THE NORM OF THE APPROXIMATE
C*       VECTOR IS COMPARED WITH A NORM OF THE SUBMATRIX TO DETERMINE
C*       WHETHER THE GROWTH IS SUFFICIENT TO ACCEPT IT AS AN
C*       EIGENVECTOR.  IF THIS VECTOR IS ACCEPTED, ITS EUCLIDIAN NORM
C*       IS MADE 1.  IF THE GROWTH IS NOT SUFFICIENT, THIS VECTOR IS
C*       USED AS THE INITIAL VECTOR IN COMPUTING THE NEXT APPROXIMATE
C*       VECTOR.  THIS ITERATION PROCESS IS REPEATED AT MOST  5  TIMES.
C*
C*       EIGENVECTORS COMPUTED IN THE ABOVE WAY CORRESPONDING TO WELL-
C*       SEPARATED EIGENVALUES OF THIS SUBMATRIX WILL BE ORTHOGONAL.
C*       HOWEVER, EIGENVECTORS CORRESPONDING TO CLOSE EIGENVALUES OF
C*       THIS SUBMATRIX MAY NOT BE SATISFACTORILY ORTHOGONAL.  HENCE,
C*       TO INSURE ORTHOGONAL EIGENVECTORS, EACH APPROXIMATE VECTOR IS
C*       MADE ORTHOGONAL TO THOSE PREVIOUSLY COMPUTED EIGENVECTORS
C*       WHOSE EIGENVALUES ARE CLOSE TO THE CURRENT EIGENVALUE, AS
C*       DETERMINED BY THE PARAMETER GRPTOL.
C*       IF THE ORTHOGONALIZATION PROCESS PRODUCES A ZERO VECTOR, A
C*       COLUMN OF THE IDENTITY MATRIX IS USED AS AN INITIAL VECTOR
C*       FOR THE NEXT ITERATION.
C*
C*       IDENTICAL EIGENVALUES ARE PERTURBED SLIGHTLY IN AN ATTEMPT TO
C*       OBTAIN INDEPENDENT EIGENVECTORS.  THESE PERTURBATIONS ARE NOT
C*       RECORDED IN THE EIGENVALUE ARRAY  W.
C*
C*       THE ABOVE STEPS ARE REPEATED ON EACH SUBMATRIX UNTIL ALL THE
C*       EIGENVECTORS ARE COMPUTED.
C*
C*
C*    ON ENTRY-
C*
C*       N      - INTEGER
C*                THE ORDER OF THE MATRIX.
C*
C*       D      - W.P. REAL (N)
C*                THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.
C*
C*       E      - W.P. REAL (N)
C*                CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C*                IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.
C*
C*       E2     - W.P. REAL (N)
C*          CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C*          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C*          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C*          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C*          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C*          0.0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0
C*          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C*          TRIDIB, OR  EIMQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C*          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE.
C*
C*       M      - INTEGER
C*                THE NUMBER OF SPECIFIED EIGENVALUES.
C*
C*       W      - W.P. REAL (M)
C*                THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER.
C*
C*       IND    - INTEGER (M)
C*           CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C*           ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C*           1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C*           TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C*
C*    ON EXIT-
C*
C*       Z      - W.P. REAL (NM,M)
C*          CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C*          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO.
C*
C*       IERR   - INTEGER
C*                SET TO
C*                ZERO   FOR NORMAL RETURN,
C*                -R     IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C*                       EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS.
C*
C*       RV1, RV2, RV3, RV4, AND RV6  -  W.P. REAL (N)
C*                TEMPORARY STORAGE ARRAYS.
C*
C*    COMPLEXITY -
C*       O(N*M)
C*
C*    NOTE -
C*    QUESTIONS AND COMMENTS ABOUT EISPACK SHOULD BE DIRECTED TO
C*    B. S. GARBOW, APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LAB.
C*
C***********************************************************************
C
C                    DECLARATIONS
C
      INTEGER IND(M)
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(N,M)
     *                ,RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      DOUBLE PRECISION ANORM,ANORM2,EPMACH,EPS2,EPS3,EPS4,ORDER,DDOT
     *                ,TEST,U,UK,V,XU,RTPRE,RTCUR ,ABS,PARM1,PARM2
C
C-----------------------------------------------------------------------
C
      IERR = 0
      IF (M .LE. 0) RETURN
C
C     ********** EPMACH IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
*ckk
*      EPMACH = 0.5D0**20   !Original
      EPMACH = 0.5D0**20
      DO 110 I=1,100
         EPMACH = EPMACH*0.5D0
         TEST = 1.D0 + EPMACH
         IF(TEST .EQ. 1.D0) GO TO 120
  110 CONTINUE
  120 CONTINUE
      EPMACH=EPMACH+EPMACH
      EPMACH=EPMACH+EPMACH
C
      ITAG = 0
      ORDER = 1.D0 - E2(1)
      IP = 1
      DO 930 IQ=1,N
         IF (IQ.NE.N .AND. E2(IQ+1) .NE. 0.D0) GO TO 930
C
C           ********** ESTABLISH AND PROCESS NEXT SUBMATRIX **********
C                           WITH INDEX RANGE  IP  TO  IQ
C
         ITAG = ITAG + 1
         IS = 0
         RTCUR=W(1)
C
         DO 920 IR = 1, M
            IF (IND(IR) .NE. ITAG) GO TO 920
               RTPRE = RTCUR
               RTCUR = W(IR)
               IPP1 = IP + 1
               IPM1 = IP - 1
               NSUB = IQ - IPM1
               IF (IS .NE. 0) GO TO 510
C
C              ********** CHECK FOR ISOLATED ROOT **********
C
                  XU = 1.D0
                  IF (IP .NE. IQ) GO TO 490
C
C                        SET ISOLATED ROOT TO UNIT VECTOR
C
                     DO 410 I=1,N
                        Z(I,IR)=0.D0
  410                CONTINUE
                     Z(IP,IR)=1.D0
                     GO TO 920
C
  490             CONTINUE
                  ANORM = DABS(D(IP))
C
                  DO 500 I = IPP1, IQ
                     ANORM = ANORM + DABS(D(I)) + DABS(E(I))
  500             CONTINUE
C
C                 ********** EPS2 IS THE CRITERION FOR GROUPING,
C                            EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                            ROOTS ARE MODIFIED BY EPS3,
C                            EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW
C
                  UK = DBLE(NSUB)
                  ANORM2 = ANORM / UK

* Original!
ckk                  EPS2 = 1.D-5 * ANORM2
ckk                  EPS3 = EPMACH * ANORM
ckk                  EPS4 = UK * EPS3

                  EPS2 = 1.D-9 * ANORM2
                  EPS3 =  5.10000001d0*( EPMACH * ANORM )
                  EPS4 = 100.d0*( UK * EPS3 )


                  UK = EPS4 / DSQRT(UK)
                  IS = IP
                  IGROUP = 0
                  GO TO 520
C
C              ********** LOOK FOR CLOSE OR COINCIDENT ROOTS **********
C
  510          CONTINUE
               IGROUP = IGROUP + 1
               IF (DABS(RTCUR-RTPRE) .GE. EPS2) IGROUP =0
               IF (ORDER * (RTCUR - RTPRE) .LE. 0.D0)
     *            RTCUR = RTPRE + ORDER * EPS3
C
C              ********** ELIMINATION WITH INTERCHANGES AND
C                         INITIALIZATION OF VECTOR **********
C
  520          CONTINUE
               V = 0.D0
               IQM1 = IQ - 1
               RV6(IP) = UK
               U = D(IP) - RTCUR
               IF(IP .EQ. IQ) GO TO 595
                  V = E(IPP1)
                  IF (IPP1 .GT. IQM1) GO TO 590
                     DO 580 I = IPP1,IQM1
                        RV6(I) = UK
                        IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
C                         **** WARNING -- A DIVIDE CHECK MAY OCCUR HERE
C                          IF E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY
                           XU = U / E(I)
                           RV4(I) = XU
                           RV1(I-1) = E(I)
                           RV2(I-1) = D(I) - RTCUR
                           RV3(I-1) = E(I+1)
                           U = V - XU * RV2(I-1)
                           V = -XU * RV3(I-1)
                           GO TO 580
C
  540                   CONTINUE
                        XU = E(I) / U
                        RV4(I) = XU
                        RV1(I-1) = U
                        RV2(I-1) = V
                        RV3(I-1) = 0.D0
                        U = D(I) - RTCUR - XU * V
                        V = E(I+1)
  580                CONTINUE
  590             CONTINUE
                  RV6(IQ) = UK
                  RV3(IQM1) = 0.D0
                  IF (DABS(E(IQ)) .LT. DABS(U)) GO TO 592
                     XU = U / E(IQ)
                     RV4(IQ) = XU
                     RV1(IQM1) = E(IQ)
                     RV2(IQM1) = D(IQ) - RTCUR
                     U = V - XU * RV2(IQM1)
                     V = 0.D0
                     GO TO 595
C
  592             CONTINUE
                  XU = E(IQ) / U
                  RV4(IQ) = XU
                  RV1(IQM1) = U
                  RV2(IQM1) = V
                  U = D(IQ) - RTCUR - XU * V
  595          CONTINUE
C
               IF (DABS(U) .LT. EPS3) U = EPS3
               RV1(IQ) = U
               RV2(IQ) = 0.D0
               RV3(IQ) = 0.D0
C
C              ********** BACK SUBSTITUTION
C                         FOR I=IQ STEP -1 UNTIL IP DO -- **********
Cckk standard 5 iter.
*               DO 830 ITERS=1,5
               DO 830 ITERS=1,25
                  RV6(IQ) = RV6(IQ)/RV1(IQ)
                  U = RV6(IQ)
                  RV6(IQM1) = (RV6(IQM1) - U * RV2(IQM1)) / RV1(IQM1)
                  V = U
                  U = RV6(IQM1)
                  ANORM = DABS(U) + DABS(V)
                  IF (IPP1 .GT. IQM1) GO TO 625
                    IQM2=IQM1-1
                    DO 620 I = IQM2,IP,-1
                      RV6(I) = (RV6(I) - U*RV2(I) - V*RV3(I)) / RV1(I)
                      V = U
                      U = RV6(I)
                      ANORM = ANORM + DABS(U)
  620               CONTINUE
  625             CONTINUE
                  IF (IGROUP .EQ. 0) GO TO 700
C
C                    ********** ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                               MEMBERS OF GROUP **********
C
                     J = IR
                     DO 680 JJ = 1,IGROUP
  630                   J = J - 1
                        IF (IND(J) .NE. ITAG) GO TO 630
                        xu = 0d0
                        nsubup = nsub+ip-1
                        do 640 i = ip, nsubup
                          xu = xu + rv6(i)*z(i,j)
  640                   continue
                        xu = -xu
                        if (xu .eq. 0d0) goto 680
                        do 650 i = ip, nsubup
                          rv6(i) = rv6(i) + xu*z(i,j)
  650                   continue
  680                CONTINUE
                     ANORM = 0.D0
                     DO 690 I = IP,IQ
                        ANORM = ANORM + DABS(RV6(I))
  690                CONTINUE
C
  700             CONTINUE
                  IF (ANORM .GE. 1.D0) GO TO 840
C
C                    ********** FORWARD SUBSTITUTION **********
C
                     IF (ANORM .NE. 0.D0) GO TO 740
                        RV6(IS) = EPS4
                        IS = IS + 1
                        IF (IS .GT. IQ) IS = IP
                        GO TO 780
C
  740                CONTINUE
                     XU = EPS4 / ANORM
C
                     DO 760 I =IP,IQ
                        RV6(I) = RV6(I) * XU
  760                CONTINUE
C
C                    ********** ELIMINATION OPERATIONS ON NEXT VECTOR
C                               ITERATE **********
C
  780                CONTINUE
                     DO 820 I = IPP1,IQ
                        U = RV6(I)
C
C                   ********** IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                              WAS PERFORMED EARLIER IN THE
C                              TRIANGULARIZATION PROCESS **********
C
                        IF (RV1(I-1) .NE. E(I)) GO TO 800
                           U = RV6(I-1)
                           RV6(I-1) = RV6(I)
  800                   CONTINUE
                        RV6(I) = U - RV4(I) * RV6(I-1)
  820                CONTINUE
  830          CONTINUE
C
C                 ********** SET ERROR -- NON-CONVERGED EIGENVECTOR
C
               IERR = -IR
               DO 835 I=1,N
                  Z(I,IR)=0.D0
  835          CONTINUE
               GO TO 920
C
C              ********** NORMALIZE SO THAT SUM OF SQUARES IS
C                         1 AND EXPAND TO FULL ORDER **********
C
  840          CONTINUE
               U = 0.D0
               DO 860 I = IP,IQ
                  U = U + RV6(I)**2
  860          CONTINUE
               XU = 1.D0 / DSQRT(U)
C
               IF(IPM1.LT.1) GO TO 880
                  DO 870 I = 1, IPM1
                     Z(I,IR) = 0.D0
  870             CONTINUE
  880          CONTINUE
               DO 890 I = IP,IQ
                  Z(I,IR) = RV6(I) * XU
  890          CONTINUE
               IQP1 = IQ + 1
               IF(IQP1 .GT. N) GO TO 905
                  DO 900 I=IQP1,N
                     Z(I,IR) = 0.D0
  900             CONTINUE
  905          CONTINUE
  920    CONTINUE
         IP = IQ + 1
  930 CONTINUE
      RETURN
      END
      subroutine btribk(n,ar,ai,tau,m,zr,zi)
c
      integer i,j,k,l,m,n
      double precision ar(n,n),ai(n,n),tau(2,n),zr(n,m),zi(n,m)
      double precision h,s,si,smin
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure trbak1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine forms the eigenvectors of a complex hermitian
c     matrix by back transforming those of the corresponding
c     real symmetric tridiagonal matrix determined by  btridi.
c
c     on input
c
c        n is the order of the matrix.
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction by  btridi  in their
c          full upper triangles except for the diagonal of ar.
c
c        tau contains further information about the transformations.
c
c        m is the number of eigenvectors to be back transformed.
c
c        zr contains the eigenvectors to be back transformed
c          in its first m columns.
c
c     on output
c
c        zr and zi contain the real and imaginary parts,
c          respectively, of the transformed eigenvectors
c          in their first m columns.
c
c     note that the last component of each returned vector
c     is real and that vector euclidean norms are preserved.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version is adapted from the eispack htridi. In contrast
c     only the upper triangle of a is used for input
c
c     -----------------------------------------------------------------
c
      if (m .eq. 0) go to 200
c     .......... transform the eigenvectors of the real symmetric
c                tridiagonal matrix to those of the hermitian
c                tridiagonal matrix. ..........

      do 50 k = 1, n
c
         do 50 j = 1, m
            zi(k,j) = -zr(k,j) * tau(2,k)
            zr(k,j) = zr(k,j) * tau(1,k)
   50 continue
c
      if (n .eq. 1) go to 200
c     .......... recover and apply the householder matrices ......... .
      do 140 i = 2, n
         l = i - 1
         h = ai(i,i)
         if (h .eq. 0.0d0) go to 140
c
         do 130 j = 1, m
            s = 0.0d0
            si = 0.0d0
c
            do 110 k = 1, l
               s = s + ar(k,i) * zr(k,j) - ai(k,i) * zi(k,j)
               si = si + ar(k,i) * zi(k,j) + ai(k,i) * zr(k,j)
  110       continue
c     .......... double divisions avoid possible underflow ..........
            s = (s / h) / h
            smin = -s
            si = - (si / h) / h
c
            do 120 k = 1, l
               zr(k,j) = zr(k,j) + smin * ar(k,i) + si * ai(k,i)
               zi(k,j) = zi(k,j) + si * ar(k,i) + s * ai(k,i)
  120       continue
c
  130    continue
c
  140 continue
c
  200 return
      end
C#elseC
C      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
Cc
C      integer i,j,k,l,n,ii,nm,jp1
C      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
C      double precision f,g,h,fi,gi,hh,si,scale,pythag
Cc
Cc     this subroutine is a translation of a complex analogue of
Cc     the algol procedure tred1, num. math. 11, 181-195(1968)
Cc     by martin, reinsch, and wilkinson.
Cc     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
Cc
Cc     this subroutine reduces a complex hermitian matrix
Cc     to a real symmetric tridiagonal matrix using
Cc     unitary similarity transformations.
Cc
Cc     on input
Cc
Cc        nm must be set to the row dimension of two-dimensional
Cc          array parameters as declared in the calling program
Cc          dimension statement.
Cc
Cc        n is the order of the matrix.
Cc
Cc        ar and ai contain the real and imaginary parts,
Cc          respectively, of the complex hermitian input matrix.
Cc          only the lower triangle of the matrix need be supplied.
Cc
Cc     on output
Cc
Cc        ar and ai contain information about the unitary trans-
Cc          formations used in the reduction in their full lower
Cc          triangles.  their strict upper triangles and the
Cc          diagonal of ar are unaltered.
Cc
Cc        d contains the diagonal elements of the the tridiagonal matrix
Cc
Cc        e contains the subdiagonal elements of the tridiagonal
Cc          matrix in its last n-1 positions.  e(1) is set to zero.
Cc
Cc        e2 contains the squares of the corresponding elements of e.
Cc          e2 may coincide with e if the squares are not needed.
Cc
Cc        tau contains further information about the transformations.
Cc
Cc     calls pythag for  dsqrt(a*a + b*b) .
Cc
Cc     questions and comments should be directed to burton s. garbow,
Cc     mathematics and computer science div, argonne national laboratory
Cc
Cc     this version dated august 1983.
Cc
Cc     -----------------------------------------------------------------
Cc
C      tau(1,n) = 1.0d0
C      tau(2,n) = 0.0d0
Cc
C      do 100 i = 1, n
C  100 d(i) = ar(i,i)
Cc     .......... for i=n step -1 until 1 do -- ..........
C      do 300 ii = 1, n
C         i = n + 1 - ii
C         l = i - 1
C         h = 0.0d0
C         scale = 0.0d0
C         if (l .lt. 1) go to 130
Cc     .......... scale row (algol tol then not needed) ..........
C         do 120 k = 1, l
C  120    scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
Cc
C         if (scale .ne. 0.0d0) go to 140
C         tau(1,l) = 1.0d0
C         tau(2,l) = 0.0d0
C  130    e(i) = 0.0d0
C         e2(i) = 0.0d0
C         go to 290
Cc
C  140    do 150 k = 1, l
C            ar(i,k) = ar(i,k) / scale
C            ai(i,k) = ai(i,k) / scale
C            h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
C  150    continue
Cc
C         e2(i) = scale * scale * h
C         g = dsqrt(h)
C         e(i) = scale * g
C         f = pythag(ar(i,l),ai(i,l))
Cc     .......... form next diagonal element of matrix t ......... .
C         if (f .eq. 0.0d0) go to 160
C         tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
C         si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
C         h = h + f * g
C         g = 1.0d0 + g / f
C         ar(i,l) = g * ar(i,l)
C         ai(i,l) = g * ai(i,l)
C         if (l .eq. 1) go to 270
C         go to 170
C  160    tau(1,l) = -tau(1,i)
C         si = tau(2,i)
C         ar(i,l) = g
C  170    f = 0.0d0
Cc
C         do 240 j = 1, l
C            g = 0.0d0
C            gi = 0.0d0
Cc     .......... form element of a*u ..........
C            do 180 k = 1, j
C               g = g + ar(j,k) * ar(i,k) + ai(j,k) * ai(i,k)
C               gi = gi - ar(j,k) * ai(i,k) + ai(j,k) * ar(i,k)
C  180       continue
Cc
C            jp1 = j + 1
C            if (l .lt. jp1) go to 220
Cc
C            do 200 k = jp1, l
C               g = g + ar(k,j) * ar(i,k) - ai(k,j) * ai(i,k)
C               gi = gi - ar(k,j) * ai(i,k) - ai(k,j) * ar(i,k)
C  200       continue
Cc     .......... form element of p ..........
C  220       e(j) = g / h
C            tau(2,j) = gi / h
C            f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
C  240    continue
Cc
C         hh = f / (h + h)
Cc     .......... form reduced a ..........
C         do 260 j = 1, l
C            f = ar(i,j)
C            g = e(j) - hh * f
C            e(j) = g
C            fi = -ai(i,j)
C            gi = tau(2,j) - hh * fi
C            tau(2,j) = -gi
Cc
C            do 260 k = 1, j
C               ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
C     x                           + fi * tau(2,k) + gi * ai(i,k)
C               ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
C     x                           - fi * e(k) - gi * ar(i,k)
C  260    continue
Cc
C  270    do 280 k = 1, l
C            ar(i,k) = scale * ar(i,k)
C            ai(i,k) = scale * ai(i,k)
C  280    continue
Cc
C         tau(2,l) = -si
C  290    hh = d(i)
C         d(i) = ar(i,i)
C         ar(i,i) = hh
C         ai(i,i) = scale * dsqrt(h)
C  300 continue
Cc
C      return
C      end
C      subroutine imtql2(nm,n,d,e,z,ierr)
Cc
C      integer i,j,k,l,m,n,ii,nm,mml,ierr
C      double precision d(n),e(n),z(nm,n)
C      double precision b,c,f,g,p,r,s,tst1,tst2,pythag
Cc
Cc     this subroutine is a translation of the algol procedure imtql2,
Cc     num. math. 12, 377-383(1968) by martin and wilkinson,
Cc     as modified in num. math. 15, 450(1970) by dubrulle.
Cc     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
Cc
Cc     this subroutine finds the eigenvalues and eigenvectors
Cc     of a symmetric tridiagonal matrix by the implicit ql method.
Cc     the eigenvectors of a full symmetric matrix can also
Cc     be found if  tred2  has been used to reduce this
Cc     full matrix to tridiagonal form.
Cc
Cc     on input
Cc
Cc        nm must be set to the row dimension of two-dimensional
Cc          array parameters as declared in the calling program
Cc          dimension statement.
Cc
Cc        n is the order of the matrix.
Cc
Cc        d contains the diagonal elements of the input matrix.
Cc
Cc        e contains the subdiagonal elements of the input matrix
Cc          in its last n-1 positions.  e(1) is arbitrary.
Cc
Cc        z contains the transformation matrix produced in the
Cc          reduction by  tred2, if performed.  if the eigenvectors
Cc          of the tridiagonal matrix are desired, z must contain
Cc          the identity matrix.
Cc
Cc      on output
Cc
Cc        d contains the eigenvalues in ascending order.  if an
Cc          error exit is made, the eigenvalues are correct but
Cc          unordered for indices 1,2,...,ierr-1.
Cc
Cc        e has been destroyed.
Cc
Cc        z contains orthonormal eigenvectors of the symmetric
Cc          tridiagonal (or full) matrix.  if an error exit is made,
Cc          z contains the eigenvectors associated with the stored
Cc          eigenvalues.
Cc
Cc        ierr is set to
Cc          zero       for normal return,
Cc          j          if the j-th eigenvalue has not been
Cc                     determined after 30 iterations.
Cc
Cc     calls pythag for  dsqrt(a*a + b*b) .
Cc
Cc     questions and comments should be directed to burton s. garbow,
Cc     mathematics and computer science div, argonne national laboratory
Cc
Cc     this version dated august 1983.
Cc
Cc     -----------------------------------------------------------------
Cc
C      ierr = 0
C      if (n .eq. 1) go to 1001
Cc
C      do 100 i = 2, n
C  100 e(i-1) = e(i)
Cc
C      e(n) = 0.0d0
Cc
C      do 240 l = 1, n
C         j = 0
Cc     .......... look for small sub-diagonal element ..........
C  105    do 110 m = l, n
C            if (m .eq. n) go to 120
C            tst1 = dabs(d(m)) + dabs(d(m+1))
C            tst2 = tst1 + dabs(e(m))
C            if (tst2 .eq. tst1) go to 120
C  110    continue
Cc
C  120    p = d(l)
C         if (m .eq. l) go to 240
C         if (j .eq. 30) go to 1000
C         j = j + 1
Cc     .......... form shift ..........
C         g = (d(l+1) - p) / (2.0d0 * e(l))
C         r = pythag(g,1.0d0)
C         g = d(m) - p + e(l) / (g + dsign(r,g))
C         s = 1.0d0
C         c = 1.0d0
C         p = 0.0d0
C         mml = m - l
Cc     .......... for i=m-1 step -1 until l do -- ..........
C         do 200 ii = 1, mml
C            i = m - ii
C            f = s * e(i)
C            b = c * e(i)
C            r = pythag(f,g)
C            e(i+1) = r
C            if (r .eq. 0.0d0) go to 210
C            s = f / r
C            c = g / r
C            g = d(i+1) - p
C            r = (d(i) - g) * s + 2.0d0 * c * b
C            p = s * r
C            d(i+1) = g + p
C            g = c * r - b
Cc     .......... form vector ..........
C            do 180 k = 1, n
C               f = z(k,i+1)
C               z(k,i+1) = s * z(k,i) + c * f
C               z(k,i) = c * z(k,i) - s * f
C  180       continue
Cc
C  200    continue
Cc
C         d(l) = d(l) - p
C         e(l) = g
C         e(m) = 0.0d0
C         go to 105
Cc     .......... recover from underflow ..........
C  210    d(i+1) = d(i+1) - p
C         e(m) = 0.0d0
C         go to 105
C  240 continue
Cc     .......... order eigenvalues and eigenvectors ..........
C      do 300 ii = 2, n
C         i = ii - 1
C         k = i
C         p = d(i)
Cc
C         do 260 j = ii, n
C            if (d(j) .ge. p) go to 260
C            k = j
C            p = d(j)
C  260    continue
Cc
C         if (k .eq. i) go to 300
C         d(k) = d(i)
C         d(i) = p
Cc
C         do 280 j = 1, n
C            p = z(j,i)
C            z(j,i) = z(j,k)
C            z(j,k) = p
C  280    continue
Cc
C  300 continue
Cc
C      go to 1001
Cc     .......... set error -- no convergence to an
Cc                eigenvalue after 30 iterations ..........
C 1000 ierr = l
C 1001 return
C      end
C      subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)
Cc
C      integer i,j,k,l,m,n,nm
C      double precision ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
C      double precision h,s,si
Cc
Cc     this subroutine is a translation of a complex analogue of
Cc     the algol procedure trbak1, num. math. 11, 181-195(1968)
Cc     by martin, reinsch, and wilkinson.
Cc     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
Cc
Cc     this subroutine forms the eigenvectors of a complex hermitian
Cc     matrix by back transforming those of the corresponding
Cc     real symmetric tridiagonal matrix determined by  htridi.
Cc
Cc     on input
Cc
Cc        nm must be set to the row dimension of two-dimensional
Cc          array parameters as declared in the calling program
Cc          dimension statement.
Cc
Cc        n is the order of the matrix.
Cc
Cc        ar and ai contain information about the unitary trans-
Cc          formations used in the reduction by  htridi  in their
Cc          full lower triangles except for the diagonal of ar.
Cc
Cc        tau contains further information about the transformations.
Cc
Cc        m is the number of eigenvectors to be back transformed.
Cc
Cc        zr contains the eigenvectors to be back transformed
Cc          in its first m columns.
Cc
Cc     on output
Cc
Cc        zr and zi contain the real and imaginary parts,
Cc          respectively, of the transformed eigenvectors
Cc          in their first m columns.
Cc
Cc     note that the last component of each returned vector
Cc     is real and that vector euclidean norms are preserved.
Cc
Cc     questions and comments should be directed to burton s. garbow,
Cc     mathematics and computer science div, argonne national laboratory
Cc
Cc     this version dated august 1983.
Cc
Cc     -----------------------------------------------------------------
Cc
C      if (m .eq. 0) go to 200
Cc     .......... transform the eigenvectors of the real symmetric
Cc                tridiagonal matrix to those of the hermitian
Cc                tridiagonal matrix. ..........
C      do 50 k = 1, n
Cc
C         do 50 j = 1, m
C            zi(k,j) = -zr(k,j) * tau(2,k)
C            zr(k,j) = zr(k,j) * tau(1,k)
C   50 continue
Cc
C      if (n .eq. 1) go to 200
Cc     .......... recover and apply the householder matrices ......... .
C      do 140 i = 2, n
C         l = i - 1
C         h = ai(i,i)
C         if (h .eq. 0.0d0) go to 140
Cc
C         do 130 j = 1, m
C            s = 0.0d0
C            si = 0.0d0
Cc
C            do 110 k = 1, l
C               s = s + ar(i,k) * zr(k,j) - ai(i,k) * zi(k,j)
C               si = si + ar(i,k) * zi(k,j) + ai(i,k) * zr(k,j)
C  110       continue
Cc     .......... double divisions avoid possible underflow ..........
C            s = (s / h) / h
C            si = (si / h) / h
Cc
C            do 120 k = 1, l
C               zr(k,j) = zr(k,j) - s * ar(i,k) - si * ai(i,k)
C               zi(k,j) = zi(k,j) - si * ar(i,k) + s * ai(i,k)
C  120       continue
Cc
C  130    continue
Cc
C  140 continue
Cc
C  200 return
C      end
C#endif
      double precision function pythag(a,b)
C
C     Finds dsqrt(a**2+b**2) without overflow or destructive underflow
C
      double precision a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
C#endif

* blas.f below.
C#ifdefC SLATEC
C      subroutine daxpy(n,da,dx,incx,dy,incy)
C      call saxpy(n,da,dx,incx,dy,incy)
C      end
C      subroutine dcopy(n,dx,incx,dy,incy)
C      call scopy(n,dx,incx,dy,incy)
C      end
C      subroutine drot (n,dx,incx,dy,incy,c,s)
C      call srot (n,dx,incx,dy,incy,c,s)
C      end
C      subroutine drotg(da,db,c,s)
C      call srotg(da,db,c,s)
C      end
C      subroutine dscal(n,da,dx,incx)
C      call sscal(n,da,dx,incx)
C      end
C      subroutine dswap (n,dx,incx,dy,incy)
C      call sswap (n,dx,incx,dy,incy)
C      end
C      double precision function dasum(n,dx,incx)
C      dasum = sasum(n,dx,incx)
C      end
C      double precision function ddot(n,dx,incx,dy,incy)
C      ddot = sdot(n,dx,incx,dy,incy)
C      end
C      double precision function dnrm2 ( n, dx, incx)
C      dnrm2 = snrm2 ( n, dx, incx)
C      end
C      integer function idamax(n,dx,incx)
C      idamax = isamax(n,dx,incx)
C      end
C      double precision function dsum(n,dx,incx)
C      dsum = ssum(n,dx,incx)
C      end
C#else
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1)
      integer i,incx,n,nincx
c
      dasum = 0d0
      if (n .le. 0) return

      nincx = n*incx
      do  10  i = 1, nincx,incx
        dasum = dasum + dabs(dx(i))
   10 continue
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,n
c
      if (da .eq. 0d0) return
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (1-n)*incx + 1
      if (incy .lt. 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      end
      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,n
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (1-n)*incx + 1
      if (incy .lt. 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,n
c
      ddot = 0d0
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (1-n)*incx + 1
      if (incy .lt. 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = ddot
      end
      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end

      subroutine drot (n,dx,incx,dy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (1-n)*incx + 1
      if(incy.lt.0)iy = (1-n)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return
      end
      subroutine drotg(da,db,c,s)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c                    modified 9/27/86.
c
      double precision da,db,c,s,roe,scale,r,z
c
      roe = db
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
      if( scale .ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         go to 20
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
   20 z = s
      if( dabs(c) .gt. 0.0d0 .and. dabs(c) .le. s ) z = 1.0d0/c
      da = r
      db = z
      return
      end
      subroutine dscal(n,da,dx,incx)
c
c     scales a vector by a constant.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision da,dx(1)
      integer i,incx,n,nincx
c
      if (n .le. 0) return
      nincx = n*incx
      do  10  i = 1, nincx,incx
        dx(i) = da*dx(i)
   10 continue
      end
      subroutine dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors. Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,n
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (1-n)*incx + 1
      if (incy.lt.0) iy = (1-n)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. abs. value.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if ( n .lt. 1 ) return
      idamax = 1
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do  10  i = 2, n
         if (dabs(dx(ix)) .le. dmax) goto 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      end
      double precision function dsum(n,dx,incx)
c
c     takes the sum of the values.  Adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1)
      integer i,incx,n,nincx
c
      dsum = 0d0
      if (n .le. 0) return

      nincx = n*incx
      do  10  i = 1, nincx,incx
        dsum = dsum + dx(i)
   10 continue
      end
C#endif

* linpack.f below.
C#ifdefC CRAY
C      subroutine dgedi(a,lda,n,ipvt,det,work,job)
C      call sgedi(a,lda,n,ipvt,det,work,job)
C      end
C      subroutine dgefa(a,lda,n,ipvt,info)
C      call sgefa(a,lda,n,ipvt,info)
C      end
C      subroutine dgesl(a,lda,n,ipvt,b,job)
C      call sgesl(a,lda,n,ipvt,b,job)
C      end
C      subroutine dsidi(a,lda,n,kpvt,det,inert,work,job)
C      call ssidi(a,lda,n,kpvt,det,inert,work,job)
C      end
C      subroutine dsifa(a,lda,n,kpvt,info)
C      call ssifa(a,lda,n,kpvt,info)
C      end
C      subroutine dsisl(a,lda,n,kpvt,b)
C      call ssisl(a,lda,n,kpvt,b)
C      end
C      subroutine dspdi(ap,n,kpvt,det,inert,work,job)
C      call sspdi(ap,n,kpvt,det,inert,work,job)
C      end
C      subroutine dspfa(ap,n,kpvt,info)
C      call sspfa(ap,n,kpvt,info)
C      end
C      subroutine dspsl(ap,n,kpvt,b)
C      call sspsl(ap,n,kpvt,b)
C      end
C#else
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) goto 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do  50  i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) goto 60
   10       if (dabs(det(1)) .ge. 1.0d0) goto 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            goto 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) goto 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            goto 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) goto 150
         do  100  k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) goto 90
            do  80  j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) goto 140
         do  130  kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do  110  i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do  120  j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) goto 70
      do  60  k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) goto 40
c
c           interchange if necessary
c
            if (l .eq. k) goto 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do  30  j = kp1, n
               t = a(l,j)
               if (l .eq. k) goto 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         goto 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) goto ...
c           do  10  j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) goto 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) goto 30
         do  20  k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) goto 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do  40  kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      goto 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do  60  k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) goto 90
         do  80  kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) goto 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine dsidi(a,lda,n,kpvt,det,inert,work,job)
      integer lda,n,job
      double precision a(lda,1),work(1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     dsidi computes the determinant, inertia and inverse
c     of a double precision symmetric matrix using the factors from
c     dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array a.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        a      contains the upper triangle of the inverse of
c               the original matrix.  the strict lower triangle
c               is never referenced.
c
c        det    double precision(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. dabs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero may occur if the inverse is requested
c        and  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab
c
c     subroutines and functions
c
c     blas daxpy,dcopy,ddot,dswap
c     fortran dabs,iabs,mod
c
c     internal variables.
c
      double precision akkp1,ddot,temp
      double precision ten,d,t,ak,akp1
      integer j,jb,k,km1,ks,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) goto 140
         if (noert) goto 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) goto 20
            det(1) = 1.0d0
            det(2) = 0.0d0
            ten = 10.0d0
   20    continue
         t = 0.0d0
         do  130  k = 1, n
            d = a(k,k)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) goto 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0d0) goto 30
                  t = dabs(a(k,k+1))
                  d = (d/t)*a(k+1,k+1) - t
               goto 40
   30          continue
                  d = t
                  t = 0.0d0
   40          continue
   50       continue
c
            if (noert) goto 60
               if (d .gt. 0.0d0) inert(1) = inert(1) + 1
               if (d .lt. 0.0d0) inert(2) = inert(2) + 1
               if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) goto 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0d0) goto 110
   70             if (dabs(det(1)) .ge. 1.0d0) goto 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0d0
                  goto 70
   80             continue
   90             if (dabs(det(1)) .lt. ten) goto 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0d0
                  goto 90
  100             continue
  110          continue
  120       continue
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) goto 270
         k = 1
  150    if (k .gt. n) goto 260
            km1 = k - 1
            if (kpvt(k) .lt. 0) goto 180
c
c              1 by 1
c
               a(k,k) = 1.0d0/a(k,k)
               if (km1 .lt. 1) goto 170
                  call dcopy(km1,a(1,k),1,work,1)
                  do  160  j = 1, km1
                     a(j,k) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  160             continue
                  a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  170          continue
               kstep = 1
            goto 220
  180       continue
c
c              2 by 2
c
               t = dabs(a(k,k+1))
               ak = a(k,k)/t
               akp1 = a(k+1,k+1)/t
               akkp1 = a(k,k+1)/t
               d = t*(ak*akp1 - 1.0d0)
               a(k,k) = akp1/d
               a(k+1,k+1) = ak/d
               a(k,k+1) = -akkp1/d
               if (km1 .lt. 1) goto 210
                  call dcopy(km1,a(1,k+1),1,work,1)
                  do  190  j = 1, km1
                     a(j,k+1) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k+1),1)
  190             continue
                  a(k+1,k+1) = a(k+1,k+1) + ddot(km1,work,1,a(1,k+1),1)
                  a(k,k+1) = a(k,k+1) + ddot(km1,a(1,k),1,a(1,k+1),1)
                  call dcopy(km1,a(1,k),1,work,1)
                  do  200  j = 1, km1
                     a(j,k) = ddot(j,a(1,j),1,work,1)
                     call daxpy(j-1,work(j),a(1,j),1,a(1,k),1)
  200             continue
                  a(k,k) = a(k,k) + ddot(km1,work,1,a(1,k),1)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) goto 250
               call dswap(ks,a(1,ks),1,a(1,k),1)
               do  230  jb = ks, k
                  j = k + ks - jb
                  temp = a(j,k)
                  a(j,k) = a(ks,j)
                  a(ks,j) = temp
  230          continue
               if (kstep .eq. 1) goto 240
                  temp = a(ks,k+1)
                  a(ks,k+1) = a(k,k+1)
                  a(k,k+1) = temp
  240          continue
  250       continue
            k = k + kstep
         goto 150
  260    continue
  270 continue
      return
      end
      subroutine dsifa(a,lda,n,kpvt,info)
      integer lda,n,kpvt(1),info
      double precision a(lda,1)
c
c     dsifa factors a double precision symmetric matrix by elimination
c     with symmetric pivoting.
c
c     to solve  a*x = b , follow dsifa by dsisl.
c     to compute  inverse(a)*c , follow dsifa by dsisl.
c     to compute  determinant(a) , follow dsifa by dsidi.
c     to compute  inertia(a) , follow dsifa by dsidi.
c     to compute  inverse(a) , follow dsifa by dsidi.
c
c     on entry
c
c        a       double precision(lda,n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that dsisl or dsidi may
c                     divide by zero if called.
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dswap,idamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
      logical swap
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) goto 200
         if (k .gt. 1) goto 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0d0) info = 1
c     ......exit
            goto 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         absakk = dabs(a(k,k))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = idamax(k-1,a(1,k),1)
         colmax = dabs(a(imax,k))
         if (absakk .lt. alpha*colmax) goto 30
            kstep = 1
            swap = .false.
         goto 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do  40  j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40       continue
            if (imax .eq. 1) goto 50
               jmax = idamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50       continue
            if (dabs(a(imax,imax)) .lt. alpha*rowmax) goto 60
               kstep = 1
               swap = .true.
            goto 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) goto 70
               kstep = 1
               swap = .false.
            goto 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) goto 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         goto 190
  100    continue
         if (kstep .eq. 2) goto 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) goto 120
c
c              perform an interchange.
c
               call dswap(imax,a(1,imax),1,a(1,k),1)
               do  110  jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
c
c           perform the elimination.
c
            do  130  jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call daxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         goto 190
  140    continue
c
c           2 x 2 pivot block.
c
            if (.not.swap) goto 160
c
c              perform an interchange.
c
               call dswap(imax,a(1,imax),1,a(1,k-1),1)
               do  150  jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) goto 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do  170  jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call daxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      goto 10
  200 continue
      return
      end
      subroutine dsisl(a,lda,n,kpvt,b)
      integer lda,n,kpvt(1)
      double precision a(lda,1),b(1)
c
c     dsisl solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dsifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) goto ...
c           do  10  j = 1, p
c              call dsisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran iabs
c
c     internal variables.
c
      double precision ak,akm1,bk,bkm1,ddot,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
   10 if (k .eq. 0) goto 80
         if (kpvt(k) .lt. 0) goto 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) goto 30
               kp = kpvt(k)
               if (kp .eq. k) goto 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call daxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         goto 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) goto 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) goto 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call daxpy(k-2,b(k),a(1,k),1,b(1),1)
               call daxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      goto 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) goto 160
         if (kpvt(k) .lt. 0) goto 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) goto 110
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) goto 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         goto 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) goto 140
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) goto 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      goto 90
  160 continue
      return
      end
      subroutine dspdi(ap,n,kpvt,det,inert,work,job)
      integer n,job
      double precision ap(1),work(1)
      double precision det(2)
      integer kpvt(1),inert(3)
c
c     dspdi computes the determinant, inertia and inverse
c     of a double precision symmetric matrix using the factors from
c     dspfa, where the matrix is stored in packed form.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the output from dspfa.
c
c        n       integer
c                the order of the matrix a.
c
c        kpvt    integer(n)
c                the pivot vector from dspfa.
c
c        work    double precision(n)
c                work vector.  contents ignored.
c
c        job     integer
c                job has the decimal expansion  abc  where
c                   if  c .ne. 0, the inverse is computed,
c                   if  b .ne. 0, the determinant is computed,
c                   if  a .ne. 0, the inertia is computed.
c
c                for example, job = 111  gives all three.
c
c     on return
c
c        variables not requested by job are not used.
c
c        ap     contains the upper triangle of the inverse of
c               the original matrix, stored in packed form.
c               the columns of the upper triangle are stored
c               sequentially in a one-dimensional array.
c
c        det    double precision(2)
c               determinant of original matrix.
c               determinant = det(1) * 10.0**det(2)
c               with 1.0 .le. dabs(det(1)) .lt. 10.0
c               or det(1) = 0.0.
c
c        inert  integer(3)
c               the inertia of the original matrix.
c               inert(1)  =  number of positive eigenvalues.
c               inert(2)  =  number of negative eigenvalues.
c               inert(3)  =  number of zero eigenvalues.
c
c     error condition
c
c        a division by zero will occur if the inverse is requested
c        and  dspco  has set rcond .eq. 0.0
c        or  dspfa  has set  info .ne. 0 .
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dcopy,ddot,dswap
c     fortran dabs,iabs,mod
c
c     internal variables.
c
      double precision akkp1,ddot,temp
      double precision ten,d,t,ak,akp1
      integer ij,ik,ikp1,iks,j,jb,jk,jkp1
      integer k,kk,kkp1,km1,ks,ksj,kskp1,kstep
      logical noinv,nodet,noert
c
      noinv = mod(job,10) .eq. 0
      nodet = mod(job,100)/10 .eq. 0
      noert = mod(job,1000)/100 .eq. 0
c
      if (nodet .and. noert) go to 140
         if (noert) go to 10
            inert(1) = 0
            inert(2) = 0
            inert(3) = 0
   10    continue
         if (nodet) go to 20
            det(1) = 1.0d0
            det(2) = 0.0d0
            ten = 10.0d0
   20    continue
         t = 0.0d0
         ik = 0
         do 130 k = 1, n
            kk = ik + k
            d = ap(kk)
c
c           check if 1 by 1
c
            if (kpvt(k) .gt. 0) go to 50
c
c              2 by 2 block
c              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s)
c                      (s  c)
c              to avoid underflow/overflow troubles.
c              take two passes through scaling.  use  t  for flag.
c
               if (t .ne. 0.0d0) go to 30
                  ikp1 = ik + k
                  kkp1 = ikp1 + k
                  t = dabs(ap(kkp1))
                  d = (d/t)*ap(kkp1+1) - t
               go to 40
   30          continue
                  d = t
                  t = 0.0d0
   40          continue
   50       continue
c
            if (noert) go to 60
               if (d .gt. 0.0d0) inert(1) = inert(1) + 1
               if (d .lt. 0.0d0) inert(2) = inert(2) + 1
               if (d .eq. 0.0d0) inert(3) = inert(3) + 1
   60       continue
c
            if (nodet) go to 120
               det(1) = d*det(1)
               if (det(1) .eq. 0.0d0) go to 110
   70             if (dabs(det(1)) .ge. 1.0d0) go to 80
                     det(1) = ten*det(1)
                     det(2) = det(2) - 1.0d0
                  go to 70
   80             continue
   90             if (dabs(det(1)) .lt. ten) go to 100
                     det(1) = det(1)/ten
                     det(2) = det(2) + 1.0d0
                  go to 90
  100             continue
  110          continue
  120       continue
            ik = ik + k
  130    continue
  140 continue
c
c     compute inverse(a)
c
      if (noinv) go to 270
         k = 1
         ik = 0
  150    if (k .gt. n) go to 260
            km1 = k - 1
            kk = ik + k
            ikp1 = ik + k
            kkp1 = ikp1 + k
            if (kpvt(k) .lt. 0) go to 180
c
c              1 by 1
c
               ap(kk) = 1.0d0/ap(kk)
               if (km1 .lt. 1) go to 170
                  call dcopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 160 j = 1, km1
                     jk = ik + j
                     ap(jk) = ddot(j,ap(ij+1),1,work,1)
                     call daxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  160             continue
                  ap(kk) = ap(kk) + ddot(km1,work,1,ap(ik+1),1)
  170          continue
               kstep = 1
            go to 220
  180       continue
c
c              2 by 2
c
               t = dabs(ap(kkp1))
               ak = ap(kk)/t
               akp1 = ap(kkp1+1)/t
               akkp1 = ap(kkp1)/t
               d = t*(ak*akp1 - 1.0d0)
               ap(kk) = akp1/d
               ap(kkp1+1) = ak/d
               ap(kkp1) = -akkp1/d
               if (km1 .lt. 1) go to 210
                  call dcopy(km1,ap(ikp1+1),1,work,1)
                  ij = 0
                  do 190 j = 1, km1
                     jkp1 = ikp1 + j
                     ap(jkp1) = ddot(j,ap(ij+1),1,work,1)
                     call daxpy(j-1,work(j),ap(ij+1),1,ap(ikp1+1),1)
                     ij = ij + j
  190             continue
                  ap(kkp1+1) = ap(kkp1+1)
     *                         + ddot(km1,work,1,ap(ikp1+1),1)
                  ap(kkp1) = ap(kkp1)
     *                       + ddot(km1,ap(ik+1),1,ap(ikp1+1),1)
                  call dcopy(km1,ap(ik+1),1,work,1)
                  ij = 0
                  do 200 j = 1, km1
                     jk = ik + j
                     ap(jk) = ddot(j,ap(ij+1),1,work,1)
                     call daxpy(j-1,work(j),ap(ij+1),1,ap(ik+1),1)
                     ij = ij + j
  200             continue
                  ap(kk) = ap(kk) + ddot(km1,work,1,ap(ik+1),1)
  210          continue
               kstep = 2
  220       continue
c
c           swap
c
            ks = iabs(kpvt(k))
            if (ks .eq. k) go to 250
               iks = (ks*(ks - 1))/2
               call dswap(ks,ap(iks+1),1,ap(ik+1),1)
               ksj = ik + ks
               do 230 jb = ks, k
                  j = k + ks - jb
                  jk = ik + j
                  temp = ap(jk)
                  ap(jk) = ap(ksj)
                  ap(ksj) = temp
                  ksj = ksj - (j - 1)
  230          continue
               if (kstep .eq. 1) go to 240
                  kskp1 = ikp1 + ks
                  temp = ap(kskp1)
                  ap(kskp1) = ap(kkp1)
                  ap(kkp1) = temp
  240          continue
  250       continue
            ik = ik + k
            if (kstep .eq. 2) ik = ik + k + 1
            k = k + kstep
         go to 150
  260    continue
  270 continue
      return
      end
      subroutine dspfa(ap,n,kpvt,info)
      integer n,kpvt(1),info
      double precision ap(1)
c
c     dspfa factors a double precision symmetric matrix stored in
c     packed form by elimination with symmetric pivoting.
c
c     to solve  a*x = b , follow dspfa by dspsl.
c     to compute  inverse(a)*c , follow dspfa by dspsl.
c     to compute  determinant(a) , follow dspfa by dspdi.
c     to compute  inertia(a) , follow dspfa by dspdi.
c     to compute  inverse(a) , follow dspfa by dspdi.
c
c     on entry
c
c        ap      double precision (n*(n+1)/2)
c                the packed form of a symmetric matrix  a .  the
c                columns of the upper triangle are stored sequentially
c                in a one-dimensional array of length  n*(n+1)/2 .
c                see comments below for details.
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        ap      a block diagonal matrix and the multipliers which
c                were used to obtain it stored in packed form.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices , trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if the k-th pivot block is singular. this is
c                     not an error condition for this subroutine,
c                     but it does indicate that dspsl or dspdi may
c                     divide by zero if called.
c
c     packed storage
c
c          the following program segment will pack the upper
c          triangle of a symmetric matrix.
c
c                k = 0
c                do 20 j = 1, n
c                   do 10 i = 1, j
c                      k = k + 1
c                      ap(k)  = a(i,j)
c             10    continue
c             20 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,dswap,idamax
c     fortran dabs,dmax1,dsqrt
c
c     internal variables
c
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer idamax,ij,ijj,ik,ikm1,im,imax,imaxp1,imim,imj,imk
      integer j,jj,jk,jkm1,jmax,jmim,k,kk,km1,km1k,km1km1,km2,kstep
      logical swap
c
c
c     initialize
c
c     alpha is used in choosing pivot block size.
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
c
      info = 0
c
c     main loop on k, which goes from n to 1.
c
      k = n
      ik = (n*(n - 1))/2
   10 continue
c
c        leave the loop if k=0 or k=1.
c
c     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (ap(1) .eq. 0.0d0) info = 1
c     ......exit
            go to 200
   20    continue
c
c        this section of code determines the kind of
c        elimination to be performed.  when it is completed,
c        kstep will be set to the size of the pivot block, and
c        swap will be set to .true. if an interchange is
c        required.
c
         km1 = k - 1
         kk = ik + k
         absakk = dabs(ap(kk))
c
c        determine the largest off-diagonal element in
c        column k.
c
         imax = idamax(k-1,ap(ik+1),1)
         imk = ik + imax
         colmax = dabs(ap(imk))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
c
c           determine the largest off-diagonal element in
c           row imax.
c
            rowmax = 0.0d0
            imaxp1 = imax + 1
            im = imax*(imax - 1)/2
            imj = im + 2*imax
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(ap(imj)))
               imj = imj + j
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,ap(im+1),1)
               jmim = jmax + im
               rowmax = dmax1(rowmax,dabs(ap(jmim)))
   50       continue
            imim = imax + im
            if (dabs(ap(imim)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
c
c           column k is zero.  set info and iterate the loop.
c
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
c
c           1 x 1 pivot block.
c
            if (.not.swap) go to 120
c
c              perform an interchange.
c
               call dswap(imax,ap(im+1),1,ap(ik+1),1)
               imj = ik + imax
               do 110 jj = imax, k
                  j = k + imax - jj
                  jk = ik + j
                  t = ap(jk)
                  ap(jk) = ap(imj)
                  ap(imj) = t
                  imj = imj - (j - 1)
  110          continue
  120       continue
c
c           perform the elimination.
c
            ij = ik - (k - 1)
            do 130 jj = 1, km1
               j = k - jj
               jk = ik + j
               mulk = -ap(jk)/ap(kk)
               t = mulk
               call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
               ijj = ij + j
               ap(jk) = mulk
               ij = ij - (j - 1)
  130       continue
c
c           set the pivot array.
c
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
c
c           2 x 2 pivot block.
c
            km1k = ik + k - 1
            ikm1 = ik - (k - 1)
            if (.not.swap) go to 160
c
c              perform an interchange.
c
               call dswap(imax,ap(im+1),1,ap(ikm1+1),1)
               imj = ikm1 + imax
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  jkm1 = ikm1 + j
                  t = ap(jkm1)
                  ap(jkm1) = ap(imj)
                  ap(imj) = t
                  imj = imj - (j - 1)
  150          continue
               t = ap(km1k)
               ap(km1k) = ap(imk)
               ap(imk) = t
  160       continue
c
c           perform the elimination.
c
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = ap(kk)/ap(km1k)
               km1km1 = ikm1 + k - 1
               akm1 = ap(km1km1)/ap(km1k)
               denom = 1.0d0 - ak*akm1
               ij = ik - (k - 1) - (k - 2)
               do 170 jj = 1, km2
                  j = km1 - jj
                  jk = ik + j
                  bk = ap(jk)/ap(km1k)
                  jkm1 = ikm1 + j
                  bkm1 = ap(jkm1)/ap(km1k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,ap(ik+1),1,ap(ij+1),1)
                  t = mulkm1
                  call daxpy(j,t,ap(ikm1+1),1,ap(ij+1),1)
                  ap(jk) = mulk
                  ap(jkm1) = mulkm1
                  ijj = ij + j
                  ij = ij - (j - 1)
  170          continue
  180       continue
c
c           set the pivot array.
c
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         ik = ik - (k - 1)
         if (kstep .eq. 2) ik = ik - (k - 2)
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine dspsl(ap,n,kpvt,b)
      integer n,kpvt(1)
      double precision ap(1),b(1)
c
c     dsisl solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dspfa.
c
c     on entry
c
c        ap      double precision(n*(n+1)/2)
c                the output from dspfa.
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from dspfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  dspco  has set rcond .eq. 0.0
c        or  dspfa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dspfa(ap,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dspsl(ap,n,kpvt,c(1,j))
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran iabs
c
c     internal variables.
c
      double precision ak,akm1,bk,bkm1,ddot,denom,temp
      integer ik,ikm1,ikp1,k,kk,km1k,km1km1,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
      k = n
      ik = (n*(n - 1))/2
   10 if (k .eq. 0) go to 80
         kk = ik + k
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call daxpy(k-1,b(k),ap(ik+1),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/ap(kk)
            k = k - 1
            ik = ik - k
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            ikm1 = ik - (k - 1)
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call daxpy(k-2,b(k),ap(ik+1),1,b(1),1)
               call daxpy(k-2,b(k-1),ap(ikm1+1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            km1k = ik + k - 1
            kk = ik + k
            ak = ap(kk)/ap(km1k)
            km1km1 = ikm1 + k - 1
            akm1 = ap(km1km1)/ap(km1k)
            bk = b(k)/ap(km1k)
            bkm1 = b(k-1)/ap(km1k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
            ik = ik - (k + 1) - k
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
      ik = 0
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,ap(ik+1),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            ik = ik + k
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,ap(ik+1),1,b(1),1)
               ikp1 = ik + k
               b(k+1) = b(k+1) + ddot(k-1,ap(ikp1+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            ik = ik + k + k + 1
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end

* extens.f below.
C#define IEEE
C#define NEWVERS
      SUBROUTINE CROSS(X1,X2,CR12)
C- Cross product CR12 = X1 cross X2
C
C CROSS PRODUCT OF TWO VECTORS IN 3 DIMENSIONS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X1(3),X2(3),CR12(3)
      CR12(1) = X1(2)*X2(3) - X1(3)*X2(2)
      CR12(2) = X1(3)*X2(1) - X1(1)*X2(3)
      CR12(3) = X1(1)*X2(2) - X1(2)*X2(1)
      RETURN
      END
C --- First line of derfc ---
C#ifdefC DOUBLE16
C      double precision function derfc (x)
CC- complement of error function, real*16 precision
CC ----------------------------------------------------------------
CCi Inputs
CCi   x
CCo Outputs
CCo   complement of error function
CCr Remarks
CCr   series for erf on the interval  0 to 1
CCr           with weighted error   1.28e-32
CCr                log weighted error  31.89
CCr      significant figures required  31.05
CCr           decimal places required  32.55
CCr   Adapted from July 1977 edition, W. Fullerton, c3,
CCr   Los Alamos Scientific Lab.
CCr   For real*8 double precision see below
CC ----------------------------------------------------------------
C      double precision x, erfcs(21), erfccs(59), erc2cs(49), sqeps,
C     .  sqrtpi, xmax, xsml, y,  d1mach, dcsevl, dexp, dlog, dsqrt
C      external d1mach, dcsevl, initds
C      save
C      data erf cs(  1) / -.4904612123 4691808039 9845440333 76 d-1 /
C      data erf cs(  2) / -.1422612051 0371364237 8247418996 31 d+0 /
C      data erf cs(  3) / +.1003558218 7599795575 7546767129 33 d-1 /
C      data erf cs(  4) / -.5768764699 7674847650 8270255091 67 d-3 /
C      data erf cs(  5) / +.2741993125 2196061034 4221607914 71 d-4 /
C      data erf cs(  6) / -.1104317550 7344507604 1353812959 05 d-5 /
C      data erf cs(  7) / +.3848875542 0345036949 9613114981 74 d-7 /
C      data erf cs(  8) / -.1180858253 3875466969 6317518015 81 d-8 /
C      data erf cs(  9) / +.3233421582 6050909646 4029309533 54 d-10/
C      data erf cs( 10) / -.7991015947 0045487581 6073747085 95 d-12/
C      data erf cs( 11) / +.1799072511 3961455611 9672454866 34 d-13/
C      data erf cs( 12) / -.3718635487 8186926382 3168282094 93 d-15/
C      data erf cs( 13) / +.7103599003 7142529711 6899083946 66 d-17/
C      data erf cs( 14) / -.1261245511 9155225832 4954248533 33 d-18/
C      data erf cs( 15) / +.2091640694 1769294369 1705002666 66 d-20/
C      data erf cs( 16) / -.3253973102 9314072982 3641600000 00 d-22/
C      data erf cs( 17) / +.4766867209 7976748332 3733333333 33 d-24/
C      data erf cs( 18) / -.6598012078 2851343155 1999999999 99 d-26/
C      data erf cs( 19) / +.8655011469 9637626197 3333333333 33 d-28/
C      data erf cs( 20) / -.1078892517 7498064213 3333333333 33 d-29/
C      data erf cs( 21) / +.1281188399 3017002666 6666666666 66 d-31/
Cc
Cc series for erc2       on the interval  2.50000e-01 to  1.00000e+00
Cc                                        with weighted error   2.67e-32
Cc                                         log weighted error  31.57
Cc                               significant figures required  30.31
Cc                                    decimal places required  32.42
Cc
C      data erc2cs(  1) / -.6960134660 2309501127 3915082619 7 d-1 /
C      data erc2cs(  2) / -.4110133936 2620893489 8221208466 6 d-1 /
C      data erc2cs(  3) / +.3914495866 6896268815 6114370524 4 d-2 /
C      data erc2cs(  4) / -.4906395650 5489791612 8093545077 4 d-3 /
C      data erc2cs(  5) / +.7157479001 3770363807 6089414182 5 d-4 /
C      data erc2cs(  6) / -.1153071634 1312328338 0823284791 2 d-4 /
C      data erc2cs(  7) / +.1994670590 2019976350 5231486770 9 d-5 /
C      data erc2cs(  8) / -.3642666471 5992228739 3611843071 1 d-6 /
C      data erc2cs(  9) / +.6944372610 0050125899 3127721463 3 d-7 /
C      data erc2cs( 10) / -.1371220902 1043660195 3460514121 0 d-7 /
C      data erc2cs( 11) / +.2788389661 0071371319 6386034808 7 d-8 /
C      data erc2cs( 12) / -.5814164724 3311615518 6479105031 6 d-9 /
C      data erc2cs( 13) / +.1238920491 7527531811 8016881795 0 d-9 /
C      data erc2cs( 14) / -.2690639145 3067434323 9042493788 9 d-10/
C      data erc2cs( 15) / +.5942614350 8479109824 4470968384 0 d-11/
C      data erc2cs( 16) / -.1332386735 7581195792 8775442057 0 d-11/
C      data erc2cs( 17) / +.3028046806 1771320171 7369724330 4 d-12/
C      data erc2cs( 18) / -.6966648814 9410325887 9586758895 4 d-13/
C      data erc2cs( 19) / +.1620854541 0539229698 1289322762 8 d-13/
C      data erc2cs( 20) / -.3809934465 2504919998 7691305772 9 d-14/
C      data erc2cs( 21) / +.9040487815 9788311493 6897101297 5 d-15/
C      data erc2cs( 22) / -.2164006195 0896073478 0981204700 3 d-15/
C      data erc2cs( 23) / +.5222102233 9958549846 0798024417 2 d-16/
C      data erc2cs( 24) / -.1269729602 3645553363 7241552778 0 d-16/
C      data erc2cs( 25) / +.3109145504 2761975838 3622741295 1 d-17/
C      data erc2cs( 26) / -.7663762920 3203855240 0956671481 1 d-18/
C      data erc2cs( 27) / +.1900819251 3627452025 3692973329 0 d-18/
C      data erc2cs( 28) / -.4742207279 0690395452 2565599996 5 d-19/
C      data erc2cs( 29) / +.1189649200 0765283828 8068307845 1 d-19/
C      data erc2cs( 30) / -.3000035590 3257802568 4527131306 6 d-20/
C      data erc2cs( 31) / +.7602993453 0432461730 1938527709 8 d-21/
C      data erc2cs( 32) / -.1935909447 6068728815 6981104913 0 d-21/
C      data erc2cs( 33) / +.4951399124 7733378810 0004238677 3 d-22/
C      data erc2cs( 34) / -.1271807481 3363718796 0862198988 8 d-22/
C      data erc2cs( 35) / +.3280049600 4695130433 1584165205 3 d-23/
C      data erc2cs( 36) / -.8492320176 8228965689 2479242239 9 d-24/
C      data erc2cs( 37) / +.2206917892 8075602235 1987998719 9 d-24/
C      data erc2cs( 38) / -.5755617245 6965284983 1281950719 9 d-25/
C      data erc2cs( 39) / +.1506191533 6392342503 5414405119 9 d-25/
C      data erc2cs( 40) / -.3954502959 0187969531 0428569599 9 d-26/
C      data erc2cs( 41) / +.1041529704 1515009799 8464505173 3 d-26/
C      data erc2cs( 42) / -.2751487795 2787650794 5017890133 3 d-27/
C      data erc2cs( 43) / +.7290058205 4975574089 9770368000 0 d-28/
C      data erc2cs( 44) / -.1936939645 9159478040 7750109866 6 d-28/
C      data erc2cs( 45) / +.5160357112 0514872983 7005482666 6 d-29/
C      data erc2cs( 46) / -.1378419322 1930940993 8964480000 0 d-29/
C      data erc2cs( 47) / +.3691326793 1070690422 5109333333 3 d-30/
C      data erc2cs( 48) / -.9909389590 6243654206 5322666666 6 d-31/
C      data erc2cs( 49) / +.2666491705 1953884133 2394666666 6 d-31/
Cc
Cc series for erfc       on the interval  0.          to  2.50000e-01
Cc                                        with weighted error   1.53e-31
Cc                                         log weighted error  30.82
Cc                               significant figures required  29.47
Cc                                    decimal places required  31.70
Cc
C      data erfccs(  1) / +.7151793102 0292477450 3697709496 d-1 /
C      data erfccs(  2) / -.2653243433 7606715755 8893386681 d-1 /
C      data erfccs(  3) / +.1711153977 9208558833 2699194606 d-2 /
C      data erfccs(  4) / -.1637516634 5851788416 3746404749 d-3 /
C      data erfccs(  5) / +.1987129350 0552036499 5974806758 d-4 /
C      data erfccs(  6) / -.2843712412 7665550875 0175183152 d-5 /
C      data erfccs(  7) / +.4606161308 9631303696 9379968464 d-6 /
C      data erfccs(  8) / -.8227753025 8792084205 7766536366 d-7 /
C      data erfccs(  9) / +.1592141872 7709011298 9358340826 d-7 /
C      data erfccs( 10) / -.3295071362 2528432148 6631665072 d-8 /
C      data erfccs( 11) / +.7223439760 4005554658 1261153890 d-9 /
C      data erfccs( 12) / -.1664855813 3987295934 4695966886 d-9 /
C      data erfccs( 13) / +.4010392588 2376648207 7671768814 d-10/
C      data erfccs( 14) / -.1004816214 4257311327 2170176283 d-10/
C      data erfccs( 15) / +.2608275913 3003338085 9341009439 d-11/
C      data erfccs( 16) / -.6991110560 4040248655 7697812476 d-12/
C      data erfccs( 17) / +.1929492333 2617070862 4205749803 d-12/
C      data erfccs( 18) / -.5470131188 7543310649 0125085271 d-13/
C      data erfccs( 19) / +.1589663309 7626974483 9084032762 d-13/
C      data erfccs( 20) / -.4726893980 1975548392 0369584290 d-14/
C      data erfccs( 21) / +.1435873376 7849847867 2873997840 d-14/
C      data erfccs( 22) / -.4449510561 8173583941 7250062829 d-15/
C      data erfccs( 23) / +.1404810884 7682334373 7305537466 d-15/
C      data erfccs( 24) / -.4513818387 7642108962 5963281623 d-16/
C      data erfccs( 25) / +.1474521541 0451330778 7018713262 d-16/
C      data erfccs( 26) / -.4892621406 9457761543 6841552532 d-17/
C      data erfccs( 27) / +.1647612141 4106467389 5301522827 d-17/
C      data erfccs( 28) / -.5626817176 3294080929 9928521323 d-18/
C      data erfccs( 29) / +.1947443382 2320785142 9197867821 d-18/
C      data erfccs( 30) / -.6826305642 9484207295 6664144723 d-19/
C      data erfccs( 31) / +.2421988887 2986492401 8301125438 d-19/
C      data erfccs( 32) / -.8693414133 5030704256 3800861857 d-20/
C      data erfccs( 33) / +.3155180346 2280855712 2363401262 d-20/
C      data erfccs( 34) / -.1157372324 0496087426 1239486742 d-20/
C      data erfccs( 35) / +.4288947161 6056539462 3737097442 d-21/
C      data erfccs( 36) / -.1605030742 0576168500 5737770964 d-21/
C      data erfccs( 37) / +.6063298757 4538026449 5069923027 d-22/
C      data erfccs( 38) / -.2311404251 6979584909 8840801367 d-22/
C      data erfccs( 39) / +.8888778540 6618855255 4702955697 d-23/
C      data erfccs( 40) / -.3447260576 6513765223 0718495566 d-23/
C      data erfccs( 41) / +.1347865460 2069650682 7582774181 d-23/
C      data erfccs( 42) / -.5311794071 1250217364 5873201807 d-24/
C      data erfccs( 43) / +.2109341058 6197831682 8954734537 d-24/
C      data erfccs( 44) / -.8438365587 9237891159 8133256738 d-25/
C      data erfccs( 45) / +.3399982524 9452089062 7359576337 d-25/
C      data erfccs( 46) / -.1379452388 0732420900 2238377110 d-25/
C      data erfccs( 47) / +.5634490311 8332526151 3392634811 d-26/
C      data erfccs( 48) / -.2316490434 4770654482 3427752700 d-26/
C      data erfccs( 49) / +.9584462844 6018101526 3158381226 d-27/
C      data erfccs( 50) / -.3990722880 3301097262 4224850193 d-27/
C      data erfccs( 51) / +.1672129225 9444773601 7228709669 d-27/
C      data erfccs( 52) / -.7045991522 7660138563 8803782587 d-28/
C      data erfccs( 53) / +.2979768402 8642063541 2357989444 d-28/
C      data erfccs( 54) / -.1262522466 4606192972 2422632994 d-28/
C      data erfccs( 55) / +.5395438704 5424879398 5299653154 d-29/
C      data erfccs( 56) / -.2380992882 5314591867 5346190062 d-29/
C      data erfccs( 57) / +.1099052830 1027615735 9726683750 d-29/
C      data erfccs( 58) / -.4867713741 6449657273 2518677435 d-30/
C      data erfccs( 59) / +.1525877264 1103575676 3200828211 d-30/
Cc
C      data sqrtpi / 1.772453850 9055160272 9816748334 115d0 /
C      data nterf, nterfc, nterc2, xsml, xmax, sqeps / 3*0, 3*0.d0 /
Cc
C      if (nterf .ne. 0) goto 10
C      eta = 0.1*sngl(d1mach(3))
C      nterf = initds (erfcs, 21, eta)
C      nterfc = initds (erfccs, 59, eta)
C      nterc2 = initds (erc2cs, 49, eta)
Cc
C      xsml = -dsqrt (-dlog(sqrtpi*d1mach(3)))
C      xmax = dsqrt (-dlog(sqrtpi*d1mach(1)) )
C      xmax = xmax - 0.5d0*dlog(xmax)/xmax - 0.01d0
C      sqeps = dsqrt (2.0d0*d1mach(3))
Cc
C   10 if (x .gt. xsml) goto 20
Cc
Cc === erfc(x) = 1.0 - erf(x)  for  x .lt. xsml ===
C      derfc = 2.0d0
C      return
Cc
C   20 if (x .gt. xmax) goto 40
C      y = dabs(x)
C      if (y .gt. 1.0d0) goto 30
Cc
Cc === erfc(x) = 1.0 - erf(x)  for abs(x) .le. 1.0 ===
C      if (y .lt. sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
C      if (y .ge. sqeps)
C     . derfc = 1.0d0 - x*(1.0d0 + dcsevl (2.d0*x*x-1.d0,erfcs, nterf))
C      return
Cc
Cc === erfc(x) = 1.0 - erf(x)  for  1.0 .lt. abs(x) .le. xmax ===
C   30 y = y*y
C      if (y .le. 4.d0) derfc = dexp(-y)/dabs(x) *
C     .  (0.5d0 + dcsevl ((8.d0/y-5.d0)/3.d0, erc2cs, nterc2))
C      if (y .gt. 4.d0) derfc = dexp(-y)/dabs(x) *
C     .   (0.5d0 + dcsevl (8.d0/y-1.d0, erfccs, nterfc))
C      if (x .lt. 0.d0) derfc = 2.0d0 - derfc
C      return
Cc
C   40 call errmsg ('DERFC: underflow', 1)
C      derfc = 0.d0
C      return
Cc
C      end
C#else
      double precision function derfc (x)
C- complement of error function, real*8 precision
C ----------------------------------------------------------------
Ci Inputs
Ci   x
Co Outputs
Co   complement of error function
Cr Remarks
Cr   erfcs: series for erf on the interval  0 to 1
Cr              with weighted error  7.10d-18
Cr                  log weighted error  17.15
Cr        significant figures required  16.31
Cr             decimal places required  17.71
Cr   erc2s: series for erc2 on the interval  .25 to 1
Cr             with weighted error   5.22d-17
Cr                  log weighted error  16.28
Cr         significant figures required  15.0
Cr             decimal places required  16.96
Cr   erfccs: series for erfc on the interval  0 to  .25
Cr              with weighted error  4.81d-17
Cr                  log weighted error  16.32
Cr         significant figures required  15.0
Cr             decimal places required  17.01
Cr   Adapted from July 1977 edition, W. Fullerton, c3,
Cr   Los Alamos Scientific Lab.
Cr   For real*16 double precision see below
C ----------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension erfcs(13), erfccs(24), erc2cs(23)
      real eta
      external dcsevl, initds, d1mach
      integer iprint
      save

      data erf cs( 1) /   -.0490461212 34691808d0 /
      data erf cs( 2) /   -.1422612051 0371364d0 /
      data erf cs( 3) /    .0100355821 87599796d0 /
      data erf cs( 4) /   -.0005768764 69976748d0 /
      data erf cs( 5) /    .0000274199 31252196d0 /
      data erf cs( 6) /   -.0000011043 17550734d0 /
      data erf cs( 7) /    .0000000384 88755420d0 /
      data erf cs( 8) /   -.0000000011 80858253d0 /
      data erf cs( 9) /    .0000000000 32334215d0 /
      data erf cs(10) /   -.0000000000 00799101d0 /
      data erf cs(11) /    .0000000000 00017990d0 /
      data erf cs(12) /   -.0000000000 00000371d0 /
      data erf cs(13) /    .0000000000 00000007d0 /
c
      data erc2cs( 1) /   -.0696013466 02309501d0 /
      data erc2cs( 2) /   -.0411013393 62620893d0 /
      data erc2cs( 3) /    .0039144958 66689626d0 /
      data erc2cs( 4) /   -.0004906395 65054897d0 /
      data erc2cs( 5) /    .0000715747 90013770d0 /
      data erc2cs( 6) /   -.0000115307 16341312d0 /
      data erc2cs( 7) /    .0000019946 70590201d0 /
      data erc2cs( 8) /   -.0000003642 66647159d0 /
      data erc2cs( 9) /    .0000000694 43726100d0 /
      data erc2cs(10) /   -.0000000137 12209021d0 /
      data erc2cs(11) /    .0000000027 88389661d0 /
      data erc2cs(12) /   -.0000000005 81416472d0 /
      data erc2cs(13) /    .0000000001 23892049d0 /
      data erc2cs(14) /   -.0000000000 26906391d0 /
      data erc2cs(15) /    .0000000000 05942614d0 /
      data erc2cs(16) /   -.0000000000 01332386d0 /
      data erc2cs(17) /    .0000000000 00302804d0 /
      data erc2cs(18) /   -.0000000000 00069666d0 /
      data erc2cs(19) /    .0000000000 00016208d0 /
      data erc2cs(20) /   -.0000000000 00003809d0 /
      data erc2cs(21) /    .0000000000 00000904d0 /
      data erc2cs(22) /   -.0000000000 00000216d0 /
      data erc2cs(23) /    .0000000000 00000052d0 /
c
      data erfccs( 1) /   0.0715179310 202925d0 /
      data erfccs( 2) /   -.0265324343 37606719d0 /
      data erfccs( 3) /    .0017111539 77920853d0 /
      data erfccs( 4) /   -.0001637516 63458512d0 /
      data erfccs( 5) /    .0000198712 93500549d0 /
      data erfccs( 6) /   -.0000028437 12412769d0 /
      data erfccs( 7) /    .0000004606 16130901d0 /
      data erfccs( 8) /   -.0000000822 77530261d0 /
      data erfccs( 9) /    .0000000159 21418724d0 /
      data erfccs(10) /   -.0000000032 95071356d0 /
      data erfccs(11) /    .0000000007 22343973d0 /
      data erfccs(12) /   -.0000000001 66485584d0 /
      data erfccs(13) /    .0000000000 40103931d0 /
      data erfccs(14) /   -.0000000000 10048164d0 /
      data erfccs(15) /    .0000000000 02608272d0 /
      data erfccs(16) /   -.0000000000 00699105d0 /
      data erfccs(17) /    .0000000000 00192946d0 /
      data erfccs(18) /   -.0000000000 00054704d0 /
      data erfccs(19) /    .0000000000 00015901d0 /
      data erfccs(20) /   -.0000000000 00004729d0 /
      data erfccs(21) /    .0000000000 00001432d0 /
      data erfccs(22) /   -.0000000000 00000439d0 /
      data erfccs(23) /    .0000000000 00000138d0 /
      data erfccs(24) /   -.0000000000 00000048d0 /
c
      data sqrtpi /1.772453850 9055160d0/
      data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0.d0/
c
      if (nterf .ne. 0d0) goto 10
      eta = 0.1*sngl(d1mach(3))
      nterf = initds (erfcs, 13, eta)
      nterfc = initds (erfccs, 24, eta)
      nterc2 = initds (erc2cs, 23, eta)

      xsml = -dsqrt (-dlog(sqrtpi*d1mach(3)))
      xmax = dsqrt (-dlog(sqrtpi*d1mach(1)))
      xmax = xmax - 0.5d0*dlog(xmax)/xmax - 0.01d0
      sqeps = dsqrt (2.0d0*d1mach(3))
c
   10 if (x .gt. xsml) goto 20
c
C --- derfc(x) = 1.0d0 - erf(x) for x .lt. xsml ---
      derfc = 2.d0
      return
c
   20 if (x .gt. xmax) goto 40
      y = dabs(x)
      if (y .gt. 1.0d0) goto 30
c
C --- derfc(x) = 1.0d0 - erf(x) for -1.d0 .le. x .le. 1.d0 ---
      if (y .lt. sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
      if (y .ge. sqeps) derfc = 1.0d0 -
     .  x*(1.0d0 + dcsevl (2.d0*x*x-1.d0, erfcs, nterf) )
      return
c
C --- derfc(x) = 1.0d0 - erf(x) for 1.d0 .lt. dabs(x) .le. xmax ---
   30 y = y*y
      if (y .le. 4.d0) derfc = dexp(-y)/dabs(x) *
     .  (0.5d0 + dcsevl ((8.d0/y-5.d0)/3.d0, erc2cs, nterc2) )
      if (y .gt. 4.d0) derfc = dexp(-y)/dabs(x) *
     .  (0.5d0 + dcsevl (8.d0/y-1.d0, erfccs, nterfc) )
      if (x .lt. 0.d0) derfc = 2.0d0 - derfc
      return
c
   40 if (iprint() .gt. 50) call errmsg ('DERFC: underflow', 1)
      derfc = 0.d0
      return
c
      end
C#endif
      double precision function derf (x)
      double precision derfc,x
      derf = 1 - derfc(x)
      end

      subroutine ishell(n,iarray)
      integer n
      integer iarray(1)
      integer lognb2,i,j,k,l,m,nn,it

      lognb2 = int(log(float(n+1))*1.4426950)
      m = n
      do  12  nn = 1, lognb2
        m = m/2
        k = n - m
        do  11  j = 1, k
          i = j
    3     continue
          l = i + m
          if (iarray(l) .lt. iarray(i)) then
            it = iarray(i)
            iarray(i) = iarray(l)
            iarray(l) = it
            i = i - m
            if (i .ge. 1) goto 3
          endif
   11   continue
   12 continue
      return
      end
      subroutine dinv33(matrix,iopt,inverse,det)
C- Inverts 3X3 matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   inverse: input matrix
Ci   iopt:  if 0, usual inverse
Ci             1, transpose of inverse
Ci             2, 2*pi*inverse
Ci             3, 2*pi*transpose of inverse
Co Outputs
Co   inverse, as modified according to iopt
Co   det:      determinant, or det/2*pi (sign ok ??)
Cr Remarks
Cr   To generate reciprocal lattice vectors, call dinv33(plat,3,rlat)
C ----------------------------------------------------------------
C     implicit none
      integer iopt,i,j
      double precision matrix(3,3),inverse(3,3),det,ddot
      double precision xx
      call cross(matrix(1,2),matrix(1,3),inverse     )
      call cross(matrix(1,3),matrix     ,inverse(1,2))
      call cross(matrix     ,matrix(1,2),inverse(1,3))
      det = ddot(3,matrix,1,inverse,1)
      if (det .eq. 0.d0) stop 'INV33: vanishing determinant'
      if (iopt .ge. 2) det = det/(8*datan(1.d0))
      if (mod(iopt,2) .eq. 0) then
        do  10  i = 1, 3
          do  10  j = i+1, 3
          xx = inverse(i,j)
          inverse(i,j) = inverse(j,i)
          inverse(j,i) = xx
   10   continue
      endif
      call dscal(9,1/det,inverse,1)
      return
      end

      SUBROUTINE DMCPY(A,NCA,NRA,B,NCB,NRB,N,M)
C- general matrix copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of columns and rows, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix copy routine,
Cr   copying a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b} is the row dimension of matrix {a,b}
Cr   and nr{a,b} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Example: Set 3-by-2 block of matrix c to constant z
Cr     call dmcpy(z,0,0,c,nc,1,3,2)
Cr   Note scalar z is represented by an array of 0 dimension
C ----------------------------------------------------------------
C
      INTEGER NCA,NRA,NCB,NRB,N,M
      DOUBLE PRECISION A(0:*), B(0:*)
      INTEGER I,J,IA,IB

      DO  200  I = N-1, 0, -1
        IA = I*NRA+M*NCA
        IB = I*NRB+M*NCB
        DO  200  J = M-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        B(IB) = A(IA)
  200 CONTINUE
      RETURN
      END
      logical function bittst(n,bit)
C- returns true when a bit is set in an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n: integer
Ci   bit: a bit, ie 1,2,4,8,16, etc
Co Outputs
Co   bittst: true when bit in n is set, false otherwise
C ----------------------------------------------------------------
      integer n,bit
      bittst = (mod(n,bit+bit) - mod(n,bit) .eq. bit)
      end
      integer function getdig(n,i,base)
C- Returns the a digit from an integer
C ----------------------------------------------------------------
Ci Inputs
Ci   n,i,base
Co Outputs
Co   getdig = ith digit from n, base "base"; eg 4=getdig(12345,1,10)
C ----------------------------------------------------------------
      integer n,i,base
      getdig = mod(n/base**i,base)
      end
C#ifdef NEWVERS
      subroutine zhmul(ndim,hk,d,nid,ck)
C- Multiplies matrix ck = hk*d*hk, where d is a diagonal matrix.
C- Only the upper triangle of ck is calculated.
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim,hk,d
Ci   nid: skip length between elements of d
Co Outputs
Co   ck
Cr Remarks
Cr   Hk and ck are dimensioned (ndim*ndim*2) (real followed by imag)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,nid
      double precision hk(ndim*ndim*2),ck(ndim*ndim*2),d(nid*ndim)
C Local parameters
      integer n2,i,j,k,jndim,j2ndim,kndim,k2ndim
      double precision xikr,xiki,xikim

      n2 = ndim**2
      do  100  j = 1, ndim
        jndim = (j-1)*ndim
        j2ndim = jndim+n2
        do  100  k = 1, ndim
          kndim = (k-1)*ndim
          k2ndim = kndim+n2
          xikr = d(k)*hk(jndim+k)
          xiki = d(k)*hk(j2ndim+k)
          xikim = -xiki
         do  110  i = 1, j
           ck(jndim+i) = ck(jndim+i) + xikr * hk(kndim+i) +
     .                                 xikim * hk(k2ndim+i)
           ck(j2ndim+i) = ck(j2ndim+i) + xiki * hk(kndim+i) +
     .                                   xikr * hk(k2ndim+i)
  110    continue
  100 continue
      end
C#elseC
C      subroutine zhmul(ndim,hk,d,nid,ck)
CC- Multiplies matrix ck = hk*d*hk, where d is a diagonal matrix
CC ----------------------------------------------------------------
CCi Inputs
CCi   ndim,hk,d
CCi   nid: skip length between elements of d
CCo Outputs
CCo   ck
CCr Remarks
CCr   Hk and ck are dimensioned (ndim,ndim,2) (real followed by imag)
CCr   The looping order is chosen so as to vectorize the inner loop
CC ----------------------------------------------------------------
CC     implicit none
CC Passed parameters
C      integer ndim,nid
C      double precision hk(ndim,ndim),ck(ndim,ndim),d(nid,1)
CC Local parameters
C      integer n2,i,j,k
C      double precision xikr,xiki
CC
C      n2 = ndim**2
C      do  100  k = 1, ndim
C        do  100  i = 1, ndim
C        xikr = d(1,k)*hk(i,k)
C        xiki = d(1,k)*hk(n2+i,k)
C#ifdefC BLAS
C        call daxpy(ndim, xikr,hk(   k,1),ndim,ck(i,1),   ndim)
C        call daxpy(ndim,-xiki,hk(n2+k,1),ndim,ck(i,1),   ndim)
C        call daxpy(ndim, xikr,hk(n2+k,1),ndim,ck(n2+i,1),ndim)
C        call daxpy(ndim, xiki,hk(   k,1),ndim,ck(n2+i,1),ndim)
C#elseC
C        do  110  j = 1, ndim
C          ck(i,j)    = ck(i,j)    + xikr*hk(k,j) - xiki*hk(n2+k,j)
C          ck(n2+i,j) = ck(n2+i,j) + xikr*hk(n2+k,j) + xiki*hk(k,j)
C  110 continue
C#endifC
C  100 continue
C      end
C#endif
      subroutine zmul0(ndim,a,b,c)
C- Multiplies complex matrix c = a*b,
C- All matrices have the same dimension
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim,a,b
Co Outputs
Co   c
Cr Remarks
Cr   a,b, and c are dimensioned (ndim*ndim*2) (real followed by imag)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim
      double precision a(ndim*ndim*2),b(ndim*ndim*2),c(ndim*ndim*2)
C Local parameters
      integer n2,i,j,k,jndim,j2ndim,kndim,k2ndim
      double precision bkj,bk2j,bk2jm

      n2 = ndim**2

C#ifdefC CRAY
C      call mxma(a,1,ndim,b,1,ndim,c,1,ndim,ndim,ndim,ndim)
C      call mxma(a,1,ndim,b(n2+1),1,ndim,c(n2+1),1,ndim,ndim,ndim,ndim)
C      call dinit(a,n2)
C      call mxma(a(n2+1),1,ndim,b(n2+1),1,ndim,a,1,ndim,ndim,ndim,ndim)
C      do 10 i = 1, n2
C        c(i) = c(i) - a(i)
C   10 continue
C      call dinit(a,n2)
C      call mxma(a(n2+1),1,ndim,b,1,ndim,a,1,ndim,ndim,ndim,ndim)
C      do 20 i = 1, n2
C        c(n2+i) = c(n2+i) + a(i)
C   20 continue
C#else
      do  100  j = 1, ndim
        jndim = (j-1)*ndim
        j2ndim = jndim+n2
        do  100  k = 1, ndim
          kndim = (k-1)*ndim
          k2ndim = kndim+n2
          bkj = b(jndim+k)
          bk2j = b(j2ndim+k)
          bk2jm = -bk2j
          do  110  i = 1, ndim
            c(jndim+i) = c(jndim+i) + bkj*a(kndim+i)
     .                                + bk2jm*a(k2ndim+i)
            c(j2ndim+i) = c(j2ndim+i) + bk2j*a(kndim+i)
     .                                  + bkj*a(k2ndim+i)
  110     continue
  100 continue
C#endif
      end
      subroutine zmul(a,nleft,b,nmid,c,nright)
C- General complex matrix copy c = a*b, cij = \sum aik*bkj
C ----------------------------------------------------------------
Ci Inputs
Ci   nleft,nmid,nright,a,b
Co Outputs
Co   c
Cr Remarks
Cr   c(nleft,nright) = a(nleft,nmid) * b(nmid,nright)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nleft,nmid,nright
      double precision a(nleft*nmid*2),b(nmid*nright*2),
     .                 c(nleft*nright*2)
C Local parameters
      integer nlm,nmr,nlr,i,j,k,jb,j2b,jc,j2c,ka,k2a
      double precision bkj,bk2j,bk2jm

      nlm = nleft * nmid
      nmr = nmid * nright
      nlr = nleft * nright

      do  100  j = 1, nright
        jc = (j-1) * nleft
        j2c = jc + nlr
        jb = (j-1) * nmid
        j2b = jb + nmr
        do  100  k = 1, nmid
          ka = (k-1) * nleft
          k2a = ka + nlm
          bkj = b(jb+k)
          bk2j = b(j2b+k)
          bk2jm = -bk2j
          do  110  i = 1, nleft
            c(jc+i)  = c(jc+i)  + bkj*a(ka+i)  + bk2jm*a(k2a+i)
            c(j2c+i) = c(j2c+i) + bk2j*a(ka+i) + bkj*a(k2a+i)
  110     continue
  100 continue
      end
      subroutine zmult(a,nleft,b,nmid,c,nright)
C- General complex matrix copy c = at*b,
C-                             cij = \sum atik*bkj = \sum aki*bkj
C- at is the transposed of a
C ----------------------------------------------------------------
Ci Inputs
Ci   nleft,nmid,nright,a,b
Co Outputs
Co   c
Cr Remarks
Cr   c(nleft,nright) = at(nleft,nmid) * b(nmid,nright)
Cr   c(nleft,nright) = a (nmid,nleft) * b(nmid,nright)
Cr   The looping order is chosen so as to vectorize the inner loop
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nleft,nmid,nright
      double precision a(nmid*nleft*2),b(nmid*nright*2),
     .                 c(nleft*nright*2)
C Local parameters
      integer nml,nmr,nlr,i,j,k,ia,i2a,jb,j2b,jc,j2c
      double precision cic,ci2c

      nml = nmid * nleft
      nmr = nmid * nright
      nlr = nleft * nright

      do  100  j = 1, nright
        jc = (j-1) * nleft
        j2c = jc + nlr
        jb = (j-1) * nmid
        j2b = jb + nmr
        do  100  i = 1, nleft
          ia = (i-1) * nmid
          i2a = ia + nml
          cic = 0
          ci2c = 0
          do  110  k = 1, nmid
            cic  = cic  + a(ia+k)*b(jb+k) - a(i2a+k)*b(j2b+k)
            ci2c = ci2c + a(i2a+k)*b(jb+k) + a(ia+k)*b(j2b+k)
  110     continue
        c(jc+i) = c(jc+i) + cic
        c(j2c+i) = c(j2c+i) + ci2c
  100 continue
      end
      SUBROUTINE DVCPY(A,NCA,B,NCB,N)
C- general vector copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca is the source vector and the number of
Ci      elements separating each element in the vector
Ci   b,ncb,nrb is the destination vector and the number of
Ci      elements separating each element in the vector
Ci   n: the number of elements to calculate
Co Outputs:
Co   result matrix stored in b
Cr Remarks:
Cr   This is a general-purpose vectore copy routine
Cr   Example: Set all elements of 3-by-2 matrix c to -1.d0
Cr     call dvcpy(-1.d0,0,c,1,3*2)
Cr   Example: Set block (n,m) of array a(p,m)=0, and a(i,i)=1, i=1,m
Cr     call dmcpy(0.d0,0,0,a,p,1,n,m)
Cr     call dvcpy(1.d0,0,a,p+1,m)
C ----------------------------------------------------------------
      INTEGER NCA,NCB,N
      DOUBLE PRECISION A(0:*), B(0:*)
      INTEGER I,IA,IB

      IA = N*NCA
      IB = N*NCB
      DO  200  I = N-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        B(IB) = A(IA)
  200 CONTINUE

      return
      end
      subroutine ivshel(m,n,iarray,iwk,lopt)
C- shell sort of a array of integer vectors
C ----------------------------------------------------------------
Ci Inputs
Ci   iarray(m,n)
Ci   iwk: a work array of dimension n
Ci   lopt:if T, iwk returned as list of indices to iarray to sort it,
Ci              while iarray is unchanged.
Ci        if F, iarray returned sorted
Co Outputs
Co   iwk a table of indices to array iarray (lopt true)
C ----------------------------------------------------------------
      integer m,n
      logical lopt
      integer iarray(m,0:n-1),iwk(n)
      integer lognb2,i,j,k,l,n2,nn,it
      lognb2 = int(log(float(n+1))*1.4426950)
      n2 = n
      do  2  i = 1, n
    2 iwk(i) = i-1
      do  12  nn = 1, lognb2
        n2 = n2/2
        k = n - n2
        do  11  j = 1, k
          i = j
    3     continue
          l = i + n2
c      print *, 'test ',i,l,iwk(i),iwk(l)
          do  15  mm = 1, m
            if (iarray(mm,iwk(l)) - iarray(mm,iwk(i))) 16,15,11
   16       continue
            if (lopt) then
              it = iwk(i)
              iwk(i) = iwk(l)
              iwk(l) = it
c      print 800, (iwk(nnn), nnn=1,n)
c800   format(' swap', 11i5)
            else
c      print 800,  (iarray(nnn,i-1), nnn=1,m)
c      print 800,  (iarray(nnn,l-1), nnn=1,m)
              do  14  mmm = 1, m
                it = iarray(mmm,i-1)
                iarray(mmm,i-1) = iarray(mmm,l-1)
                iarray(mmm,l-1) = it
   14         continue
            endif
            i = i - n2
            if (i .ge. 1) goto 3
            goto 11
   15     continue
   11   continue
   12 continue

      return
      end



      SUBROUTINE DMPY(A,NCA,NRA,B,NCB,NRB,C,NCC,NRC,N,M,L)
C- matrix multiplication
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   c,ncc,nrc is the product matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of (n,l) subblock of a into (l,m) subblock of b
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
C ----------------------------------------------------------------
      INTEGER NCA,NRA,NCB,NRB,NCC,NRC,N,M,L
      DOUBLE PRECISION A(0:*), B(0:*), C(0:*)
      DOUBLE PRECISION SUM
      INTEGER I,J,K,NAKPI,NBJPK
C
C#ifdefC CRAY
C      CALL MXMA(A,NRA,NCA,B,NRB,NCB,C,NRC,NCC,N,L,M)
C#else
      DO  200  I = N-1, 0, -1
        DO  200  J = M-1, 0, -1
        SUM = 0
        NAKPI = NRA*I
        NBJPK = NCB*J
        DO  210  K = L-1, 0, -1
          SUM = SUM + A(NAKPI)*B(NBJPK)
          NAKPI = NAKPI + NCA
          NBJPK = NBJPK + NRB
  210   CONTINUE
        C(I*NRC+J*NCC) = SUM
  200 CONTINUE
C#endif
      END
      double precision function dcsevl (x, a, n)
C- Evaluate the n-term Chebyshev series a at x.
C ----------------------------------------------------------------
Ci Inputs
Ci   x:  dble prec value at which the series is to be evaluated.
Ci   a:  dble prec array of n terms of a chebyshev series.
Ci       In evaluating a, only half the first coef is summed.
Ci   n:  number of terms in array a.
Co Outputs
Co   dcsevl
Cr Remarks
Cr   Adapted from R. Broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
C ----------------------------------------------------------------
      double precision a(n), x, twox, b0, b1, b2
c
      if (dabs(x) .gt. 1.1d0) call errmsg('DCSEVL:  x outside (-1,1)',2)
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do  10  i = 1, n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
   10 continue
      dcsevl = .5d0*(b0-b2)

      return
      end
      subroutine errmsg (messg, iopt)
C- Write error message to standard error device
C ----------------------------------------------------------------
Ci Inputs
Ci   iopt: 0, return without message printed
Ci         1, return with message printed
Ci         2, stop with message printed
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
      character*(*) messg

      if (iopt .ne. 0) write(i1mach(4),*) messg
      if (iopt .lt. 2) return
      end
      function initds (dos, nos, eta)
C- Initialize things for Chebychev series
C ----------------------------------------------------------------
Ci Inputs
Ci   dos: dble prec array of nos coefficients in an orthogonal series.
Ci   nos: number of coefficients in dos.
Ci   eta: requested accuracy of series (real).
Co Outputs
Co   Returns number of terms necessary in series for spec'd eta
Cr Remarks
Cr   Initialize the double precision orthogonal series dos so that
Cr   initds is the number of terms needed to insure the error is no
Cr   larger than eta.  ordinarily eta will be chosen to be one-tenth
Cr   machine precision.
Cr   Adapted from June 1977 edition W. Fullerton,
Cr   c3, los alamos scientific lab.
C ----------------------------------------------------------------
      double precision dos(nos)
      real eta
c
      err = 0.
      do  10  ii = 1, nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err .gt. eta) goto 20
   10 continue
c
   20 continue
c     if (i .eq. nos) call errmsg('INITDS: eta may be too small',1)
      initds = i
c
      return
      end
      subroutine zampy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- complex matrix multiplication: matrix incremented, not overwritten
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of a (in real words)
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of b (in real words)
Ci   c,ncc,nrc is the product matrix and respectively the number of
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of c (in real words)
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored added into c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of complex matrix c = a*b  (arrays b,c
Cr     dimensioned complex*16; a real*8 with imag following real)
Cr     call zampy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b
Cr     call zampy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   This version suitable for Cray
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nrci,nccj,ncbj,nrcicj
      double precision ar,ai

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
        ar = a(      nra*i + nca*k)
        ai = a(nia + nra*i + nca*k)
        nrci = nrc*i
C#ifdefC BLAS
C        nrbk = nrb*k
C        call daxpy(m, ar,b(    nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m,-ai,b(nib+nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m, ar,b(nib+nrbk),ncb,c(nic+nrci),ncc)
C        call daxpy(m, ai,b(    nrbk),ncb,c(nic+nrci),ncc)
C#else
        nccj = -ncc
        ncbj = -ncb + nrb*k
        do  20  j = m-1, 0, -1
        nccj = nccj + ncc
        ncbj = ncbj + ncb
        nrcicj = nrci + nccj
        c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
        c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
C#endif
   20 continue
      end
      subroutine zmpy(a,nca,nra,nia,b,ncb,nrb,nib,c,ncc,nrc,nic,n,m,l)
C- Double precision complex matrix multiplication
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of a (in real words)
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of b (in real words)
Ci   c,ncc,nrc is the product matrix and respectively the number of
Ci      between column elements and row elements and between the
Ci      real and imaginary parts of c (in real words)
Ci   n,m: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed n*m*l times.  No attempt is made to optimize
Cr   the outer loops, executed n*m times.
Cr     Examples: product of complex matrix c = a*b  (arrays b,c
Cr     dimensioned complex*16; a real*8 with imag following real)
Cr     call zmpy(a,n,1,ndim**2,b,2*n,2,1,c,2*n,2,1,n,n,n)
Cr     To generate c = a*b
Cr     call zmpy(a,2*n,2,1,b,2,2*n,1,c,2*n,2,1,n,n,n)
Cr   Warning: this routine has not been thoroughly checked!
Cr   This version suitable for Cray
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,nia,ncb,nrb,nib,ncc,nrc,nic,n,m,l
      double precision a(0:*), b(0:*), c(0:*)
      integer i,j,k,nrci,nccj,ncbj,nrcicj
      double precision ar,ai

C --- Initialize array to zero ---
      do  10  i = n-1, 0, -1
        nrci = nrc*i
        nccj = -ncc
        do  10  j = m-1, 0, -1
        nccj = nccj + ncc
        nrcicj = nrci + nccj
        c(nrcicj)     = 0
        c(nic+nrcicj) = 0
   10 continue

C --- Do multiplication ---
      do  20  k = l-1, 0, -1
        do  20  i = n-1, 0, -1
        ar = a(      nra*i + nca*k)
        ai = a(nia + nra*i + nca*k)
        nrci = nrc*i
C#ifdefC BLAS
C        nrbk = nrb*k
C        call daxpy(m, ar,b(    nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m,-ai,b(nib+nrbk),ncb,c(nrci),    ncc)
C        call daxpy(m, ar,b(nib+nrbk),ncb,c(nic+nrci),ncc)
C        call daxpy(m, ai,b(    nrbk),ncb,c(nic+nrci),ncc)
C#else
        nccj = -ncc
        ncbj = -ncb + nrb*k
        do  20  j = m-1, 0, -1
        nccj = nccj + ncc
        ncbj = ncbj + ncb
        nrcicj = nrci + nccj
        c(nrcicj)     = c(nrcicj)     + ar*b(ncbj) - ai*b(nib+ncbj)
        c(nic+nrcicj) = c(nic+nrcicj) + ar*b(nib+ncbj) + ai*b(ncbj)
C#endif
   20 continue
      end

      function ran1()
C- Return a random deviate between 0.0 and 1.0.
C ----------------------------------------------------------------
Ci Inputs
Co Outputs
Co   ran1
Cr Remarks
Cr   Algorithm from Knuth; adapted here from Numerical Recipes, chapter
Cr   Uses three linear
Cr   congruential generators (two for high and low order, a third
Cr   to shuffle).  Use ran1in to initialize.
C ----------------------------------------------------------------
C     implicit none
C Local
      real ran1
      integer j
      integer m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
      parameter (m1=259200,ia1=7141,ic1=54773)
      parameter (m2=134456,ia2=8121,ic2=28411)
      parameter (m3=243000,ia3=4561,ic3=51349)
C Static
      integer ix1,ix2,ix3
      real r(97),rm1,rm2
      common /dran1/ ix1,ix2,ix3,rm1,rm2,r

C Generate next number for each sequence;
c use third to generate random integer between 1 and 97
      ix1 = mod(ia1*ix1+ic1,m1)
      ix2 = mod(ia2*ix2+ic2,m2)
      ix3 = mod(ia3*ix3+ic3,m3)
      j = 1 + (97*ix3)/m3
C#ifdefC TEST
C      if (j .gt. 97 .or. j .lt. 1) pause
C#endif
C Return the table entry ...
      ran1 = r(j)
C And refill it.
      r(j) = (float(ix1) + float(ix2)*rm2)*rm1
      return
      end
      subroutine ran1in(iseed)
C- A simple one-parameter initializer for ran1
C ----------------------------------------------------------------
Ci Inputs
Ci   iseed
Co Outputs
Co   ran1 is set up
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
C Passed
      integer iseed
C Local
      integer j
      integer m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
      parameter (m1=259200,ia1=7141,ic1=54773)
      parameter (m2=134456,ia2=8121,ic2=28411)
      parameter (m3=243000,ia3=4561,ic3=51349)
C To preserve
      integer ix1,ix2,ix3
      real r(97),rm1,rm2
      common /dran1/ ix1,ix2,ix3,rm1,rm2,r

      rm1 = 1./m1
      rm2 = 1./m2
C Seed the first, second and third sequences
      ix1 = mod(ic1-iseed,m1)
      ix1 = mod(ia1*ix1+ic1,m1)
      ix2 = mod(ix1,m2)
      ix1 = mod(ia1*ix1+ic1,m1)
      ix3 = mod(ix1,m3)
C Fill table with sequential uniform deviates generated by first two
      do  11  j = 1, 97
        ix1 = mod(ia1*ix1+ic1,m1)
        ix2 = mod(ia2*ix2+ic2,m2)
        r(j) = (float(ix1) + float(ix2)*rm2)*rm1
   11 continue
      return
      end
C#ifdefC TEST
C      PROGRAM D7R2
CC     Driver for routine RAN1
CC     Calculates pi statistically using volume of unit n-sphere
C      parameter(pi=3.1415926)
C      dimension iy(3),yprob(3)
C      fnc(x1,x2,x3,x4) = sngl(dsqrt(
C     .  dble(x1)**2+dble(x2)**2+dble(x3)**2+dble(x4)**2))
CC
C      call ran1in(-1)
C      do  11  i = 1, 3
C        iy(i) = 0
C   11 continue
C      write(*,'(1x,/,t15,a)') 'Volume of unit n-sphere, n=2,3,4'
C      write(*,'(1x,/,t3,a,t17,a,t26,a,t37,a)')
C     .  '# points','pi','(4/3)*pi','(1/2)*pi^2'
C      do  14  j = 1, 15
C        do  12  k=2**(j-1),2**j
C          x1 = ran1()
C          x2 = ran1()
C          x3 = ran1()
C          x4 = ran1()
C          if (fnc(x1,x2,0.0,0.0) .lt. 1.0) iy(1) = iy(1)+1
C          if (fnc(x1,x2,x3,0.0) .lt. 1.0) iy(2) = iy(2)+1
C          if (fnc(x1,x2,x3,x4) .lt. 1.0) iy(3) = iy(3)+1
C   12   continue
C        do  13  i = 1, 3
C          yprob(i) = 1.0*(2**(i+1))*iy(i)/(2**j)
C   13   continue
C        write(*,'(1x,i8,3f12.6)') 2**j,(yprob(i),i=1,3)
C   14 continue
C      write(*,'(1x,/,t4,a,3f12.6,/)') 'actual',pi,4.0*pi/3.0,0.5*(pi**
C      end
C#endif
      integer function amix(wk,nelts,nmix,mmix,beta,norm,t,kpvt,ipr,tm,
     .                      rmsdel)
C- Anderson mixing of a vector
C ----------------------------------------------------------------
Ci Inputs
Ci   nmix: number of previous iterations to fold into mix
Ci         nmix = 0 => linear mixing (x* = x0)
Ci         Practical experience for LMTO programs shows that nmix>2
Ci         does not work well.
Ci   mmix: maximum number of previous iterations to fold into mix
Ci         (used for dimensioning of work array)
Ci   nelts:number of elements to mix
Ci   wk:   array of dimension (nelts,1+mmix,2) where:
Ci         wk(*,i,1) holds f(xi) (see remarks)
Ci         wk(*,i,2) holds   xi  (see remarks)
Ci   beta: new x is beta f(x*) + (1-beta) x*
Ci   tm:   upper limit to any tj:  if any tj exceeds tm, effective
Ci         nmix is set to zero.
Co Outputs
Co   f(x_i) => f(x_i+1); x_i => x_i+1
Co   wk(*,1,2): new x = beta f(x*) + (1 - beta) x*
Co   rmsdel:rms difference between x_0 and f(x_0)
Co   amix:  returns effective nmix (see input tm)
Cr Remarks
Cr   Given a vector function f(x), where x is some
Cr   vector, we want to find x* such that f(x*) = x*.  We want to find
Cr   x* with the minimum number of computations of f(x).  Supposing
Cr   that we have x0,f(x0); x1,f(x1); x2,f(x2); ...; x_n+1,f(x_n+1).
Cr   (Usually x_j corresponds to x at the jth previous iteration.)
Cr   We take a linear combination x* of x_0, x_1, x_2, ... x_n that
Cr   minimizes <(x* - f(x*))^2>.  We then seek t_1, t_2, ... t_n in
Cr     x* = x_0 -
Cr   To evaluate f(x*) we linearize d(x) = f(x)-x as
Cr     f(x*)-x*  =  d*  =  d_0  -
Cr   Then
Cr     < (d_0 -
Cr   constitute a set of n simultaneous equations in the n unknowns t_k.
Cr   Note that d's enter into these equations, not the f's.
Cr   The d's are stored on disk, unit n1; the x's on unit n2.
Cr   Given the t_k's, x* can be estimated from (2).  To dampen
Cr   instablities, a linear combination of (1) and (2) is taken as
Cr       beta f(x*)  +  (1-beta) x*  = beta d*  +  x*           (4)
Cr   beta is an input parameter, which can usually be taken to be 1.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nelts,mmix,nmix,kpvt(1),ipr
      double precision norm(mmix,mmix),wk(nelts,0:mmix+1,2),t(mmix),
     .  beta,tm,rmsdel
C Local parameters
      integer i,j,nwk,jmix
      character*1 mixis0
      double precision dsqrt,ddot

      nwk = nelts*(mmix+2)

C --- d_0 = f-x  =>  wk(*,0,1) ---
      call daxpy(nelts,-1d0,wk(1,0,2),1,wk,1)

C --- Copy x_0 and d_0 to end of wk:  x*, d* constructed here ---
      call dmcpy(wk,nwk,1,wk(1,mmix+1,1),nwk,1,nelts,2)

C --- Make < (d_0 - d_j) (d_0 - d_k) > and  < d_0 (d_0 - d_j) >  ---
      do  20  i = 1, nmix
        t(i) = ddot(nelts,wk(1,0,1),1,wk(1,0,1),1) -
     .         ddot(nelts,wk(1,0,1),1,wk(1,i,1),1)
        do  20  j = 1, nmix
        norm(i,j) =  ddot(nelts,wk(1,0,1),1,wk(1,0,1),1)
     .             - ddot(nelts,wk(1,0,1),1,wk(1,j,1),1)
     .             - ddot(nelts,wk(1,i,1),1,wk(1,0,1),1)
     .             + ddot(nelts,wk(1,i,1),1,wk(1,j,1),1)
   20 continue

C --- Solve the simultaneous equations ---
      call dsifa(norm,mmix,nmix,kpvt,i)
      if (i .ne. 0) then
        if (ipr .ge. 20) print *,  'AMIX: normal eqns singular'
      else
        call dsisl(norm,mmix,nmix,kpvt,t)
      endif

C --- Set jmix = effective nmix to zero if any t_j exceeds tm ---
      jmix = nmix
      mixis0 = ' '
      do  30  j = 1, nmix
   30 if (dabs(t(j)) .gt. dabs(tm)) jmix = 0
      if (jmix .ne. nmix) mixis0 = '*'
      amix = jmix

C --- Make (d,x)* = (d,x)_0 -
      do  40  j = 1, jmix
        call dscal(nelts,1-t(j),wk(1,mmix+1,1),1)
        call daxpy(nelts,t(j),wk(1,j,1),1,wk(1,mmix+1,1),1)
        call dscal(nelts,1-t(j),wk(1,mmix+1,2),1)
        call daxpy(nelts,t(j),wk(1,j,2),1,wk(1,mmix+1,2),1)
   40 continue

C --- Copy arrays to new positions ---
      do  50  i = mmix, 1, -1
   50 call dmcpy(wk(1,i-1,1),nwk,1,wk(1,i,1),nwk,1,nelts,2)

C --- x* + beta d* ---
   10 call daxpy(nelts,beta,wk(1,mmix+1,1),1,wk(1,mmix+1,2),1)

C --- Calculate rms change ---
      rmsdel = dsqrt(ddot(nelts,wk,1,wk,1)/nelts)

C --- Printout ---
      if (ipr .lt. 20) goto 60
ca    write(*,133) nmix,mixis0,mmix,nelts,beta,tm,rmsdel
      write(*,133) nmix,mixis0,mmix,nelts,beta,rmsdel
ca133 format(/' AMIX:  nmix=',i1,a1,' mmix=',i1,'  nelts=',i3,
ca   .        '  beta=',f7.5,'  tm=',f8.5,'  rmsdel=',1pd8.2)
  133 format(/' AMIX:  nmix=',i1,a1,' mmix=',i1,'  nelts=',i3,
     .        '  beta=',f7.5,'  rmsdel=',1pd8.2)
      if (nmix .gt. 0) write(*,134) (t(j), j=1,nmix)
  134 format(3x,'tj:',7(f8.5,2x))

      if (ipr .le. 30) goto 60
      write(*,110)
      do  12  i = 1, nelts
        if (dabs(wk(i,0,1)) + dabs(wk(i,mmix+1,2)-wk(i,0,2)) .ge. 5.d-7)
     .  write(*,111) i,wk(i,0,2),wk(i,0,2)+wk(i,0,1),
     .                 wk(i,0,1),wk(i,mmix+1,2)
   12 continue

C --- Restore d* and x* + beta d* from end of wk ---
   60 call dmcpy(wk(1,mmix+1,1),nwk,1,wk,nwk,1,nelts,2)

  104 format(1p,4d18.11)
  111 format(i5,4f14.6)
  110 format(14x,'OLD',11X,' NEW',9X,'DIFF',10X,'MIXED')
      end
C --- Simple FORTRAN dynamic memory allocation ---
C Memory is allocated from a single heap.  The heap is declared as
C an integer array in common block W by the main program.  Its size
C is fixed by the size of W as declared by the main program.  It
C is possible to allocated memory from the top of the heap, and to
C free memory already allocated from the top of the heap.  It is
C not possible to free memory allocated from the middle of the heap.
C
C To use these routines, first call WKINIT with NSIZE = size of W.
C
C Character, integer, real, double precision, complex or double
C complex arrays can then be allocated by calling DEFCH, DEFI,
C DEFR, DEFDR, DEFC or DEFDC with an index IPNAME and LENGTH the
C desired length.  IPNAME is returned as the effective offset of
C the integer array W.  If it is desired that the array be
C initialized to zero, pass LENGTH as the negative of the desired
C length.
C
C Allocated memory can be freed by calling RLSE(IPNAME) where
C IPNAME must one of the indices from a previous allocation.  All
C blocks of memory allocated after IPNAME are also released.
C
C WKPRNT turns on or off a trace that prints a message each time
C        memory is allocated or released.  Also when turned on,
C        runs through links each time RLSE is called.
C DEFASK returns the amount of memory remaining in the pool
C WKSAV  sets an internal switch so that a heap allocation request
C        returns the the negative of length of the request in IPNAME
C        when not enough memory is available, rather than aborting.
C        In this case, no memory is allocated and no other internal
C        variables are changed.
C WKINFO prints information about the arrays allocated
C WKCHK  runs through links to see if any array was overwritten
C        call wkchk(string), where string terminated by '$'
C ----------------------------------------------------------------
      SUBROUTINE WKINIT(NSIZE)
C- Initialize conditions for the 'common pool' data W
ca !!!!!!!!!!!
ca    INTEGER W(2)
      INTEGER W(1)
      CHARACTER*1 STRING(60),STR(60)
      logical lopt,lerr
      SAVE
      COMMON /W/ W

      IRND(I) = (I+499)/1000

C ----- DEFINE STORAGE SIZE ------
C  START OF FIRST ARRAY AND MAX NUMBER TO BE DEFINED:
      IP0 = 5
      NDEFMX = 100
      LIMIT = NSIZE
      IPMAX = 0
      NDFDMX = 0
      NDFSUM = 0
      JPR = 0
      IPFREE = 5
      lerr = .false.
      WRITE(*,*) 'WKINIT:  size=',IRND(NSIZE),'K'
      WRITE(*,*) ' '
      RETURN
      ENTRY WKPRNT(JPRINT)
C- Set debug switch for heap management
      if (jprint .eq. 2) then
        jpr = 1-jpr
      else
        JPR = JPRINT
      endif
      RETURN
      entry wksav(lopt)
      lerr = lopt
      return
C ------ SUBROUTINES TO DEFINE ARRAYS OF VARIOUS TYPES -----
      ENTRY DEFCH(IPNAME,LENG)
C- Allocate character array
      LENGTH = (LENG+3)/4
      GOTO 10
      ENTRY DEFI(IPNAME,LENG)
C- Allocate integer array
      LENGTH = LENG
      JOPT = 1
      GOTO 10
      ENTRY DEFR(IPNAME,LENG)
C- Allocate single precision real array
      LENGTH = LENG*I1MACH(17)
      JOPT = 2
      GOTO 10
      ENTRY DEFRR(IPNAME,LENG)
      ENTRY DEFDR(IPNAME,LENG)
C- Allocate double precision real array
      LENGTH = LENG*I1MACH(18)
      JOPT = 3
      GOTO 10
      ENTRY DEFC(IPNAME,LENG)
C- Allocate single precision complex array
      LENGTH = LENG*2*I1MACH(17)
      JOPT = 4
      GOTO 10
      ENTRY DEFCC(IPNAME,LENG)
      ENTRY DEFDC(IPNAME,LENG)
C- Allocate double precision complex array
      LENGTH = LENG*2*I1MACH(18)
      JOPT = 5
   10 IOPT = 0
      IF (LENGTH .LT. 0) THEN
        IOPT = 1
        LENGTH = -LENGTH
      ENDIF
      IF (LENGTH .EQ. 0) LENGTH = 1
      IMOD = 1
      GOTO 83
   84 IPNAME = IPFREE
      if (lerr .and. ipfree+length+2 .gt. limit) then
        ipname = -LENGTH
        if (jpr .gt. 0)  print *,
     .    'ALLOC: heap storage exceeded; returning -LENGTH=',-LENGTH
        return
      endif
      IPFREE = IPFREE + LENGTH + 1
      IPFREE = 4*((IPFREE+2)/4)+1
      IPMAX = MAX0(IPMAX,IPFREE)
      W(IPNAME-1)=IPFREE
      NDEFD = NDEFD + 1
      NDFDMX = MAX0(NDFDMX,NDEFD)
      NDFSUM = NDFSUM + LENGTH
      IF (JPR .GT. 0) WRITE(*,100) NDEFD,LENG,LENGTH,IPNAME,IPFREE-1
  100 FORMAT(' define array',I4,':   els=',I8,'   length=',I8,',',
     .   I8,'  to',I8)
      IF (IPFREE .LE. LIMIT) THEN
        IF (IOPT .NE. 0) GOTO (201,202,203,204,205) JOPT
        RETURN
  201   CALL IINIT(W(IPNAME),-LENG)
        RETURN
  202   CALL AINIT(W(IPNAME),-LENG)
        RETURN
  203   CALL DINIT(W(IPNAME),-LENG)
        RETURN
  204   CALL CINIT(W(IPNAME),-LENG)
        RETURN
  205   CALL ZINIT(W(IPNAME),-LENG)
        RETURN
      ENDIF
      WRITE(*,101) IPFREE
  101 FORMAT(' ALLOC: WORKSPACE OVERFLOW, NEED AT LEAST',I8)
      STOP
C- Release data up to pointer
      ENTRY RLSE(IPNAME)
      IF (IPNAME .GT. LIMIT) STOP 'RLSE: release pointer exceeds limit'
      IF (IPNAME .LT. 3) STOP 'RLSE: release pointer less than 3'
      IF (JPR .eq. 0) goto 82
      imod = 3
      goto 83
   87 WRITE(*,*) 'RLSE from: ',IPNAME
   82 IPFREE = IPNAME
      return
      ENTRY DEFASK(LREST)
C- Return number of words left in common pool
      LREST = LIMIT - IPFREE - 2
      IF (JPR .GT. 0) WRITE(*,*) 'SPACE LEFT=',LREST,'  SINGLE WORDS'
      RETURN
      ENTRY WKINFO
C- Output workspace information
      IMOD = 2
      GOTO 83
  81  WRITE(*,601) IRND(LIMIT),IRND(NDFSUM),IRND(IPMAX-1),
     .             IRND(IPFREE-1),NDFDMX,NDEFD
  601 FORMAT(
     .  /'  total workspace size =',I5,' K',
     .  /'  total space allocated=',I5,' K',
     .  /'  workspace used:    max',I5,' K   now',I4,' K',
     .  /'  arrays defined:    max',I7,  '   now',I6)
      IF (IPFREE .EQ. IP0) RETURN
      if (jpr .gt. 0) WRITE(*,602)
  602 FORMAT(/'  array',6X,'begin',7X,'end',7X,'length')
      IPX = IP0
      DO  30  I = 1, NDEFMX
        IPY = W(IPX-1)
        IF (IPX .EQ. IPFREE) RETURN
        IF (IPY .LT. IP0 .OR. IPY .GT. LIMIT) WRITE(*,*) '   . . . . . '
        IF (IPY .LT. IP0 .OR. IPY .GT. LIMIT) RETURN
        IF (JPR .GT. 0) WRITE(*,603) I,IPX,IPY,IPY-IPX
  603   FORMAT(4(I6,5X))
        IPX = IPY
   30 CONTINUE
      RETURN
      ENTRY WKCHK(STRING)
C- Run through links to see if any dynamic array was overwritten
      IMOD = 0
      DO  88  I = 1, 60
        STR(I)=STRING(I)
        NSTR = I-1
   88 IF (STRING(I) .EQ. '$') GOTO 89
   89 WRITE(*,*) 'WKCHK: ',(STR(I),I = 1,NSTR)
   83 NDEFD = 0
      IPX = IP0
      IPPLOC = -999
      DO  35  I = 1, NDEFMX
        IF (IPX .LT. IP0 .OR. IPX .GT. LIMIT) THEN
          WRITE(*,888) NDEFD,IPX,IPPLOC
  888     FORMAT(' ALLOC: LINK DESTROYED AT START OF ARRAY',I3,
     .     ',  PTR=',I8,' AT',I8)
          STOP
        ENDIF
        IF (IPX .EQ. IPFREE) GOTO 86
        NDEFD = NDEFD + 1
        IPPLOC = IPX - 1
  35  IPX = W(IPPLOC)
  86  CONTINUE
      GOTO (84,81,87), imod
      WRITE(*,360) NDEFD,IPFREE-1
  360 FORMAT('     LINKS OK   NDEFD=',I3,'   SPACE USED=',I7)
      RETURN
      END
      SUBROUTINE CINIT(ARRAY,LENG)
C- Initializes complex array to zero
      INTEGER LENG
      REAL ARRAY(2*LENG)
      CALL AINIT(ARRAY,LENG+LENG)
      RETURN
      END
      SUBROUTINE ZINIT(ARRAY,LENG)
C- Initializes complex array to zero
      INTEGER LENG
      DOUBLE PRECISION ARRAY(2*LENG)
      CALL DINIT(ARRAY,LENG+LENG)
      RETURN
      END
      SUBROUTINE IINIT(ARRAY,LENG)
C- Initializes integer array to zero
      INTEGER LENG
      INTEGER ARRAY(LENG)
      DO 10 I=1,LENG
        ARRAY(I)=0
   10 CONTINUE
      RETURN
      END
      SUBROUTINE AINIT(ARRAY,LENG)
C- Initializes real array to zero
      INTEGER LENG
      REAL ARRAY(LENG)
      DO 10 I=1,LENG
        ARRAY(I)=0.
   10 CONTINUE
      RETURN
      END
      SUBROUTINE DINIT(ARRAY,LENG)
C- Initializes double precision array to zero
      INTEGER LENG
      DOUBLE PRECISION ARRAY(LENG)
      DO 10 I=1,LENG
        ARRAY(I)=0.D0
   10 CONTINUE
      RETURN
      END
      FUNCTION RVAL(ARRAY,INDEX)
C- Returns the real value of ARRAY(INDEX)
      REAL ARRAY(INDEX)
      RVAL=ARRAY(INDEX)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DRVAL(ARRAY,INDEX)
C- Returns the double precision value of ARRAY(INDEX)
      DOUBLE PRECISION ARRAY(INDEX)
      DRVAL=ARRAY(INDEX)
      RETURN
      END
      INTEGER FUNCTION IVAL(ARRAY,INDEX)
C- Returns the integer value of ARRAY(INDEX)
      INTEGER ARRAY(INDEX)
      IVAL=ARRAY(INDEX)
      RETURN
      END
      COMPLEX FUNCTION CVAL(ARRAY,INDEX)
C- Returns the complex value of ARRAY(INDEX)
      COMPLEX ARRAY(INDEX)
      CVAL=ARRAY(INDEX)
      RETURN
      END



      SUBROUTINE DMADD(A,NCA,NRA,SCALEA,B,NCB,NRB,SCALEB,C,NCC,NRC,N,M)
C- general matrix addition
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   b,ncb,nrb is the right matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   c,ncc,nrc is the result matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   n,m: the number of rows and columns, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix linear combination routine,
Cr   adding a subblock of matrix a to a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop.
Cr
Cr   Destination matrix c can coincide with either a or b, provided that
Cr   the transpose of the coincident matrix is not taken.
Cr   Example: Add 3-by-2 block of (transpose of a - .5*b) into c
Cr     call dmadd(a,1,na,1.d0,b,nb,0,-.5d0,c,nc,1,3,2)
Cr     OLD call dmadd(a,na,1,b,nb,0,-.5d0,c,nc,3,2)
C ----------------------------------------------------------------
C
      INTEGER NCA,NRA,NCB,NRB,NCC,NRC,N,M
      DOUBLE PRECISION A(0:1), B(0:1), C(0:1), SCALEA, SCALEB
      INTEGER I,J,IA,IB,IC

      DO  200  I = N-1, 0, -1
        IA = I*NRA+M*NCA
        IB = I*NRB+M*NCB
        IC = I*NRC+M*NCC
        DO  200  J = M-1, 0, -1
        IA = IA-NCA
        IB = IB-NCB
        IC = IC-NCC
        C(IC) = A(IA)*SCALEA + B(IB)*SCALEB
  200 CONTINUE
      RETURN
      END
      subroutine dpmpy(a,b,nscb,nsrb,c,nscc,nsrc,nr,nc,l)
C- matrix multiplication, (packed) (normal) -> (normal)
C ----------------------------------------------------------------
Ci Inputs:
Ci   a is the left matrix (packed)
Ci   b,nscb,nsrb is the right matrix and respectively the spacing
Ci      between column elements and row elements.
Ci   c,nscc,nsrc is the product matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   nr,nc: the number of rows and columns, respectively, to calculate
Ci   l:   length of vector for matrix multiply
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This is a general-purpose matrix multiplication routine,
Cr   multiplying a subblock of matrix a by a subblock of matrix b.
Cr   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
Cr   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
Cr   Arrays are locally one-dimensional so as to optimize inner loop,
Cr   which is executed nr*nc*l times.  No attempt is made to optimize
Cr   the outer loops, executed nr*nc times.
Cr     Examples: product of (nr,l) subblock of a into (l,nc) subblock of
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,nr,nc,l)
Cr     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
Cr     To generate the tranpose of that product, use:
Cr   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,nr,nc,l)
C ----------------------------------------------------------------
C Passed Parameters
      integer nscb,nsrb,nscc,nsrc,nr,nc,l
      double precision a(0:*), b(0:*), c(0:*)
C Local parameters
      double precision sum
      integer i,j,k,offa,offb

      do  20  i = 0, nr-1
        do  20  j = 0, nc-1
        sum = 0
        offa = (i*(i+1))/2
        offb = nscb*j
        do  21  k = 0, i-1
          sum = sum + a(offa)*b(offb)
          offa = offa + 1
          offb = offb + nsrb
   21   continue
        do  22  k = i, l-1
          sum = sum + a(offa)*b(offb)
          offa = offa + k+1
          offb = offb + nsrb
   22   continue
        c(i*nsrc+j*nscc) = sum
   20 continue
      end
      REAL FUNCTION R1MACH(I)
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
C  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
C  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C#ifdef IEEE | RIEEE
      DATA SMALL(1) /     8388608 /
      DATA LARGE(1) /  2139095039 /
      DATA RIGHT(1) /   864026624 /
      DATA DIVER(1) /   872415232 /
      DATA LOG10(1) /  1050288283 /
C#endif
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00014000000000000000B /
C      DATA RMACH(2) / 37767777777777777777B /
C      DATA RMACH(3) / 16404000000000000000B /
C      DATA RMACH(4) / 16414000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C#ifdefC CRAY
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /
C#endif
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C#ifdefC IBM_VM | IBM370 | IBM3080 | IBM3090
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /
C#endif
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C#ifdefC VMS
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /
C#endif
C
c aek
      write(i1mach(2),123)
 123  format('WARNING: R1MACH IS USED !!')
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      R1MACH = RMACH(I)
      RETURN
  999 WRITE(I1MACH(2),1999) I
 1999 FORMAT(' R1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
c aek this version is only for use of d1mach(1)--d1mach(3)
c aek dmach(1): smallest possible value
c aek dmach(2): largest possible value
c aek dmach(3): smallest eps so that 1.d0+eps > 1.d0
      double precision dmach(3)
      integer i
      if (i .lt. 1 .or. i .gt. 3) go to 999
      dmach(1) = 1.d-99
      dmach(2) = 1.d+99
      dmach(3) = 1.d-15
      d1mach = dmach(i)
      return
 999  write(i1mach(2),123) i
 123  format(' D1MACH: I=',i1,' OUT OF BOUNDS !')
      stop
      end
C#define EXTENDED
C This adaptation of I1MACH, is identical to the public version, except
C that I1MACH has additional elements (17,18) specifying the length of a
C real word and of a double precision word.  BEWARE that these numbers
C have so far only been included for machines this program has been
C tested on.
C
      INTEGER FUNCTION I1MACH(I)
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
C                 FOR  FORTRAN 77, THIS IS ALWAYS 1.  FOR FORTRAN 66,
C                 CHARACTER STORAGE UNIT = INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C *** Extensions: *** (only implemented for IEEE machines)
C
C    I1MACH(17) = number of integer words that fit into a real word
C
C    I1MACH(18) = number of integer words that fit into a real*8 word
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
C  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
C  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
C  FOR IMACH(1) - IMACH(4).
C
C#ifdef EXTENDED
      INTEGER IMACH(18),OUTPUT,SANITY
C#elseC
C      INTEGER IMACH(16),OUTPUT,SANITY
C#endif
C
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
C#ifdef IEEE | RIEEE
c aek      DATA IMACH( 1) /    5 /
      DATA IMACH( 2) /    6 /
c aek      DATA IMACH( 3) /    6 /
      DATA IMACH( 4) /    6 /
c aek      DATA IMACH( 5) /   32 /
      DATA IMACH( 6) /    4 /
c aek      DATA IMACH( 7) /    2 /
c aek      DATA IMACH( 8) /   31 /
c aek      DATA IMACH( 9) / 2147483647 /
c aek      DATA IMACH(10) /    2 /
c aek      DATA IMACH(11) /   24 /
c aek      DATA IMACH(12) / -125 /
c aek      DATA IMACH(13) /  128 /
c aek      DATA IMACH(14) /   53 /
c aek      DATA IMACH(15) / -1021 /
C#ifdef EXTENDED
c aek      DATA IMACH(16) /  1024 /
      DATA IMACH(17) /     1 /
      DATA IMACH(18) /     2 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  1024 /, SANITY/987/
C#endif
C#endif
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C#ifdefC CRAY
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C#ifdefC EXTENDED
C      DATA IMACH(16) /  8190 /
C      DATA IMACH(17) /     1 /
C      DATA IMACH(18) /     1 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  8190 /, SANITY/987/
C#endifC
C#endif
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C#ifdefC IBM_VM | IBM370 | IBM3080 | IBM3090
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C#ifdefC EXTENDED
C      DATA IMACH(16) /  63 /
C      DATA IMACH(17) /   1 /
C      DATA IMACH(18) /   2 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  63 /, SANITY/987/
C#endifC   EXTENDED
C#endif   IBM_VM
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA IMACH( 1) /            1 /
C      DATA IMACH( 2) /            1 /
C      DATA IMACH( 3) /            2 /
C      DATA IMACH( 4) /            1 /
C      DATA IMACH( 5) /           32 /
C      DATA IMACH( 6) /            4 /
C      DATA IMACH( 7) /            2 /
C      DATA IMACH( 8) /           31 /
C      DATA IMACH( 9) / :17777777777 /
C      DATA IMACH(10) /            2 /
C      DATA IMACH(11) /           23 /
C      DATA IMACH(12) /         -127 /
C      DATA IMACH(13) /         +127 /
C      DATA IMACH(14) /           47 /
C      DATA IMACH(15) /       -32895 /
C      DATA IMACH(16) /       +32637 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR VAX.
C
C#ifdefC VAX | VMS
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
CC
C#ifdefC EXTENDED
C      DATA IMACH(16) /  127 /
C      DATA IMACH(17) /    1 /
C      DATA IMACH(18) /    2 /, SANITY/987/
C#elseC
C      DATA IMACH(16) /  127 /, SANITY/987/
C#endifC
C#endif
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) STOP 777
C#ifdef EXTENDED
      IF (I .LT. 1  .OR.  I .GT. 18) GO TO 999
c aek
      if (i.ne.2.and.i.ne.4.and.i.ne.6.and.i.ne.17.and.i.ne.18) goto 99
      go to 98
 99   write(output,123) i
 123  format('WARNING: i1mach used with i=',i2)
C#elseC
C      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 999
C#endif
 98   I1MACH=IMACH(I)

C#ifdefC FTN77
C      IF(I.EQ.6) I1MACH=1
C#endif
C#ifdefC GFLOAT_VAX
CC machine constants for g_float optional compilation on the VAX
CC -------------------------------------------------------------
C      double precision dmach(5)
C      data dmach(1) / 1.112536929253601E-308 /
C      data dmach(2) / 4.494232837155789D+307 /
C      data dmach(3) / 1.110223024625157D-016 /
C      data dmach(4) / 2.220446049250313D-016 /
C      data dmach(5) / 0.301029995663981      /
C      d1mach = dmach(i)
C#endif
      RETURN
  999 WRITE(OUTPUT,1999) I
 1999 FORMAT(' I1MACH:  I OUT OF BOUNDS',I10)
      STOP
      END

* iolib.cmp.f below.
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
C#define FORMATTED_READ
C --- FORTRAN INPUT LIBRARY ---
C
C The i/o library provides fortran support for unformatted input
C resembling a NAMELIST statement, though it is more general.  As a
C guiding philosophy, the execution a set of fortran programs is
C intended to be governed by a single, small control file, that serves
C the dual purpose of governing program execution and document what is
C intended to calculate.  It is intended that this file is to be only
C read.  Only portions of a control file might be used by a single
C program; several programs can share a common control file.
C
C As with the fortran NAMELIST statement a quantity is identified by a
C token, such as 'NSPIN='.  The "contents" of the token are a vector (of
C length 0,1 or more) elements of cast char, logical, integer, real or
C double: 'NSPIN= 1' is a typical input for the control file
C
C Input tokens are actually separated into "categories" and data is read
C category by category.  A token always belongs to a category and tokens
C within different categories are inequivalent.  A new category begins
C (and any previous category ends) as soon as a non-blank character
C appears in the first column of a control file.  In the control file
C
C HEADER  Example of the beginning of a typical
C         control file
C VERS    1
C IO      SHOW=F HELP=F VERBOS=3 WKP=F
C
C the category 'HEADER' begins with the 'H' and continues until the 'V'
C in 'VERS', and so on.  Apart from this, input is essentially
C free-format.  Categories need not be in any particular order, but only
C the first occurence of two categories with equivalent names will be
C read.  This is handy for keeping alternative input within one file.
C Categories can either be optional or mandatory; similarly tokens
C within categories can be either optional or mandatory.  It is up to
C the calling routine to handle optional input that is not found.
C
C Names of both tokens and categories end by specifying the last
C character of the name, such as '=' in 'NSPIN='.
C
C
C --- Use of the routines ---
C
C In practice the routines are called as followed: (1) routine rdfile is
C called to load the entire control file into RAM; (2) a category is
C sought and delimited by calling routine getcat; (3) a token is sought
C and contents loaded by routine partok.  Step 3 is repeated for as many
C elements within a category as are desired, and a new category may be
C sought after the first is finished.
C
C Routines getcat and partok can be called in one of three "modes" (the
C mode in effect is governed by variable optio).  In the first mode
C getcat and partok do not attempt to search for anything; instead they
C display what they would have searched for had the next mode been
C selected.  In the next mode, they carry out the requested searches; in
C the next they print out the contents of the variables corresponding to
C the tokens (essentially an all-purpose output to verify what was
C input).  Again, it is best to look at a practical implementation to
C see how this, and many other small options can be used.
C
C
C --- Subcategories ---
C
C There is a facility to further restrict the range of search for tokens
C within a category.  Routine partok search between offsets iostart and
C subsiz from the beginning of a category.  Unless otherwise specified,
C subsiz defaults to the size of the whole category, but it is possible
C to shrink subsiz by a call to routine subcat.  This is particularly
C useful when looking for several occurences of the same token (and
C related information) within a category
C
C
C --- Other routines ---
C
C Routine iprint is intended as a function that returns an integer used
C to control how verbose the output is.  There is a stack of
C levels; iprint returns the top of the stack.  A new value is pushed
C onto the stack by calling pshprt; the previous value is popped by
C calling popprt.  The top two levels can be toggled by calling togprt.
C
C Routine fopna is the file opening routine; there are associated
C entries fopn, fopnn, fopno, fxst, fhndl which when used in conjunction
C with entry fadd, that file logical units with names, so that files
C can be subsequently referenced by name only.
C Routine fclose closes files; dfclose closes and deletes a file.
C All files should be closed by a call to fclose.
C
C Subroutine Query provides a simple "interactive mode".  When query is
C called, nothing happens unless the "query" switch has been set by a
C call to entry initqu.  But if this switch has been set, the user has
C an option to abort program execution, change the verbosity, turn on
C the work array debug switch, or to change a value of a single number
C passed to query.  There is also the option to reset the query switch
C
C There are several string-manipulation routines (copy, concatenate,
C equality, display).
C
C Routine finits performs machine-dependent initialization.
C Routine fexit performs machine-dependent program termination.
C
C --- Variables with common meaning ---
C In common block iolib:
C   recoff: offset to recrd corresponding to first character in file
C           (set by calling program before rdfile is called)
C   reclen: length of file records
C           (can be set in function recln before rdfile is called)
C   nrecs:  number of file records (set by rdfile)
C   maxlen: the maximum allowed number of characters in a category or
C           token (set by calling program before rdfile is called)
C   catbeg: offset to recrd demarcating beginning of a category
C           (set whenever getcat is called)
C   catsiz: the size of a category
C           (set whenever getcat is called)
C   subsiz: the size of a collection of tokens within a subcategory
C           (defaults to catsiz when getcat is called;
C            shrunk when subcat called)
C   iend:   points to the first character after the last match
C           (set whenever getcat or partok is called)
C   ichoos,nchoos: index to current choice of possible inputs and
C           total number of choices, respectively (see description)
C   noerr:  error on return from getcat or partok; see getcat and partok
C           (set whenever getcat or partok is called)
C   optio:  0, 1, or 2 (see description)
C ==============
      SUBROUTINE FINITS
C- Machine and compiler-dependent initialization
C#ifdefC SVS | unix
C      logical lsequ
C      integer j,fext,iargc,iarg,n
C      character*4 extns
C#endif
 
       EXTERNAL DIOLIB
 
       CALL INITQU (.FALSE.)
 
C --- For Lahey F77L, open standard output as list-directed ---
C#ifdefC F77LAHEY
C      open(unit=*,carriage control='LIST')
C#endif
 
C --- Handle floating point exceptions in the IBM VM environment ---
C#ifdefC IBM_VM
C      call errset(208,999,-1)
C#endif
 
C --- Read in extension from command line if SVS or unix compiler ---
C#ifdefC SVS | unix
C      iarg = 1
C   10 iarg = iarg+1
C        if (iargc() .ge. iarg) then
C          call getarg(iarg,extns)
C          if (lsequ(extns,'-',1,' ',n)) goto 10
C          j = fext(extns)
C        endif
C#endif
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE FEXIT ( RETVAL, STRNG )
C- Machine and compiler-dependent program termination
       INTEGER RETVAL
       INTEGER FOPN, FHNDL, LGUNIT
       DOUBLEPRECISION CPUSEC
C#ifdefC APOLLO_BUG
C      character*(10) strng
C#else
       CHARACTER*(*) STRNG
       INTEGER II3, II2, II1
       PARAMETER (II3 = 2, II2 = 1, II1 = 10)
C#endif
       PRINT *
       PRINT *, 'Stop:  ', STRNG
       IF (IPRINT( ) .GE. II1 .AND. CPUSEC( ) .NE. 0) THEN
        WRITE (LGUNIT (II2), 10) CPUSEC( )
        WRITE (LGUNIT (II3), 10) CPUSEC( )
   10   FORMAT(' CPU time, sec:',F12.3)
       END IF
       CALL WKINFO
       IF (FHNDL ('TMP') .GE. 0) CALL DFCLOS (FOPN ('TMP'))
C#ifdefC unix
C      call set_ret_val(retval,1)
C#endif
       STOP 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      DOUBLEPRECISION FUNCTION CPUSEC (  )
C- process cputime, in seconds
C ----------------------------------------------------------------------
Ci Inputs:
Ci   none
Co Outputs:
Co   returns cpu time, in seconds
Cr  Remarks
Cr    On the Apollo: time (in 4 microsecond units) is
Cr    (time(1)*65536 + time(2))*65536 + time(3)
C ----------------------------------------------------------------------
C#ifdefC APOLLO
C      integer time(3)
C      integer*2 t2(3)
C      call proc1_$get_cput(t2)
C      time(1) = t2(1)
C      time(2) = t2(2)
C      time(3) = t2(3)
C      if (time(1) .lt. 0) time(1) = time(1) + 65536
C      if (time(2) .lt. 0) time(2) = time(2) + 65536
C      if (time(3) .lt. 0) time(3) = time(3) + 65536
C      cpusec = 4d-6 *
C     .  ((dble(time(1))*65536 + dble(time(2)))*65536 + dble(time(3)))
C#else
       CPUSEC = 0.D0
C#endif
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      INTEGER FUNCTION RECLN ( I )
C- Returns (and optionally sets) the record length of ascii input files
       INTEGER I
C For io routines ...
       INTEGER RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, SUBSIZ, 
     X   IEND, ICHOOS, NCHOOS, OPTIO
       LOGICAL NOERR
       COMMON /IOLIB/ RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, 
     X   SUBSIZ, IEND, ICHOOS, NCHOOS, NOERR, OPTIO
       IF (I .GT. 0) RECLEN = I
       RECLN = RECLEN
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      INTEGER FUNCTION IPRINT (  )
C- get last integer off print priority stack
C     implicit none
       INTEGER NSTACK
       PARAMETER (NSTACK = 5)
       INTEGER VSTACK(0:4), STACKP
       COMMON /IPRNT/ VSTACK, STACKP
       IPRINT = VSTACK(STACKP)
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE PSHPRT ( VB )
C     implicit none
       INTEGER VB
       INTEGER NSTACK
       PARAMETER (NSTACK = 5)
       INTEGER VSTACK(0:4), STACKP
       INTEGER II1
       PARAMETER (II1 = 5)
       COMMON /IPRNT/ VSTACK, STACKP
       STACKP = MOD (STACKP + 1, II1)
       VSTACK(STACKP) = VB
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE POPPRT
C     implicit none
       INTEGER NSTACK
       PARAMETER (NSTACK = 5)
       INTEGER VSTACK(0:4), STACKP
       INTEGER II2, II1
       PARAMETER (II2 = 5, II1 = 4)
       COMMON /IPRNT/ VSTACK, STACKP
       STACKP = MOD (STACKP + II1, II2)
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE TOGPRT
C     implicit none
       INTEGER NSTACK
       PARAMETER (NSTACK = 5)
       INTEGER VSTACK(0:4), STACKP
       COMMON /IPRNT/ VSTACK, STACKP
       INTEGER ITMP, JTMP
 
       ITMP = VSTACK(STACKP)
       CALL POPPRT
       JTMP = VSTACK(STACKP)
       VSTACK(STACKP) = ITMP
       CALL PSHPRT (JTMP)
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
C --- INTERACTIVE MODE ROUTINES ---
      SUBROUTINE QUERY ( INSTR, CAST, VAR )
C- interactive flow control
C ----------------------------------------------------------------
Ci Inputs
Ci   strng: prompt string
Ci   cast:  <0, if nothing to change, otherwise
Ci          cast is  0=,logical, 2=int, 3=real 4=double
Co Outputs
Co   var:   query will change if requested
Cr Remarks
Cr   At the prompt, user enters either nothing, or one of
Cr     'Snnn', where nnn is number (or T or F for logical variable);
Cr     'Vnnn', where nnn is the new verbosity;
Cr     'W' to toggle printing of work array;
Cr     'I' to turn off interactive prompt;
Cr     'A' to abort program execution
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
       CHARACTER*(*) INSTR
       INTEGER CAST, VAR
       LOGICAL LSET
C Local parameters
       INTEGER LPRMPT
       PARAMETER (LPRMPT = 62)
       CHARACTER*62 PRMPT, RDSTR
       LOGICAL LQUERY, LSEQU, CNVT
       INTEGER I, IMIN, IVBS, IPRINT, J
       CHARACTER*1 PRMPT2(1), RDSTR2(2)
       INTEGER II5, II4, II3, II2, II1
       PARAMETER (II5 = 2, II4 = 0, II3 = 62, II2 = 14, II1 = 1)
       EQUIVALENCE (PRMPT, PRMPT2), (RDSTR, RDSTR2)
       EXTERNAL LSEQU, CNVT, IPRINT
       SAVE IVBS, LQUERY
 
       DATA PRMPT/
     X   '(S)et parm   (V)erbos   (W)kp toggle   (I)act toggle   (A)bort
     X'/ 
 
C#ifndef CRAY | BATCH
       IF (.NOT.LQUERY) RETURN 
 
       PRINT *
       IMIN = II1
       GO TO (1, 2, 1, 1, 1), CAST + 1
C Case nothing to change
C#ifdefC APOLLO_BUG
C      print *, 'QUERY:  '
C#else
       PRINT *, 'QUERY:  ', INSTR
C#endif
       IMIN = II2
       GO TO 10
    2  STOP 'QUERY:  cannot change char or real'
    1  CONTINUE
       CALL DISPLY (INSTR,'=',VAR,VAR,VAR,VAR,' ',II1,CAST)
   10  CONTINUE
       PRINT 333, (PRMPT2(I), I=IMIN,II3)
  333  FORMAT(' ',62A1)
  334  FORMAT(A62)
       PRINT 335, IPRINT( )
C#ifdefC SVS
C  335 format(' V=',i3,' ? '\)
C#elseifC APOLLO
C  335 format(' V=',i3,' ? '$)
C#else
  335  FORMAT(' V=',I3,' ?')
C#endif
       RDSTR = ' '
       READ (*, 334) RDSTR
 
       I = II4
       IF (LSEQU (RDSTR,'S',II1,' ',J) .OR. LSEQU (RDSTR,'s',II1,' ',J)
     X   ) THEN
        IF (.NOT.CNVT (RDSTR2(2),VAR,VAR,VAR,VAR,CAST,II4,' ',I)) PRINT 
     X    *, 'conversion error'
        GO TO 10
       ELSE IF (LSEQU (RDSTR,'V',II1,' ',J) .OR. LSEQU (RDSTR,'v',II1,
     X   ' ',J)) THEN
        IF (.NOT.CNVT (RDSTR2(2),II4,IVBS,II4,II4,II5,II4,' ',I)) THEN
         PRINT *, 'conversion error'
        ELSE
         CALL POPPRT
         CALL PSHPRT (IVBS)
        END IF
        GO TO 10
       ELSE IF (LSEQU (RDSTR,'W',II1,' ',J) .OR. LSEQU (RDSTR,'w',II1,
     X   ' ',J)) THEN
        CALL WKPRNT (II5)
        GO TO 10
       ELSE IF (LSEQU (RDSTR,'I',II1,' ',J) .OR. LSEQU (RDSTR,'i',II1,
     X   ' ',J)) THEN
        LQUERY = .FALSE.
       ELSE IF (LSEQU (RDSTR,'A',II1,' ',J) .OR. LSEQU (RDSTR,'a',II1,
     X   ' ',J)) THEN
        CALL FEXIT (II1,'QUERY')
       END IF
 
C#endif
       RETURN 
 
      ENTRY INITQU ( LSET )
       LQUERY = LSET
 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
C --- FREE-FORMAT INPUT PROCEDURES ---
      SUBROUTINE RDFILE ( UNIT, RECRD, MXRECS, A, RECLEN, NRECS )
C- Read entire file into RAM, up to a maximum allowed number of records
C ----------------------------------------------------------------
Ci Inputs
Ci   unit,recrd,mxrecs,a,reclen,nrecs
Co Outputs
Co   nrecs
Cr Remarks
Cr   Reads in an entire file from unit.
Cr   Error if nrecs exceeds mxrecs
C ----------------------------------------------------------------
C     implicit none
       INTEGER UNIT, MXRECS, RECLEN, NRECS
       LOGICAL RDSTRN
       INTEGER J
       CHARACTER*1 RECRD(0:*), A(*)
       INTEGER II1
       PARAMETER (II1 = 0)
 
       NRECS = II1
   10  CONTINUE
       IF (.NOT.RDSTRN (UNIT,A,RECLEN,.FALSE.)) RETURN 
       CALL STRCOP (RECRD(RECLEN*NRECS),A,RECLEN,'',J)
       NRECS = NRECS + 1
       IF (NRECS .GE. MXRECS) RETURN 
 
c      print *, '------------------------'
c      print 20, (recrd(j), j= 0, reclen*nrecs-1)
c   20 format(1x,72a1)
 
       GO TO 10
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE GETCAT ( RECRD, CAT, TERM, LOPT )
C- find a category in a record
C ----------------------------------------------------------------
Ci Inputs
Ci   recrd,cat,term
Ci   lopt: if error, aborts if unmatched category
Ci   optio: 0: show category that would have been sought
Ci          1: attempt to find category
Ci          2: print out category sought
Co Outputs
Co   catbeg,catsiz,noerr (see remarks)
Cr Remarks
Cr   on return, noerr is always true unless a category was actually
Cr   sought (optio=1) and none was found.
C ----------------------------------------------------------------
C     implicit none
       INTEGER RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, SUBSIZ
       CHARACTER*1 RECRD(0:*), TERM
       CHARACTER*(*) CAT
       INTEGER IRECS, I, J, IEND, ICHOOS, NCHOOS, OPTIO
       LOGICAL LOPT, NOERR
       LOGICAL LSEQU
       INTEGER II4, II3, II2, II1
       PARAMETER (II4 = 3, II3 = 0, II2 = 1, II1 = 5)
       COMMON /IOLIB/ RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, 
     X   SUBSIZ, IEND, ICHOOS, NCHOOS, NOERR, OPTIO
 
       NOERR = .TRUE.
   55  FORMAT(70A1)
       GO TO (1, 2, 3), OPTIO + 1
 
C --- print message showing category sought ---
    1  CONTINUE
       J = II1
       IF (LOPT) J = II2
       CALL MSGIO (RECRD(NRECS*RECLEN),CAT,' ',' ',II2,J,.TRUE.,I)
       RETURN 
 
C --- print out category ---
    3  CONTINUE
C#ifdefC APOLLO_BUG
C      j = 0
C      call chrpos(cat,term,maxlen-1,j)
C      print 11, cat(1:j)
C#else
       PRINT 11, CAT
C#endif
   11  FORMAT(1X,60A)
       RETURN 
 
C --- find a category ---
    2  CONTINUE
       CATBEG = RECOFF
       CATSIZ = II3
       SUBSIZ = II3
       IRECS = NRECS
   10  CONTINUE
       IRECS = IRECS - 1
       IF (.NOT.LSEQU (RECRD(CATBEG),CAT,MAXLEN,TERM,IEND)) THEN
        CATBEG = CATBEG + RECLEN
        IF (IRECS .EQ. 0) THEN
         NOERR = .FALSE.
         IF (LOPT) THEN
          CALL MSGIO (RECRD(NRECS*RECLEN),CAT,' ',' ',II4,II2,.TRUE.,I)
          STOP 
         END IF
         RETURN 
        END IF
        GO TO 10
       END IF
C Found a category.  Next determine the length of the input string
       CATSIZ = CATBEG
   20  CONTINUE
       IRECS = IRECS - 1
       CATSIZ = CATSIZ + RECLEN
       IF (RECRD(CATSIZ) .EQ. TERM .AND. IRECS .GE. 0) GO TO 20
C Clean up and quit
       CATSIZ = CATSIZ - CATBEG
       SUBSIZ = CATSIZ
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      INTEGER FUNCTION PARTOK ( CATEG, TOKEN, TERM, RESULT, CHR, COUNT, 
     X  CAST, ISTART, LOPT )
C- reads a vector of numbers found after a token
C ----------------------------------------------------------------
Ci Inputs
Ci   categ,token,term
Ci   count: maximum number of elements to read (if count passed as <0
Ci          partok will read exactly -count elements or abort w/ error)
Ci   cast:  0=,logical, 1=char, 2=int, 3=real, 4=double
Ci   istart:offset to first character of category where search begins
Ci   In the case of a character field, result should be the max size
Ci   optio: 0: show token that would have been sought
Ci          1: attempt to find token and read contents
Ci          2: print out token and contents
Co Outputs
Co   result:array into which elements are read
Co   noerr: see remarks
Co    iend: offset to first character past match
Co  partok: number of elements actually read (-1 if no attempt to read
Co          an element -- if optio is 0 or 2.)
Cr Remarks
Cr   on return, noerr is always false unless a token was actually
Cr   sought (optio=1) and found (note that this differs from getcat).
Cr   Passing common-block variable iend in place of istart has the
Cr   effect of updating istart to point beyond the end of the match.
Cr
Cr   ichoos is automatically incremented on each call if nchoos .ne. 0
Cr   ichoos and nchoos are automatically reset to 0 if a match is found
Cr   or if ichoos=nchoos
C ----------------------------------------------------------------
C     implicit none
       INTEGER CAST
       CHARACTER*1 CATEG(0:*), TERM
       CHARACTER*(*) CHR, TOKEN
       INTEGER COUNT, RESULT, ISTART
       LOGICAL PARSTR, CNVT, ERRFLG, LOPT
 
C local variables
       INTEGER J, K
       CHARACTER*1 STRN(72)
       CHARACTER*8 NMCAST(0:4)
       CHARACTER*11 OPTION
 
C common block
       INTEGER RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, SUBSIZ, 
     X   IEND, ICHOOS, NCHOOS, OPTIO
       LOGICAL NOERR
       INTEGER II8, II7, II6, II5, II4, II3, II2, II1
       PARAMETER (II8 = 4, II7 = 2, II6 = 0, II5 = 100, II4 = 9)
       PARAMETER (II3 = 1, II2 = 6, II1 = -1)
       COMMON /IOLIB/ RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, 
     X   SUBSIZ, IEND, ICHOOS, NCHOOS, NOERR, OPTIO
       SAVE K, J
 
       DATA NMCAST/'logical$','char$','integer$','real$','double$'/ 
 
       ERRFLG = COUNT .LT. 0
       PARTOK = II1
       NOERR = .FALSE.
       IF (NCHOOS .NE. 0) ICHOOS = ICHOOS + 1
       GO TO (1, 2, 3), OPTIO + 1
 
C --- Print message showing token sought, its cast and length ---
    1  CONTINUE
       OPTION = ' (optional)'
       IF (LOPT) OPTION = ' '
       IF (NCHOOS .NE. 0 .AND. ICHOOS .NE. NCHOOS) OPTION = 
     X   '    --- OR:'
       J = II2
       IF (COUNT .EQ. 0) J = II3
       CALL MSGIO (STRN,' ',TOKEN,TERM,II3,J,.FALSE.,K)
       IF (COUNT .EQ. 0) THEN
        J = K
       ELSE
        CALL STRCOP (STRN(K+1),NMCAST(CAST),II4,'$',J)
        J = J + K - 1
       END IF
       IF (IABS (COUNT) .LE. II3) THEN
        PRINT *, (STRN(K), K=1,J), OPTION
       ELSE
        CALL STRCOP (STRN(J+1),' and length $',II5,'$',K)
        J = J + K - 1
        PRINT *, (STRN(K), K=1,J), IABS (COUNT), OPTION
       END IF
       GO TO 42
 
C --- print out the contents of token ---
    3  IF (NCHOOS .EQ. 0 .OR. ICHOOS .EQ. II3 .AND. COUNT .NE. 0) CALL 
     X   DISPLY (TOKEN,TERM,RESULT,RESULT,RESULT,RESULT,CHR,IABS (COUNT)
     X   ,CAST)
       GO TO 42
 
C --- seek token within a category ---
    2  IEND = ISTART
       PARTOK = II6
       IF (.NOT.PARSTR (CATEG,TOKEN,SUBSIZ,MAXLEN,TERM,IEND,J)) THEN
        IF (ICHOOS .LT. NCHOOS .OR. .NOT.LOPT) GO TO 40
        CALL MSGIO (STRN,CATEG,TOKEN,TERM,II7,II8,.TRUE.,J)
        STOP 
       END IF
       NOERR = .TRUE.
 
       PARTOK = II1
   10  PARTOK = PARTOK + 1
       IF (PARTOK .NE. IABS (COUNT)) THEN
        CALL SKIPBL (CATEG,SUBSIZ,J)
        IF (CAST .EQ. II3) THEN
         CHR = ' '
         CALL STRCOP (CHR,CATEG(J),RESULT,' ',K)
         J = J + K
         GO TO 10
        END IF
        IF (CNVT (CATEG,RESULT,RESULT,RESULT,RESULT,CAST,PARTOK,' ',J)) 
     X    GO TO 10
 
        IF (PARTOK .LT. IABS (COUNT) .AND. ERRFLG) THEN
         CALL MSGIO (STRN,CATEG,TOKEN,TERM,II7,II3,.TRUE.,K)
         PRINT 22, IABS (COUNT), PARTOK
   22    FORMAT(' sought',I3,' elements but found only',I3)
         STOP 
        END IF
       END IF
 
   40  IEND = J
   42  IF (NOERR .OR. ICHOOS .EQ. NCHOOS) THEN
        NCHOOS = II6
        ICHOOS = II6
       END IF
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      LOGICAL FUNCTION SCAT ( IFI, CATEG, TERM, LREWND )
C- scan file for a category
C ----------------------------------------------------------------
Ci Inputs
Ci   ifi,categ,term
Ci   lrewnd: true, rewind file before searching; false, do not
Co Outputs
Co   true if category found, false if not
Cr Remarks
Cr   This is a file version of subroutine getcat
C ----------------------------------------------------------------
       INTEGER IFI
       LOGICAL LREWND
       CHARACTER CATEG, TERM
       CHARACTER*72 A
C Local variables:
       INTEGER CATL0, RECL0
       PARAMETER (CATL0 = 7, RECL0 = 72)
       INTEGER I
       LOGICAL RDSTRN, LSEQU
       INTEGER II2, II1
       PARAMETER (II2 = 7, II1 = 72)
       EXTERNAL RDSTRN, LSEQU
       IF (LREWND) REWIND IFI
       SCAT = .FALSE.
   10  IF (.NOT.RDSTRN (IFI,A,II1,.FALSE.)) RETURN 
       IF (.NOT.LSEQU (CATEG,A,II2,TERM,I)) GO TO 10
       SCAT = .TRUE.
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE SUBCAT ( CATEG, TOKEN, TERM, I )
C- find a subcategory within a category
C ----------------------------------------------------------------
Ci Inputs
Ci   categ,token,term
Ci   i:     offset to first character of category where search begins
Co Outputs
Co   subsiz set to less than range of new token
Cr Remarks
Cr   noerr always returns unchanged.  Does nothing unless optio=1.
C ----------------------------------------------------------------
C passed variables
       INTEGER I
       CHARACTER*1 CATEG(0:*), TERM
       CHARACTER*(*) TOKEN
       LOGICAL PARSTR
 
C local variables
       INTEGER I2, J
 
C Common block
       INTEGER RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, SUBSIZ, 
     X   IEND, ICHOOS, NCHOOS, OPTIO
       LOGICAL NOERR
       COMMON /IOLIB/ RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, 
     X   SUBSIZ, IEND, ICHOOS, NCHOOS, NOERR, OPTIO
 
       GO TO (1, 2, 1), OPTIO + 1
 
C --- seek token within a category ---
    2  SUBSIZ = CATSIZ
       I2 = I
       IF (PARSTR (CATEG,TOKEN,SUBSIZ,MAXLEN,TERM,I2,J)) SUBSIZ = I2
 
    1  RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      LOGICAL FUNCTION PARSTR ( S1, S2, RECLEN, LEN, TERM, I, J )
C- find a substring within a given string
C ----------------------------------------------------------------
Ci Inputs
Ci   i: character where search within string should begin
Co Outputs
Co   i: index to first position of token
Co   j: index to first position after token
Cr Remarks
Cr    seeks match at string(i), string(i+1), ... string(i+reclen) until
Cr    match is found.  returns false if no match found.
C ----------------------------------------------------------------
       INTEGER RECLEN, LEN, I, J
       CHARACTER*1 S1(0:*), S2(0:*), TERM
       LOGICAL LSEQU
 
       PARSTR = .FALSE.
       I = I - 1
   10  I = I + 1
       IF (I .EQ. RECLEN) RETURN 
       IF (.NOT.LSEQU (S1(I),S2,LEN,TERM,J)) GO TO 10
       PARSTR = .TRUE.
       J = J + I
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE CHRPOS ( S, CH, MAXCH, ICH )
C- Finds position of character in string
C ----------------------------------------------------------------
Ci Inputs
Ci   s:   string (declared as s(0:*)
Ci   ch:  character sought
Ci   ich: start search at s(ich)
ci   maxch: see ich
Co Outputs
Co   ich: position of character ch, not to exceed maxch
Cr Remarks
Cr    seeks match at string(i0), string(i0+1) until ch is found or until
Cr    ich = maxch.
C ----------------------------------------------------------------
       INTEGER ICH, MAXCH
       CHARACTER*1 CH, S(0:*)
 
   10  IF (ICH .EQ. MAXCH .OR. S(ICH) .EQ. CH) RETURN 
       ICH = ICH + 1
       GO TO 10
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE SKIPBL ( T, NT, I )
C- Parses string T(I) for blanks
       INTEGER NT, I
       CHARACTER*1 T(0:NT)
   99  IF (T(I) .NE. ' ') RETURN 
       I = I + 1
       IF (I .GE. NT) RETURN 
       GO TO 99
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE TOKMAT ( STRING, TOKEN, N, LEN1, TERM, ITOKEN, LOPT )
C- compare a string to a list of strings
C ----------------------------------------------------------------
Ci Inputs
Ci   string: test string
Ci   token: vector of strings to compare
Ci   n,len: number and length of strings in token
Ci   lopt:  if true, tokmat stops with error message when no match found
Co Outputs
Co   itoken: index to string tokened, -1 if none
C ----------------------------------------------------------------
       INTEGER N, LEN1, ITOKEN
       LOGICAL LSEQU, LOPT
       CHARACTER*1 STRING(*), TERM
       CHARACTER*(*) TOKEN(0:*)
       INTEGER I
       INTEGER II3, II2
       PARAMETER (II3 = -1, II2 = 0)
       INTEGER II1
       DO II1=0,N-1
        ITOKEN = II1
        IF (LSEQU (STRING,TOKEN(II1),LEN1,TERM,I)) RETURN 
       END DO
       ITOKEN = MAX0 (N, II2)
       ITOKEN = II3
       IF (LOPT) THEN
        PRINT *, 'TOKMAT: unmatched ', (STRING(I), I=1,LEN1)
        STOP 
       END IF
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE MSGIO ( STRNG, CATEG, TOKEN, TERM, M1, M2, LOPT, I )
C- Messages for getcat and partok
C ----------------------------------------------------------------
Ci Inputs
Ci   strng,categ,token,term,m1,m2,lopt,i
Co Outputs
Co   string is printed to std out
Co   i: length of printed string
Cr Remarks
Cr   String is composed as: msg(m1) msg(m2)
Cr   If the first character in <categ> is not blank, the string
Cr   'category <categ>' is appended to msg(m1)
Cr   If the first character in <token> is not blank, the string
Cr   '  token <token>' is prepended to msg(m2)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
       CHARACTER*1 STRNG(72)
       CHARACTER*1 CATEG(*), TOKEN(*), TERM
       INTEGER I, M1, M2
       LOGICAL LOPT
C Local parameters
       INTEGER IMSG(11), K
       CHARACTER*1 MSGS(1)
       CHARACTER*70 MSG2
       INTEGER II2, II1
       PARAMETER (II2 = 100, II1 = 72)
       EQUIVALENCE (MSG2, MSGS)
       SAVE K
       DATA IMSG/1,2,17,28,41,52,1,1,1,1,1/ 
       DATA MSG2/
     X   '$read error in $unmatched $  is missing$(optional)$  of cast $
     X'/ 
       CALL STRCOP (STRNG,MSGS(IMSG(M1)),II1,'$',I)
 
       IF (CATEG(1) .NE. ' ') THEN
        CALL STRCAT (STRNG,II2,'$','category: $',II2,'$',I)
        CALL STRCAT (STRNG,II2,'$',CATEG,II2,' ',I)
       END IF
       IF (TOKEN(1) .NE. ' ') THEN
        STRNG(I) = '$'
        CALL STRCAT (STRNG,II2,'$','  token  $',II2,'$',I)
        CALL STRCAT (STRNG,II2,'$',TOKEN,II2,TERM,I)
       END IF
       CALL STRCOP (STRNG(I+1),MSGS(IMSG(M2)),II1,'$',K)
       I = I + K - 1
       IF (LOPT) PRINT 10, (STRNG(K), K=1,I)
   10  FORMAT(' ',72A1)
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      LOGICAL FUNCTION CNVT ( INSTR, RESL, RESI, RESR, RESD, CAST, COUNT
     X  , TERM, J )
C- ASCII conversion to logical, integer, real, double
C ----------------------------------------------------------------
Ci Inputs
Ci   cast:  0=,logical, 1=char, 2=int, 3=real, 4=double
Ci   j:     offset to first character in string to read
Co Outputs
Co   count'th element of resL,resI,resR, or resD is converted
Co   j is put past last character read
Cr Remarks
Cr   An ascii string is converted to logical, integer, real or double
Cr   using the FORTRAN internal read facility.  An unformatted read
Cr   is preferable, but does not conform to the ANSII 77 standard.
Cr
Cr   IOLIB only uses cnvt as a nucleus for partok.
C ----------------------------------------------------------------
       CHARACTER*1 INSTR(0:1), TERM
       INTEGER CAST, COUNT
C these are all equivalent address spaces:
       DOUBLEPRECISION RESD(0:1)
       REAL RESR(0:1)
       INTEGER RESI(0:1)
       LOGICAL RESL(0:1)
 
       INTEGER K, MAXSIZ, J
       PARAMETER (MAXSIZ = 15)
       CHARACTER*15 STRN
       CHARACTER*1 STRN2(15), CJ
       INTEGER II3, II2, II1
       PARAMETER (II3 = 17, II2 = 15, II1 = 100)
       EQUIVALENCE (STRN, STRN2)
       SAVE K
 
C --- Early error checking to avoid problems on some compilers ---
       CNVT = .FALSE.
       IF (INSTR(J) .EQ. ' ') CALL SKIPBL (INSTR(J),II1,J)
       CJ = INSTR(J)
       IF (CAST .NE. 0 .AND. .NOT.(CJ .EQ. '+' .OR. CJ .EQ. '-' .OR. CJ
     X    .EQ. '.' .OR. CJ .GE. '0' .AND. CJ .LE. '9')) RETURN 
 
       STRN = ' '
       CALL STRCOP (STRN,INSTR(J),II2,TERM,K)
C Right justify input (required for some compilers)
       IF (K .LE. II2) THEN
        STRN = ' '
        CALL STRCOP (STRN2(II3-K),INSTR(J),K - 1,TERM,K)
       END IF
       J = J + K + 1
 
       GO TO (1, 2, 3, 4, 5), CAST + 1
C#ifdef FORMATTED_READ
    1  READ (STRN, '(L15)', ERR=20) RESL(COUNT)
C     print *, 'logical string ', strn, 'converted to ', resL(count)
    2  GO TO 10
    3  READ (STRN, '(I15)', ERR=20) RESI(COUNT)
C     print *, 'integer string ', strn, 'converted to ', resI(count)
       GO TO 10
    4  READ (STRN, '(E15.0)', ERR=20) RESR(COUNT)
C     print *, 'real string ', strn, 'converted to ', resR(count)
       GO TO 10
    5  CONTINUE
       READ (STRN, '(E15.0)', ERR=20) RESD(COUNT)
C     print *, 'real*8 string ', strn, 'converted to ', resD(count)
C#elseC
C    1 read(strn,*,err=20) resL(count)
CC     print *, 'logical string ', strn, 'converted to ', resL(count)
C    2 goto 10
C    3 read(strn,*,err=20) resI(count)
CC     print *, 'integer string ', strn, 'converted to ', resI(count)
C      goto 10
C    4 read(strn,*,err=20) resR(count)
CC     print *, 'real string ', strn, 'converted to ', resR(count)
C      goto 10
C    5 continue
C      read(strn,*,err=20) resD(count)
CC     print *, 'real*8 string ', strn, 'converted to ', resD(count)
C      goto 10
C#endif
   10  CNVT = .TRUE.
   20  RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE DISPLY ( TOKEN, TERM, RESL, RESI, RESR, RESD, CHR, 
     X  COUNT, CAST )
C- Called as a nucleus only by partok
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr   Used because fortran cannot treat same space as different casts.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
       INTEGER COUNT, CAST
C#ifdefC APOLLO_BUG
C      character*(*) token
C      character*4 chr(0:*)
C#else
       CHARACTER*(*) TOKEN, CHR(0:*)
C#endif
       CHARACTER TERM
C these are all equivalent address spaces:
       DOUBLEPRECISION RESD(0:COUNT)
       REAL RESR(0:COUNT)
       INTEGER RESI(0:COUNT)
       LOGICAL RESL(0:COUNT)
C Local variables:
       INTEGER CATL0
       PARAMETER (CATL0 = 7)
       INTEGER I
       CHARACTER*7 STRN
       INTEGER II1
       PARAMETER (II1 = 7)
 
       STRN = ' '
       CALL STRCOP (STRN,TOKEN,II1,TERM,I)
 
       GO TO (1, 2, 3, 4, 5), CAST + 1
       PRINT 10, STRN
       RETURN 
    1  PRINT 10, STRN, (RESL(I), I=0,COUNT-1)
   10  FORMAT(3X,1P,A8,1X,20(20(L1,1X)/12X))
       RETURN 
    2  PRINT *, STRN, (CHR(I), I=0,COUNT-1)
       RETURN 
    3  PRINT 30, STRN, (RESI(I), I=0,COUNT-1)
   30  FORMAT(3X,1P,A8,1X,20(10I5/12X))
       RETURN 
    4  PRINT 40, STRN, (RESR(I), I=0,COUNT-1)
   40  FORMAT(3X,1P,A8,1X,20(5G14.7/12X))
       RETURN 
    5  PRINT 50, STRN, (RESD(I), I=0,COUNT-1)
   50  FORMAT(3X,1P,A8,1X,20(3G15.8/12X))
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      LOGICAL FUNCTION RDSTRN ( UNIT, A, LEN, LOPT )
       INTEGER UNIT, LEN
       LOGICAL LOPT
       CHARACTER*1 A(LEN)
 
       RDSTRN = .TRUE.
       READ (UNIT, 10, END=20) A
       IF (LOPT) PRINT 10, A
       RETURN 
   10  FORMAT(72A1)
   20  RDSTRN = .FALSE.
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      LOGICAL FUNCTION LSEQU ( S1, S2, LEN1, TERM, I )
C- Determine whether two strings are equal or not
C ----------------------------------------------------------------
Ci Inputs
Ci   s1,s2: strings to compare
Ci   len:   maximum length of string
Ci   term:  terminator
Co Outputs
Co   lsequ: returned true or false
Co   i:     number of characters tokened (including terminator)
Cr Remarks
Cr   string comparison continues until terminator encountered or
Cr   len characters are checked.
C ----------------------------------------------------------------
       INTEGER I, LEN1
       CHARACTER*1 S1(LEN1), S2(LEN1), TERM
       INTEGER II2
       PARAMETER (II2 = 0)
       INTEGER II1
c     lsequ = .true.
c     do  10  i = 1, len
c       lsequ = (lsequ .and. s1(i) .eq. s2(i))
c       if (s1(i) .eq. term) return
c  10 continue
       LSEQU = .FALSE.
       DO II1=1,LEN1
        I = II1
        IF (S1(II1) .NE. S2(II1)) RETURN 
        IF (S1(II1) .EQ. TERM) GO TO 15
       END DO
       I = MAX0 (LEN1, II2) + 1
   15  LSEQU = .TRUE.
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE STRCOP ( DEST, SOURCE, LEN, TERM, I )
C- copy one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   dest,source: source and destination strings
Ci   len:   maximum number of characters to copy
Ci   term:  terminator
Co   i:     number of characters copied (including terminator)
Cr Remarks
Cr   string copy continues until term encountered or
Cr   len characters are checked.
C ----------------------------------------------------------------
       INTEGER LEN
       CHARACTER*1 DEST(LEN), SOURCE(LEN), TERM
       INTEGER I
       INTEGER II1
       PARAMETER (II1 = 0)
       IF (LEN .EQ. 0) RETURN 
       I = II1
   10  I = I + 1
       DEST(I) = SOURCE(I)
       IF (DEST(I) .NE. TERM .AND. I .LT. LEN) GO TO 10
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE STRCAT ( S1, LEN1, TERM1, S2, LEN2, TERM2, I )
C- concatenate one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   s1,s2: source string and string to concatenate
Ci   len1:  maximum length of s1
Ci   term1: terminator for s1
Ci   len2:  maximum length of s2
Ci   term2: terminator for s2
Co Outputs
Co   i:     number of characters in s1 (including terminator)
Cr Remarks
Cr   concatenation continues until term encountered or
Cr   len2 characters are concatenated.
C ----------------------------------------------------------------
       INTEGER LEN1, LEN2
       CHARACTER*1 S1(LEN1), S2(LEN2), TERM1, TERM2
       INTEGER I, J
       INTEGER II1
       PARAMETER (II1 = 0)
 
       I = II1
   10  I = I + 1
       IF (S1(I) .NE. TERM1 .AND. I .LT. LEN1) GO TO 10
       J = II1
       IF (S1(I) .EQ. TERM1) I = I - 1
   20  J = J + 1
       I = I + 1
       S1(I) = S2(J)
       IF (S2(J) .NE. TERM2 .AND. J .LT. LEN2) GO TO 20
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE STRCAT2 ( S1, S2, LEN, I )
C- concatenate one string to another
C ----------------------------------------------------------------
Ci Inputs/Outputs
Ci   s2,s1: source string and string to concatenate
Ci   len :  number of characters in s1 and s2
Ci   i    : last len - i cahracters of s2 copied
Ci          into first len - i char. in s1
Co Outputs
Cr Remarks
C ----------------------------------------------------------------
       INTEGER LEN
       CHARACTER*1 S1(LEN), S2(LEN)
       INTEGER I, J
 
       DO J=1,LEN
        IF (J + I .GT. LEN) THEN
         S1(J) = ' '
        ELSE
         S1(J) = S2(J+I)
        END IF
       END DO
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
C --- FILE HANDLING ---
      INTEGER FUNCTION FOPNA ( NAM, UNIT, SWITCH )
C- File opening
C ----------------------------------------------------------------
Ci Inputs
Ci    nam:  LMTO file name
Ci   unit:  logical unit
Ci   switch:switch governing mode of file opening
Co Outputs
Co   fopna returns logical unit number for file name
Co     (file is opened if it is not already open)
Cr Remarks
Cr   fopna(name,integer unit,integer switch): general file opening:
Cr      logical unit is passed as part of the call.
Cr      Returns logical unit.
Cr   fadd(name,unit,switch): sends a new name to append to a
Cr     list of known file names and associated logical units.
Cr     Initializer for fopn, fopno and fopno, routines which return
Cr     logical units for files from the name.
Cr     Returns 1 + logical unit.
Cr   fopn(name) file opening, for which a logical unit and switch has
Cr     already been given by a call to fadd.  Same as fopna, except
Cr     that logical unit and switches have been previously stored.
Cr     Returns logical unit.
Cr   fopnn(name) same as fopn, but open with status='NEW'
Cr   fopno(name) same as fopn, but open with status='OLD'.
Cr   fxst(name) same as fopn, but inquires as to whether file exists
Cr   fext(ext) changes the extension from the default ".dat".  Valid
Cr     only for operating systems which permit extensions.
Cr   fhndl(name) returns with logical unit associated with name,
Cr     and -1 if none exists.
Cr
Cr   Switch is a composite of integers a, b, c, d, stored as digits
Cr   abcd, base 2.
Cr     Bits 0,1: = 0, open the file as 'UNKNOWN'
Cr               = 1, open the file as 'OLD'
Cr               = 2, open the file as 'NEW'
Cr     Bit 2:    if set, open file as unformatted.
Cr     Thus switch=5 opens file unformatted, status='old'
C ----------------------------------------------------------------
C Passed parameters
       CHARACTER*(*) NAM
       INTEGER SWITCH, UNIT
C Local parameters
       INTEGER FEXT, FOPN, FOPNO, FOPNN, FHNDL, FADD, FXST, I, ISTA, 
     X   IUNIT, ISW, N, BIT, IPRINT, I1MACH
       CHARACTER*11 FTNFMT
       CHARACTER*8 FNAM, FTNSTA
       LOGICAL BITTST, JSOPEN, LDUM
       INTEGER MXNAM
       PARAMETER (MXNAM = 12)
       CHARACTER*4 PASNAM(12)
       INTEGER PASSW(12), PASUNI(12), NNAM
       INTEGER II10, II9, II8, II7, II6, II5, II4, II3, II2
       PARAMETER (II10 = 100, II9 = 8, II8 = 12, II7 = -2, II6 = -1)
       PARAMETER (II5 = 2, II4 = 1, II3 = 4, II2 = 0)
       SAVE NNAM, PASUNI, PASSW, PASNAM
       INTEGER II1
 
C --- If an extension is allowed or desired ... ---
C#ifdefC EXTDOTDAT
C      integer extlen
C      parameter (extlen=3)
C      character*(extlen+1) ext
C      save ext
C#endif
 
       BITTST(N, BIT) = MOD (N, BIT * 2) .EQ. BIT + MOD (N, BIT)
 
       DATA NNAM/0/ 
C --- Extension for all file names, if desired ---
C#ifdefC EXTDOTDAT
C      data ext /'.dat'/
C#endif
 
C --- Purge any leading blanks ---
       FNAM = NAM
       I = II2
       CALL SKIPBL (FNAM,II3,I)
C#ifdefC CRAY
C        if ( i .ne. 0 ) then
C        call strcat2(fnam,fnam,4,i)
C        endif
C#else
       FNAM = FNAM(I + 1:4)
C#endif
 
C#ifdefC unix
C      call str_lo_case(fnam)
C#endif
       IUNIT = UNIT
       ISW = SWITCH
       ISTA = II2
       GO TO 20
 
      ENTRY FOPN ( NAM )
       ISTA = II2
       GO TO 2
 
      ENTRY FOPNO ( NAM )
       ISTA = II4
       GO TO 2
 
      ENTRY FOPNN ( NAM )
       ISTA = II5
       GO TO 2
 
      ENTRY FXST ( NAM )
       ISTA = II6
       GO TO 2
 
      ENTRY FHNDL ( NAM )
       ISTA = II7
       FHNDL = II6
       GO TO 2
 
      ENTRY FADD ( NAM, UNIT, SWITCH )
       FADD = UNIT + 1
       NNAM = NNAM + 1
       IF (NNAM .GT. II8) STOP 'iolib: too many file names'
       PASNAM(NNAM) = NAM
C#ifdefC unix
C      call str_lo_case(pasnam(nnam))
C#endif
       PASUNI(NNAM) = UNIT
       PASSW(NNAM) = SWITCH
       RETURN 
 
      ENTRY FEXT ( NAM )
C#ifdefC EXTDOTDAT
C      ext = '.'//nam(1:extlen)
C      fext = 0
C#endif
       RETURN 
 
    2  CONTINUE
       FNAM = NAM
C#ifdefC unix
C      call str_lo_case(fnam)
C#endif
       DO II1=1,NNAM
        IF (FNAM .EQ. PASNAM(II1)) THEN
c          fnam = pasnam(i)
         IUNIT = PASUNI(II1)
         FHNDL = IUNIT
         ISW = PASSW(II1) - MOD (PASSW(II1), II3) + ISTA
         GO TO 20
        END IF
c        PRINT *, 'dbg fopn name check', i, nam, pasnam(i)
       END DO
       IF (ISTA .EQ. -2) RETURN 
       WRITE (I1MACH (II3), *) 'fopn: name mismatch: file ', FNAM
       STOP 
 
   20  CONTINUE
       IF (ISTA .EQ. -2) RETURN 
 
       FTNFMT = 'FORMATTED'
       IF (MOD (ISW, II9) .EQ. MOD (ISW, II3) + II3) FTNFMT = 
     X   'UNFORMATTED'
 
C --- Attach an extension ---
C#ifdefC EXTDOTDAT
C      call strcat(fnam,4,' ',ext,4,' ',i)
C#endif
 
C --- Handle INQUIRE statements ---
       IF (ISTA .LT. 0) THEN
        FXST = II2
C#ifdefC PRECONNECTED_UNITS | IBM_VM
C        INQUIRE(UNIT=iunit,EXIST=ldum)
C#else
        INQUIRE (FILE=FNAM, EXIST=LDUM) 
C#endif
        IF (LDUM) FXST = II4
        RETURN 
       END IF
 
       FTNSTA = 'UNKNOWN'
       IF (MOD (ISW, II3) .EQ. II4) FTNSTA = 'OLD'
       IF (MOD (ISW, II3) .EQ. II5) THEN
        FTNSTA = 'NEW'
        CALL FCLOSE (IUNIT)
C#ifdefC PRECONNECTED_UNITS | IBM_VM
C        open(iunit,FORM=ftnfmt,STATUS='UNKNOWN')
C#else
        OPEN (IUNIT, FILE=FNAM, FORM=FTNFMT, STATUS='UNKNOWN') 
C#endif
        CLOSE (UNIT=IUNIT, STATUS='DELETE') 
       END IF
 
c      print *, 'dbg: fopn',fnam,iunit,isw,ista,ftnfmt,ftnsta
 
       FOPNA = IUNIT
       FOPN = IUNIT
       FOPNO = IUNIT
       FOPNN = IUNIT
 
       IF (.NOT.JSOPEN (IUNIT)) THEN
        IF (IPRINT( ) .GE. II10) PRINT 300, FTNFMT, FNAM, FTNSTA, IUNIT
  300   FORMAT(/' FOPEN: opening',A12,' file ',A9,' status=',A7,
     X    ', unit=',I2)
C#ifdefC PRECONNECTED_UNITS | IBM_VM
C        open(iunit,FORM=ftnfmt,STATUS=ftnsta)
C#else
C#ifdefC APOLLO_PRE_OS10.1
C        if (ftnfmt .eq. 'FORMATTED') then
C          open(iunit,FILE=fnam,FORM=ftnfmt,STATUS=ftnsta)
C        else
C          open(iunit,FILE=fnam,FORM=ftnfmt,STATUS=ftnsta,RECL=9999999)
C        endif
C#else
        OPEN (IUNIT, FILE=FNAM, FORM=FTNFMT, STATUS=FTNSTA) 
C#endif
C#endif
       END IF
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      LOGICAL FUNCTION JSOPEN ( UNIT )
C- adds unit to list of open files, returns whether or not aleady open
C ----------------------------------------------------------------
Ci Inputs
Ci   unit
Co Outputs
Co   jsopen
Cr Remarks
Cr   maxfil is maximum number of files that may be open at one time
Cr   minuni is the lowest logical unit number
Cr   All file closings should be made through entry fclose
C ----------------------------------------------------------------
C     implicit none
       INTEGER UNIT
       INTEGER I, IPRINT
       INTEGER MAXFIL, MINUNI
       PARAMETER (MAXFIL = 15, MINUNI = 10)
       INTEGER UNITAB(0:14), NOPEN
       INTEGER II7, II6, II5, II4
       PARAMETER (II7 = 15, II6 = -1, II5 = 0, II4 = 102)
       COMMON /FUNITS/ UNITAB, NOPEN
       INTEGER II3, II2, II1
 
       IF (IPRINT( ) .GE. II4) PRINT 20, UNIT, (UNITAB(I), I=0,NOPEN-1)
   20  FORMAT(/' JSOPEN: check logical unit',I3,' among open units: ',15
     X   I3)
 
       JSOPEN = .TRUE.
       II2 = II5
       II3 = II6
       DO II1=NOPEN-1,II2,II3
        IF (UNITAB(II1) .EQ. UNIT) RETURN 
       END DO
 
       JSOPEN = .FALSE.
       UNITAB(NOPEN) = UNIT
       NOPEN = NOPEN + 1
       IF (NOPEN .GT. II7) STOP 'JSOPEN: too many files'
       RETURN 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE FCLOSE ( UNIT )
C- closes an open file, removing unit from the stack
C ----------------------------------------------------------------
Ci Inputs
Ci   unit
Co Outputs
Co   none
Cr Remarks
Cr   use in conjunction with jsopen
C ----------------------------------------------------------------
c      implicit none
       INTEGER UNIT
       INTEGER I, IPRINT
       CHARACTER*6 CLSTAT
       INTEGER MAXFIL, MINUNI
       PARAMETER (MAXFIL = 15, MINUNI = 10)
       INTEGER UNITAB(0:14), NOPEN
       INTEGER II3, II2, II1
       PARAMETER (II3 = 0, II2 = 999, II1 = 100)
       COMMON /FUNITS/ UNITAB, NOPEN
 
       CLSTAT = 'KEEP'
 
   10  CONTINUE
       IF (IPRINT( ) .GE. II1) PRINT 20, UNIT
   20  FORMAT(/' FCLOSE: closing',I3)
 
       DO I=NOPEN-1,0,-1
        IF (UNITAB(I) .EQ. UNIT) THEN
         UNITAB(I) = II2
        END IF
       END DO
       CALL ISHELL (NOPEN,UNITAB)
       IF (NOPEN .EQ. 0 .OR. UNITAB(MAX (NOPEN-1, II3)) .NE. II2) THEN
        IF (IPRINT( ) .GE. II1) PRINT *, 
     X    'FCLOSE: attempt to close unopened file'
       ELSE
        CLOSE (UNIT=UNIT, STATUS=CLSTAT) 
        NOPEN = NOPEN - 1
       END IF
       RETURN 
      ENTRY DFCLOS ( UNIT )
       CLSTAT = 'DELETE'
       GO TO 10
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      INTEGER FUNCTION LGUNIT ( I )
C- Returns stdout for i=1, log for i=2
       INTEGER I, FOPNO, I1MACH
       INTEGER II1
       PARAMETER (II1 = 2)
       EXTERNAL FOPNO, I1MACH
 
       LGUNIT = I1MACH (II1)
       IF (I .EQ. II1) LGUNIT = FOPNO ('LOG')
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE HEADL2 ( NAME, IFI )
C-  Puts a heading line into file, unit ifi
       INTEGER IFI
       CHARACTER*8 NAME
       WRITE (IFI, 300) NAME
  300  FORMAT(' -------------------------  START ',A8,
     X   ' -------------------------')
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE POSEOF ( IUNIT )
C- Positions file at end-of-file
C Passed parameters
       INTEGER IUNIT
C Local parameters
       INTEGER I, NREC
       INTEGER II1
       PARAMETER (II1 = 0)
 
       NREC = II1
       REWIND IUNIT
       DO I=1,10000
        READ (IUNIT, 100, END=90, ERR=91) 
        NREC = I
  100   FORMAT(A1)
       END DO
       WRITE (*, 200) IUNIT
  200  FORMAT(' POSEOF: no EOF found for file',I3)
       RETURN 
   90  CONTINUE
       REWIND IUNIT
       DO I=1,NREC
        READ (IUNIT, 100) 
       END DO
   91  CONTINUE
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
C --- Virtual memory routines ---
      SUBROUTINE VMEM ( O1, O2 )
C- Virtual memory routines
C ----------------------------------------------------------------
Ci Inputs
Ci   o1,o2
Co Outputs
Co   Nothing
Cr Remarks
Cr
C ----------------------------------------------------------------
C Passed parameters
       INTEGER O1, O2
C Local parameters
       LOGICAL LSAVE
       INTEGER LEN, OFFST, FOPN, FOPNO, IFI, I, II
C heap:
       INTEGER W(1)
       COMMON /W/ W
       SAVE OFFST, LEN, LSAVE
       INTEGER II1
 
       OFFST = O1
       LEN = -O2
       LSAVE = O2 .LT. 0
       IF (O2 .LT. 0) O2 = O1
       RETURN 
 
      ENTRY VMEMS ( I )
       IF (.NOT.LSAVE) RETURN 
       IFI = FOPN ('TMP')
       REWIND IFI
       II1 = -IFI
       DO II=1,I
        CALL VMEM2 (II1,W(OFFST),LEN)
       END DO
       RETURN 
 
      ENTRY VMEMG ( I )
       IF (.NOT.LSAVE) RETURN 
       IFI = FOPNO ('TMP')
       REWIND IFI
       DO II=1,I
        CALL VMEM2 (IFI,W(OFFST),LEN)
       END DO
       CALL FCLOSE (IFI)
 
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
      SUBROUTINE VMEM2 ( IFI, W, LEN )
       INTEGER IFI, LEN, W(LEN)
 
       IF (IFI .GT. 0) READ (IFI) W
       IF (IFI .LT. 0) WRITE (-IFI) W
      END
C     KAP/Digital_UA_F      3.1a k280615 970519 o5r3so3  30-Jan-1998 11:21:07
C#ifdefC FUNIT
C      subroutine funit(unit)
CC- returns next available logical unit for file opening
CC ----------------------------------------------------------------
CCo Outputs
CCo   unit
CCr Remarks
CCr   maxfil is maximum number of files that may be open at one time
CCr   minuni is the lowest logical unit number
CCr   All file closings should be made through entry fclose
CC ----------------------------------------------------------------
CC      implicit none
C      integer unit
C      integer maxfil,minuni
C      parameter (maxfil=15,minuni=10)
C      integer unitab(0:maxfil-1),nfiles
C      integer i
C      save
C      data nfiles /0/
C
C      do  10  i = 0, nfiles-1
C   10 if (unitab(i) .ne. i+minuni) goto 20
C      i = nfiles
C   20 unit = i + minuni
C      unitab(nfiles) = unit
C      nfiles = nfiles+1
C      call ishell(nfiles,unitab)
C      return
C
C      entry fclose(unit)
C      do  30  i = 0, nfiles-1
C   30 if (unitab(i) .eq. unit) goto 40
C      i = nfiles
C   40 continue
C      unitab(i) = 1000
C      call ishell(nfiles,unitab)
C      nfiles = nfiles-1
C
C      endfile unit
C      close(unit=unit)
C
C      return
C      end
C#endif  FUNIT
C      subroutine prtstr(unit,string,term,len)
C      integer unit,len
C      character*1 string(*),term
C      i = 0
C   88 i = i+1
C      if (string(i) .ne. term .and. i .le. len) goto 88
C      write(unit,891) (string(j), j=1,i-1)
C  891 format(72a1)
C      return
C      end
C      character*1 function chrint(i)
C      integer i
C      character*1 csym(10)
C      data csym/'0','1','2','3','4','5','6','7','8','9'/
C      chrint = csym(i)
C      return
C      end
C      logical function ischar(ch)
C      character*1 ch
C      ischar = (ch .ge. 'a' .and. ch .lt. 'z' .or.
C     .          ch .ge. 'A' .and. ch .lt. 'Z' .or.
C     .          ch .ge. '0' .and. ch .le. '9')
C      return
C      end
C --- Default initialization ---
      BLOCK DATA DIOLIB
C     implicit none
       INTEGER CATL0, RECL0
       PARAMETER (CATL0 = 7, RECL0 = 72)
 
C for iprint...
       INTEGER NSTACK
       PARAMETER (NSTACK = 5)
       INTEGER VSTACK(0:NSTACK-1), STACKP
       COMMON /IPRNT/ VSTACK, STACKP
 
C For io routines ...
       INTEGER RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, SUBSIZ, 
     X   IEND, ICHOOS, NCHOOS, OPTIO
       LOGICAL NOERR
       COMMON /IOLIB/ RECOFF, RECLEN, NRECS, MAXLEN, CATBEG, CATSIZ, 
     X   SUBSIZ, IEND, ICHOOS, NCHOOS, NOERR, OPTIO
 
C for logical unit routines ...
       INTEGER MAXFIL, MINUNI
       PARAMETER (MAXFIL = 15, MINUNI = 10)
       INTEGER UNITAB(0:MAXFIL-1), NOPEN
       COMMON /FUNITS/ UNITAB, NOPEN
 
       DATA RECOFF/0/ RECLEN/RECL0/ MAXLEN/CATL0/ ICHOOS/0/ NCHOOS/0/ 
       DATA VSTACK/30,30,30,30,30/ STACKP/0/ 
       DATA NOPEN/0/ 
 
      END
