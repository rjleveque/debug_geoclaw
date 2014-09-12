
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,
     &                  ql,qr,auxl,auxr,wave,s,amdq,apdq)
c     =====================================================
c
c     # Riemann solver for the the shallow water equations

c     
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1 
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #
c     #            amdq = A^- Delta q, 
c     #            apdq = A^+ Delta q, 
c     #                   the decomposition of the flux difference
c     #                       f(qr(i-1)) - f(ql(i))
c     #                   into leftgoing and rightgoing parts respectively.
c     #               
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c     # maux=0 and aux arrays are unused in this example.
c
c
      implicit double precision (a-h,o-z)
c     common /comrp/ grav                       !hard-wired below
c
      dimension wave( meqn, mwaves, 1-mbc:maxm+mbc)
      dimension    s(mwaves, 1-mbc:maxm+mbc)
      dimension   ql(meqn, 1-mbc:maxm+mbc)
      dimension   qr(meqn, 1-mbc:maxm+mbc)
      dimension  apdq(meqn, 1-mbc:maxm+mbc)
      dimension  amdq(meqn, 1-mbc:maxm+mbc)
      dimension  auxl(maux, 1-mbc:maxm+mbc)
      dimension  auxr(maux, 1-mbc:maxm+mbc)
c
c
c     #internal variables
      dimension W1(meqn)
      dimension W2(meqn)
      dimension W3(meqn)

      grav = 9.81
      abs_tol=1.d-10
      rel_tol=1.d-4



      do 30 i = 2-mbc, mx+mbc-1
         
         
         hl=qr(1,i-1)
         hr=ql(1,i)
         hul=qr(2,i-1)
         hur=ql(2,i)
         hvl=qr(3,i-1)
         hvr=ql(3,i)
         bl=auxr(1,i-1)
         br=auxl(1,i)

c     #deal with near dry states

       if (hl.lt.abs_tol) then
          hl=0.d0
          ul=0.d0
          vl=0.d0
       else
          ul=hul/hl
          vl=hvl/hl
       endif

       if (hr.lt.abs_tol) then
          hr=0.d0
          ur=0.d0
          vr=0.d0
       else
          ur=hur/hr
          vr=hvr/hr
       endif

c     #determine direction of the Riemann problem

       if (ixy.eq.1) then !x is the normal direction

c         #determine dryness in area

         if ((hl.lt.abs_tol).and.(hr.lt.abs_tol)) then !totally dry area
            s1=0.d0
            s2=0.d0
            s3=0.d0
            W1(1)=0.d0
            W1(2)=0.d0
            W1(3)=0.d0
            W2(1)=0.d0
            W2(2)=0.d0
            W2(3)=0.d0
            W3(1)=0.d0
            W3(2)=0.d0
            W3(3)=0.d0            
            go to 25

         elseif (hr.lt.rel_tol*hl) then ! dry on the right
            s1=ul-sqrt(grav*hl)
            s3=ul+2*sqrt(grav*hl)
            s2=s3
            if ( (br-bl.gt.hl).and.(ul.lt.abs_tol) ) then
               deltab=hl
            else
               deltab=br-bl
            endif

         elseif (hl.lt.rel_tol*hr) then ! dry on the left
            s1=ur-2*sqrt(grav*hl)
            s3=ur+sqrt(grav*hl)
            s2=s1
            if ( (bl-br.gt.hr).and.(ur.gt.-abs_tol) ) then
               deltab=-hr
            else
               deltab=br-bl
            endif

         else  !wet on both sides
            !compute Roe avgs
            hhat=0.5d0*(hr+hl)
            chat=sqrt(grav*hhat)
            uhat=(sqrt(hl)*ul + sqrt(hr)*ur)/(sqrt(hl)+sqrt(hr))
            !determine traditional HLL speeds
            s1=min(ul-sqrt(grav*hl),uhat-chat)
            s3=max(ur+sqrt(grav*hr),uhat+chat)
            s2=uhat 
            
            deltab=br-bl
         endif

c     #calculate single middle state

         psidx=-0.5*(hl+hr)*grav*deltab
         fl1=hl*ul
         fl2=hl*ul**2 + 0.5d0*grav*hl**2
         fr1=hr*ur
         fr2=hr*ur**2 + 0.5d0*grav*hr**2

         hm = (s3*hr-s1*hl-(fr1-fl1))/(s3-s1)

         if (hm.lt.0.d0) then
            write(*,*) 'problem with hm: hm=',hm
         elseif(hm.lt.abs_tol) then
            hm=0.d0
            hum=0.d0
         else
            hum=(s3*hr*ur-s1*hl*ul-(fr2-fl2)+psidx)/(s3-s1)
         endif

c     #if middle state straddles interface with topography, divide into 2

         if ((s1.gt.-abs_tol).or.(s3.lt.abs_tol)
     &                  .or.(abs(psidx).lt.abs_tol)) then !only 1 middle state

            hrstar=hm
            hlstar=hm
            W1(1)=hm-hl
            W1(2)=hum-hl*ul
            W1(3)=hm*vl-hl*vl
            W2(1)=0.d0
            W2(2)=0.0d0
            W2(3)=hm*vr-hm*vl
            W3(1)=hr-hm
            W3(2)=hr*ur-hum
            W3(3)=hr*vr-hm*vr

         else!divide into 2 states w/fluxes that approx balance delta source

            alpha=s1/s3
            beta=1-s1/s3
            discrim=(alpha**2-1.d0)*(2.*psidx/grav)+(beta*hm)**2
            if (discrim.ge.0.d0) then
               hlstar=(-alpha*beta*hm-sqrt(discrim))/(alpha**2-1.d0)
            elseif (psidx.gt.0.d0) then
               hlstar=0.d0
            else !no root and psidx is negative 
               hlstar=(1-s3/s1)*hm
            endif
            hrstar=alpha*hlstar+beta*hm

c           #prevent large velocities resulting from division            
            if (hlstar.lt.abs_tol) then
               hlstar=0.d0
               hulstar=0.d0
               hrstar=(1-s1/s3)*hm
            else
               hulstar=hum
            endif

            if (hrstar.lt.abs_tol) then
               hrstar=0.d0
               hurstar=0.d0
               hlstar=(1-s3/s1)*hm
            else
               hurstar=hum
            endif

c           #determine waves
            if (s2.ge.0.d0) then
               W2(1)=0.d0
               W2(2)=0.d0
               W2(3)=hrstar*vr-hrstar*vl
            else
               W2(1)=0.d0
               W2(2)=0.d0
               W2(3)=hlstar*vr-hlstar*vl
            endif
            W1(1)=hlstar-hl
            W1(2)=hulstar-hl*ul
            W1(3)=hlstar*vl-hl*vl
            W3(1)=hr-hrstar
            W3(2)=hr*ur-hurstar
            W3(3)=hr*vr-hrstar*vr
         endif
               
      else !ixy=2, y is the normal direction

c         #determine dryness in area

         if ((hl.lt.abs_tol).and.(hr.lt.abs_tol)) then !totally dry area
            s1=0.d0
            s2=0.d0
            s3=0.d0
            W1(1)=0.d0
            W1(2)=0.d0
            W1(3)=0.d0
            W2(1)=0.d0
            W2(2)=0.d0
            W2(3)=0.d0
            W3(1)=0.d0
            W3(2)=0.d0
            W3(3)=0.d0            
            go to 25

         elseif (hr.lt.rel_tol*hl) then ! dry on the right
            s1=vl-sqrt(grav*hl)
            s3=vl+2*sqrt(grav*hl)
            s2=s3
            if ( (br-bl.gt.hl).and.(vl.lt.abs_tol) ) then
               deltab=hl
            else
               deltab=br-bl
            endif

         elseif (hl.lt.rel_tol*hr) then ! dry on the left
            s1=vr-2*sqrt(grav*hl)
            s3=vr+sqrt(grav*hl)
            s2=s1
            if ( (bl-br.gt.hr).and.(vr.gt.-abs_tol) ) then
               deltab=-hr
            else
               deltab=br-bl
            endif

         else  !wet on both sides
            !compute Roe avgs
            hhat=0.5d0*(hr+hl)
            chat=sqrt(grav*hhat)
            vhat=(sqrt(hl)*vl + sqrt(hr)*vr)/(sqrt(hl)+sqrt(hr))
            !determine traditional HLL speeds
            s1=min(vl-sqrt(grav*hl),vhat-chat)
            s3=max(vr+sqrt(grav*hr),vhat+chat)
            s2=vhat !seems to work for now
            
            deltab=br-bl
         endif

c     #calculate single middle state

         psidx=-0.5*(hl+hr)*grav*deltab
         fl1=hl*vl
         fl2=hl*vl**2 + 0.5d0*grav*hl**2
         fr1=hr*vr
         fr2=hr*vr**2 + 0.5d0*grav*hr**2

         hm = (s3*hr-s1*hl-(fr1-fl1))/(s3-s1)

         if (hm.lt.0.d0) then
            write(*,*) 'problem with hm: hm=',hm
         elseif(hm.lt.abs_tol) then
            hm=0.d0
            hvm=0.d0
         else
            hvm=(s3*hr*vr-s1*hl*vl-(fr2-fl2)+psidx)/(s3-s1)
         endif

c     #if middle state straddles interface with topography, divide into 2

         if ((s1.gt.-abs_tol).or.(s3.lt.abs_tol)
     &                  .or.(abs(psidx).lt.abs_tol)) then !only 1 middle state

            hrstar=hm
            hlstar=hm
            W1(1)=hm-hl
            W1(2)=hm*ul-hl*ul
            W1(3)=hvm-hl*vl
            W2(1)=0.d0
            W2(2)=hm*ur-hm*ul
            W2(3)=0.d0
            W3(1)=hr-hm
            W3(2)=hr*ur-hm*ur
            W3(3)=hr*vr-hvm

         else!divide into 2 states w/fluxes that approx balance delta source

            alpha=s1/s3
            beta=1-s1/s3
            discrim=(alpha**2-1.d0)*(2.*psidx/grav)+(beta*hm)**2
            if (discrim.ge.0.d0) then
               hlstar=(-alpha*beta*hm-sqrt(discrim))/(alpha**2-1.d0)
            elseif (psidx.gt.0.d0) then
               hlstar=0.d0
            else !no root and psidx is negative 
               hlstar=(1-s3/s1)*hm
            endif
            hrstar=alpha*hlstar+beta*hm

c           #prevent large velocities resulting from division            
            if (hlstar.lt.abs_tol) then
               hlstar=0.d0
               hvlstar=0.d0
               hrstar=(1-s1/s3)*hm
            else
               hvlstar=hvm
            endif

            if (hrstar.lt.abs_tol) then
               hrstar=0.d0
               hvrstar=0.d0
               hlstar=(1-s3/s1)*hm
            else
               hvrstar=hvm
            endif

c           #determine waves
            if (s2.ge.0.d0) then
               W2(1)=0.d0
               W2(2)=hrstar*ur-hrstar*ul
               W2(3)=0.d0
            else
               W2(1)=0.d0
               W2(2)=hlstar*ur-hlstar*ul
               W2(3)=0.d0
            endif
            W1(1)=hlstar-hl
            W1(2)=hlstar*ul-hl*ul
            W1(3)=hvlstar-hl*vl
            W3(1)=hr-hrstar
            W3(2)=hr*ur-hrstar*ur
            W3(3)=hr*vr-hvrstar
         endif

      endif

c     #set values for outside loop

 25      s(1,i)=s1
         s(2,i)=s2
         s(3,i)=s3

         wave(1,1,i)=W1(1)
         wave(2,1,i)=W1(2)
         wave(3,1,i)=W1(3)
         wave(1,2,i)=W2(1)
         wave(2,2,i)=W2(2)
         wave(3,2,i)=W2(3)
         wave(1,3,i)=W3(1)
         wave(2,3,i)=W3(2)
         wave(3,3,i)=W3(3)


c
   30    continue
c
         
         do 90 i=2-mbc,mx+mbc-1
            do 95 m=1,meqn
               amdq(m,i)=0.0d0
               apdq(m,i)=0.0d0
               do 100 n=1,mwaves
                  amdq(m,i)=amdq(m,i)+min(0.0d0,s(n,i))*wave(m,n,i)
                  apdq(m,i)=apdq(m,i)+max(0.0d0,s(n,i))*wave(m,n,i)
 100           continue
 95         continue
 90      continue


      return
      end
