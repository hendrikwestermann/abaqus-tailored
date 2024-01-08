C**C*DECK ABPHAS
c
C*****************************************************************************
C
C     Interface to Abaqus/Umat for phase transformation
C                                                                   
C     Chair of Engineering Mechanics, University of Paderborn 
c
c     imate = 65:   (fbphas.for)
c
c         Hendrik Westermann                                 05.08.2020
c         last modified (  )                                 05.08.2020
c
c**********************************************************************
c
C     VERZEICHNIS:        FEAP/FEABA
*USER SUBROUTINE
      SUBROUTINE UMAT     (stress,statev,ddsdde,sse,spd,scd,
c      SUBROUTINE UMAT_PHAS (stress,statev,ddsdde,sse,spd,scd,
     &                 rpl,ddsddt,drplde,drpldt,
     &                 stran,dstran,time,dtim,temp,dtemp,
     &                 predef,dpred,cmname,
     &                 ndi,nshr,ntens,nstatv,props,nprops,
     &                 coords,drot,pnewdt,
     & CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
c
c     Variables passed to this program:
c
c     stress(6):   stress tensor at beginning of increment. 
c     statev(16):  solution dependent state variables. 
c     ddsdde(6,6): jacobian of constitutive model 
c     sse:         specific elastic strain energy. (unused here) 
c     spd:         specific plastic dissipation. (unused here) 
c     scd:         specific creep dissipation. (unused here) 
c     rpl:         volumetric heat generation. (unused) 
c     ddsddt(6):   variation of stress increments w.r.t. temperature 
c     drplde:      variation of RPL w.r.t. temperature 
c     drpldt:      variation of RPL w.r.t. temperature increment 
c     stran(6):    total strains at beginning of increment 
c     dstran(6):   array of mechanical strain increments 
c     time(2):     array of step time and total time 
c     dtim:        time increment 
c     temp:        temperature at start of increment (different to fb.glue) 
c     dtemp:       increment of temperature 
c     predef(20):  interpolated predefined state variables (unused) 
c     dpred(20):   array of increments of predefined state variables 
c     cmname(8):   name of the material option (unused) 
c     ndi=3:       number of direct stress components at this point 
c     nshr=3:      numb of engineering shear stresses 
c     ntens=6:     size of the stress or strain component array 
c     nstatv=nhis: number of solution state variables in this model 
c     nprops:      number of material constants 
c     coords(3):   current coordinates of this point 
c     drot(3):     rotation increment matrix 
c     pnewdt:      ratio of suggested new time increment to current one 
c     celent:      characteristic element length 
c
c     additional variables:
c     ipoint               number of integration point
c     nelement             number of element
c     imsg1                printout in cuser/umat
c     imsg2                printout in odeint
c     ddt                  time step subincrement
c     nvar,x1,x2,h1,       parameters of integration method
c     hmin,epsilon
c
c     contents of array props (material property definition)
c
c     header information:
c     props(1)=  integer, not used here
c     props(2)=  Young's modulus
c     props(3)=  Poisson Ratio
c     props(*)=  ltes (for printing)
c---------------------------------------------------------------------------
      include 'aba_param.inc'
      parameter (nstrh=6,nhish=50,nparah=500,nplotvh=100,isen=1)
      character*8 cmname    

      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     2 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     3 props(*),coords(3),drot(3,3),TempP(40)
 
c...dimension for internal variables
      dimension  eps(nstrh), epsn(nstrh), sig(nstrh), cm(nstrh*nstrh),
     &           his1(nhish),his2(nhish), plotv(nplotvh),para(nparah),
     &           energy(10)
      dimension  kpar(isen),his3(isen),temp_par(isen),                                   
     &           dksig(isen),dkeps(isen),dkepsn(isen),dkth(isen),
     &           dkthn(isen)
      dimension  dksigvol(10),verbsgp(100)

      character*10 tpvpf(nplotvh)                                                     
      common /kiofile/ ifile,jfile
      common /kqqqeleinf3d/ ipoint,nelement

c....Copy external variables from ABAQUS
      nstr     = ntens       !! No of stresses
      do i = 1,nstr
         eps(i)  = stran(i) + dstran(i)   !! total strains at end of time step
         epsn(i) = stran(i)               !! total strains at beginning of time step
      enddo
      
c....convert engineering strain to tensor strain
      do i = 4,nstr
         eps(i)  = eps(i)/2d0
         epsn(i) = epsn(i)/2d0
         if (abs(eps(i)) .le. 1e-16) then
             eps(i) = 0d0
         end if
         if (abs(epsn(i)) .le. 1e-16) then
             epsn(i) = 0d0
         end if         
      enddo

      tempn     = temp
      temp      = temp + dtemp  !! temperature
      dtime     = dtim          !! time increment

c...Copy material block
       npar  = nprops
       do i = 1,nprops
          para(i) = props(i) 
       enddo

      ndm = 1
      lpv = 1

c....Internal variables
      nq       = 0           !! No of internal variables
      nstrf    = 0           !! for coupled problems, not used here
      icoup    = 0           !!               -- " --
      ncoup    = 0           !!               -- " --
      npar_act = 0           !! for identification, not used here
c    jfile    = 6           !! output parameter for common kiofile
      if (ntens .eq. 3) iptyp = 1  !!For ESZ(KS)
      if (ntens .eq. 4) iptyp = 3
      if (ntens .eq. 6) iptyp = 4                    

c....Call FEAP/INA material subroutine
      do i = 1,3
         if (i.eq. 1) jsw = 3      !! Prepare for plot variables
         if (i.eq. 2) jsw = 2      !! Prepare for iteration 
         if (i.eq. 3) jsw = 4      !! Compute stresses and tangent

c      write(jfile,fmt='(''++++ jsw = '',I4 )') jsw

      CALL FBPHAS (jsw,imate,incr,isym,nhis,nhis3,iptyp,
     &                 ntens,nq,npvpf5,tpvpf,
     &                 para,npar,eps,epsn,
     &                 his1,his2,plotv,nplotv,
     &                 sig,cm,ifail,dtime, temp, tempn, 
     &                 npar_act,kpar,his3,dksig,dkeps,dkepsn)

c...test
      if (i.eq. 1)  then
              nplotv = npvpf5
      endif

c....Check again npar
         if (i .eq. 3) then 
            if(nprops .ne. npar) then
                write(jfile,fmt='(''**** phas nprops = '',i4,
     &            '' ne npar = '',i8 )') nprops,npar                 
                stop'*** umatphas: npar .ne. nprops, see *.dat file'
            endif
         endif

c...get material parameter and history variables for i=2
         if (i.eq. 2) then 
           if(ntens .ne. nstr) then
               write(jfile,fmt='(''**** umatphas: ntens = '',i4,
     &           '' ne nstr = '',i8 )') ntens,nstr                 
               stop'*** umat: ntens .ne. nstr'
           endif
             nhis = max(1,nhis)     
             if(nstatv .ne. nhis) then
                 write(jfile,fmt='(''**** umatphas: nstatv = '',i4,
     &                             '' ne nhis = '',i4 )')
     &                                nstatv,nhis                 
                 stop'*** umat: nhis .ne. nstatv, see *.dat file'
             endif     
c...test
c           do ipar = 1,npar
c              para(ipar) = props(ipar)
c              write(jfile,fmt='(''para('',i3,'')  = '',e12.6)')
c     &                              ipar,para(ipar)
c           enddo
             do ihis = 1,nhis
                his1(ihis) = statev(ihis) 
                his2(ihis) = statev(ihis)
             enddo
         endif                   
      enddo

      ltes = 0  !! test
c....for ifail=1 write warning and reduce stepsize
      if ( ifail .gt. 0) then
           write(6,9992)
           write(6,9994) nelement,ipoint, nok,nbad,epsilon
           pnewdt=0.25
      else

c...store  stresses, tangent, history variables
           do itens=1,ntens
                stress(itens)=sig(itens)
           enddo
           if (ltes .eq. 5)
     &     write(jfile,1000)(sig(i), i = 1,ntens)
 1000      format ('sig = ',6e14.6)

           do itens=1,ntens
             do jtens=1,ntens
              ddsdde(itens,jtens)=cm((jtens-1)*ntens+itens)
             enddo
              if (ltes .eq. 5)
     &        write(jfile,1001)(ddsdde(itens,j), j = 1,ntens)
 1001         format ('ddsdde = ',6e14.6)
           enddo

c           do itens=4,ntens
c             do jtens=1,ntens
c              ddsdde(itens,jtens)=ddsdde(itens,jtens)/2.d0
c             enddo
c           enddo

           if (ltes .eq. 5)
     &         write(jfile,1002)(statev(ihis), ihis = 1,nhis)
 1002      format ('UMAT-FINAL his-vals.  = ',8e14.6)

           do ihis = 1,nhis
              statev(ihis) = his2(ihis)
           enddo
      endif

 9991 format('    **********SUBROUTINE UMATPHAS  ***********')
 9992 format('stop (SUBROUTINE CUSER/umat): istop=1')
 9993 format('element,ipoint= ',2i5)
 9994 format('element,ipoint,nok,nbad,epsilon=',4i10,f10.3)

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS (jsw     ,imate, incr ,isym  ,nhis  ,
     &                   nhis3   ,iptyp, nstrn,nq    ,npvpf5,
     &                   tpvpf   ,para , npar ,eps   ,epsn  ,
     &                   his1    ,his2 , plotv,nplotv,sig   ,
     &                   cm      ,ifail, dtime,temp  ,tempn ,
     &                   npar_act,kpar , his3,dksig ,dk_eps,
     &                   dk_epsn )
c***********************************************************************
      implicit double precision (a-h, o-z)
      common /iofile/ ifile, jfile
      parameter (nparh = 57)
      dimension eps(*), epsn(*), his1(*), his2(*),
     &          para(*), cm(*), plotv(*), sig(*)
      dimension kpar(*), his3(*), dksig(*), dk_eps(*), dk_epsn(*)
      character*10 tpvpf(100)

c**********************
c....read input data
c**********************
      if (jsw .eq. 1 .or. jsw .eq. 2) then
        CALL FBPHAS_IN(jsw, imate, para, npar, nhis, incr, nhis3,
     &                 iptyp, isym)
        return
c*********************************
c....prepare for plot variables
c*********************************
      elseif (jsw .eq. 3) then
        CALL FBPHAS_PP(npvpf5, tpvpf, iptyp, para)
        return
c**********************************
c....compute stresses and tangent
c**********************************
      elseif (jsw .eq. 4) then
        ifail = 0
        CALL FBPHAS_COM (para  ,nparh,his1 ,his2 , nhis   ,
     &                   nhis3 ,eps  ,epsn ,plotv, nplotv ,
     &                   sig   ,cm   ,ifail,nstrn, dtime  ,
     &                   jfile ,temp ,tempn,npvpf5)
c        write(6,fmt='(2x,''eps(i)= '',12e12.4)') (eps(i),i=1,4)
c        write(6,fmt='(2x,''sig(i)= '',12e12.4)') (sig(i),i=1,4)
        return
c*********************************
c....compute sensitivity
c*********************************
      elseif (jsw .eq. 5 .or. jsw .eq. 6) then
        stop '**** FBPHAS: sensitivity not implemented'
        return
      else
        stop 'FBPHAS:  wrong jsw'
      endif
      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_COM (para, npar, his1, his2, nhis,
     &                       nhis3 ,eps  ,epsn ,plotv, nplotv ,
     &                       sig   ,cm   ,ifail,nstrn, dtime  ,
     &                       jfile ,temp ,tempn,npvpf5)
c***********************************************************************
c....compute stresses and tangent
c***********************************************************************
      implicit double precision (a-h, o-z)
      parameter (nstrh = 6, neq = 1)
      dimension eps(nstrn), epsn(nstrn), his1(nhis), his2(nhis),
     &       para(*), cm(nstrn, nstrn), plotv(nplotv), sig(nstrn),
     &       deveps_Tr(nstrh), depsp(nstrh)

c....prepare input data for iteration
      ltes  = 0 !int(para(npar))

c...Temperature change from temp [Celsius] to theta [Kelvin]

      theta  = temp  + 273.15
      thetan = tempn + 273.15

c....get variables from history array
      CALL FBPHAS_HIS1 (eps, his1, nhis, nstrn,
     &       neq, deveps_Tr, tr_eps, epsv_n,
     &       thetamax, theta, thetan, sigvn, zmn,
     &       zbn, rink_n, zan, zfn, zpn, zamax, dsigvn,
     &       ZAF, ZAP, ZAB , ZAM)

c....iterative  calculation of stresses and internal variables
      CALL FBPHAS_STRE(para, npar, sig, dtime, cm,
     &       depsp, tr_eps, nstrn, deveps_Tr, ifail, sigv,
     &       epsv_n, epsv, ltes, jfile,
     &       theta, thetan, thetamax,  sigvn, zmn, zm, dsigvn,
     &       rink_n, rink, zbn, zb,
     &       za, zan,zf, zfn, zp, zpn, zamax, dsigv,
     &       ZAF, ZAP, ZAB , ZAM)

c....insert variables into history array
      CALL FBPHAS_HIS2(0, his1, his2, nhis, nstrn,
     &       neq, depsp, epsv, thetamax, sigv, zm, dsigv,
     &       zb, rink, za, zf, zp, zamax ,
     &       ZAF, ZAP, ZAB , ZAM)

c....insert plot variables
      CALL FBPHAS_IP(plotv, nplotv, para, sig,  sigv,
     &       nstrn, epsv, zm, theta, zb, za, zf, zp,
     &       ZAF, ZAP, ZAB , ZAM)

      return
      end

c**********************************************************************
      SUBROUTINE FBPHAS_STRE (para, npar, sig, dtime, cm,
     &   depsp, tr_eps, nstrn, deveps_Tr, ifail, sigv,
     &   epsv_n, epsv, ltes, jfile, theta, thetan, thetamax,
     &   sigvn, zmn, zm, dsigvn, rink_n, rink, zbn, zb,
     &   za, zan, zf, zfn, zp, zpn, zamax, dsigv , ZAF, ZAP, ZAB , ZAM)
c***********************************************************************
c     compute stresses and tangent moduli
c***********************************************************************
      implicit double precision (a-h, o-z)
      parameter (neq = 1, nstrh = 6)
      dimension para(*), sig(nstrn), cm(nstrn, nstrn),
     &          deveps_Tr(nstrn), depsp(nstrn)
      dimension sig_tr(nstrh), vn(nstrh), sig_dev(nstrh), sig_vol(nstrh)
      dimension xjac(neq, neq), res(neq), y(neq), jact(neq),
     &          deps_res(neq*(nstrh+1)), deps_Y(neq*(nstrh+1))

      itype = para(1)
      y0  = para(5)

      par_a = para(55)
      par_b = para(56)

      za = 0d0
      zb = 0d0
      zm = 0d0
      zf = 0d0
      zp = 0d0

c....phase fraction austenite
      CALL FBPHAS_AUST (para, za, zan, dza, thetan,
     &                  theta, thetamax,zf,zfn)

c....phase fraction ferrite
      CALL FBPHAS_FERRIT (para, dtime, ZAF, zfn, zf, dzf, za, zan,
     &                    theta, thetan, zamax, dza)

c....phase fraction pearlite
      CALL FBPHAS_PEARLITE (para, dtime, ZAP, zpn, zp, dzp,
     &                      za, zan, theta, dza, zf, dzf,thetamax)

c....phase fraction bainite
      CALL FBPHAS_BAIN_IMP (para, ZAB, zb, zbn, dzb,
     &       za, zan, rink_n, rink, dtime,
     &       theta, epsv_n, sigvn, thetamax)

c....phase fraction martensite
      CALL FBPHAS_MART (para, ZAM, zmn, zm, dzm, za, zan, theta)

c....initializing unknowns y and plastic strain increment
      Y(1) = 0d0
      CALL FBPHAS_ZEROS(depsp, nstrn)

c....check predictor
      lsw = 0
      CALL FBPHAS_RES (lsw, res, xjac, Y, deps_res,
     &       xmu, compar, para, dtime, cm, depsp,
     &       tr_eps, nstrn, deveps_Tr, ifail, sigv,
     &       epsv_n, epsv, ltes, jfile, theta,
     &       resTr, res_norm, sig_tr, vn, phi_tr, sigvn,zm)

      if (ifail .gt. 0) return

c....purely elastic material behaviour
      if (phi_tr .lt. 1.d-12)  then

c....deviatoric part of elastic stress
        do i = 1, nstrn
            sig_dev(i) = 2d0*xmu*deveps_Tr(i)
        enddo

c...volumetric part of elastic stress
        CALL FBPHAS_ZEROS(sig_vol,nstrn)
        do i = 1,3
            sig_vol(i) = compar*tr_eps
        enddo

c....total elastic stress
        do i = 1,nstrn
            sig(i) = sig_dev(i) + sig_vol(i)
        enddo

c....elasticity matrix
        CALL FBPHAS_ZEROS(cm,nstrn*nstrn)
        if (nstrn .eq. 1) then
            cm(1,1) = 2d0 * xmu
        else
            do i = 1,3
                cm(i,i) = cm(i,i) + 2d0*xmu
            enddo
            do i = 4,nstrn
                cm(i,i) = cm(i,i) + 2d0*xmu/2d0
            enddo

            xlm = compar - 2d0*xmu/3d0
            do i = 1,3
                do j = 1,3
                    cm(i,j) = cm(i,j) + xlm
                enddo
            enddo
        endif

        dsigv = sigv - sigvn
        return
      endif

c      if (y0 .lt. 500) then

c....perform nonlinear iteration loop
      CALL FBPHAS_ITE (para    ,npar  ,sig      ,dtime  ,cm    ,
     &                 nstrn   ,ifail ,deveps_Tr,tr_eps ,resTr ,
     &                 deps_res,deps_Y,depsp    ,sigv   ,xjac  ,
     &                 jact    ,Y     ,res      ,neq    ,epsv_n,
     &                 epsv    ,ltes  ,jfile    ,theta  ,thetan,
     &                 dzm   ,sigvn ,zmn     ,zm  ,zb,
     &                 dzb ,za,zan )

c      endif

      dsigv = sigv - sigvn
      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_ITE (para, npar, sig, dtime, cm,
     &   nstrn, ifail, deveps_Tr, tr_eps, resTr, deps_res, deps_Y,
     &   depsp, sigv, xjac, jact, Y, res, neq,
     &   epsv_n, epsv, ltes, jfile, theta, thetan, dzm, sigvn,
     &   zmn, zm, zb, dzb, za, zan)
c***********************************************************************
c     compute stresses and tangent moduli
c**********************************************************************
      implicit double precision (a-h, o-z)
      parameter (neqh = 1, nparh = 57, nstrh = 6)
      double precision FBPHAS_norm
      dimension para(*), sig(nstrn), cm(nstrn, nstrn),
     &          deveps_Tr(nstrn), depsp(nstrn), deps_Y(neq, nstrn+1)
      dimension xjac(neq, neq), res(neq), Y(neq), jact(neq)
      dimension deps_res(neqh*(6+1)),
     &          sig_dev(nstrh), sig_vol(nstrh), sig_tr(nstrh),
     &          vn(nstrh), s(neqh)
      dimension res_store(neqh,2), y_store(neqh,2), i_store(neqh,2)

      logical done, tol_ch

      tol_ch(resr, resa, su, tolr, tols) =
     &       (((resr .lt. tolr   .or. resa .lt. tolr * 10000)
     &               .and. resr .lt.  0.98 * resa) .or.
     &               (resr .lt. tolr * 10000 .and. su .lt. tols))

 555  continue

c....initialize working set
      do i = 1, neq
        jact(i) = 1
      enddo

      nact = neq
      tolr = 1.d-8
      tols = 1.d-10
      maxite = 500
      resa = 1.d25
      s_norm = 1.d25
      rstore = 1.d12 + 1.d0

c....determine inital value for
      dY = 0.0001d0 !0.01d0
      facts = 1.5d0    !3d0

      do ite = 1, maxite
        y(1) = dY*facts**(ite-1.d0) - dY + 1.d-16

        lsw = 1
        CALL FBPHAS_RES (lsw, res, xjac, Y, deps_res,
     &           xmu, compar, para, dtime, cm,
     &           depsp, tr_eps, nstrn, deveps_Tr, ifail, sigv,
     &           epsv_n, epsv, ltes, jfile, theta,
     &           resTr, res_norm, sig_tr, vn, phi_tr, sigvn,zm)

c....initialize working set
        ilevl = 1
        if (ite .eq. 1) then
            i_store(ilevl, 1) = 1
            res_store(ilevl, 1) = res(ilevl)
            y_store(ilevl, 1) = y(ilevl)
        else if (res(ilevl) .lt. 0d0) then
            i_store(ilevl, 2) = 1
            res_store(ilevl, 2) = res(ilevl)
            y_store(ilevl, 2) = y(ilevl)
            goto 301
        endif
        if (ite .ge. maxite) stop '**** FBPHAS_ite: 400 ite = maxite'
      enddo
  301 continue

c....mean starting value
      Y(ilevl) = (y_store(ilevl,1) + y_store(ilevl,2))/2d0

      if (res(ilevl) .lt. 0d0) then
        Y(ilevl) = y_store(ilevl,1)
      end if

c....unknowns y and search direction
      done = .false.
      CALL FBPHAS_ZEROS(s,neq)

      if (ltes .gt. 0)
     &   write(jfile, fmt = '(/,3x,''--- Iteration results ---- '',//,
     &               4x,''ite'','' ilin'',4x,''y'',10x,''res'',10x,
     &               ''res-str'',10x,
     &               ''nact'',2x,''ire '',2x,''nabh'')')

c....local iteration start
      do ite = 1, maxite
c....residual and jacobian
        lsw = 1
        CALL FBPHAS_RES(lsw, res, xjac, Y, deps_res,
     &           xmu, compar, para, dtime, cm, depsp,
     &           tr_eps, nstrn, deveps_Tr, ifail, sigv,
     &           epsv_n, epsv, ltes, jfile, theta, resTr, res_norm,
     &           sig_tr, vn, phi_tr, sigvn,zm)

c....write  results
c        write(jfile,fmt='(/,4x,'' Results '')')
        if (ltes .gt. 3) then
            write(jfile, fmt = '(/,4x,'' vector Y '')')
            write(jfile, fmt = '(2x,4(1x,e15.6))') (Y(j), j = 1, neq)
            write(jfile, fmt = '(4x,'' vector res '')')
            write(jfile, fmt = '(2x,4(1x,e15.6))') (res(j), j = 1, neq)
            write(jfile, fmt = '(4x,'' working set jact '')')
            write(jfile, fmt = '(2x,4(1x,i4))') (jact(j), j = 1, neq)
            write(jfile, fmt = '(/,4x,'' Jacobian '')')
            do i = 1, neq
                write(jfile, fmt = '(3x,i3,4(1x,e15.6))')
     &                   i, (xjac(i, j), j = 1, neq)
            enddo
        endif

        if (ltes .gt. 0)
     &           write(jfile, fmt = '(''**'',i3,2x,i3,2x,e12.5,
     &                   2x,e12.5,2x,e12.5,
     &                   2x,i3,3x,i3,
     &              3x,i3,3x,i3)') ite, ilin, Y(1), res(1), xjac(1, 1),
     &                   nact, ire, nabh

c....check for convergence
        if (tol_ch(res_norm, res_norm, abs(s(1)), tolr, tols)) then
            do i = 1, neq
                jact(i) = 1
            enddo
            ilevl = 1
            done = .true.
            goto 77
        else
c            resa2 = res_norm
c            resa1 = 10d4
        endif

c....update boundaries if necessary
        if (res(ilevl) .gt. 0d0) then
c     &  res(ilevl) .lt. res_store(ilevl,1)) then
            i_store(ilevl, 1) = 1
            res_store(ilevl, 1) = res(ilevl)
            y_store(ilevl, 1) = y(ilevl)
        elseif (res(ilevl) .lt. 0d0) then
c     &  res(ilevl) .gt. res_store(ilevl,2)) then
            i_store(ilevl, 2) = 1
            res_store(ilevl, 2) = res(ilevl)
            y_store(ilevl, 2) = y(ilevl)
        endif

c....solve for search direction
        CALL FBPHAS_GAUSS (xjac, res, s, jact, neq, 1, nabh,nact,s_norm)
        y(ilevl) = y(ilevl) - s(ilevl)

c....bisection
        if (y(ilevl) .lt. y_store(ilevl,1) .or.
     &           y(ilevl) .gt. y_store(ilevl,2)) then
            y(ilevl) = (y_store(ilevl,1) + y_store(ilevl,2))/2d0
        endif

        resa = res_norm
      if (ite .ge. maxite) then
        write(*,*) 'ite =',ite
        stop 'error: ite: maxite reached'
      endif
      enddo
c....end of local iteration

  77  continue

      if (.not. done)  then
        if (ifail .eq. 0) ifail = 1
        ites = 0
        if (ites .gt. 1) goto 555
        return
      endif

c....perform postiteration steps
      dlam = Y(1)

c....deviatoric stress component
      do i = 1,nstrn
        sig_dev(i) = sig_tr(i) - 2d0*xmu*3d0*dlam*vn(i)/2d0
      enddo
      sig_dev_norm = FBPHAS_norm(sig_dev,nstrn)

c....volumetric stress component
      CALL FBPHAS_ZEROS(sig_vol,nstrn)
      do i = 1,3
        sig_vol(i) = compar*tr_eps
      enddo

c....total stress
      do i = 1,nstrn
        sig(i) = sig_dev(i) + sig_vol(i)
      enddo

c....plastic strain increment
      do i = 1,nstrn
        depsp(i) = dlam*3d0*vn(i)/2d0
      enddo

c....equivalent von mises stress
      sigv = sqrt(3d0/2d0)*sig_dev_norm
      epsv = epsv_n + dlam

c....tangent stiffnes matrix
      lsw = 2
      CALL FBPHAS_RES (lsw, res, xjac, Y, deps_res,
     &       xmu, compar, para, dtime, cm, depsp,
     &       tr_eps, nstrn, deveps_Tr, ifail, sigv,
     &       epsv_n, epsv, ltes, jfile, theta, resTr, res_norm,
     &       sig_tr, vn, phi_tr, sigvn,zm)

c....calculate deps_Y from deps_res
      CALL FBPHAS_GAUSS (xjac, deps_res, deps_Y, jact, neq, nstrn+1,
     &                   nabh, nact, g0)

c....tangent
      sig_tr_norm = max(FBPHAS_norm(sig_tr,nstrn),1.d-12)
      CALL FBPHAS_ZEROS (cm,nstrn*nstrn)

      cm1 = 2d0*xmu
     &       -dlam*sqrt(2d0/3d0)*3d0*(2d0*xmu)**2/(2d0*sig_tr_norm)
      cm2 = dlam*sqrt(3d0/2d0)*3d0*(2d0*xmu)**2/(2d0*sig_tr_norm)
      cm3 = 3d0*2d0*xmu/2d0

      do i = 1,3
        cm(i,i) = cm(i,i) + cm1
      enddo
      do i = 4,nstrn
        cm(i,i) = cm(i,i) + cm1/2d0
      enddo
      do i = 1,3
        do j = 1,3
            cm(i,j) = cm(i,j) - cm1/3d0
        enddo
      enddo

      do i = 1,nstrn
        do j = 1,nstrn
            cm(i,j) = cm(i,j) + cm2*vn(i)*vn(j)
        enddo
      enddo

      do i = 1,nstrn
        do j = 1,nstrn
            cm(i,j) = cm(i,j) + cm3*vn(i)*deps_Y(1,j)
        enddo
      enddo

c....volumetric part of tangent moduli
      do i = 1,3
        do j = 1,3
            cm(i,j) = cm(i,j) + compar
        enddo
      enddo

c....test
c      ktes = 0
c      if (ktes .eq. 1) then
c        do i = 1,1
c            sig(i) =  vn(i)
c            do j = 1,nstrn
c                cm(i,j) = -deps_Y(i,j)
c            enddo
c        enddo
c        do i = 1,nstrn
c            sig(i) = vn(i)
c            do j = 1,nstrn
c               cm(i,j) = -deps_Y(i,j)
c            enddo
c        enddo
c      endif

      return
      end

c**********************************************************************
      SUBROUTINE FBPHAS_RES  (lsw   ,res   ,xjac    ,Y        ,deps_res,
     &                        xmu   ,compar,para    ,dtime    ,cm      ,
     &                        depsp ,tr_eps,nstrn   ,deveps_Tr,ifail   ,
     &                        sigv  ,epsv_n,epsv    ,ltes     ,jfile   ,
     &                        theta ,resTr ,res_norm,sig_tr   ,vn      ,
     &                        phi_tr,sigvn ,zm)
c***********************************************************************
c....compute residual and its derivatives
c***********************************************************************
      implicit double precision (a-h, o-z)
      parameter (neq = 1, nstrh = 6)
      dimension xjac(neq, neq), Y(neq), res(neq),
     &          deps_res(neq,nstrn+1)
      dimension para(*), cm(nstrn,nstrn), deveps_Tr(nstrn),
     &          depsp(nstrn), sig_tr(nstrn), vn(nstrh)
      double precision FBPHAS_norm

c....parameters
      iYtype = para(1)
      E0  = para(2)
      cE  = para(3)
      xnu = para(4)
      y0  = para(5)
      cy0 = para(6)
      harda = para(7)
      hardx = para(8)
      bm = para(9)
      qm = para(10)
      pKp = para(11)
      xm = para(12)

      if (iYtype .gt. 10)   iYtype = iYtype - (iYtype/10)*10
      if (iYtype .gt. 100)  iYtype = iYtype - (iYtype/100)*100

c....elastic material parameters
      Emod   = max(E0+cE*theta,400d0)
      compar = Emod/(3d0*(1d0-2d0*para(4)))
      xmu    = Emod/(2*(1+xnu))

c....plastic material parameters
      yt = y0 - cy0*theta
      hard = harda * exp(hardx*zm)

c....trial stress
      do i = 1, nstrn
        sig_tr(i) = 2d0*xmu*deveps_Tr(i)
      enddo
      sig_tr_norm = max(FBPHAS_norm(sig_tr,nstrn),1.d-12)

      if (iYtype .eq. 1) goto 888

c....unknown
      dlam = Y(1)

c....projection tensor
      do i = 1,nstrn
        vn(i) = sqrt(2d0/3d0)*sig_tr(i)/sig_tr_norm
      enddo

c....equivalent plastic strain
      epsv = epsv_n + dlam

c....yield condition
      phi = (sqrt(3d0/2d0)*sig_tr_norm-3d0*xmu*dlam) - yt
     &      - qm*(1d0-exp(-bm*epsv)) - hard*epsv
     &        - pKp*(dlam/dtime)**(1d0/xm)
c    write(jfile, fmt = '(4x,'' phi = '',12e12.4)') phi

c....residual
        res(1) = phi

c....norm of residual
      res_norm = 0d0
      do i = 1, neq
        res_norm = res_norm + res(i)*res(i)
      end do
      res_norm = dsqrt(res_norm)

c....predictor
      if (lsw .eq. 0) then
        phi_tr = phi
        resTr = res(1)
        epsv = epsv_n
        sigv = sqrt(3d0/2d0)*sig_tr_norm
        return
      endif

c....upper boundary for dlam
      if (lsw .eq. -5) then
        res(1) = phi
        epsv = epsv_nin
        sigv = sqrt(3d0/2d0)*sig_tr_norm
        return
      endif

c....derivative of residual w.r.t deltalambda (jacobian) (dr/dlam)
        xjac(1,1) = -3d0 * xmu  - bm*qm*exp(-bm*epsv) - hard
     &              -pKp/(xm*dtime)*(dlam/dtime)**(1d0/xm-1d0)

c....derivative of residual w.r.t eps (dr/deps)
      if (abs(lsw) .eq. 2 .or. lsw .eq. 6) then
        if (dlam .lt. 0d0) then
       write(jfile, fmt = '('' **** FBPHAS_res: dlam= '',e20.10)') dlam
        endif

        CALL FBPHAS_ZEROS(deps_res,(nstrn+1)*neq)

        do i = 1, nstrn
           deps_res(1,i) = sqrt(3d0/2d0)*2d0*xmu*sig_tr(i)/sig_tr_norm
        enddo
      endif

      return

  888 continue
      epsv = epsv_n
      sigv = sqrt(3d0/2d0)*sig_tr_norm

      end
c***********************************************************************
      SUBROUTINE FBPHAS_AUST (para, za, zan, dza, thetan,
     &                        theta, thetamax,zf ,zfn)
c***********************************************************************
c     compute austenite phase fraction based on koistinen-marbuger
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)

      pAc3 = para(21)
      pktheta_A = para(22)

      if (theta .GE. thetamax .AND. theta .GT. pAc3) then
        za = max(1d0-exp(-pktheta_A*(theta-pAc3)),0d0)
        dza = za-zan
      else
        dza = 0d0
        za = zan
      endif
      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_FERRIT (para, dtime, ZAF, zfn, zf, dzf, za, zan,
     &                          theta, thetan, zamax, dza)
c***********************************************************************
c     compute ferrite phase fraction based on koistinen-marbuger
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)

      pPs     = para(54)
      pFs     = para(53)
      pktheta = para(22)
      AC2     = para(31)

      cn = para(43)
      ck = para(44)
      t0 = para(45)
      isw = 2

      if (theta .gt. thetan) then
        zf = 1 - za
      elseif (theta .LE. pFs .and. theta .GT. pPs .and. za .NE. 0) then

       if (ZAF .EQ. 0d0) then
         ZAF = za
       endif

            zfn = max(zfn,1d-6)
            if (isw .eq. 1) then
                dzf =  dtime*(1d0-zfn)*ck*cn/t0*
     &                  (-dlog(1d0-zfn)/ck)**((cn-1)/cn) *ZAF
            elseif (isw .eq. 2) then
                dzf =  dtime*(ZAF-zfn)*ck*cn/t0*
     &                  (-dlog(1d0-zfn/ZAF)/ck)**((cn-1)/cn)
            endif
            zf =  zfn + dzf
            za = zan - dzf
      else
        zf = zfn
        dzf = 0d0
      endif

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_PEARLITE (para, dtime, ZAP, zpn, zp,
     &                            dzp, za, zan, theta, dza, zf, dzf ,
     &                            thetamax)
c***********************************************************************
c     compute pearlite phase fraction based on koistinen-marbuger
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)

      AC2 = para(31)
      pPs     = para(54)
      pktheta = para(22)

      cn = para(46)
      ck = para(47)
      t0 = para(48)

      theta_star = para(49)
      theta_m = para(50)

      isw = 2
      if (theta .LE. pPs .and. theta .GT. AC2 .and. za .NE. 0) then

       if (ZAP .EQ. 0d0) then
         ZAP = za
       endif
            zpn = max(zpn,1d-6)
            if (isw .eq. 1) then
                dzp =  dtime*(1d0-zpn)*ck*cn/t0*
     &                 (-dlog(1d0-zpn)/ck)**((cn-1)/cn)* ZAP
            elseif (isw .eq. 2) then

                dzp =  dtime*(ZAP-zpn)*ck*cn/t0*
     &                 (-dlog(max(1d0-zpn/ZAP,1d-6))/ck)**((cn-1)/cn)
     &                 *(theta_star-thetamax)/theta_m
            endif
            zp =  min(zpn + dzp,ZAP)
            za = max(zan - dzp,0d0)
        else
            zp = zpn
            dzp = 0d0
      endif

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_BAIN_IMP (para, ZAB, zb, zbn, dzb,
     &   za, zan, rink_n, rink, dtime,
     &   theta, epsv_n, sigvn, thetamax)
c***********************************************************************
c     compute bainite phase fraction with implicit integration scheme
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)

      pMs = para(20)
      alpha1 = para(25)
      beta1 = para(26)
      alpha2 = para(27)
      beta2 = para(28)
      A1 = para(29)
      A2 = para(30)
      AC2 = para(31)
      BU = para(32)
      BL = para(33)
      A3 = para(34)
      A4 = para(35)
      A5 = para(36)
      B1 = para(37)
      gama = para(38)
      pn = para(39)
      Gc = para(40)
      Q = para(41)

      theta_star = para(51)
      theta_m = para(52)

c....nucleation of bainite
      DT = AC2-theta
      if (theta .LE. AC2 .and. theta .GT. pMs .and. za .NE. 0
     &       .and. theta .lt. thetamax) then

        if (ZAB .EQ. 0d0) then
            ZAB = za
        endif

        r_str = A1*AC2*DT**(0.5)/(theta*DT*Q+A2*AC2*epsv_n*sigvn)
        Gibs = A3*DT**(1.5)*AC2**2d0/((DT*Q+A4*epsv_n*sigvn*AC2)**2)

        theta_str = alpha1*exp(beta1*epsv_n)
        A_rink = alpha2*exp(-beta2*epsv_n)

        absTemp = abs(theta - theta_str)
        if (theta .ge. theta_str) then
            rink = rink_n + dtime*A_rink*exp((-absTemp)/BU)
        else
            rink = rink_n + dtime*A_rink*exp((-absTemp)/BL)
        endif

        r_norm = (rink - r_str)/rink

        if (r_norm .GT. 0d0) then
            if (rink_n .LE. 0d0) then
                rink_n = rink
            endif

            r_norm_n = max((rink_n - r_str) / rink_n,0d0)
            fun = A5*exp((theta_str-theta)/B1)
     &          *(r_norm**pn + r_norm_n**pn)
     &          *exp(-Gibs/theta*Gc)*(gama-1d0)
            zb = (1d0 - (dtime*fun/2d0 + (1d0-zbn)
     &          **(1d0-gama))**(1d0/(1d0-gama)))

            zb = zb * (theta_star-thetamax)/theta_m
            zb =  min(max(zb*ZAB,0d0),ZAB)
            dzb = zb - zbn
            za = max(zan - dzb,0d0)
        else
            r_norm = 0d0
            zb = zbn
            dzb = 0d0
        endif

      else
        zb = zbn
        rink = rink_n
        r_norm = 0d0
        dzb = 0d0
      endif

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_MART (para, ZAM, zmn, zm, dzm, za, zan, theta)
c***********************************************************************
c     compute martensite phase fraction (based on koistinen marburger)
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)

      pMs     = para(20)
      pktheta = para(22)

      if (theta .LE. pMs .and. za .NE. 0) then

        if (ZAM .EQ. 0d0 ) then
            ZAM = za
        endif
c        write(*,*) ZAM

        zm =  max((1d0 - exp(-1.8*pktheta * (pMs - theta)))*ZAM,0d0)
        dzm = zm - zmn
        za = zan - dzm
        else
            zm = zmn
            dzm = 0d0
      endif

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_HIS1 (eps, his1, nhis, nstrn,
     &   neq, deveps_Tr, tr_eps, epsv_n,
     &   thetamax, theta, thetan, sigvn, zmn,
     &   zbn, rink_n, zan, zfn, zpn, zamax, dsigvn ,
     &   ZAF, ZAP, ZAB , ZAM)
c***********************************************************************
c....compute trial strains and extract history variables
c***********************************************************************
      implicit double precision (a-h, o-z)
      parameter (nstrh = 6)
      dimension eps(nstrn), deveps_Tr(nstrn), his1(nhis),
     &          epsn(nstrh), eps_tr(nstrh)

c....deviatoric trial strains
      do i = 1, nstrn
        eps_tr(i) = eps(i) - his1(i)
      enddo

      if (nstrn .eq. 1) then
        deveps_Tr(1) = eps_tr(1)
      else
        tr_eps = eps_tr(1) + eps_tr(2) + eps_tr(3)
c....deviatoric elastic trail strain
        do i = 1, 3
            deveps_Tr(i) = eps_tr(i) - tr_eps/3d0
        enddo
        do i = 4, nstrn
            deveps_Tr(i) = eps_tr(i)
        enddo
      endif

c....history variables
      epsv_n   = his1(nstrn + 1)
      sigvn    = his1(nstrn + 2)
      dsigvn   = his1(nstrn + 3)
      thetamax = his1(nstrn + 4)
      rink_n   = his1(nstrn + 5)
      zamax    = his1(nstrn + 6)
      zfn      = his1(nstrn + 7)
      zan      = his1(nstrn + 8)
      zbn      = his1(nstrn + 9)
      zmn      = his1(nstrn + 10)
      zpn      = his1(nstrn + 11)
      ZAF      = his1(nstrn + 12)
      ZAP      = his1(nstrn + 13)
      ZAB      = his1(nstrn + 14)
      ZAM      = his1(nstrn + 15)

      zamax = max(zamax,zan)
      thetamax = max(thetamax, theta, thetan)

      return
      end

c**********************************************************************
      SUBROUTINE FBPHAS_HIS2(ksw, his1, his2, nhis, nstrn,
     &   neq, depsp,  epsv, thetamax, sigv, zm, dsigv,
     &   zb, rink, za ,zf, zp, zamax , ZAF, ZAP, ZAB , ZAM)
c***********************************************************************
c....insert results into history arrray
c***********************************************************************
      implicit double precision (a-h,o-z)
      parameter (nstrh = 6)
      dimension his1(nhis), his2(nhis), depsp(nstrn)

      do i = 1, nstrn
        his2(i) = his1(i) + depsp(i)
c        his2(nstrn + i) = his1(nstrn + i)
      enddo

      his2(nstrn + 1) = epsv

c....only for direct problem
      if (ksw .gt. 0) return

      his2(nstrn + 2) = sigv
      his2(nstrn + 3) = dsigv
      his2(nstrn + 4) = thetamax
      his2(nstrn + 5) = rink
      his2(nstrn + 6) = zamax
      his2(nstrn + 7) = zf
      his2(nstrn + 8) = za
      his2(nstrn + 9) = zb
      his2(nstrn + 10) = zm
      his2(nstrn + 11) = zp
      his2(nstrn + 12) = ZAF
      his2(nstrn + 13) = ZAP
      his2(nstrn + 14) = ZAB
      his2(nstrn + 15) = ZAM

      return
      end
c***********************************************************************
      SUBROUTINE FBPHAS_IN (jsw ,imate ,para , npar, nhis,
     &                      incr,nhis3 ,iptyp, isym)
c***********************************************************************
c     read input data and prepare for feap routine
c***********************************************************************
      parameter (nparh = 57)
      implicit double precision (a-h, o-z)
      common /iofile/ ifile, jfile
      dimension para(*)

      if (jsw .eq. 1) then
        npar = nparh
        read (ifile,*) (para(i), i = 1,npar)
        itype = int(para(1))
        if (itype .gt. 10) itype = itype - (itype/10)*10
      endif

      incr = 0
      isym = 0

      if (iptyp .eq. 1) then
        nstrn = 1
      else if (iptyp .lt. 4) then
        nstrn = 4
      else
        nstrn = 6
      endif

      nhis  = nstrn + 16
      nhis3 = nhis

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_PP (npvpf5, tpvpf, iptyp, para)
c***********************************************************************
c     prepare for plot variables
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)
      character*10 tpvpf(100)

c....stress
      if (iptyp .eq. 1) then
        tpvpf(1) = 'SIG '
        tpvpf(2) = 'EPSV'
        npvpf5 = 2
        ns = 1
        return
      else if (iptyp .eq. 2) then
        tpvpf(1) = 'SXX'
        tpvpf(2) = 'SYY'
        tpvpf(3) = 'SZZ'
        tpvpf(4) = 'SXY'
        ns = 4
      else if (iptyp .eq. 3) then
        tpvpf(1) = 'SRR'
        tpvpf(2) = 'SZZ'
        tpvpf(3) = 'STT'
        tpvpf(4) = 'SRZ'
        ns = 4
       else if (iptyp .eq. 4) then
        tpvpf(1) = 'SXX'
        tpvpf(2) = 'SYY'
        tpvpf(3) = 'SZZ'
        tpvpf(4) = 'SXY'
        tpvpf(5) = 'SYZ'
        tpvpf(6) = 'SXZ'
        ns = 6
      else
        tpvpf(1) = 'SXX'
        tpvpf(2) = 'SYY'
        tpvpf(3) = 'SZZ'
        tpvpf(4) = 'SXY'
        ns = 4
      endif

        tpvpf(ns +  1) = 'TEMP'
        tpvpf(ns +  2) = 'SIGV'
        tpvpf(ns +  3) = 'EPSV'
        tpvpf(ns +  4) = 'ZF'
        tpvpf(ns +  5) = 'ZA'
        tpvpf(ns +  6) = 'ZB'
        tpvpf(ns +  7) = 'ZM'
        tpvpf(ns +  8) = 'ZP'

		npvpf5 = ns + 8

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_IP (plotv ,nplotv, para , sig  ,sigv  ,nstr  ,
     &                      epsv  ,zm ,theta,zb,za,zf, zp,
     &       ZAF, ZAP, ZAB , ZAM)
c***********************************************************************
c     insert plot variables
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension para(*)
      dimension sig(nstr), plotv(nplotv)

      do i = 1, nstr
        plotv(i) = sig(i)
      enddo

        plotv(nstr + 1) = theta-273.15
        plotv(nstr + 2) = sigv
        plotv(nstr + 3) = epsv
        plotv(nstr + 4) = zf
        plotv(nstr + 5) = za
        plotv(nstr + 6) = zb
        plotv(nstr + 7) = zm
        plotv(nstr + 8) = zp

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_GAUSS (a, c, x, jact,nres,nrs,nabh,nact,s_norm)
c***********************************************************************
c     Gauss equation solver for singular matrices a
c     (also negative definite a)
c***********************************************************************
      implicit double precision (a-h, o-z)
      parameter (nresh = 4)
      dimension a(nres, nres), c(nres, nrs), x(nres, nrs), jact(nres)
      dimension ip(nresh), is(nresh), ah(nresh, nresh), del_res(nresh)

      nact = 0
      if (nrs .eq. 1) then
        do i = 1, nres
            do j = 1, nres
               if (jact(i) .gt. 0 .and. jact(j) .gt. 0) ah(i,j) = a(i,j)
            enddo
        enddo
c        CALL FBPHAS_COPY(a,nres,nres,ah)
c        CALL FBPHAS_COPY(c,nres,1,ch)
      endif

      eps = 1.d-14

c....maxium diagonal element
      gam = -1.d30
      do  i = 1, nres
        if (jact(i) .gt. 0) gam = max(gam, abs(a(i, i)))
      enddo
      eps = gam * eps

      CALL FBPHAS_ZEROS (x, nres * nrs)
      do i = 1, nres
        ip(i) = i
        is(i) = i !!is not used here
      enddo

      nabh = 0

c....factorization
      nun = 0   !! independent columns
      mS = 0   !! max row/column element

      do 30 k = 1, nres
        kh = is(k)
        if (jact(kh) .gt. 0) then
c....pivoting
            pivo = 0d0
            mrow = 0
            diag = a(ip(k), kh)
            do i = k, nres
                if (jact(ip(i)) .gt. 0) then
                    if (abs(a(ip(i), kh)) .gt. pivo) then
                        pivo = abs(a(ip(i), kh))
                        mrow = i
                    endif
                endif
            enddo
            if (pivo .gt. eps)  then
                nun = nun + 1
                mS = k
            else
                jact(k) = 0
                nabh = nabh + 1
                goto 30
            endif
            if (abs(diag) .le. eps)  then
                ldum = ip(mrow)
                ip(mrow) = ip(k)
                ip(k) = ldum
            endif
            do i = k + 1, nres
                if (jact(ip(i)) .gt. 0) then
                    r = a(ip(i), kh) / a(ip(k), kh)
                    do irs = 1, nrs
                        c(ip(i), irs) = c(ip(i), irs) - r *c(ip(k),irs)
                    enddo
                    do j = k, nres
                        if (jact(is(j)) .gt. 0)  then
                            a(ip(i), is(j)) = a(ip(i), is(j)) -
     &                               r * a(ip(k), is(j))
                        endif
                    enddo
                endif
            enddo
        endif

c....nresh independent columns
        if (nun .eq. nresh) goto 35
   30     continue
  35   continue

c....delete dependent rows (not necessary here ?!)
      do j = k + 1, nres
        if (jact(is(j)) .gt. 0) then
            jact(is(j)) = 0
            nabh = nabh + 1
        endif
      enddo
c      nact = nact - nabh

c....backsubstitution
      do irs = 1, nrs
        do k = mS, 1, -1
            if (jact(is(k)) .gt. 0) then
                sum = 0d0
                do j = k + 1, mS
                    if (jact(is(j)) .gt. 0) then
                        sum = sum + a(ip(k), is(j)) * x(is(j), irs)
                    endif
                enddo
                x(is(k), irs) = (c(ip(k), irs) - sum) / a(ip(k), is(k))
            endif
        enddo
      enddo

c....check for increasing
      if (nrs .eq. 1)   then
        s_norm = 0d0
        do i = 1, nres
            if (jact(i) .gt. 0)  then
                s_norm = s_norm + x(i, 1) * x(i, 1)
                del_res(i) = 0d0
                do j = 1, nres
                if (jact(j).gt.0) del_res(i) = del_res(i)+ah(i,j)*x(j,1)
                enddo
            else if (jact(i) .eq. 0)  then
                jact(i) = 1
            endif
        enddo
        s_norm = sqrt(s_norm)
      endif
      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_COPY(a, m, n, r)
c***********************************************************************
c     copys a matrix r(m,n) = a(m,n)
c***********************************************************************
      implicit double precision (a-h, o-z)
      dimension a(m * n), r(m * n)

      mn = m * n
      do i = 1, mn
        r(i) = a(i)
      enddo

      return
      end

c***********************************************************************
      SUBROUTINE FBPHAS_ZEROS(v, nn)
c***********************************************************************
c     fills arrays with zeros
c***********************************************************************
      double precision v(nn)

      do n = 1, nn
        v(n) = 0d0
      enddo

      return
      end

c***********************************************************************
      FUNCTION FBPHAS_NORM(a, nstrn)
c***********************************************************************
c     calcutates norm of deviatoric tensor a
c***********************************************************************
      implicit double precision (a-h, o-z)
      double precision FBPHAS_NORM
      dimension a(nstrn)

      sum = 0d0
      do i = 1, nstrn
        sum = sum + a(i) * a(i)
      enddo

      if (nstrn .gt. 3) then
        sum = sum + a(4) * a(4)
      endif
      if (nstrn .gt. 4) then
        sum = sum + a(5) * a(5)
        sum = sum + a(6) * a(6)
      endif

      FBPHAS_NORM = sqrt(sum)

      return
      end
c***********************************************************************
