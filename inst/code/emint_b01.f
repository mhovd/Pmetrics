      subroutine emint(psi,ldpsi,theta,ldtheta,npoint,nsub,ijob,
     &                 x,dx,y,dy,fobj,gap,nvar,keep,IHESS,ERRFIL,
     &                 isupres)
      use npag_utils, only: pause, verifyval
      implicit real*8 (a-h,o-z)
      integer ldpsi,ldtheta,npoint,nsub,ijob,nvar,IHESS
      dimension psi(:,:),theta(:,:),x(:),dx(:),y(:),dy(:)
      CHARACTER(len=20) ERRFIL
      integer isupres
      real*8 mu, s
      parameter (MAXSUBem=999,MAXACTem=10000000)
      dimension w(MAXSUBem),dw(MAXSUBem),Ptx(MAXSUBem),
     &          hess(MAXSUBem,2*MAXSUBem)
      dimension psisum(MAXSUBem),XVERIFY(100)
      integer kpvt(MAXSUBem), ipivot(MAXACTem), list(MAXACTem)
       keep = 0
       s = 0.d0
      write (*,*) "Enterred emint(); ijob =",ijob
      if(nsub.gt.MAXSUBem) then
      write(6,*) 'nsub =',nsub, ' is greater than MAXSUBem=',MAXSUBem
      write(6,*) 'MAXSUBem needs to be reset as large as nsub'
      write(6,*) 'in PARAMETER statement in subroutine emint'
        OPEN(42,FILE=ERRFIL)
      write(42,*) 'nsub =',nsub, ' is greater than MAXSUBem=',MAXSUBem
      write(42,*) 'MAXSUBem needs to be reset as large as nsub'
      write(42,*) 'in PARAMETER statement in subroutine emint'
        CLOSE(42)
      CALL PAUSE
      stop
      endif
      if(npoint.gt.MAXACTem) then
      write(6,*) 'npoint=',npoint,' is larger than MAXACTem=',MAXACTem
      write(6,*) 'MAXACTem needs to be reset as large as npoint'
      write(6,*) 'in PARAMETER statement in subroutine emint'
        OPEN(42,FILE=ERRFIL)
      write(42,*) 'npoint=',npoint,' is larger than MAXACTem=',MAXACTem
      write(42,*) 'MAXACTem needs to be reset as large as npoint'
      write(42,*) 'in PARAMETER statement in subroutine emint'
        CLOSE(42)
      CALL PAUSE
      stop
      endif
      psimin=0.
      do j=1,nsub
      do i=1,npoint
      if(psi(j,i).le.psimin) psimin=psi(j,i)
      enddo
      enddo
      if(psimin.lt.0) then
        write(6,*) 'Psi matrix not non-negative -stop'
        OPEN(42,FILE=ERRFIL)
        write(42,*) 'Psi matrix not non-negative -stop'
        CLOSE(42)
        CALL PAUSE
        stop
      endif
      colsummin=1.e10
      do j=1,nsub
        s=0.
        do i=1,npoint
           x(i)=1.d0
           s=s+psi(j,i)
        enddo
        psisum(j) = s
        Ptx(j)=s
        if(s.le.colsummin) colsummin=s
        if(s.le.0) then
           write(6,*) 'psi has a zero row -stop'
        OPEN(42,FILE=ERRFIL)
        write(42,*) 'psi has a zero row -stop'
        CLOSE(42)
        CALL PAUSE
        stop
        endif
        w(j)=1./s
      enddo
      shrink=0.
      do i=1,npoint
        sum=0.d0
        do j=1,nsub
           sum=sum+psi(j,i)*w(j)
        enddo
        y(i)=sum
        if(sum.gt.shrink) shrink=sum
      enddo
      shrink=2.d0*shrink
      if(s.le.0) then
        write(6,*) 'Psi has a zero column -stop'
        OPEN(42,FILE=ERRFIL)
        write(42,*) 'Psi has a zero column -stop'
        CLOSE(42)
        CALL PAUSE
        stop
      endif
      eps=1.d-10
      sig=0.d0
      mu=0.d0
      do i=1,npoint
        x(i)=1.d0*shrink
        y(i)=y(i)/shrink
        y(i)=1.d0-y(i)
        mu=mu+x(i)*y(i)
      enddo
      mu=mu/npoint
      rmax = -1.e38
      do j=1,nsub
        w(j)=w(j)/shrink
        Ptx(j)=Ptx(j)*shrink
        if(dabs(1.d0-w(j)*Ptx(j)).ge.rmax) rmax =
     &  dabs(1.d0-w(j)*Ptx(j))
      enddo
      gap=1.d0
      iter=0
100   continue
      conval = mu
      if(conval .lt. rmax) conval = rmax
      if(conval .lt. gap) conval = gap
      convcrit = eps/conval
      XVERIFY(1) = convcrit
      CALL VERIFYVAL(1,XVERIFY)
       WRITE(*,123) iter,XVERIFY(1)
  123 FORMAT(' Iteration ',I9,' CONV. CRIT = ',G15.2,' (1 OR HIGHER FOR
     1CONVERGENCE)')
      if(mu.le.eps.and.rmax.le.eps.and.gap.le.eps) go to 9000
      iter=iter+1
        ILOOP = ILOOP + 1
      tbuilda=0
      smu=sig*mu
      do j=1,nsub
        do k=1,nsub
           hess(j,k)=0.
        enddo
      enddo
      do i=1,npoint
        scale=x(i)/y(i)
       do j=1,nsub
         fact=scale*psi(j,i)
         do k=j,nsub
           hess(k,j)=hess(k,j)+fact*psi(k,i)
         enddo
       enddo
      enddo
      do j=1,nsub-1
      do k=j+1,nsub
      hess(j,k)=hess(k,j)
      enddo
      enddo
      do j=1,nsub
        hess(j,j)=hess(j,j)+Ptx(j)/w(j)
      enddo
      tbuildb=0
      tbuild=tbuildb-tbuilda
      IF(ISUPRES .EQ. 0) write(6,*) 'tbuild=',tbuild
      call dsifa(hess,MAXSUBem,nsub,kpvt,info)
      tbuildc=0
      tfactor=tbuildc-tbuildb
      IF(ISUPRES .EQ. 0) write(6,*) 'tfactor=',tfactor
      write(6,*) 'gap,info,tfactor=',gap,info,tfactor
      IHESS = 0
      if(info .ne. 0) then
       IHESS = -1
       WRITE(25,163)
       WRITE(*,163)
  163  FORMAT(//' Hessian matrix in interior point EM algorithm'/
     1' is singular.  Possibly number of grid points is too small,'/
     2' or assay coefficients are too large. '//
     3' Try again with a new assay polynomial or larger grid.'//
     4' Suggested quick fix: rerun and select error model 2)'/
     5' in response to the initial question; then enter a'/
     6' initial value gamma = 10.0 in response to the prompt for'/
     7' that value.'//
     8' THE PROGRAM WILL CREATE OUTPUT FILES BEFORE STOPPING. '//)
      WRITE(25,164) info
      WRITE(*,164) info
  164 FORMAT(//' NOTE THAT IN SUBROUTINE emint, THE VALUE OF INFO'/
     1' IS ',i6,//
     2' IF THIS VALUE IS POSTIVE, IT IS LIKELY THE NO. OF THE SUBJECT'/
     3' (OR AT LEAST THE FIRST SUBJECT) WHICH CAUSED THE HESSIAN '/
     4' ERROR. SO IN THIS CASE, YOU MIGHT ALSO WANT TO EXAMINE THE'/
     5' DATA IN THIS SUBJECT TO VERIFY THEY ARE CORRECT.'//)
        OPEN(42,FILE=ERRFIL)
         WRITE(42,163)
         WRITE(42,164) info
        CLOSE(42)
      return
      endif
      do j=1,nsub
        sum=0.d0
        do i=1,npoint
          sum=sum+psi(j,i)*smu/y(i)
        enddo
        dw(j)=1.d0/w(j)-sum
      enddo
       call dsisl(hess,MAXSUBem,nsub,kpvt,dw)
      do i=1,npoint
        sum=0.
        do j=1,nsub
          sum=sum+psi(j,i)*dw(j)
        enddo
        dy(i)=-sum
        dx(i)=smu/y(i)-x(i)-dy(i)*x(i)/y(i)
      enddo
      alfpri=-.5
      do i=1,npoint
        if(dx(i)/x(i).le.alfpri) alfpri=dx(i)/x(i)
      enddo
      alfpri=-1.d0/alfpri
      alfpri=min(1.d0,0.99995*alfpri)
      alfdual=-0.5d0
      do i=1,npoint
        if(dy(i)/y(i).le.alfdual) alfdual=dy(i)/y(i)
      enddo
      alfdual=-1.d0/alfdual
      alfdual=min(1.d0,0.99995*alfdual)
      mu=0.d0
      do i=1,npoint
        x(i)=x(i)+alfpri*dx(i)
        y(i)=y(i)+alfdual*dy(i)
        mu=mu+x(i)*y(i)
      enddo
      mu=mu/npoint
      do j=1,nsub
        sum=0.d0
        do i=1,npoint
          sum=sum+psi(j,i)*x(i)
        enddo
        Ptx(j)=sum
      enddo
      do j=1,nsub
        w(j)=w(j)+alfdual*dw(j)
      enddo
      rmax=0.
      do j=1,nsub
        rtest=1.d0-w(j)*Ptx(j)
        if(dabs(rtest).gt.rmax) rmax=dabs(rtest)
      enddo
      sumlogw=0.d0
      sumlgPtx=0.d0
      do j=1,nsub
        sumlogw=sumlogw+dlog(w(j))
        sumlgPtx=sumlgPtx+dlog(Ptx(j))
      enddo
      gap = dabs(sumlogw+sumlgPtx)/(1.d0+dabs(sumlgPtx))
      if(mu.lt.eps.and.rmax.gt.eps) then
        sig=1.d0
      else
        c2=1.d2
        term1=(1.d0-alfpri)**2
        term2=(1.d0-alfdual)**2
        term3=(rmax-mu)/(rmax+c2*mu)
        term=max(term1,term2)
        term=max(term,term3)
        sig=min(0.3d0,term)
      endif
      sumx=0.d0
      do i=1,npoint
        sumx=sumx+x(i)
      enddo
      fobj=0.
      do j=1,nsub
        fobj=fobj+dlog(Ptx(j)/sumx)
      enddo
      go to 100
9000  continue
      sumx=0.
      do i=1,npoint
        sumx=sumx+x(i)
      enddo
      do i=1,npoint
        x(i)=x(i)/sumx
      enddo
      if(ijob.eq.0) return
      isum=0
      xlim=0.
      do i=1,npoint
        if(x(i).gt.xlim) xlim=x(i)
      enddo
      xlim=xlim*1.d-3
      isum = 0
      do i=1,npoint
        if(x(i).gt.xlim) then
          isum = isum + 1
          list(isum) = i
          do j=1,nsub
            psi(j,isum) = psi(j,i)
          enddo
          do j=1,nvar
            theta(isum,j)=theta(i,j)
          enddo
          x(isum)=x(i)
        endif
      enddo
      job=1
      do k=1,npoint
        ipivot(k)=0
      enddo
      do i=1,isum
        do j=1,nsub
          psi(j,i+isum)=psi(j,i)
        enddo
      enddo
      do i=1,isum
        do j=1,nsub
          psi(j,i) = psi(j,i)/psisum(j)
        enddo
      enddo
      call dqrdc(psi,ldpsi,nsub,isum,y,ipivot,dy,job)
      keep = 0
      limloop = nsub
      if(isum.lt.nsub) then
        limloop = isum
      else
        write (*,*) "in emint(ijob=1,...) isum > nsub",isum,nsub
      end if
      do i=1,limloop
        test=dnrm2(i,psi(1,i:),1)
        if(dabs(psi(i,i)/test).ge.1.d-8) keep=keep+1
      enddo
      if(isum.gt.1) then
      do i=1,keep-1
      do j=i,keep
      if(ipivot(i)*ipivot(j).ne.0.and.ipivot(i).gt.ipivot(j)) then
         itemp=ipivot(i)
         ipivot(i)=ipivot(j)
         ipivot(j)=itemp
      endif
      enddo
      enddo
      endif
      do i=1,isum
        do j=1,nsub
          psi(j,i)=psi(j,i+isum)
        enddo
      enddo
      do k=1,npoint
        dx(k)=0
      enddo
      sumkeep = 0.
      do k=1,keep
        j=ipivot(k)
        if(j.ne.0) then
          do jj=1,nsub
            psi(jj,k)=psi(jj,j+isum)
          enddo
          do jvar=1,nvar
            theta(k,jvar) = theta(j,jvar)
          enddo
        endif
        if(j.gt.0) dx(list(j))=1.
        if(j.gt.0) sumkeep = sumkeep + x(list(j))
        if(j.gt.0) w(k)=x(list(j))
      enddo
      write (*,*)
      write (*,*) "(isum,keep,sum(x(list))=", isum, keep, sumkeep
      write (*,*)
      return
      end
      subroutine dpoco(a,lda,n,rcond,z,info)
      integer lda,n,info
      double precision a(lda,1),z(1)
      double precision rcond
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer i,j,jm1,k,kb,kp1
      do 30 j = 1, n
         z(j) = dasum(j,a(1,j),1)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = z(i) + dabs(a(i,j))
   10    continue
   20    continue
   30 continue
      anorm = 0.0d0
      do 40 j = 1, n
         anorm = dmax1(anorm,z(j))
   40 continue
      call dpofa(a,lda,n,info)
      if (info .ne. 0) go to 180
         ek = 1.0d0
         do 50 j = 1, n
            z(j) = 0.0d0
   50    continue
         do 110 k = 1, n
            if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
            if (dabs(ek-z(k)) .le. a(k,k)) go to 60
               s = a(k,k)/dabs(ek-z(k))
               call dscal(n,s,z,1)
               ek = s*ek
   60       continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = dabs(wk)
            sm = dabs(wkm)
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
            kp1 = k + 1
            if (kp1 .gt. n) go to 100
               do 70 j = kp1, n
                  sm = sm + dabs(z(j)+wkm*a(k,j))
                  z(j) = z(j) + wk*a(k,j)
                  s = s + dabs(z(j))
   70          continue
               if (s .ge. sm) go to 90
                  t = wkm - wk
                  wk = wkm
                  do 80 j = kp1, n
                     z(j) = z(j) + t*a(k,j)
   80             continue
   90          continue
  100       continue
            z(k) = wk
  110    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         do 130 kb = 1, n
            k = n + 1 - kb
            if (dabs(z(k)) .le. a(k,k)) go to 120
               s = a(k,k)/dabs(z(k))
               call dscal(n,s,z,1)
  120       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call daxpy(k-1,t,a(1,k),1,z(1),1)
  130    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = 1.0d0
         do 150 k = 1, n
            z(k) = z(k) - ddot(k-1,a(1,k),1,z(1),1)
            if (dabs(z(k)) .le. a(k,k)) go to 140
               s = a(k,k)/dabs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  140       continue
            z(k) = z(k)/a(k,k)
  150    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = s*ynorm
         do 170 kb = 1, n
            k = n + 1 - kb
            if (dabs(z(k)) .le. a(k,k)) go to 160
               s = a(k,k)/dabs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  160       continue
            z(k) = z(k)/a(k,k)
            t = -z(k)
            call daxpy(k-1,t,a(1,k),1,z(1),1)
  170    continue
         s = 1.0d0/dasum(n,z,1)
         call dscal(n,s,z,1)
         ynorm = s*ynorm
         if (anorm .ne. 0.0d0) rcond = ynorm/anorm
         if (anorm .eq. 0.0d0) rcond = 0.0d0
  180 continue
      return
      end
      subroutine dpofa(a,lda,n,info)
      integer lda,n,info
      double precision a(lda,1)
      double precision ddot,t
      double precision s
      integer j,jm1,k
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
            if (s .le. 0.0d0) go to 40
            a(j,j) = dsqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
      subroutine dposl(a,lda,n,b)
      integer lda,n
      double precision a(lda,1),b(1)
      double precision ddot,t
      integer k,kb
      do 10 k = 1, n
         t = ddot(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k) - t)/a(k,k)
   10 continue
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call daxpy(k-1,t,a(1,k),1,b(1),1)
   20 continue
      return
      end
      subroutine dsifa(a,lda,n,kpvt,info)
      integer lda,n,kpvt(1),info
      double precision a(lda,1)
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
      logical swap
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
      info = 0
      k = n
   10 continue
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0d0) info = 1
            go to 200
   20    continue
         km1 = k - 1
         absakk = dabs(a(k,k))
         imax = idamax(k-1,a(1,k),1)
         colmax = dabs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50       continue
            if (dabs(a(imax,imax)) .lt. alpha*rowmax) go to 60
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
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
            if (.not.swap) go to 120
               call dswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call daxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
            if (.not.swap) go to 160
               call dswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
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
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      subroutine dsisl(a,lda,n,kpvt,b)
      integer lda,n,kpvt(1)
      double precision a(lda,1),b(1)
      double precision ak,akm1,bk,bkm1,ddot,denom,temp
      integer k,kp
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
               call daxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
               call daxpy(k-2,b(k),a(1,k),1,b(1),1)
               call daxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
            if (k .eq. 1) go to 110
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
            if (k .eq. 1) go to 140
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      integer ldx,n,p,job
      integer jpvt(1)
      double precision x(ldx,1),qraux(1),work(1)
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dnrm2,tt
      double precision ddot,nrmxl,t
      logical negj,swapj
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0d0
         if (l .eq. n) go to 190
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
      SUBROUTINE DPOTRF( UPLO, N, A, LDA, INFO, ERRFIL )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      CHARACTER*20       ERRFIL
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      LOGICAL            UPPER
      INTEGER            J, JB, NB
      LOGICAL            LSAME
      EXTERNAL           LSAME
      EXTERNAL           DGEMM, DPOTF2, DSYRK, DTRSM, XERBLA
      INTRINSIC          MAX, MIN
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRF', -INFO, ERRFIL )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      nb = 16
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
         CALL DPOTF2( UPLO, N, A, LDA, INFO, ERRFIL )
      ELSE
         IF( UPPER ) THEN
            DO 10 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Upper', JB, A( J, J ), LDA, INFO, ERRFIL )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
                  CALL DGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
         ELSE
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
               CALL DSYRK( 'Lower', 'No transpose', JB, J-1, -ONE,
     $                     A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL DPOTF2( 'Lower', JB, A( J, J ), LDA, INFO, ERRFIL )
               IF( INFO.NE.0 )
     $            GO TO 30
               IF( J+JB.LE.N ) THEN
                  CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1, JB,
     $                        J-1, -ONE, A( J+JB, 1 ), LDA, A( J, 1 ),
     $                        LDA, ONE, A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'Transpose', 'Non-unit',
     $                        N-J-JB+1, JB, ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
   30 CONTINUE
      INFO = INFO + J - 1
   40 CONTINUE
      RETURN
      END
      SUBROUTINE DPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO, ERRFIL )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
      CHARACTER*20       ERRFIL
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      LOGICAL            UPPER
      LOGICAL            LSAME
      EXTERNAL           LSAME
      EXTERNAL           DTRSM, XERBLA
      INTRINSIC          MAX
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTRS', -INFO, ERRFIL )
         RETURN
      END IF
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
      IF( UPPER ) THEN
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
      END IF
      RETURN
      END
      SUBROUTINE DPOTF2( UPLO, N, A, LDA, INFO, ERRFIL )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   A( LDA, * )
      CHARACTER*20       ERRFIL
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
      LOGICAL            UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
      LOGICAL            LSAME
      DOUBLE PRECISION   DDOT
      EXTERNAL           LSAME, DDOT
      EXTERNAL           DGEMV, DSCAL, XERBLA
      INTRINSIC          MAX, SQRT
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DPOTF2', -INFO, ERRFIL )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( UPPER ) THEN
         DO 10 J = 1, N
            AJJ = A( J, J ) - DDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
            IF( J.LT.N ) THEN
               CALL DGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
               CALL DSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
         DO 20 J = 1, N
            AJJ = A( J, J ) - DDOT( J-1, A( J, 1 ), LDA, A( J, 1 ),
     $            LDA )
            IF( AJJ.LE.ZERO ) THEN
               A( J, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ
            IF( J.LT.N ) THEN
               CALL DGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ),
     $                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
               CALL DSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
      GO TO 40
   30 CONTINUE
      INFO = J
   40 CONTINUE
      RETURN
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
      CHARACTER          CA, CB
      INTRINSIC          ICHAR
      INTEGER            INTA, INTB, ZCODE
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
      ZCODE = ICHAR( 'Z' )
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
      END
      SUBROUTINE XERBLA( SRNAME, INFO, ERRFIL )
      use npag_utils, only : pause
      CHARACTER*6        SRNAME
      CHARACTER*20       ERRFIL
      INTEGER            INFO
      WRITE( *, FMT = 9999 )SRNAME, INFO
        OPEN(42,FILE=ERRFIL)
      WRITE( 42, FMT = 9999 )SRNAME, INFO
        CLOSE(42)
      CALL PAUSE
      STOP
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
      END
