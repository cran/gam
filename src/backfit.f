C Output from Public domain Ratfor, version 1.0
      subroutine bakfit(x,npetc,y,w,which,spar,dof,match,nef, etal,s,eta
     *,beta,var,tol, qr,qraux,qpivot,effect,work)
      implicit double precision(a-h,o-z)
      logical ifvar
      integer npetc(7),iter
      integer n,p,q,which(1),match(1),nef(1),nit,maxit,qrank,qpivot(1)
      double precision x(1),y(1),w(1),spar(1),dof(1), etal(1),s(1),eta(1
     *),beta(1),var(1),tol, qr(1),qraux(1),effect(1),work(1)
      n=npetc(1)
      p=npetc(2)
      q=npetc(3)
      ifvar=.false.
      if(npetc(4).eq.1)then
      ifvar=.true.
      endif
      maxit=npetc(6)
      qrank=npetc(7)
      do23002 i=1,q
      work(i)=dof(i)
23002 continue
23003 continue
      call backf1(x,n,p,y,w,q,which,spar,dof,match,nef, etal,s,eta,beta,
     *var,ifvar,tol,nit,maxit, qr,qraux,qrank,qpivot,effect,work(q+1),wo
     *rk(q+n+1), work(q+2*n+1),work(q+3*n+1),work(q+4*n+1))
      npetc(7)=qrank
      return
      end
      subroutine backf1(x,n,p,y,w,q,which,spar,dof,match,nef, etal,s,eta
     *,beta,var,ifvar,tol,nit,maxit, qr,qraux,qrank,qpivot,effect,z,old,
     *sqwt,sqwti,work)
      implicit double precision(a-h,o-z)
      logical ifvar
      integer n,p,q,which(q),match(n,q),nef(q),nit,maxit,qrank,qpivot(p)
      double precision x(n,p),y(n),w(n),spar(q),dof(q), etal(n),s(n,q),e
     *ta(n),beta(p),var(n,q),tol, qr(n,p),qraux(p),effect(n),work(1)
      double precision z(1),old(1),dwrss,ratio
      double precision sqwt(n),sqwti(n)
      logical anyzwt
      double precision deltaf, normf,onedm7
      integer job,info
      onedm7=1d-7
      job=1101
      info=1
      if(q.eq.0)then
      maxit=1
      endif
      ratio=1d0
      anyzwt=.false.
      do23006 i=1,n
      if(w(i).gt.0d0)then
      sqwt(i)=dsqrt(w(i))
      sqwti(i)=1d0/sqwt(i)
      else
      sqwt(i)=0d0
      sqwti(i)=0d0
      anyzwt=.true.
      endif
23006 continue
23007 continue
      if(qrank.eq.0)then
      do23012 i=1,n
      do23014 j=1,p
      qr(i,j)=x(i,j)*sqwt(i)
23014 continue
23015 continue
23012 continue
23013 continue
      do23016 j=1,p
      qpivot(j)=j
23016 continue
23017 continue
      call dqrdca(qr,n,n,p,qraux,qpivot,work,qrank,onedm7)
      endif
      do23018 i=1,n
      eta(i)=0d0
      j=1
23020 if(.not.(j.le.q))goto 23022
      eta(i)=eta(i)+s(i,j)
23021 j=j+1
      goto 23020
23022 continue
23018 continue
23019 continue
      nit=0
23023 if((ratio .gt. tol ).and.(nit .lt. maxit))then
      deltaf=0d0
      nit=nit+1
      do23025 i=1,n
      z(i)=(y(i)-eta(i))*sqwt(i)
      old(i)=etal(i)
23025 continue
23026 continue
      call dqrsl(qr,n,n,qrank,qraux,z,work(1),effect(1),beta, work(1),et
     *al,job,info)
      do23027 i=1,n
      etal(i)=etal(i)*sqwti(i)
23027 continue
23028 continue
      k=1
23029 if(.not.(k.le.q))goto 23031
      j=which(k)
      do23032 i=1,n
      old(i)=s(i,k)
      z(i)=y(i)-etal(i)-eta(i)+old(i)
23032 continue
23033 continue
      if(nit.gt.1)then
      dof(k)=0d0
      endif
      call splsm(x(1,j),z,w,n,match(1,k),nef(k),spar(k), dof(k),s(1,k),s
     *0,var(1,k),ifvar,work)
      do23036 i=1,n
      eta(i)=eta(i)+s(i,k)-old(i)
      etal(i)=etal(i)+s0
23036 continue
23037 continue
      deltaf=deltaf+dwrss(n,old,s(1,k),w)
23030 k=k+1
      goto 23029
23031 continue
      normf=0d0
      do23038 i=1,n
      normf=normf+w(i)*eta(i)*eta(i)
23038 continue
23039 continue
      if(normf.gt.0d0)then
      ratio=dsqrt(deltaf/normf)
      else
      ratio = 0d0
      endif
      goto 23023
      endif
23024 continue
      do23042 j=1,p 
      work(j)=beta(j)
23042 continue
23043 continue
      do23044 j=1,p 
      beta(qpivot(j))=work(j)
23044 continue
23045 continue
      if(anyzwt)then
      do23048 i=1,n 
      if(w(i) .le. 0d0)then
      etal(i)=0d0
      do23052 j=1,p
      etal(i)=etal(i)+beta(j)*x(i,j)
23052 continue
23053 continue
      endif
23048 continue
23049 continue
      endif
      do23054 i=1,n
      eta(i)=eta(i)+etal(i)
23054 continue
23055 continue
      do23056 j=1,q 
      call unpck(n,nef(j),match(1,j),var(1,j),old)
      do23058 i=1,n 
      var(i,j)=old(i)
23058 continue
23059 continue
23056 continue
23057 continue
      return
      end
