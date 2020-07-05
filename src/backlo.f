C Output from Public domain Ratfor, version 1.0
      subroutine baklo(x,y,w,npetc,wddnfl,spatol,match, etal,s,eta,beta,
     *var,dof, qr,qraux,qpivot,effect,iv,v,iwork,work)
      implicit double precision(a-h,o-z)
      integer n,p,q,nit,maxit,qrank
      integer npetc(7),wddnfl(*),match(*),qpivot(*),iv(*),iwork(*)
      double precision x(*),y(*),w(*),spatol(*), etal(*),s(*),eta(*),bet
     *a(*),var(*),dof(*), qr(*),qraux(*),v(*),effect(*),work(*)
      n=npetc(1)
      p=npetc(2)
      q=npetc(3)
      maxit=npetc(5)
      qrank=npetc(6)
      call baklo0(x,n,p,y,w,q,wddnfl(1),wddnfl(q+1),wddnfl(2*q+1), spato
     *l(1),wddnfl(3*q+1),dof,match,wddnfl(4*q+1), etal,s,eta,beta,var,sp
     *atol(q+1), nit,maxit,qr,qraux,qrank,qpivot,effect, work(1),work(n+
     *1),work(2*n+1),work(3*n+1), iv,wddnfl(5*q+1),wddnfl(6*q+1),v,wddnf
     *l(7*q+1), iwork(1),work(4*n+1))
      npetc(4)=nit
      npetc(6)=qrank
      return
      end
      subroutine baklo0(x,n,p,y,w,q,which,dwhich,pwhich,span,degree,dof,
     *match,nef, etal,s,eta,beta,var,tol,nit,maxit, qr,qraux,qrank,qpivo
     *t,effect,z,old,sqwt,sqwti, iv,liv,lv,v,nvmax,iwork,work)
      implicit double precision(a-h,o-z)
      integer n,p,q,which(q),dwhich(q),pwhich(q),degree(q),match(n,q),ne
     *f(q),nit, maxit,qrank,qpivot(p),iv(*),liv(q),lv(q),nvmax(q),iwork(
     *q)
      double precision x(n,p),y(n),w(n),span(q),dof(q), etal(n),s(n,q),e
     *ta(n),beta(p),var(n,q),tol, qr(n,p),qraux(p),v(*),effect(n),work(*
     *)
      double precision z(*),old(*),dwrss,ratio
      double precision sqwt(n),sqwti(n)
      logical anyzwt
      double precision deltaf, normf,onedm7
      integer job,info,slv,sliv,iw,j,dj,pj
      onedm7=1d-7
      job=1101
      info=1
      if(q.eq.0)then
      maxit=1
      endif
      ratio=1d0
      anyzwt=.false.
      do23002 i=1,n
      if(w(i).gt.0d0)then
      sqwt(i)=dsqrt(w(i))
      sqwti(i)=1d0/sqwt(i)
      else
      sqwt(i)=0d0
      sqwti(i)=0d0
      anyzwt=.true.
      endif
23002 continue
      continue
      if(qrank.eq.0)then
      do23008 i=1,n
      do23010 j=1,p
      qr(i,j)=x(i,j)*sqwt(i)
23010 continue
      continue
23008 continue
      continue
      do23012 j=1,p
      qpivot(j)=j
23012 continue
      continue
      call dqrdca(qr,n,n,p,qraux,qpivot,work,qrank,onedm7)
      endif
      do23014 i=1,n
      eta(i)=0d0
      j=1
23016 if(.not.(j.le.q))goto 23018
      eta(i)=eta(i)+s(i,j)
      j=j+1
      goto 23016
23018 continue
23014 continue
      continue
      nit=0
23019 if((ratio .gt. tol ).and.(nit .lt. maxit))then
      deltaf=0d0
      nit=nit+1
      do23021 i=1,n
      z(i)=(y(i)-eta(i))*sqwt(i)
      old(i)=etal(i)
23021 continue
      continue
      call dqrsl(qr,n,n,qrank,qraux,z,work(1),effect(1),beta, work(1),et
     *al,job,info)
      do23023 i=1,n
      etal(i)=etal(i)*sqwti(i)
23023 continue
      continue
      sliv=1
      slv=1
      iw=5*n+1
      k=1
23025 if(.not.(k.le.q))goto 23027
      j=which(k)
      dj=dwhich(k)
      pj=pwhich(k)
      do23028 i=1,n
      old(i)=s(i,k)
      z(i)=y(i)-etal(i)-eta(i)+old(i)
23028 continue
      continue
C Trevor edited this 06/28/2020      
      call lo1(x(1,j),z,w,n,dj,pj,nvmax(k),span(k),degree(k),match(1,k),
     * nef(k),nit,dof(k),s(1,k),var(1,k),work(iw), work(iw+pj+1),work(iw
     *+nef(k)*dj+pj+1), work(iw+nef(k)*(dj+1)+pj+2),work(iw + nef(k)*(dj
     *+2)+pj+2), work(iw+nef(k)*(dj+3)+pj+2),work(iw+nef(k)*(pj+dj+4)+pj
     *+2), iwork(1),work(iw+nef(k)*(pj+dj+4)+4+2*pj), iv(sliv),liv(k),lv
     *(k),v(slv), work(1) )
      sliv=sliv+liv(k)
      slv=slv+lv(k)
      iw=iw+nef(k)*(pj+dj+4)+5+3*pj
      do23030 i=1,n
      eta(i)=eta(i)+s(i,k)-old(i)
23030 continue
      continue
      deltaf=deltaf+dwrss(n,old,s(1,k),w)
      k=k+1
      goto 23025
23027 continue
      normf=0d0
      do23032 i=1,n
      normf=normf+w(i)*eta(i)*eta(i)
23032 continue
      continue
      if(normf.gt.0d0)then
      ratio=dsqrt(deltaf/normf)
      else
      ratio = 0d0
      endif
      goto 23019
      endif
      continue
      do23036 j=1,p 
      work(j)=beta(j)
23036 continue
      continue
      do23038 j=1,p 
      beta(qpivot(j))=work(j)
23038 continue
      continue
      if(anyzwt)then
      do23042 i=1,n 
      if(w(i) .le. 0d0)then
      etal(i)=0d0
      do23046 j=1,p
      etal(i)=etal(i)+beta(j)*x(i,j)
23046 continue
      continue
      endif
23042 continue
      continue
      endif
      do23048 i=1,n
      eta(i)=eta(i)+etal(i)
23048 continue
      continue
      return
      end
