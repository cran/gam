C Output from Public domain Ratfor, version 1.0
      subroutine baklo(bchat,x,y,w,npetc,wddnfl,spatol,match, etal,s,eta
     *,beta,var,dof, qr,qraux,qpivot,iv,v,work)
      implicit double precision(a-h,o-z)
      integer n,p,q,nit,maxit,qrank
      integer npetc(7),wddnfl(1),match(1),qpivot(1),iv(1)
      double precision x(1),y(1),w(1),spatol(1), etal(1),s(1),eta(1),bet
     *a(1),var(1),dof(1), qr(1),qraux(1),v(1),work(1)
      logical bchat,onemor
      onemor=.true.
      n=npetc(1)
      p=npetc(2)
      q=npetc(3)
      maxit=npetc(5)
      qrank=npetc(6)
      call baklo0(x,n,p,y,w,q,wddnfl(1),wddnfl(q+1),wddnfl(2*q+1), spato
     *l(1),wddnfl(3*q+1),dof,match,wddnfl(4*q+1), etal,s,eta,beta,var,sp
     *atol(q+1), nit,maxit,qr,qraux,qrank,qpivot, work(1),work(n+1),work
     *(2*n+1),work(3*n+1), iv,wddnfl(5*q+1),wddnfl(6*q+1),v,wddnfl(7*q+1
     *), work(4*n+1))
      npetc(4)=nit
      npetc(6)=qrank
      return
      end
      subroutine baklo0(x,n,p,y,w,q,which,dwhich,pwhich,span,degree,dof,
     *match,nef, etal,s,eta,beta,var,tol,nit,maxit, qr,qraux,qrank,qpivo
     *t,z,old,sqwt,sqwti, iv,liv,lv,v,nvmax,work)
      implicit double precision(a-h,o-z)
      integer n,p,q,which(q),dwhich(q),pwhich(q),degree(q),match(n,q),ne
     *f(q),nit, maxit,qrank,qpivot(p),iv(1),liv(q),lv(q),nvmax(q)
      double precision x(n,p),y(n),w(n),span(q),dof(q), etal(n),s(n,q),e
     *ta(n),beta(p),var(n,q),tol, qr(n,p),qraux(p),v(1),work(1)
      double precision z(1),old(1),dwrss,ratio
      double precision sqwt(n),sqwti(n)
      logical anyzwt
      double precision deltaf, normf,onedm7
      integer job,info,slv,sliv,iw,j,dj,pj
      onedm7=1d-7
      job=101
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
23003 continue
      if(qrank.eq.0)then
      do23008 i=1,n
      do23010 j=1,p
      qr(i,j)=x(i,j)*sqwt(i)
23010 continue
23011 continue
23008 continue
23009 continue
      do23012 j=1,p
      qpivot(j)=j
23012 continue
23013 continue
      call dqrdca(qr,n,n,p,qraux,qpivot,work,qrank,onedm7)
      endif
      do23014 i=1,n
      eta(i)=0d0
      j=1
23016 if(.not.(j.le.q))goto 23018
      eta(i)=eta(i)+s(i,j)
23017 j=j+1
      goto 23016
23018 continue
23014 continue
23015 continue
      nit=0
23019 if((ratio .gt. tol ).and.(nit .lt. maxit))then
      deltaf=0d0
      nit=nit+1
      do23021 i=1,n
      z(i)=(y(i)-eta(i))*sqwt(i)
      old(i)=etal(i)
23021 continue
23022 continue
      call dqrsl(qr,n,n,qrank,qraux,z,work(1),work(1),beta, work(1),etal
     *,job,info)
      do23023 i=1,n
      etal(i)=etal(i)*sqwti(i)
23023 continue
23024 continue
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
23029 continue
      call lo1(x(1,j),z,w,n,dj,pj,nvmax(k),span(k),degree(k),match(1,k),
     * nef(k),nit,dof(k),s(1,k),var(1,k),work(iw), work(iw+pj+1),work(iw
     *+nef(k)*dj+pj+1), work(iw+nef(k)*(dj+1)+pj+2),work(iw + nef(k)*(dj
     *+2)+pj+2), work(iw+nef(k)*(dj+3)+pj+2),work(iw+nef(k)*(pj+dj+4)+pj
     *+2), work(iw+nef(k)*(pj+dj+4)+pj+3),work(iw+nef(k)*(pj+dj+4)+4+2*p
     *j), iv(sliv),liv(k),lv(k),v(slv), work(1) )
      sliv=sliv+liv(k)
      slv=slv+lv(k)
      iw=iw+nef(k)*(pj+dj+4)+5+3*pj
      do23030 i=1,n
      eta(i)=eta(i)+s(i,k)-old(i)
23030 continue
23031 continue
      deltaf=deltaf+dwrss(n,old,s(1,k),w)
23026 k=k+1
      goto 23025
23027 continue
      normf=0d0
      do23032 i=1,n
      normf=normf+w(i)*eta(i)*eta(i)
23032 continue
23033 continue
      if(normf.gt.0d0)then
      ratio=dsqrt(deltaf/normf)
      else
      ratio = 0d0
      endif
      goto 23019
      endif
23020 continue
      do23036 j=1,p 
      work(j)=beta(j)
23036 continue
23037 continue
      do23038 j=1,p 
      beta(qpivot(j))=work(j)
23038 continue
23039 continue
      if(anyzwt)then
      do23042 i=1,n 
      if(w(i) .le. 0d0)then
      etal(i)=0d0
      do23046 j=1,p
      etal(i)=etal(i)+beta(j)*x(i,j)
23046 continue
23047 continue
      endif
23042 continue
23043 continue
      endif
      do23048 i=1,n
      eta(i)=eta(i)+etal(i)
23048 continue
23049 continue
      return
      end
      subroutine lo0(x,y,w,n,d,p,nvmax,span,degree,match,nef,dof,s,var, 
     *beta,iv,liv,lv,v,work)
      integer n,d,p,nvmax,degree,match(1),nef,liv,lv,iv(liv)
      double precision x(n,d),y(n),w(n),span,dof,s(n),var(n),v(lv),work(
     *1)
      double precision beta(p+1)
      integer qrank
      call lo1(x,y,w,n,d,p,nvmax,span,degree,match,nef,0,dof,s,var,beta,
     * work(1),work(nef*d+1),work(nef*(d+1)+2),work(nef*(d+2)+2), work(n
     *ef*(d+3)+2),qrank,work(nef*(p+d+4)+2),work(nef*(p+d+4)+3+p), iv,li
     *v,lv,v, work(nef*(p+d+4)+4+2*p) )
      return
      end
      subroutine lo1(x,y,w,n,d,p,nvmax,span,degree,match,nef,nit,dof,s,v
     *ar,beta, xin,win,sqwin,sqwini,xqr,qrank,qpivot,qraux, iv,liv,lv,v,
     * work)
      integer n,d,p,nvmax,degree,match(1),nef,nit,qrank,qpivot(p+1)
      integer iv(liv),liv,lv
      double precision x(n,d),y(n),w(n),span,dof,s(n),var(n),beta(p+1), 
     *xin(nef,d),win(nef+1),sqwin(nef),sqwini(nef),xqr(nef,p+1), qraux(p
     *+1),v(lv), work(1)
      call lo2(x,y,w,n,d,p,nvmax,span,degree,match,nef,nit,dof,s,var,bet
     *a, xin,win,sqwin,sqwini,xqr,qrank,qpivot,qraux, iv,liv,lv,v, work(
     *1),work(nef+2),work(2*nef+3),work(3*nef+4))
      return
      end
      subroutine lo2(x,y,w,n,d,p,nvmax,span,degree,match,nef,nit,dof,s,v
     *ar,beta, xin,win,sqwin,sqwini,xqr,qrank,qpivot,qraux, iv,liv,lv,v,
     * levout,sout,yin,work)
      integer n,d,p,nvmax,degree,match(1),nef,nit,qrank,qpivot(p+1)
      integer iv(liv),liv,lv
      double precision x(n,d),y(n),w(n),span,dof,s(n),var(n),beta(p+1), 
     *xin(nef,d),win(nef+1),sqwin(nef),sqwini(nef),xqr(nef,p+1), qraux(p
     *+1),v(lv), levout(nef+1), sout(nef+1),yin(nef+1),work(1)
      double precision junk, onedm7
      integer job, info
      logical setlf, ifvar
      job=110
      info=1
      ifvar=.true.
      onedm7=1d-7
      if(nit.le.1)then
      call pack(n,nef,match,w,win)
      do23052 i=1,nef
      if(win(i).gt.0d0)then
      sqwin(i)=dsqrt(win(i))
      sqwini(i)=1d0/sqwin(i)
      else
      sqwin(i)=1d-5
      sqwini(i)=1d5
      endif
23052 continue
23053 continue
      do23056 i=1,n
      k=match(i)
      if(k.le.nef)then
      do23060 j=1,d
      xin(k,j)=x(i,j)
23060 continue
23061 continue
      j=d+1
23062 if(.not.(j.le.p))goto 23064
      xqr(k,j+1)=x(i,j)
23063 j=j+1
      goto 23062
23064 continue
      endif
23056 continue
23057 continue
      do23065 i=1,nef
      xqr(i,1)=sqwin(i)
      do23067 j=1,d
      xqr(i,j+1)=xin(i,j)*sqwin(i)
23067 continue
23068 continue
      j=d+2
23069 if(.not.(j.le.p+1))goto 23071
      xqr(i,j)=xqr(i,j)*sqwin(i)
23070 j=j+1
      goto 23069
23071 continue
23065 continue
23066 continue
      j=1
23072 if(.not.(j.le.p+1))goto 23074
      qpivot(j)=j
23073 j=j+1
      goto 23072
23074 continue
      call dqrdca(xqr,nef,nef,p+1,qraux,qpivot,work,qrank,onedm7)
      setlf = (nit.eq.1)
      call lowesd(106,iv,liv,lv,v,d,nef,span,degree,nvmax,setlf)
      v(2)=span/5d0
      endif
      do23075 i=1,n
      work(i)=y(i)*w(i)
23075 continue
23076 continue
      call pack(n,nef,match,work,yin)
      do23077 i=1,nef
      yin(i)=yin(i)*sqwini(i)*sqwini(i)
23077 continue
23078 continue
      if(nit.le.1)then
      call lowesb(xin,yin,win,levout,ifvar,iv,liv,lv,v)
      else
      call lowesr(yin,iv,liv,lv,v)
      endif
      call lowese(iv,liv,lv,v,nef,xin,sout)
      do23081 i=1,nef
      sout(i)=sout(i)*sqwin(i)
23081 continue
23082 continue
      call dqrsl(xqr,nef,nef,qrank,qraux,sout,work(1),work(1),beta, sout
     *,work(1),job,info)
      do23083 i=1,nef
      sout(i)=sout(i)*sqwini(i)
23083 continue
23084 continue
      if(nit.le.1)then
      job=10000
      j=1
23087 if(.not.(j.le.p+1))goto 23089
      do23090 i=1,nef
      work(i)=0d0
23090 continue
23091 continue
      work(j)=1d0
      call dqrsl(xqr,nef,nef,qrank,qraux,work,var,junk,junk, junk,junk,j
     *ob,info)
      do23092 i=1,nef
      levout(i)=levout(i) - var(i)**2
23092 continue
23093 continue
23088 j=j+1
      goto 23087
23089 continue
      dof=0d0
      do23094 i=1,nef 
      if(win(i).gt.0d0)then
      levout(i)=levout(i)/win(i)
      else
      levout(i)=0d0
      endif
23094 continue
23095 continue
      do23098 i=1,nef 
      dof=dof+levout(i)*win(i)
23098 continue
23099 continue
      call unpack(n,nef,match,levout,var)
      j=1
23100 if(.not.(j.le.p+1))goto 23102
      work(j)=beta(j)
23101 j=j+1
      goto 23100
23102 continue
      j=1
23103 if(.not.(j.le.p+1))goto 23105
      beta(qpivot(j))=work(j)
23104 j=j+1
      goto 23103
23105 continue
      endif
      call unpack(n,nef,match,sout,s)
      return
      end
      subroutine pack(n,p,match,x,xbar)
      integer match(n),p,n
      double precision x(n),xbar(n)
      do23106 i=1,p
      xbar(i)=0d0
23106 continue
23107 continue
      do23108 i=1,n
      xbar(match(i))=xbar(match(i))+x(i)
23108 continue
23109 continue
      return
      end
      subroutine suff(n,p,match,x,y,w,xbar,ybar,wbar,work)
      integer match(n),p,n
      double precision x(n),xbar(n),y(n),ybar(n),w(n),wbar(n),work(n)
      call pack(n,p,match,w,wbar)
      do23110 i=1,n
      xbar(match(i))=x(i)
23110 continue
23111 continue
      do23112 i=1,n
      work(i)=y(i)*w(i)
23112 continue
23113 continue
      call pack(n,p,match,work,ybar)
      do23114 i=1,p
      if(wbar(i).gt.0d0)then
      ybar(i)=ybar(i)/wbar(i)
      else
      ybar(i)=0d0
      endif
23114 continue
23115 continue
      return
      end
      subroutine unpack(n,p,match,xbar,x)
      integer match(n),p,n
      double precision x(n),xbar(p+1)
      if(p.lt.n)then
      xbar(p+1)=0d0
      endif
      do23120 i = 1,n
      x(i)=xbar(match(i))
23120 continue
23121 continue
      return
      end
      double precision function dwrss(n,y,eta,w)
      integer n
      double precision y(n),w(n),wtot,wsum,work,eta(n)
      wsum=0d0
      wtot=0d0
      do23122 i = 1,n
      work=y(i)-eta(i)
      wsum=wsum+w(i)*work*work
      wtot=wtot+w(i)
23122 continue
23123 continue
      if(wtot .gt. 0d0)then
      dwrss=wsum/wtot
      else
      dwrss=0d0
      endif
      return
      end
