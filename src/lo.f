C Output from Public domain Ratfor, version 1.0
      subroutine lo0(x,y,w,n,d,p,nvmax,span,degree,match,nef,dof,s,var, 
     *beta,iv,liv,lv,v,iwork,work)
      integer n,d,p,nvmax,degree,match(1),nef,liv,lv,iv(liv),iwork(1)
      double precision x(n,d),y(n),w(n),span,dof,s(n),var(n),v(lv),work(
     *1)
      double precision beta(p+1)
      integer qrank
      call lo1(x,y,w,n,d,p,nvmax,span,degree,match,nef,0,dof,s,var,beta,
     * work(1),work(nef*d+1),work(nef*(d+1)+2),work(nef*(d+2)+2), work(n
     *ef*(d+3)+2),qrank,iwork(1),work(nef*(p+d+4)+3+p), iv,liv,lv,v, wor
     *k(nef*(p+d+4)+4+2*p) )
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
      call pck(n,nef,match,w,win)
      do23002 i=1,nef
      if(win(i).gt.0d0)then
      sqwin(i)=dsqrt(win(i))
      sqwini(i)=1d0/sqwin(i)
      else
      sqwin(i)=1d-5
      sqwini(i)=1d5
      endif
23002 continue
23003 continue
      do23006 i=1,n
      k=match(i)
      if(k.le.nef)then
      do23010 j=1,d
      xin(k,j)=x(i,j)
23010 continue
23011 continue
      j=d+1
23012 if(.not.(j.le.p))goto 23014
      xqr(k,j+1)=x(i,j)
23013 j=j+1
      goto 23012
23014 continue
      endif
23006 continue
23007 continue
      do23015 i=1,nef
      xqr(i,1)=sqwin(i)
      do23017 j=1,d
      xqr(i,j+1)=xin(i,j)*sqwin(i)
23017 continue
23018 continue
      j=d+2
23019 if(.not.(j.le.p+1))goto 23021
      xqr(i,j)=xqr(i,j)*sqwin(i)
23020 j=j+1
      goto 23019
23021 continue
23015 continue
23016 continue
      j=1
23022 if(.not.(j.le.p+1))goto 23024
      qpivot(j)=j
23023 j=j+1
      goto 23022
23024 continue
      call dqrdca(xqr,nef,nef,p+1,qraux,qpivot,work,qrank,onedm7)
      setlf = (nit.eq.1)
      call lowesd(106,iv,liv,lv,v,d,nef,span,degree,nvmax,setlf)
      v(2)=span/5d0
      endif
      do23025 i=1,n
      work(i)=y(i)*w(i)
23025 continue
23026 continue
      call pck(n,nef,match,work,yin)
      do23027 i=1,nef
      yin(i)=yin(i)*sqwini(i)*sqwini(i)
23027 continue
23028 continue
      if(nit.le.1)then
      call lowesb(xin,yin,win,levout,ifvar,iv,liv,lv,v)
      else
      call lowesr(yin,iv,liv,lv,v)
      endif
      call lowese(iv,liv,lv,v,nef,xin,sout)
      do23031 i=1,nef
      sout(i)=sout(i)*sqwin(i)
23031 continue
23032 continue
      call dqrsl(xqr,nef,nef,qrank,qraux,sout,work(1),work(1),beta, sout
     *,work(1),job,info)
      do23033 i=1,nef
      sout(i)=sout(i)*sqwini(i)
23033 continue
23034 continue
      if(nit.le.1)then
      job=10000
      j=1
23037 if(.not.(j.le.p+1))goto 23039
      do23040 i=1,nef
      work(i)=0d0
23040 continue
23041 continue
      work(j)=1d0
      call dqrsl(xqr,nef,nef,qrank,qraux,work,var,junk,junk, junk,junk,j
     *ob,info)
      do23042 i=1,nef
      levout(i)=levout(i) - var(i)**2
23042 continue
23043 continue
23038 j=j+1
      goto 23037
23039 continue
      dof=0d0
      do23044 i=1,nef 
      if(win(i).gt.0d0)then
      levout(i)=levout(i)/win(i)
      else
      levout(i)=0d0
      endif
23044 continue
23045 continue
      do23048 i=1,nef 
      dof=dof+levout(i)*win(i)
23048 continue
23049 continue
      call unpck(n,nef,match,levout,var)
      j=1
23050 if(.not.(j.le.p+1))goto 23052
      work(j)=beta(j)
23051 j=j+1
      goto 23050
23052 continue
      j=1
23053 if(.not.(j.le.p+1))goto 23055
      beta(qpivot(j))=work(j)
23054 j=j+1
      goto 23053
23055 continue
      endif
      call unpck(n,nef,match,sout,s)
      return
      end
      subroutine pck(n,p,match,x,xbar)
      integer match(n),p,n
      double precision x(n),xbar(n)
      do23056 i=1,p
      xbar(i)=0d0
23056 continue
23057 continue
      do23058 i=1,n
      xbar(match(i))=xbar(match(i))+x(i)
23058 continue
23059 continue
      return
      end
      subroutine suff(n,p,match,x,y,w,xbar,ybar,wbar,work)
      integer match(n),p,n
      double precision x(n),xbar(n),y(n),ybar(n),w(n),wbar(n),work(n)
      call pck(n,p,match,w,wbar)
      do23060 i=1,n
      xbar(match(i))=x(i)
23060 continue
23061 continue
      do23062 i=1,n
      work(i)=y(i)*w(i)
23062 continue
23063 continue
      call pck(n,p,match,work,ybar)
      do23064 i=1,p
      if(wbar(i).gt.0d0)then
      ybar(i)=ybar(i)/wbar(i)
      else
      ybar(i)=0d0
      endif
23064 continue
23065 continue
      return
      end
      subroutine unpck(n,p,match,xbar,x)
      integer match(n),p,n
      double precision x(n),xbar(p+1)
      if(p.lt.n)then
      xbar(p+1)=0d0
      endif
      do23070 i = 1,n
      x(i)=xbar(match(i))
23070 continue
23071 continue
      return
      end
      double precision function dwrss(n,y,eta,w)
      integer n
      double precision y(n),w(n),wtot,wsum,work,eta(n)
      wsum=0d0
      wtot=0d0
      do23072 i = 1,n
      work=y(i)-eta(i)
      wsum=wsum+w(i)*work*work
      wtot=wtot+w(i)
23072 continue
23073 continue
      if(wtot .gt. 0d0)then
      dwrss=wsum/wtot
      else
      dwrss=0d0
      endif
      return
      end
