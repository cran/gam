C Output from Public domain Ratfor, version 1.0
      subroutine sknotl(x,n,knot,k)
      implicit double precision(a-h,o-z)
      double precision x(n),knot(n+6),a1,a2,a3,a4
      integer n,k,ndk,j
      a1 = log(50d0)/log(2d0) 
      a2 = log(100d0)/log(2d0)
      a3 = log(140d0)/log(2d0) 
      a4 = log(200d0)/log(2d0)
      if(n.lt.50)then
      ndk = n 
      else
      if(n.ge.50 .and. n.lt.200)then
      ndk = 2.**(a1+(a2-a1)*(n-50.)/150.) 
      else
      if(n.ge.200 .and. n.lt.800)then
      ndk = 2.**(a2+(a3-a2)*(n-200.)/600.) 
      else
      if(n.ge.800 .and. n.lt.3200)then
      ndk = 2.**(a3+(a4-a3)*(n-800.)/2400.) 
      else
      if(n.ge.3200)then
      ndk = 200. + float(n-3200)**.2 
      endif
      endif
      endif
      endif
      endif
      k = ndk + 6
      do23010 j=1,3 
      knot(j) = x(1) 
23010 continue
23011 continue
      do23012 j=1,ndk 
      knot(j+3) = x( 1 + (j-1)*(n-1)/(ndk-1) ) 
23012 continue
23013 continue
      do23014 j=1,3 
      knot(ndk+3+j) = x(n) 
23014 continue
23015 continue
      return
      end
      subroutine splsm(x,y,w,n,match,nef,spar,dof,smo,s0,cov,ifcov,work)
      implicit double precision(a-h,o-z)
      double precision x(*),y(*),w(*),spar,dof,smo(*),s0,cov(*),work(*)
      integer n,match(*),nef
      integer ifcov
      call splsm1(x,y,w,n,match,nef,spar,dof,smo,s0,cov,ifcov, work(1), 
     *work(nef+2),work(2*nef+3),work(3*nef+4), work(3*nef+n+10))
      return
      end
      subroutine splsm1(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov, xin
     *,yin,win,knot, work)
      implicit double precision(a-h,o-z)
      double precision x(*),y(*),w(*),spar,dof,smo(*),s0,lev(*),work(*)
      integer n,match(*),nef
      integer ifcov
      double precision xin(nef+1),yin(nef+1),win(nef+1),knot(nef+6)
      integer nk,ldnk,ld4,k
      double precision xmin,xrange
      call suff(n,nef,match,x,y,w,xin,yin,win,work(1))
      xmin=xin(1)
      xrange=xin(nef)-xin(1)
      do23016 i=1,nef 
      xin(i)=(xin(i)-xmin)/xrange
23016 continue
23017 continue
      call sknotl(xin,nef,knot,k)
      nk=k-4
      ld4=4
      ldnk=1
      call splsm2(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov, xin,yin,w
     *in,knot, work(1), work(nk+1), work(nk+nef+2),work(nk+2*nef+3), wor
     *k(2*nk+2*nef+3),work(3*nk+2*nef+3),work(4*nk+2*nef+3), work(5*nk+2
     **nef+3), work(6*nk+2*nef+3),work(7*nk+2*nef+3),work(8*nk+2*nef+3),
     * work(9*nk+2*nef+3), work(10*nk+2*nef+3),work((10+ld4)*nk+2*nef+3)
     *, work((10+2*ld4)*nk+2*nef+3), ld4,ldnk,nk)
      return
      end
      subroutine splsm2(x,y,w,n,match,nef,spar,dof,smo,s0,lev,ifcov, xin
     *,yin,win,knot, coef,sout,levout,xwy, hs0,hs1,hs2,hs3, sg0,sg1,sg2,
     *sg3, abd,p1ip,p2ip,ld4,ldnk,nk)
      implicit double precision(a-h,o-z)
      double precision x(*),y(*),w(*),spar,dof,smo(*),s0,lev(*)
      integer n,match(*),nef
      integer nk,ldnk,ld4
      integer ifcov
      double precision xin(nef+1),yin(nef+1),win(nef+1),knot(nk+4)
      double precision coef(nk),sout(nef+1),levout(nef+1),xwy(nk), hs0(n
     *k),hs1(nk),hs2(nk),hs3(nk), sg0(nk),sg1(nk),sg2(nk),sg3(nk), abd(l
     *d4,nk),p1ip(ld4,nk),p2ip(ldnk,*)
      integer ispar,icrit,isetup,ier
      double precision lspar,uspar,tol,penalt, sumwin,dofoff,crit,xbar,d
     *sum,xsbar
      double precision yssw, eps
      integer maxit
      double precision wmean
      crit=0d0
      if(dof .le. 0d0)then
      ispar=1
      icrit=3
      dofoff=0d0
      else
      if( dof .lt. 1d0 )then
      dof=1d0
      endif
      ispar=0
      icrit=3
      dofoff=dof+1d0
      endif
      isetup=0
      ier=1
      penalt=1d0
      lspar= -1.5
      uspar= 2.0
      tol=1d-4
      eps=2d-8
      maxit=200
      do23022 i=1,nef
      sout(i)=yin(i)*yin(i)
23022 continue
23023 continue
      sumwin=0d0
      do23024 i=1,nef
      sumwin=sumwin+win(i)
23024 continue
23025 continue
      yssw=wmean(nef,sout,win)
      s0=wmean(n,y,w)
      yssw=yssw*(sumwin-s0*s0)
      call sbart(penalt,dofoff,xin,yin,win,yssw,nef,knot,nk, coef,sout,l
     *evout,crit, icrit,spar,ispar,maxit, lspar,uspar,tol,eps, isetup, x
     *wy, hs0,hs1,hs2,hs3, sg0,sg1,sg2,sg3, abd,p1ip,p2ip,ld4,ldnk,ier)
      do23026 i=1,nef 
      win(i)=win(i)*win(i)
23026 continue
23027 continue
      sbar=wmean(nef,sout,win)
      xbar=wmean(nef,xin,win)
      do23028 i=1,nef 
      lev(i)=(xin(i)-xbar)*sout(i) 
23028 continue
23029 continue
      xsbar=wmean(nef,lev,win)
      do23030 i=1,nef 
      lev(i)=(xin(i)-xbar)**2 
23030 continue
23031 continue
      dsum=wmean(nef,lev,win)
      do23032 i=1,nef 
      if(win(i).gt.0d0)then
      lev(i)=levout(i)/win(i)-1d0/sumwin -lev(i)/(sumwin*dsum)
      else
      lev(i)=0d0
      endif
23032 continue
23033 continue
      dof=0d0
      do23036 i=1,nef 
      dof=dof+lev(i)*win(i)
23036 continue
23037 continue
      dof=dof+1d0
      do23038 i=1,nef
      sout(i)=sout(i)-sbar -(xin(i)-xbar)*xsbar/dsum
23038 continue
23039 continue
      call unpck(n,nef,match,sout,smo)
      return
      end
      double precision function wmean(n,y,w)
      integer n
      double precision y(n),w(n),wtot,wsum
      wtot=0d0
      wsum=0d0
      do23040 i=1,n
      wsum=wsum+y(i)*w(i)
      wtot=wtot+w(i)
23040 continue
23041 continue
      if(wtot .gt. 0d0)then
      wmean=wsum/wtot
      else
      wmean=0d0
      endif
      return
      end
