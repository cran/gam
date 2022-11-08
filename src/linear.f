C Output from Public domain Ratfor, version 1.0
      subroutine dqrls(x,dx,pivot,qraux,y,dy,beta,res,qt,tol,scrtch,rank
     *)
      integer pivot(*),dx(2),dy(2),rank
      double precision x(*), qraux(*), y(*), beta(*),res(*),qt(*),tol(*)
     *, scrtch(*)
      integer n,p,q,kn,kp,k,info
      n=dx(1)
      p=dx(2)
      q=dy(2)
      call dqrdca(x,n,n,p,qraux,pivot,scrtch,rank,tol(1))
      kn=1
      kp=1
      if(rank.gt.0)then
      k=1
23002 if(.not.(k.le.q))goto 23004
      call dqrsl(x,n,n,rank,qraux,y(kn),scrtch,qt(kn),beta(kp), res(kn),
     *scrtch,00110,info)
      kn = kn+n
      kp=kp+p
      k=k+1
      goto 23002
23004 continue
      endif
      return
      end
      subroutine dqrsl1(qr,dq,qra,rank,y,k,qy,qb,job,info)
      double precision qr(*),qra(*),y(*),qy(*),qb(*)
      integer dq(2),job,k,rank
      integer n,kn,kb,j
      double precision ourqty(1), ourqy(1), ourb(1), ourrsd(1), ourxb(1)
      ourqty(1) = 0d0
      ourqy(1) = 0d0
      ourb(1) = 0d0
      ourrsd(1) = 0d0
      ourxb(1) = 0d0
      n = dq(1)
      kn = 1
      kb = 1
      I23005=(job)
      goto 23005
23007 continue
      j=0
23008 if(.not.(j.lt.k))goto 23010
      call dqrsl(qr,dq(1),dq(1),rank,qra,y(kn),qy(kn),ourqty,ourb,ourrsd
     *,ourxb,job,info)
      kn = kn +n
      j = j+1
      goto 23008
23010 continue
      goto 23006
23011 continue
      j=0
23012 if(.not.(j.lt.k))goto 23014
      call dqrsl(qr,dq(1),dq(1),rank,qra,y(kn),ourqy,qy(kn),ourb,ourrsd,
     *ourxb,job,info)
      kn = kn +n
      j = j+1
      goto 23012
23014 continue
      goto 23006
23015 continue
      j=0
23016 if(.not.(j.lt.k))goto 23018
      call dqrsl(qr,dq(1),dq(1),rank,qra,y(kn),ourqy,qy(kn),qb(kb),ourrs
     *d,ourxb,job,info)
      kn = kn +n
      kb = kb +rank
      j = j+1
      goto 23016
23018 continue
      goto 23006
23019 continue
      j=0
23020 if(.not.(j.lt.k))goto 23022
      call dqrsl(qr,dq(1),dq(1),rank,qra,y(kn),ourqy,qy(kn),ourb,qb(kn),
     *ourxb,job,info)
      kn = kn +n
      j = j+1
      goto 23020
23022 continue
      goto 23006
23023 continue
      j=0
23024 if(.not.(j.lt.k))goto 23026
      call dqrsl(qr,dq(1),dq(1),rank,qra,y(kn),ourqy,qy(kn),ourb,ourrsd,
     *qb(kn),job,info)
      kn = kn +n
      j = j+1
      goto 23024
23026 continue
      goto 23006
23027 continue
      info = -1
      goto 23006
23005 continue
      if (I23005.eq.1)goto 23023
      if (I23005.eq.10)goto 23019
      if (I23005.eq.100)goto 23015
      if (I23005.eq.1000)goto 23011
      if (I23005.eq.10000)goto 23007
      goto 23027
23006 continue
      return
      end
      subroutine dqr(x,dx,pivot,qraux,tol,scrtch,rank)
      integer pivot(*),dx(2),rank
      double precision x(*), qraux(*), tol(*), scrtch(*)
      integer n,p
      n=dx(1)
      p=dx(2)
      call dqrdca(x,n,n,p,qraux,pivot,scrtch,rank,tol(1))
      return
      end
      subroutine dqrdca(x,ldx,n,p,qraux,jpvt,work,rank,eps)
      integer ldx,n,p,rank
      integer jpvt(*)
      double precision x(ldx,*),qraux(*),work(*),eps
      integer j,jj,jp,l,lup,curpvt
      double precision dnrm2,tt
      double precision ddot,nrmxl,t,ww
      do23028 j=1,p 
      qraux(j) = dnrm2(n,x(1,j),1)
      work(j) = qraux(j)
      work(j+p) = qraux(j)
23028 continue
      continue
      l=1
      lup = min0(n,p)
      curpvt = p
23030 if(l.le.lup)then
      qraux(l) = 0.0d0
      nrmxl = dnrm2(n-l+1,x(l,l),1)
      t = work(l+p)
      if(t .gt. 0.)then
      t = nrmxl/t
      endif
      if(t .lt. eps)then
      call dshift(x,ldx,n,l,curpvt)
      jp = jpvt(l)
      t=qraux(l)
      tt=work(l)
      ww = work(l+p)
      j=l+1
23036 if(.not.(j.le.curpvt))goto 23038
      jj=j-1
      jpvt(jj)=jpvt(j)
      qraux(jj)=qraux(j)
      work(jj)=work(j)
      work(jj+p) = work(j+p)
      j=j+1
      goto 23036
23038 continue
      jpvt(curpvt)=jp
      qraux(curpvt)=t
      work(curpvt)=tt
      work(curpvt+p) = ww
      curpvt=curpvt-1
      if(lup.gt.curpvt)then
      lup=curpvt
      endif
      else
      if(l.eq.n)then
      goto 23031
      endif
      if(x(l,l).ne.0.0d0)then
      nrmxl = dsign(nrmxl,x(l,l))
      endif
      call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
      x(l,l) = 1.0d0+x(l,l)
      j=l+1
23045 if(.not.(j.le.curpvt))goto 23047
      t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
      call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
      if(qraux(j).ne.0.0d0)then
      tt = 1.0d0-(dabs(x(l,j))/qraux(j))**2
      tt = dmax1(tt,0.0d0)
      t = tt
      tt = 1.0d0+0.05d0*tt*(qraux(j)/work(j))**2
      if(tt.ne.1.0d0)then
      qraux(j) = qraux(j)*dsqrt(t)
      else
      qraux(j) = dnrm2(n-l,x(l+1,j),1)
      work(j) = qraux(j)
      endif
      endif
      j=j+1
      goto 23045
23047 continue
      qraux(l) = x(l,l)
      x(l,l) = -nrmxl
      l=l+1
      endif
      goto 23030
      endif
23031 continue
      rank = lup
      return
      end
      subroutine dchdc(a,lda,p,work,jpvt,job,info)
      integer lda,p,jpvt(p),job,info
      double precision a(lda,p),work(p)
      integer pu,pl,plp1,j,jp,jt,k,kb,km1,kp1,l,maxl
      double precision temp
      double precision maxdia
      logical swapk,negk
      pl = 1
      pu = 0
      info = p
      if(job.ne.0)then
      do23054 k = 1,p 
      swapk = jpvt(k).gt.0
      negk = jpvt(k).lt.0
      jpvt(k) = k
      if(negk)then
      jpvt(k) = -jpvt(k)
      endif
      if(swapk)then
      if(k.ne.pl)then
      call dswap(pl-1,a(1,k),1,a(1,pl),1)
      temp = a(k,k)
      a(k,k) = a(pl,pl)
      a(pl,pl) = temp
      plp1 = pl+1
      if(p.ge.plp1)then
      do23064 j = plp1,p
      if(j.lt.k)then
      temp = a(pl,j)
      a(pl,j) = a(j,k)
      a(j,k) = temp
      else
      if(j.ne.k)then
      temp = a(k,j)
      a(k,j) = a(pl,j)
      a(pl,j) = temp
      endif
      endif
23064 continue
      continue
      endif
      jpvt(k) = jpvt(pl)
      jpvt(pl) = k
      endif
      pl = pl+1
      endif
23054 continue
      continue
      pu = p
      if(p.ge.pl)then
      do23072 kb = pl,p 
      k = p-kb+pl
      if(jpvt(k).lt.0)then
      jpvt(k) = -jpvt(k)
      if(pu.ne.k)then
      call dswap(k-1,a(1,k),1,a(1,pu),1)
      temp = a(k,k)
      a(k,k) = a(pu,pu)
      a(pu,pu) = temp
      kp1 = k+1
      if(p.ge.kp1)then
      do23080 j = kp1,p
      if(j.lt.pu)then
      temp = a(k,j)
      a(k,j) = a(j,pu)
      a(j,pu) = temp
      else
      if(j.ne.pu)then
      temp = a(k,j)
      a(k,j) = a(pu,j)
      a(pu,j) = temp
      endif
      endif
23080 continue
      continue
      endif
      jt = jpvt(k)
      jpvt(k) = jpvt(pu)
      jpvt(pu) = jt
      endif
      pu = pu-1
      endif
23072 continue
      continue
      endif
      endif
      do23086 k = 1,p 
      maxdia = a(k,k)
      kp1 = k+1
      maxl = k
      if(k.ge.pl.and.k.lt.pu)then
      do23090 l = kp1,pu
      if(a(l,l).gt.maxdia)then
      maxdia = a(l,l)
      maxl = l
      endif
23090 continue
      continue
      endif
      if(maxdia.le.0.0d0)then
      go to 10
      endif
      if(k.ne.maxl)then
      km1 = k-1
      call dswap(km1,a(1,k),1,a(1,maxl),1)
      a(maxl,maxl) = a(k,k)
      a(k,k) = maxdia
      jp = jpvt(maxl)
      jpvt(maxl) = jpvt(k)
      jpvt(k) = jp
      endif
      work(k) = dsqrt(a(k,k))
      a(k,k) = work(k)
      if(p.ge.kp1)then
      do23100 j = kp1,p 
      if(k.ne.maxl)then
      if(j.lt.maxl)then
      temp = a(k,j)
      a(k,j) = a(j,maxl)
      a(j,maxl) = temp
      else
      if(j.ne.maxl)then
      temp = a(k,j)
      a(k,j) = a(maxl,j)
      a(maxl,j) = temp
      endif
      endif
      endif
      a(k,j) = a(k,j)/work(k)
      work(j) = a(k,j)
      temp = -a(k,j)
      call daxpy(j-k,temp,work(kp1),1,a(kp1,j),1)
23100 continue
      continue
      endif
23086 continue
      continue
      return
10    info = k-1
      return
      end
      double precision function epslon(x)
      double precision x
      double precision a,b,c,eps
      a = 4.0d0/3.0d0
23108 continue
      b = a-1.0d0
      c = b+b+b
      eps = dabs(c-1.0d0)
      if(.not.(eps.ne.0.0d0))goto 23108
      continue
      epslon = eps*dabs(x)
      return
      end
      double precision function pythag(a,b)
      double precision a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if(p.ne.0.0d0)then
      r = (dmin1(dabs(a),dabs(b))/p)**2
23113 continue
      t = 4.0d0+r
      if(t.eq.4.0d0)then
      goto 23115
      endif
      s = r/t
      u = 1.0d0+2.0d0*s
      p = u*p
      r = (s/u)**2*r
      goto 23113
23115 continue
      endif
      pythag = p
      return
      end
      subroutine rg(nm,n,a,wr,wi,matz,z,iv1,fv1,ierr)
      integer n,nm,is1,is2,ierr,matz
      double precision a(nm,n),wr(n),wi(n),z(nm,n),fv1(n)
      integer iv1(n)
      if(n.gt.nm)then
      ierr = 10*n
      else
      call balanc(nm,n,a,is1,is2,fv1)
      call elmhes(nm,n,is1,is2,a,iv1)
      if(matz.eq.0)then
      call hqr(nm,n,is1,is2,a,wr,wi,ierr)
      else
      call eltran(nm,n,is1,is2,a,iv1,z)
      call hqr2(nm,n,is1,is2,a,wr,wi,z,ierr)
      if(ierr.eq.0)then
      call balbak(nm,n,is1,is2,fv1,n,z)
      endif
      endif
      endif
      return
      end
      subroutine chol(a,p,work,jpvt,job,info)
      integer p,jpvt(*),job,info(*)
      double precision a(p,*),work(*)
      integer i,j
      j =2
23124 if(.not.(j.le.p))goto 23126
      i=1
23127 if(.not.(i.lt.j))goto 23129
      if(a(i,j).ne.a(j,i))then
      info(1) = -1 
      return
      endif
      i = i+1
      goto 23127
23129 continue
      j = j+1
      goto 23124
23126 continue
      call dchdc(a,p,p,work,jpvt,job,info(1))
      j =2
23132 if(.not.(j.le.p))goto 23134
      i=1
23135 if(.not.(i.lt.j))goto 23137
      a(j,i) = 0.
      i = i+1
      goto 23135
23137 continue
      j = j+1
      goto 23132
23134 continue
      return
      end
      subroutine crs(x,dmx,matz,w,z,fv1,fv2,ierr)
      double precision x(*),w(*),z(*),fv1(*),fv2(*)
      integer dmx(2),nx,nv,ierr,matz
      nx=dmx(1)
      nv=dmx(2)
      call rs(nx,nv,x,w,matz,z,fv1,fv2,ierr)
      return
      end
      subroutine dqrls2(x,dx,pivot,qraux,y,dy,beta,res,qt,scrtch,eps)
      integer pivot(*),dx(2),dy(2)
      double precision x(*), qraux(*), y(*), beta(*),res(*),qt(*), scrtc
     *h(*),eps
      integer n,p,q,kn,kp,k,info,rank
      n=dx(1)
      p=dx(2)
      q=dy(2)
      call dqrdca(x,n,n,p,qraux,pivot,scrtch,rank,eps)
      kn=1
      kp=1
      k=1
23138 if(.not.(k.le.q))goto 23140
      call dqrsl(x,n,n,p,qraux,y(kn),scrtch,qt(kn),beta(kp), res(kn),scr
     *tch,00110,info)
      kn = kn+n
      kp=kp+p
      k=k+1
      goto 23138
23140 continue
      return
      end
      subroutine dsvdc1(x,dmx,job,work,e,s,u,v,info)
      double precision x(*),work(*),s(*),e(*),u(*),v(*)
      integer dmx(2),nx,nv,job,info
      nx=dmx(1)
      nv=dmx(2)
      call dsvdc(x,nx,nx,nv,s,e,u,nx,v,nv,work,job,info)
      return
      end
      subroutine balanc(nm,n,a,low,igh,scale)
      integer i,j,k,l,m,n,nm,igh,low,iexc
      double precision a(nm,n),scale(n)
      double precision c,f,g,r,s,b2,radix
      logical noconv
      radix = 16.0d0
      b2 = radix*radix
      k = 1
      l = n
23141 continue
      j=l
23144 if(.not.(j.gt.0))goto 23146
      do23147 i = 1,l
      if(i.ne.j)then
      if(a(j,i).ne.0.0d0)then
      goto 23145
      endif
      endif
23147 continue
      continue
      go to 10
23145 j=j-1 
      goto 23144
23146 continue
      go to 20
10    m = l
      iexc = 1
23153 continue
      scale(m) = j
      if(j.ne.m)then
      do23158 i = 1,l 
      f = a(i,j)
      a(i,j) = a(i,m)
      a(i,m) = f
23158 continue
      continue
      do23160 i = k,n 
      f = a(j,i)
      a(j,i) = a(m,i)
      a(m,i) = f
23160 continue
      continue
      endif
      I23162=(iexc)
      goto 23162
23164 continue
      if(l.eq.1)then
      go to 40
      endif
      l = l-1
      goto 23155
      goto 23163
23167 continue
      k = k+1
20    do23168 j = k,l 
      do23170 i = k,l
      if(i.ne.j)then
      if(a(i,j).ne.0.0d0)then
      goto 23168
      endif
      endif
23170 continue
      continue
      go to 30
23168 continue
      continue
      goto 23143
30    m = k
      iexc = 2
      goto 23163
23162 continue
      if (I23162.eq.1)goto 23164
      if (I23162.eq.2)goto 23167
23163 continue
      goto 23153
23155 continue
      goto 23141
23143 continue
      do23176 i = k,l
      scale(i) = 1.0d0
23176 continue
      continue
23178 continue
      noconv = .false.
      do23181 i = k,l 
      c = 0.0d0
      r = 0.0d0
      do23183 j = k,l
      if(j.ne.i)then
      c = c+dabs(a(j,i))
      r = r+dabs(a(i,j))
      endif
23183 continue
      continue
      if(c.ne.0.0d0.and.r.ne.0.0d0)then
      g = r/radix
      f = 1.0d0
      s = c+r
23189 if(c.lt.g)then
      f = f*radix
      c = c*b2
      goto 23189
      endif
      continue
      g = r*radix
23191 if(c.ge.g)then
      f = f/radix
      c = c/b2
      goto 23191
      endif
      continue
      if((c+r)/f.lt.0.95d0*s)then
      g = 1.0d0/f
      scale(i) = scale(i)*f
      noconv = .true.
      do23195 j = k,n
      a(i,j) = a(i,j)*g
23195 continue
      continue
      do23197 j = 1,l
      a(j,i) = a(j,i)*f
23197 continue
      continue
      endif
      endif
23181 continue
      continue
      if(.not.(.not.noconv))goto 23178
      continue
40    low = k
      igh = l
      return
      end
      subroutine balbak(nm,n,low,igh,scale,m,z)
      integer i,j,k,m,n,ii,nm,igh,low
      double precision scale(n),z(nm,m)
      double precision s
      if(m.ne.0)then
      if(igh.ne.low)then
      do23203 i = low,igh 
      s = scale(i)
      do23205 j = 1,m
      z(i,j) = z(i,j)*s
23205 continue
      continue
23203 continue
      continue
      endif
      do23207 ii = 1,n 
      i = ii
      if(i.lt.low.or.i.gt.igh)then
      if(i.lt.low)then
      i = low-ii
      endif
      k = int(scale(i))
      if(k.ne.i)then
      do23215 j = 1,m 
      s = z(i,j)
      z(i,j) = z(k,j)
      z(k,j) = s
23215 continue
      continue
      endif
      endif
23207 continue
      continue
      endif
      return
      end
      subroutine elmhes(nm,n,low,igh,a,int)
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      double precision a(nm,n)
      double precision x,y
      integer int(igh)
      la = igh-1
      kp1 = low+1
      if(la.ge.kp1)then
      do23219 m = kp1,la 
      mm1 = m-1
      x = 0.0d0
      i = m
      do23221 j = m,igh
      if(dabs(a(j,mm1)).gt.dabs(x))then
      x = a(j,mm1)
      i = j
      endif
23221 continue
      continue
      int(m) = i
      if(i.ne.m)then
      do23227 j = mm1,n 
      y = a(i,j)
      a(i,j) = a(m,j)
      a(m,j) = y
23227 continue
      continue
      do23229 j = 1,igh 
      y = a(j,i)
      a(j,i) = a(j,m)
      a(j,m) = y
23229 continue
      continue
      endif
      if(x.ne.0.0d0)then
      mp1 = m+1
      do23233 i = mp1,igh 
      y = a(i,mm1)
      if(y.ne.0.0d0)then
      y = y/x
      a(i,mm1) = y
      do23237 j = m,n
      a(i,j) = a(i,j)-y*a(m,j)
23237 continue
      continue
      do23239 j = 1,igh
      a(j,m) = a(j,m)+y*a(j,i)
23239 continue
      continue
      endif
23233 continue
      continue
      endif
23219 continue
      continue
      endif
      return
      end
      subroutine eltran(nm,n,low,igh,a,int,z)
      integer i,j,n,kl,mp,nm,igh,low,mp1
      double precision a(nm,igh),z(nm,n)
      integer int(igh)
      do23241 j = 1,n 
      do23243 i = 1,n
      z(i,j) = 0.0d0
23243 continue
      continue
      z(j,j) = 1.0d0
23241 continue
      continue
      kl = igh-low-1
      if(kl.ge.1)then
      mp = igh-1
23247 if(.not.(mp .gt. low))goto 23249
      mp1 = mp+1
      do23250 i = mp1,igh
      z(i,mp) = a(i,mp-1)
23250 continue
      continue
      i = int(mp)
      if(i.ne.mp)then
      do23254 j = mp,igh 
      z(mp,j) = z(i,j)
      z(i,j) = 0.0d0
23254 continue
      continue
      z(i,mp) = 1.0d0
      endif
      mp = mp -1
      goto 23247
23249 continue
      endif
      return
      end
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
      integer i,j,k,l,m,n,en,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      double precision h(nm,n),wr(n),wi(n)
      double precision p,q,r,s,t,w,x,y,zz,norm,tst1,tst2
      logical notlas
c     Initializing zz, b, s, r, q, m, p, zrank to avoid warnings
      zz = 0d0
      b = 0d0
      s = 0d0
      r = 0d0
      q = 0d0
      m = 0
      p = 0d0
      zrank = 0d0
      
      ierr = 0
      norm = 0.0d0
      k = 1
      do23256 i = 1,n 
      do23258 j = k,n
      norm = norm+dabs(h(i,j))
23258 continue
      continue
      k = i
      if(i.lt.low.or.i.gt.igh)then
      wr(i) = h(i,i)
      wi(i) = 0.0d0
      endif
23256 continue
      continue
      en = igh
      t = 0.0d0
      itn = 30*n
23262 continue
      if(en.lt.low)then
      return
      endif
      its = 0
      na = en-1
      enm2 = na-1
23267 continue
      l=en
23270 if(.not.(l .gt. low))goto 23272
      s = dabs(h(l-1,l-1))+dabs(h(l,l))
      if(s.eq.0.0d0)then
      s = norm
      endif
      tst1 = s
      tst2 = tst1+dabs(h(l,l-1))
      if(tst2.eq.tst1)then
      goto 23272
      endif
      l = l-1
      goto 23270
23272 continue
      x = h(en,en)
      if(l.eq.en)then
      go to 50
      endif
      y = h(na,na)
      w = h(en,na)*h(na,en)
      if(l.eq.na)then
      goto 23269
      endif
      if(itn.eq.0)then
      goto 23264
      endif
      if(its.eq.10.or.its.eq.20)then
      t = t+x
      do23285 i = low,en
      h(i,i) = h(i,i)-x
23285 continue
      continue
      s = dabs(h(en,na))+dabs(h(na,enm2))
      x = 0.75d0*s
      y = x
      w = -0.4375d0*s*s
      endif
      its = its+1
      itn = itn-1
      do23287 mm = l,enm2 
      m = enm2+l-mm
      zz = h(m,m)
      r = x-zz
      s = y-zz
      p = (r*s-w)/h(m+1,m)+h(m,m+1)
      q = h(m+1,m+1)-zz-r-s
      r = h(m+2,m+1)
      s = dabs(p)+dabs(q)+dabs(r)
      p = p/s
      q = q/s
      r = r/s
      if(m.eq.l)then
      goto 23288
      endif
      tst1 = dabs(p)*(dabs(h(m-1,m-1))+dabs(zz)+dabs(h(m+1,m+1)))
      tst2 = tst1+dabs(h(m,m-1))*(dabs(q)+dabs(r))
      if(tst2.eq.tst1)then
      goto 23288
      endif
23287 continue
23288 continue
      mp2 = m+2
      do23293 i = mp2,en 
      h(i,i-2) = 0.0d0
      if(i.ne.mp2)then
      h(i,i-3) = 0.0d0
      endif
23293 continue
      continue
      do23297 k = m,na 
      notlas = k.ne.na
      if(k.ne.m)then
      p = h(k,k-1)
      q = h(k+1,k-1)
      r = 0.0d0
      if(notlas)then
      r = h(k+2,k-1)
      endif
      x = dabs(p)+dabs(q)+dabs(r)
      if(x.eq.0.0d0)then
      goto 23297
      endif
      p = p/x
      q = q/x
      r = r/x
      endif
      s = dsign(dsqrt(p*p+q*q+r*r),p)
      if(k.ne.m)then
      h(k,k-1) = -s*x
      else
      if(l.ne.m)then
      h(k,k-1) = -h(k,k-1)
      endif
      endif
      p = p+s
      x = p/s
      y = q/s
      zz = r/s
      q = q/p
      r = r/p
      if(.not.notlas)then
      do23311 j = k,n 
      p = h(k,j)+q*h(k+1,j)
      h(k,j) = h(k,j)-p*x
      h(k+1,j) = h(k+1,j)-p*y
23311 continue
      continue
      j = min0(en,k+3)
      do23313 i = 1,j 
      p = x*h(i,k)+y*h(i,k+1)
      h(i,k) = h(i,k)-p
      h(i,k+1) = h(i,k+1)-p*q
23313 continue
      continue
      else
      do23315 j = k,n 
      p = h(k,j)+q*h(k+1,j)+r*h(k+2,j)
      h(k,j) = h(k,j)-p*x
      h(k+1,j) = h(k+1,j)-p*y
      h(k+2,j) = h(k+2,j)-p*zz
23315 continue
      continue
      j = min0(en,k+3)
      do23317 i = 1,j 
      p = x*h(i,k)+y*h(i,k+1)+zz*h(i,k+2)
      h(i,k) = h(i,k)-p
      h(i,k+1) = h(i,k+1)-p*q
      h(i,k+2) = h(i,k+2)-p*r
23317 continue
      continue
      endif
23297 continue
      continue
      goto 23267
23269 continue
      p = (y-x)/2.0d0
      q = p*p+w
      zz = dsqrt(dabs(q))
      x = x+t
      if(q.lt.0.0d0)then
      wr(na) = x+p
      wr(en) = x+p
      wi(na) = zz
      wi(en) = -zz
      else
      zz = p+dsign(zz,p)
      wr(na) = x+zz
      wr(en) = wr(na)
      if(zz.ne.0.0d0)then
      wr(en) = x-w/zz
      endif
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      endif
      en = enm2
      goto 23263
50    wr(en) = x+t
      wi(en) = 0.0d0
      en = na
23263 goto 23262
23264 continue
      ierr = en
      return
      end
      subroutine hqr2(nm,n,low,igh,h,wr,wi,z,ierr)
      integer i,j,k,l,m,n,en,ii,jj,ll,mm,na,nm,nn,igh,itn,its,low,mp2,en
     *m2,ierr
      double precision h(nm,n),wr(n),wi(n),z(nm,n)
      double precision p,q,r,s,t,w,x,y,ra,sa,vi,vr,zz,norm,tst1,tst2
      logical notlas
c     Setting zz, r, s, q, p, m, l to zero to avoid warnings
      zz = 0d0
      r = 0d0
      s = 0d0
      q = 0d0
      p = 0d0
      m = 0
      l = 0
      ierr = 0
      norm = 0.0d0
      k = 1
      do23323 i = 1,n 
      do23325 j = k,n
      norm = norm+dabs(h(i,j))
23325 continue
      continue
      k = i
      if(i.lt.low.or.i.gt.igh)then
      wr(i) = h(i,i)
      wi(i) = 0.0d0
      endif
23323 continue
      continue
      en = igh
      t = 0.0d0
      itn = 30*n
23329 continue
      if(en.lt.low)then
      go to 70
      endif
      its = 0
      na = en-1
      enm2 = na-1
23334 continue
      do23337 ll = low,en 
      l = en+low-ll
      if(l.eq.low)then
      goto 23338
      endif
      s = dabs(h(l-1,l-1))+dabs(h(l,l))
      if(s.eq.0.0d0)then
      s = norm
      endif
      tst1 = s
      tst2 = tst1+dabs(h(l,l-1))
      if(tst2.eq.tst1)then
      goto 23338
      endif
23337 continue
23338 continue
      x = h(en,en)
      if(l.eq.en)then
      go to 60
      endif
      y = h(na,na)
      w = h(en,na)*h(na,en)
      if(l.eq.na)then
      goto 23336
      endif
      if(itn.eq.0)then
      goto 23331
      endif
      if(its.eq.10.or.its.eq.20)then
      t = t+x
      do23353 i = low,en
      h(i,i) = h(i,i)-x
23353 continue
      continue
      s = dabs(h(en,na))+dabs(h(na,enm2))
      x = 0.75d0*s
      y = x
      w = -0.4375d0*s*s
      endif
      its = its+1
      itn = itn-1
      do23355 mm = l,enm2 
      m = enm2+l-mm
      zz = h(m,m)
      r = x-zz
      s = y-zz
      p = (r*s-w)/h(m+1,m)+h(m,m+1)
      q = h(m+1,m+1)-zz-r-s
      r = h(m+2,m+1)
      s = dabs(p)+dabs(q)+dabs(r)
      p = p/s
      q = q/s
      r = r/s
      if(m.eq.l)then
      goto 23356
      endif
      tst1 = dabs(p)*(dabs(h(m-1,m-1))+dabs(zz)+dabs(h(m+1,m+1)))
      tst2 = tst1+dabs(h(m,m-1))*(dabs(q)+dabs(r))
      if(tst2.eq.tst1)then
      goto 23356
      endif
23355 continue
23356 continue
      mp2 = m+2
      do23361 i = mp2,en 
      h(i,i-2) = 0.0d0
      if(i.ne.mp2)then
      h(i,i-3) = 0.0d0
      endif
23361 continue
      continue
      do23365 k = m,na 
      notlas = k.ne.na
      if(k.ne.m)then
      p = h(k,k-1)
      q = h(k+1,k-1)
      r = 0.0d0
      if(notlas)then
      r = h(k+2,k-1)
      endif
      x = dabs(p)+dabs(q)+dabs(r)
      if(x.eq.0.0d0)then
      goto 23365
      endif
      p = p/x
      q = q/x
      r = r/x
      endif
      s = dsign(dsqrt(p*p+q*q+r*r),p)
      if(k.ne.m)then
      h(k,k-1) = -s*x
      else
      if(l.ne.m)then
      h(k,k-1) = -h(k,k-1)
      endif
      endif
      p = p+s
      x = p/s
      y = q/s
      zz = r/s
      q = q/p
      r = r/p
      if(.not.notlas)then
      do23379 j = k,n 
      p = h(k,j)+q*h(k+1,j)
      h(k,j) = h(k,j)-p*x
      h(k+1,j) = h(k+1,j)-p*y
23379 continue
      continue
      j = min0(en,k+3)
      do23381 i = 1,j 
      p = x*h(i,k)+y*h(i,k+1)
      h(i,k) = h(i,k)-p
      h(i,k+1) = h(i,k+1)-p*q
23381 continue
      continue
      do23383 i = low,igh 
      p = x*z(i,k)+y*z(i,k+1)
      z(i,k) = z(i,k)-p
      z(i,k+1) = z(i,k+1)-p*q
23383 continue
      continue
      else
      do23385 j = k,n 
      p = h(k,j)+q*h(k+1,j)+r*h(k+2,j)
      h(k,j) = h(k,j)-p*x
      h(k+1,j) = h(k+1,j)-p*y
      h(k+2,j) = h(k+2,j)-p*zz
23385 continue
      continue
      j = min0(en,k+3)
      do23387 i = 1,j 
      p = x*h(i,k)+y*h(i,k+1)+zz*h(i,k+2)
      h(i,k) = h(i,k)-p
      h(i,k+1) = h(i,k+1)-p*q
      h(i,k+2) = h(i,k+2)-p*r
23387 continue
      continue
      do23389 i = low,igh 
      p = x*z(i,k)+y*z(i,k+1)+zz*z(i,k+2)
      z(i,k) = z(i,k)-p
      z(i,k+1) = z(i,k+1)-p*q
      z(i,k+2) = z(i,k+2)-p*r
23389 continue
      continue
      endif
23365 continue
      continue
      goto 23334
23336 continue
      p = (y-x)/2.0d0
      q = p*p+w
      zz = dsqrt(dabs(q))
      h(en,en) = x+t
      x = h(en,en)
      h(na,na) = y+t
      if(q.lt.0.0d0)then
      wr(na) = x+p
      wr(en) = x+p
      wi(na) = zz
      wi(en) = -zz
      else
      zz = p+dsign(zz,p)
      wr(na) = x+zz
      wr(en) = wr(na)
      if(zz.ne.0.0d0)then
      wr(en) = x-w/zz
      endif
      wi(na) = 0.0d0
      wi(en) = 0.0d0
      x = h(en,na)
      s = dabs(x)+dabs(zz)
      p = x/s
      q = zz/s
      r = dsqrt(p*p+q*q)
      p = p/r
      q = q/r
      do23395 j = na,n 
      zz = h(na,j)
      h(na,j) = q*zz+p*h(en,j)
      h(en,j) = q*h(en,j)-p*zz
23395 continue
      continue
      do23397 i = 1,en 
      zz = h(i,na)
      h(i,na) = q*zz+p*h(i,en)
      h(i,en) = q*h(i,en)-p*zz
23397 continue
      continue
      do23399 i = low,igh 
      zz = z(i,na)
      z(i,na) = q*zz+p*z(i,en)
      z(i,en) = q*z(i,en)-p*zz
23399 continue
      continue
      endif
      en = enm2
      goto 23330
60    h(en,en) = x+t
      wr(en) = h(en,en)
      wi(en) = 0.0d0
      en = na
23330 goto 23329
23331 continue
      ierr = en
      return
70    if(norm.ne.0.0d0)then
      do23403 nn = 1,n 
      en = n+1-nn
      p = wr(en)
      q = wi(en)
      na = en-1
      if(q.lt.0)then
      m = na
      if(dabs(h(en,na)).le.dabs(h(na,en)))then
      call cdiv(0.0d0,-h(na,en),h(na,na)-p,q,h(na,na),h(na,en))
      else
      h(na,na) = q/h(en,na)
      h(na,en) = -(h(en,en)-p)/h(en,na)
      endif
      h(en,na) = 0.0d0
      h(en,en) = 1.0d0
      enm2 = na-1
      if(enm2.ne.0)then
      do23411 ii = 1,enm2 
      i = na-ii
      w = h(i,i)-p
      ra = 0.0d0
      sa = 0.0d0
      do23413 j = m,en 
      ra = ra+h(i,j)*h(j,na)
      sa = sa+h(i,j)*h(j,en)
23413 continue
      continue
      if(wi(i).lt.0.0d0)then
      zz = w
      r = ra
      s = sa
      else
      m = i
      if(wi(i).eq.0.0d0)then
      call cdiv(-ra,-sa,w,q,h(i,na),h(i,en))
      else
      x = h(i,i+1)
      y = h(i+1,i)
      vr = (wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)-q*q
      vi = (wr(i)-p)*2.0d0*q
      if(vr.eq.0.0d0.and.vi.eq.0.0d0)then
      tst1 = norm*(dabs(w)+dabs(q)+dabs(x)+dabs(y)+dabs(zz))
      vr = tst1
23421 continue
      vr = 0.01d0*vr
      tst2 = tst1+vr
      if(.not.(tst2.le.tst1))goto 23421
      continue
      endif
      call cdiv(x*r-zz*ra+q*sa,x*s-zz*sa-q*ra,vr,vi,h(i,na),h(i,en))
      if(dabs(x).le.dabs(zz)+dabs(q))then
      call cdiv(-r-y*h(i,na),-s-y*h(i,en),zz,q,h(i+1,na),h(i+1,en))
      else
      h(i+1,na) = (-ra-w*h(i,na)+q*h(i,en))/x
      h(i+1,en) = (-sa-w*h(i,en)-q*h(i,na))/x
      endif
      endif
      t = dmax1(dabs(h(i,na)),dabs(h(i,en)))
      if(t.ne.0.0d0)then
      tst1 = t
      tst2 = tst1+1.0d0/tst1
      if(tst2.le.tst1)then
      do23430 j = i,en 
      h(j,na) = h(j,na)/t
      h(j,en) = h(j,en)/t
23430 continue
      continue
      endif
      endif
      endif
23411 continue
      continue
      endif
      else
      if(q.eq.0)then
      m = en
      h(en,en) = 1.0d0
      if(na.ne.0)then
      do23436 ii = 1,na 
      i = en-ii
      w = h(i,i)-p
      r = 0.0d0
      do23438 j = m,en
      r = r+h(i,j)*h(j,en)
23438 continue
      continue
      if(wi(i).lt.0.0d0)then
      zz = w
      s = r
      else
      m = i
      if(wi(i).ne.0.0d0)then
      x = h(i,i+1)
      y = h(i+1,i)
      q = (wr(i)-p)*(wr(i)-p)+wi(i)*wi(i)
      t = (x*s-zz*r)/q
      h(i,en) = t
      if(dabs(x).le.dabs(zz))then
      h(i+1,en) = (-s-y*t)/zz
      else
      h(i+1,en) = (-r-w*t)/x
      endif
      else
      t = w
      if(t.eq.0.0d0)then
      tst1 = norm
      t = tst1
23448 continue
      t = 0.01d0*t
      tst2 = norm+t
      if(.not.(tst2.le.tst1))goto 23448
      continue
      endif
      h(i,en) = -r/t
      endif
      t = dabs(h(i,en))
      if(t.ne.0.0d0)then
      tst1 = t
      tst2 = tst1+1.0d0/tst1
      if(tst2.le.tst1)then
      do23455 j = i,en
      h(j,en) = h(j,en)/t
23455 continue
      continue
      endif
      endif
      endif
23436 continue
      continue
      endif
      endif
      endif
23403 continue
      continue
      do23457 i = 1,n
      if(i.lt.low.or.i.gt.igh)then
      do23461 j = i,n
      z(i,j) = h(i,j)
23461 continue
      continue
      endif
23457 continue
      continue
      do23463 jj = low,n 
      j = n+low-jj
      m = min0(j,igh)
      do23465 i = low,igh 
      zz = 0.0d0
      do23467 k = low,m
      zz = zz+z(i,k)*h(k,j)
23467 continue
      continue
      z(i,j) = zz
23465 continue
      continue
23463 continue
      continue
      endif
      return
      end
      subroutine cdiv(ar,ai,br,bi,cr,ci)
      double precision ar,ai,br,bi,cr,ci
      double precision s,ars,ais,brs,bis
      s = dabs(br)+dabs(bi)
      ars = ar/s
      ais = ai/s
      brs = br/s
      bis = bi/s
      s = brs**2+bis**2
      cr = (ars*brs+ais*bis)/s
      ci = (ais*brs-ars*bis)/s
      return
      end
      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
      integer n,nm,ierr,matz
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
      if(n.gt.nm)then
      ierr = 10*n
      else
      if(matz.ne.0)then
      call tred2(nm,n,a,w,fv1,z)
      call tql2(nm,n,w,fv1,z,ierr)
      else
      call tred1(nm,n,a,w,fv1,fv2)
      call tqlrat(n,w,fv2,ierr)
      endif
      endif
      return
      end
      subroutine tql2(nm,n,d,e,z,ierr)
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
C     Initializing s2, c3 and ensuring ierr is not clobbered      
      s2 = 0.0
      c3 = 0.0
      ierr = 0
      if(n.ne.1)then
      do23475 i = 2,n
      e(i-1) = e(i)
23475 continue
      continue
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
      do23477 l = 1,n 
      j = 0
      h = dabs(d(l))+dabs(e(l))
      if(tst1.lt.h)then
      tst1 = h
      endif
      do23481 m = l,n 
      tst2 = tst1+dabs(e(m))
      if(tst2.eq.tst1)then
      goto 23482
      endif
23481 continue
23482 continue
      if(m.ne.l)then
23487 continue
      if(j.eq.30)then
      go to 10
      endif
      j = j+1
      l1 = l+1
      l2 = l1+1
      g = d(l)
      p = (d(l1)-g)/(2.0d0*e(l))
      r = pythag(p,1.0d0)
      d(l) = e(l)/(p+dsign(r,p))
      d(l1) = e(l)*(p+dsign(r,p))
      dl1 = d(l1)
      h = g-d(l)
      if(l2.le.n)then
      do23494 i = l2,n
      d(i) = d(i)-h
23494 continue
      continue
      endif
      f = f+h
      p = d(m)
      c = 1.0d0
      c2 = c
      el1 = e(l1)
      s = 0.0d0
      mml = m-l
      do23496 ii = 1,mml 
      c3 = c2
      c2 = c
      s2 = s
      i = m-ii
      g = c*e(i)
      h = c*p
      r = pythag(p,e(i))
      e(i+1) = s*r
      s = e(i)/r
      c = p/r
      p = c*d(i)-s*g
      d(i+1) = h+s*(c*g+s*d(i))
      do23498 k = 1,n 
      h = z(k,i+1)
      z(k,i+1) = s*z(k,i)+c*h
      z(k,i) = c*z(k,i)-s*h
23498 continue
      continue
23496 continue
      continue
      p = -s*s2*c3*el1*e(l)/dl1
      e(l) = s*p
      d(l) = c*p
      tst2 = tst1+dabs(e(l))
      if(.not.(tst2.le.tst1))goto 23487
      continue
      endif
      d(l) = d(l)+f
23477 continue
      continue
      do23500 ii = 2,n 
      i = ii-1
      k = i
      p = d(i)
      do23502 j = ii,n
      if(d(j).lt.p)then
      k = j
      p = d(j)
      endif
23502 continue
      continue
      if(k.ne.i)then
      d(k) = d(i)
      d(i) = p
      do23508 j = 1,n 
      p = z(j,i)
      z(j,i) = z(j,k)
      z(j,k) = p
23508 continue
      continue
      endif
23500 continue
      continue
      return
10    ierr = l
      endif
      return
      end
      subroutine tqlrat(n,d,e2,ierr)
      integer i,j,l,m,n,ii,l1,mml,ierr
      double precision d(n),e2(n)
      double precision b,c,f,g,h,p,r,s,t,epslon,pythag
C     Initializing b, c and making sure ierr is not globbered
      c = 0.0
      b = 0.0
      ierr = 0
      if(n.ne.1)then
      do23512 i = 2,n
      e2(i-1) = e2(i)
23512 continue
      continue
      f = 0.0d0
      t = 0.0d0
      e2(n) = 0.0d0
      do23514 l = 1,n 
      j = 0
      h = dabs(d(l))+dsqrt(e2(l))
      if(t.le.h)then
      t = h
      b = epslon(t)
      c = b*b
      endif
      do23518 m = l,n
      if(e2(m).le.c)then
      goto 23519
      endif
23518 continue
23519 continue
      if(m.ne.l)then
23524 continue
      if(j.eq.30)then
      go to 20
      endif
      j = j+1
      l1 = l+1
      s = dsqrt(e2(l))
      g = d(l)
      p = (d(l1)-g)/(2.0d0*s)
      r = pythag(p,1.0d0)
      d(l) = s/(p+dsign(r,p))
      h = g-d(l)
      do23529 i = l1,n
      d(i) = d(i)-h
23529 continue
      continue
      f = f+h
      g = d(m)
      if(g.eq.0.0d0)then
      g = b
      endif
      h = g
      s = 0.0d0
      mml = m-l
      do23533 ii = 1,mml 
      i = m-ii
      p = g*h
      r = p+e2(i)
      e2(i+1) = s*r
      s = e2(i)/r
      d(i+1) = h+s*(h+d(i))
      g = d(i)-e2(i)/g
      if(g.eq.0.0d0)then
      g = b
      endif
      h = g*p/r
23533 continue
      continue
      e2(l) = s*g
      d(l) = h
      if(h.eq.0.0d0)then
      goto 23526
      endif
      if(dabs(e2(l)).le.dabs(c/h))then
      goto 23526
      endif
      e2(l) = h*e2(l)
      if(.not.(e2(l).eq.0.0d0))goto 23524
23526 continue
      endif
      p = d(l)+f
      if(l.ne.1)then
      do23543 ii = 2,l 
      i = l+2-ii
      if(p.ge.d(i-1))then
      go to 10
      endif
      d(i) = d(i-1)
23543 continue
      continue
      endif
      i = 1
10    d(i) = p
23514 continue
      continue
      return
20    ierr = l
      endif
      return
      end
      subroutine tred1(nm,n,a,d,e,e2)
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),e2(n)
      double precision f,g,h,scale
      do23547 i = 1,n 
      d(i) = a(n,i)
      a(n,i) = a(i,i)
23547 continue
      continue
      do23549 ii = 1,n 
      i = n+1-ii
      l = i-1
      h = 0.0d0
      scale = 0.0d0
      if(l.ge.1)then
      do23553 k = 1,l
      scale = scale+dabs(d(k))
23553 continue
      continue
      if(scale.eq.0.0d0)then
      do23557 j = 1,l 
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = 0.0d0
23557 continue
      continue
      else
      do23559 k = 1,l 
      d(k) = d(k)/scale
      h = h+d(k)*d(k)
23559 continue
      continue
      e2(i) = scale*scale*h
      f = d(l)
      g = -dsign(dsqrt(h),f)
      e(i) = scale*g
      h = h-f*g
      d(l) = f-g
      if(l.ne.1)then
      do23563 j = 1,l
      e(j) = 0.0d0
23563 continue
      continue
      do23565 j = 1,l 
      f = d(j)
      g = e(j)+a(j,j)*f
      jp1 = j+1
      if(l.ge.jp1)then
      do23569 k = jp1,l 
      g = g+a(k,j)*d(k)
      e(k) = e(k)+a(k,j)*f
23569 continue
      continue
      endif
      e(j) = g
23565 continue
      continue
      f = 0.0d0
      do23571 j = 1,l 
      e(j) = e(j)/h
      f = f+e(j)*d(j)
23571 continue
      continue
      h = f/(h+h)
      do23573 j = 1,l
      e(j) = e(j)-h*d(j)
23573 continue
      continue
      do23575 j = 1,l 
      f = d(j)
      g = e(j)
      do23577 k = j,l
      a(k,j) = a(k,j)-f*e(k)-g*d(k)
23577 continue
      continue
23575 continue
      continue
      endif
      do23579 j = 1,l 
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f*scale
23579 continue
      continue
      goto 23549
      endif
      endif
      e(i) = 0.0d0
      e2(i) = 0.0d0
23549 continue
      continue
      return
      end
      subroutine tred2(nm,n,a,d,e,z)
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale
      do23581 i = 1,n 
      do23583 j = i,n
      z(j,i) = a(j,i)
23583 continue
      continue
      d(i) = a(n,i)
23581 continue
      continue
      if(n.ne.1)then
      do23587 ii = 2,n 
      i = n+2-ii
      l = i-1
      h = 0.0d0
      scale = 0.0d0
      if(l.ge.2)then
      do23591 k = 1,l
      scale = scale+dabs(d(k))
23591 continue
      continue
      if(scale.ne.0.0d0)then
      do23595 k = 1,l 
      d(k) = d(k)/scale
      h = h+d(k)*d(k)
23595 continue
      continue
      f = d(l)
      g = -dsign(dsqrt(h),f)
      e(i) = scale*g
      h = h-f*g
      d(l) = f-g
      do23597 j = 1,l
      e(j) = 0.0d0
23597 continue
      continue
      do23599 j = 1,l 
      f = d(j)
      z(j,i) = f
      g = e(j)+z(j,j)*f
      jp1 = j+1
      if(l.ge.jp1)then
      do23603 k = jp1,l 
      g = g+z(k,j)*d(k)
      e(k) = e(k)+z(k,j)*f
23603 continue
      continue
      endif
      e(j) = g
23599 continue
      continue
      f = 0.0d0
      do23605 j = 1,l 
      e(j) = e(j)/h
      f = f+e(j)*d(j)
23605 continue
      continue
      hh = f/(h+h)
      do23607 j = 1,l
      e(j) = e(j)-hh*d(j)
23607 continue
      continue
      do23609 j = 1,l 
      f = d(j)
      g = e(j)
      do23611 k = j,l
      z(k,j) = z(k,j)-f*e(k)-g*d(k)
23611 continue
      continue
      d(j) = z(l,j)
      z(i,j) = 0.0d0
23609 continue
      continue
      go to 10
      endif
      endif
      e(i) = d(l)
      do23613 j = 1,l 
      d(j) = z(l,j)
      z(i,j) = 0.0d0
      z(j,i) = 0.0d0
23613 continue
      continue
10    d(i) = h
23587 continue
      continue
      do23615 i = 2,n 
      l = i-1
      z(n,l) = z(l,l)
      z(l,l) = 1.0d0
      h = d(i)
      if(h.ne.0.0d0)then
      do23619 k = 1,l
      d(k) = z(k,i)/h
23619 continue
      continue
      do23621 j = 1,l 
      g = 0.0d0
      do23623 k = 1,l
      g = g+z(k,i)*z(k,j)
23623 continue
      continue
      do23625 k = 1,l
      z(k,j) = z(k,j)-g*d(k)
23625 continue
      continue
23621 continue
      continue
      endif
      do23627 k = 1,l
      z(k,i) = 0.0d0
23627 continue
      continue
23615 continue
      continue
      endif
      do23629 i = 1,n 
      d(i) = z(n,i)
      z(n,i) = 0.0d0
23629 continue
      continue
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
      subroutine dmatp(x,dx,y,dy,z)
      integer dx(2),dy(2)
      double precision x(*), y(*),z(*),ddot
      integer n,p,q,i,j
      n=dx(1)
      p=dx(2)
      q=dy(2)
      do23631 i = 1,n 
      jj = 1
      ij = i
      do23633 j = 1, q 
      z(ij) = ddot(p,x(i),n,y(jj),1)
      if(j.lt.q)then
      jj = jj + p
      ij = ij + n
      endif
23633 continue
      continue
23631 continue
      continue
      return
      end
      subroutine dmatpt(x,dx,y,dy,z)
      integer dx(2),dy(2)
      double precision x(*), y(*),z(*),ddot
      integer n,p,q,i,j,ii
      n=dx(1)
      p=dx(2)
      q=dy(2)
      ii=1
      do23637 i = 1,p 
      jj = 1
      ij = i
      do23639 j = 1, q 
      z(ij) = ddot(n,x(ii),1,y(jj),1)
      if(j.lt.q)then
      jj = jj + n
      ij = ij + p
      endif
23639 continue
      continue
      ii = ii +n
23637 continue
      continue
      return
      end
      subroutine matpm(x,dx,mmx,mx,y,dy,mmy,my,z)
      integer dx(2),dy(2)
      integer mmx(*), mmy(*)
      integer mx(*), my(*)
      double precision x(*), y(*),z(*),ddot
      integer n,p,q,i,j
      n=dx(1)
      p=dx(2)
      q=dy(2)
      call rowmis(mmx,dx(1),dx(2),mx)
      call colmis(mmy,dy(1),dy(2),my)
      do23643 i = 1,n 
      jj = 1
      ij = i
      do23645 j = 1, q 
      if(.not.(mx(i).ne.0 .or. my(j).ne.0))then
      z(ij) = ddot(p,x(i),n,y(jj),1)
      endif
      if(j.lt.q)then
      jj = jj + p
      ij = ij + n
      endif
23645 continue
      continue
23643 continue
      continue
      return
      end
      subroutine matptm(x,dx,mmx,mx,y,dy,mmy,my,z)
      integer dx(2),dy(2)
      integer mmx(*), mmy(*)
      integer mx(*), my(*)
      double precision x(*), y(*),z(*),ddot
      integer n,p,q,i,j
      call colmis(mmx,dx(1),dx(2),mx)
      call colmis(mmy,dy(1),dy(2),my)
      n=dx(1)
      p=dx(2)
      q=dy(2)
      ii=1
      do23651 i = 1,p 
      jj = 1
      ij = i
      do23653 j = 1, q 
      if(.not.(mx(i).ne.0 .or. my(j).ne.0))then
      z(ij) = ddot(n,x(ii),1,y(jj),1)
      endif
      if(j.lt.q)then
      jj = jj + n
      ij = ij + p
      endif
23653 continue
      continue
      ii = ii +n
23651 continue
      continue
      return
      end
      subroutine rowmis(m,n,p,vec)
      integer n,p
      integer m(n,p)
      integer vec(*)
      do23659 i = 1,n 
      vec(i)=0
      do23661 j = 1,p 
      if(m(i,j).ne.0)then
      vec(i) = 1
      endif
23661 continue
      continue
23659 continue
      continue
      return
      end
      subroutine colmis(m,n,p,vec)
      integer n,p
      integer m(n,p)
      integer vec(*)
      do23665 j = 1,p 
      vec(j)=0
      do23667 i = 1,n 
      if(m(i,j).ne.0)then
      vec(j) = 1
      endif
23667 continue
      continue
23665 continue
      continue
      return
      end
      subroutine dshift(x,ldx,n,j,k)
      integer ldx,n,j,k
      double precision x(ldx,k),tt
      integer i,jj
      if(k.gt.j)then
      do23790 i = 1,n 
      tt = x(i,j)
      do23792 jj = j+1,k
      x(i,jj-1) = x(i,jj)
23792 continue
      continue
      x(i,k) = tt
23790 continue
      continue
      endif
      return
      end
      subroutine rtod(dx,dy,n)
      real dx(*)
      double precision dy(*)
      integer i,m,mp1,n
      if(n.gt.0)then
      m = mod(n,7)
      if(m.ne.0)then
      do23798 i = 1,m
      dy(i) = dx(i)
23798 continue
      continue
      if(n.lt.7)then
      return
      endif
      endif
      mp1 = m+1
      do23802 i = mp1,n,7 
      dy(i) = dx(i)
      dy(i+1) = dx(i+1)
      dy(i+2) = dx(i+2)
      dy(i+3) = dx(i+3)
      dy(i+4) = dx(i+4)
      dy(i+5) = dx(i+5)
      dy(i+6) = dx(i+6)
23802 continue
      continue
      endif
      return
      end
      subroutine dtor(dx,dy,n)
      double precision dx(*)
      real dy(*)
      integer i,m,mp1,n
      if(n.gt.0)then
      m = mod(n,7)
      if(m.ne.0)then
      do23808 i = 1,m
      dy(i) = real(dx(i))
23808 continue
      continue
      if(n.lt.7)then
      return
      endif
      endif
      mp1 = m+1
      do23812 i = mp1,n,7 
      dy(i) = real(dx(i))
      dy(i+1) = real(dx(i+1))
      dy(i+2) = real(dx(i+2))
      dy(i+3) = real(dx(i+3))
      dy(i+4) = real(dx(i+4))
      dy(i+5) = real(dx(i+5))
      dy(i+6) = real(dx(i+6))
23812 continue
      continue
      endif
      return
      end
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      double precision x(ldx,*),qraux(*),y(*),qy(*),qty(*),b(*),rsd(*),x
     *b(*)
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
      info = 0
      cqy = job/10000.ne.0
      cqty = mod(job,10000).ne.0
      cb = mod(job,1000)/100.ne.0
      cr = mod(job,100)/10.ne.0
      cxb = mod(job,10).ne.0
      ju = min0(k,n-1)
      if(ju.eq.0)then
      if(cqy)then
      qy(1) = y(1)
      endif
      if(cqty)then
      qty(1) = y(1)
      endif
      if(cxb)then
      xb(1) = y(1)
      endif
      if(cb)then
      if(x(1,1).ne.0.0d0)then
      b(1) = y(1)/x(1,1)
      else
      info = 1
      endif
      endif
      if(cr)then
      rsd(1) = 0.0d0
      endif
      else
      if(cqy)then
      call dcopy(n,y,1,qy,1)
      endif
      if(cqty)then
      call dcopy(n,y,1,qty,1)
      endif
      if(cqy)then
      do23854 jj = 1,ju 
      j = ju-jj+1
      if(qraux(j).ne.0.0d0)then
      temp = x(j,j)
      x(j,j) = qraux(j)
      t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
      call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
      x(j,j) = temp
      endif
23854 continue
      continue
      endif
      if(cqty)then
      do23860 j = 1,ju
      if(qraux(j).ne.0.0d0)then
      temp = x(j,j)
      x(j,j) = qraux(j)
      t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
      call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
      x(j,j) = temp
      endif
23860 continue
      continue
      endif
      if(cb)then
      call dcopy(k,qty,1,b,1)
      endif
      kp1 = k+1
      if(cxb)then
      call dcopy(k,qty,1,xb,1)
      endif
      if(cr.and.k.lt.n)then
      call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
      endif
      if(cxb.and.kp1.le.n)then
      do23872 i = kp1,n
      xb(i) = 0.0d0
23872 continue
      continue
      endif
      if(cr)then
      do23876 i = 1,k
      rsd(i) = 0.0d0
23876 continue
      continue
      endif
      if(cb)then
      do23880 jj = 1,k 
      j = k-jj+1
      if(x(j,j).eq.0.0d0)then
      go to 130
      endif
      b(j) = b(j)/x(j,j)
      if(j.ne.1)then
      t = -b(j)
      call daxpy(j-1,t,x(1,j),1,b,1)
      endif
23880 continue
      continue
      go to 140
130   info = j
      endif
140   if(cr.or.cxb)then
      do23888 jj = 1,ju 
      j = ju-jj+1
      if(qraux(j).ne.0.0d0)then
      temp = x(j,j)
      x(j,j) = qraux(j)
      if(cr)then
      t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
      call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
      endif
      if(cxb)then
      t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
      call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
      endif
      x(j,j) = temp
      endif
23888 continue
      continue
      endif
      endif
      return
      end
      subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
      integer ldx,n,p,ldu,ldv,job,info
      double precision x(ldx,*),s(*),e(*),u(ldu,*),v(ldv,*),work(*)
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,mm,
     *mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      double precision ddot,t
      double precision b,c,cs,el,emm1,f,g,dnrm2,scale,shift,sl,sm,sn,smm
     *1,t1,test,ztest
      logical wantu,wantv
c     Setting l, ls, c, g, b to zero just to initialize and avoid warnings
      l = 0
      ls = 0
      c = 0d0
      g = 0d0
      b = 0d0
      maxit = 30
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if(jobu.gt.1)then
      ncu = min0(n,p)
      endif
      if(jobu.ne.0)then
      wantu = .true.
      endif
      if(mod(job,10).ne.0)then
      wantv = .true.
      endif
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if(lu.ge.1)then
      do23904 l = 1,lu 
      lp1 = l+1
      if(l.le.nct)then
      s(l) = dnrm2(n-l+1,x(l,l),1)
      if(s(l).ne.0.0d0)then
      if(x(l,l).ne.0.0d0)then
      s(l) = dsign(s(l),x(l,l))
      endif
      call dscal(n-l+1,1.0d0/s(l),x(l,l),1)
      x(l,l) = 1.0d0+x(l,l)
      endif
      s(l) = -s(l)
      endif
      if(p.ge.lp1)then
      do23914 j = lp1,p 
      if(l.le.nct)then
      if(s(l).ne.0.0d0)then
      t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
      call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
      endif
      endif
      e(j) = x(l,j)
23914 continue
      continue
      endif
      if(wantu.and.l.le.nct)then
      do23922 i = l,n
      u(i,l) = x(i,l)
23922 continue
      continue
      endif
      if(l.le.nrt)then
      e(l) = dnrm2(p-l,e(lp1),1)
      if(e(l).ne.0.0d0)then
      if(e(lp1).ne.0.0d0)then
      e(l) = dsign(e(l),e(lp1))
      endif
      call dscal(p-l,1.0d0/e(l),e(lp1),1)
      e(lp1) = 1.0d0+e(lp1)
      endif
      e(l) = -e(l)
      if(lp1.le.n.and.e(l).ne.0.0d0)then
      do23932 i = lp1,n
      work(i) = 0.0d0
23932 continue
      continue
      do23934 j = lp1,p
      call daxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
23934 continue
      continue
      do23936 j = lp1,p
      call daxpy(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
23936 continue
      continue
      endif
      if(wantv)then
      do23940 i = lp1,p
      v(i,l) = e(i)
23940 continue
      continue
      endif
      endif
23904 continue
      continue
      endif
      m = min0(p,n+1)
      nctp1 = nct+1
      nrtp1 = nrt+1
      if(nct.lt.p)then
      s(nctp1) = x(nctp1,nctp1)
      endif
      if(n.lt.m)then
      s(m) = 0.0d0
      endif
      if(nrtp1.lt.m)then
      e(nrtp1) = x(nrtp1,m)
      endif
      e(m) = 0.0d0
      if(wantu)then
      if(ncu.ge.nctp1)then
      do23952 j = nctp1,ncu 
      do23954 i = 1,n
      u(i,j) = 0.0d0
23954 continue
      continue
      u(j,j) = 1.0d0
23952 continue
      continue
      endif
      if(nct.ge.1)then
      do23958 ll = 1,nct 
      l = nct-ll+1
      if(s(l).eq.0.0d0)then
      do23962 i = 1,n
      u(i,l) = 0.0d0
23962 continue
      continue
      u(l,l) = 1.0d0
      else
      lp1 = l+1
      if(ncu.ge.lp1)then
      do23966 j = lp1,ncu 
      t = -ddot(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
      call daxpy(n-l+1,t,u(l,l),1,u(l,j),1)
23966 continue
      continue
      endif
      call dscal(n-l+1,-1.0d0,u(l,l),1)
      u(l,l) = 1.0d0+u(l,l)
      lm1 = l-1
      if(lm1.ge.1)then
      do23970 i = 1,lm1
      u(i,l) = 0.0d0
23970 continue
      continue
      endif
      endif
23958 continue
      continue
      endif
      endif
      if(wantv)then
      do23974 ll = 1,p 
      l = p-ll+1
      lp1 = l+1
      if(l.le.nrt)then
      if(e(l).ne.0.0d0)then
      do23980 j = lp1,p 
      t = -ddot(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
      call daxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
23980 continue
      continue
      endif
      endif
      do23982 i = 1,p
      v(i,l) = 0.0d0
23982 continue
      continue
      v(l,l) = 1.0d0
23974 continue
      continue
      endif
      mm = m
      iter = 0
23984 continue
      if(m.eq.0)then
      return
      endif
      if(iter.ge.maxit)then
      goto 23986
      endif
      do23991 ll = 1,m 
      l = m-ll
      if(l.eq.0)then
      goto 23992
      endif
      test = dabs(s(l))+dabs(s(l+1))
      ztest = test+dabs(e(l))
      if(ztest.eq.test)then
      go to 150
      endif
23991 continue
23992 continue
      go to 160
150   e(l) = 0.0d0
160   if(l.eq.m-1)then
      kase = 4
      else
      lp1 = l+1
      mp1 = m+1
      do23999 lls = lp1,mp1 
      ls = m-lls+lp1
      if(ls.eq.l)then
      goto 24000
      endif
      test = 0.0d0
      if(ls.ne.m)then
      test = test+dabs(e(ls))
      endif
      if(ls.ne.l+1)then
      test = test+dabs(e(ls-1))
      endif
      ztest = test+dabs(s(ls))
      if(ztest.eq.test)then
      go to 170
      endif
23999 continue
24000 continue
      go to 180
170   s(ls) = 0.0d0
180   if(ls.eq.l)then
      kase = 3
      else
      if(ls.eq.m)then
      kase = 1
      else
      kase = 2
      l = ls
      endif
      endif
      endif
      l = l+1
      I24013=(kase)
      goto 24013
24015 continue
      mm1 = m-1
      f = e(m-1)
      e(m-1) = 0.0d0
      do24016 kk = l,mm1 
      k = mm1-kk+l
      t1 = s(k)
      call drotg(t1,f,cs,sn)
      s(k) = t1
      if(k.ne.l)then
      f = -sn*e(k-1)
      e(k-1) = cs*e(k-1)
      endif
      if(wantv)then
      call drot(p,v(1,k),1,v(1,m),1,cs,sn)
      endif
24016 continue
      continue
      goto 24014
24022 continue
      f = e(l-1)
      e(l-1) = 0.0d0
      do24023 k = l,m 
      t1 = s(k)
      call drotg(t1,f,cs,sn)
      s(k) = t1
      f = -sn*e(k)
      e(k) = cs*e(k)
      if(wantu)then
      call drot(n,u(1,k),1,u(1,l-1),1,cs,sn)
      endif
24023 continue
      continue
      goto 24014
24027 continue
      scale = dmax1(dabs(s(m)),dabs(s(m-1)),dabs(e(m-1)),dabs(s(l)),dabs
     *(e(l)))
      sm = s(m)/scale
      smm1 = s(m-1)/scale
      emm1 = e(m-1)/scale
      sl = s(l)/scale
      el = e(l)/scale
      b = ((smm1+sm)*(smm1-sm)+emm1**2)/2.0d0
      c = (sm*emm1)**2
      shift = 0.0d0
      if(b.ne.0.0d0.or.c.ne.0.0d0)then
      shift = dsqrt(b**2+c)
      if(b.lt.0.0d0)then
      shift = -shift
      endif
      shift = c/(b+shift)
      endif
      f = (sl+sm)*(sl-sm)+shift
      g = sl*el
      mm1 = m-1
      do24032 k = l,mm1 
      call drotg(f,g,cs,sn)
      if(k.ne.l)then
      e(k-1) = f
      endif
      f = cs*s(k)+sn*e(k)
      e(k) = cs*e(k)-sn*s(k)
      g = sn*s(k+1)
      s(k+1) = cs*s(k+1)
      if(wantv)then
      call drot(p,v(1,k),1,v(1,k+1),1,cs,sn)
      endif
      call drotg(f,g,cs,sn)
      s(k) = f
      f = cs*e(k)+sn*s(k+1)
      s(k+1) = -sn*e(k)+cs*s(k+1)
      g = sn*e(k+1)
      e(k+1) = cs*e(k+1)
      if(wantu.and.k.lt.n)then
      call drot(n,u(1,k),1,u(1,k+1),1,cs,sn)
      endif
24032 continue
      continue
      e(m-1) = f
      iter = iter+1
      goto 24014
24040 continue
      if(s(l).lt.0.0d0)then
      s(l) = -s(l)
      if(wantv)then
      call dscal(p,-1.0d0,v(1,l),1)
      endif
      endif
24045 if(l.ne.mm)then
      if(s(l).ge.s(l+1))then
      goto 24046
      endif
      t = s(l)
      s(l) = s(l+1)
      s(l+1) = t
      if(wantv.and.l.lt.p)then
      call dswap(p,v(1,l),1,v(1,l+1),1)
      endif
      if(wantu.and.l.lt.n)then
      call dswap(n,u(1,l),1,u(1,l+1),1)
      endif
      l = l+1
      goto 24045
      endif
24046 continue
      iter = 0
      m = m-1
      goto 24014
24013 continue
      if (I24013.eq.1)goto 24015
      if (I24013.eq.2)goto 24022
      if (I24013.eq.3)goto 24027
      if (I24013.eq.4)goto 24040
24014 continue
      goto 23984
23986 continue
      info = m
      return
      end
      subroutine dbksl(x,p,k,b,q,info)
      integer p,k,q,info
      double precision x(p,p),b(p,q)
      double precision t
      integer j,l
      info = 0
      j=k
24053 if(.not.(j.gt.0))goto 24055
      if(x(j,j).eq.0.0d0)then
      info = j
      goto 24055
      endif
      l=1
24058 if(.not.(l.le.q))goto 24060
      b(j,l) = b(j,l)/x(j,j)
      if(j.ne.1)then
      t = -b(j,l)
      call daxpy(j-1,t,x(1,j),1,b(1,l),1)
      endif
      l = l+1
      goto 24058
24060 continue
      j = j-1
      goto 24053
24055 continue
      return
      end
      subroutine dtrsl(t,ldt,n,b,job,info)
      integer ldt,n,job,info
      double precision t(ldt,*),b(*)
      double precision ddot,temp
      integer which,j,jj
      do24063 info = 1,n
      if(t(info,info).eq.0.0d0)then
      return
      endif
24063 continue
      continue
      info = 0
      which = 1
      if(mod(job,10).ne.0)then
      which = 2
      endif
      if(mod(job,100)/10.ne.0)then
      which = which+2
      endif
      I24071=(which)
      goto 24071
24073 continue
      b(1) = b(1)/t(1,1)
      if(n.ge.2)then
      do24076 j = 2,n 
      temp = -b(j-1)
      call daxpy(n-j+1,temp,t(j,j-1),1,b(j),1)
      b(j) = b(j)/t(j,j)
24076 continue
      continue
      endif
      goto 24072
24078 continue
      b(n) = b(n)/t(n,n)
      if(n.ge.2)then
      do24081 jj = 2,n 
      j = n-jj+1
      temp = -b(j+1)
      call daxpy(j,temp,t(1,j+1),1,b(1),1)
      b(j) = b(j)/t(j,j)
24081 continue
      continue
      endif
      goto 24072
24083 continue
      b(n) = b(n)/t(n,n)
      if(n.ge.2)then
      do24086 jj = 2,n 
      j = n-jj+1
      b(j) = b(j)-ddot(jj-1,t(j+1,j),1,b(j+1),1)
      b(j) = b(j)/t(j,j)
24086 continue
      continue
      endif
      goto 24072
24088 continue
      b(1) = b(1)/t(1,1)
      if(n.ge.2)then
      do24091 j = 2,n 
      b(j) = b(j)-ddot(j-1,t(1,j),1,b(1),1)
      b(j) = b(j)/t(j,j)
24091 continue
      continue
      endif
      goto 24072
24071 continue
      if (I24071.eq.1)goto 24073
      if (I24071.eq.2)goto 24078
      if (I24071.eq.3)goto 24083
      if (I24071.eq.4)goto 24088
24072 continue
      return
      end
