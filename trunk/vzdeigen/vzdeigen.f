ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c vzdeigen - dominant eigenvalue and vector of a linear operator
c
c Original Code by Vlad Zheligovsky
c
c Interface improvements by Manuel Baptista (2005)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     Find the dominant eigenvalue for the operator
c     v1 - Array containing the initial condition and the
c       dominant eigenvector
c     xp - real part of the dominant eigenvalue
c     eim - imaginary part of the dominant eigenvalue
c     ep - tolerance (stopping criterium)
c     M - size of v1
c     v2..v7 - Set of arrays of size M
c     mp - 
c     sc - 
c     nseq - 
c     name - filename for intermediate outputs
      subroutine vzdeigen(v1,xp,eim,ep,thr,M,v2,v3,v4,v5,v6,v7,
     *mp_,sc_,nseq_,name)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter(ncomplex=20,small1=1d-9,thri=.5,
     * nmd=31,ms=3*(2*nmd+1)*(nmd+1)*(nmd+1),
     * small2=1d-12,iter0=1500000,itst=500,init=1,kopt=1)
      DIMENSION V1(M),V2(M),V3(M),V4(M),V5(M),V6(M),V7(M)
      character*100 name
c      common/cc/visc,sc,CF(2,ncomplex),ell,gamma,du(3,12),mp,nseq,nf
      common/cc/sc,CF(2,ncomplex),mp,nseq
      common/uu/qq,z1,z2
c      print *,"M=",M," ",xp,eim,ep,thr,name
      sc=sc_
      mp=mp_
      nseq=nseq_
      call EIGEN(v1,xp,eim,ep,thr,M,name,v2,v3,v4,v5,v6,v7) 
      end

c     Subroutine that interfaces with any externally supplied routine
c     to evaluate the product of the linear operator by the vector
c     NIT - current iteration
c     VI - Array with input vector
c     VO - Array with input vector
c     M - Size of the arrays
      SUBROUTINE PRODX(NIT,VI,VO,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter(ncomplex=20)
      dimension vi(M),vo(M)
      common/cc/sc,CF(2,ncomplex),mp,nseq
      call external_prodx(NIT,VI,VO,M)
      vo=vi+vo/sc
c      print *,sum(VI),sum(VO)
      END

c     Subroutine to load a fortran binary file of doubles
c     into a vector
      subroutine vzdeigen_load_ffile(v,m,name) 
      IMPLICIT REAL*8 (A-H,O-Z)
      character*100 name
      real*8 v(m)
      OPEN(3,STATUS='OLD',FILE=name,
     *     FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*m)
      READ(3,REC=1) v

      X1=0D0
      DO 4 L=M,1,-1
 4    X1=X1+V(L)*V(L)
      print *,"Norm V1: " , X1
      print *,"new conf v1: ",(v(i),i=15,30)



      CLOSE(3)
      end

c     Subroutine to save a vector of doubles into a fortran
c     binary file
      subroutine vzdeigen_save_ffile(v,m,name) 
      IMPLICIT REAL*8 (A-H,O-Z)
      character*100 name
      real*8 v(m)
      character*200 cmd 
c      cmd=ADJUSTL("if [ ! -f ")//TRIM(ADJUSTL(name))//" ]; then  mv "
c     *//TRIM(ADJUSTL(name))//" "//TRIM(ADJUSTL(name))//".old; fi"
c      cmd=TRIM(ADJUSTL(cmd))
c      print *,cmd
c      call system(cmd)
      OPEN(3,STATUS='UNKNOWN',FILE=name,
     *     FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*m)
      WRITE(3,REC=1) v
      CLOSE(3)
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Original Code by Vlad Zheligovsky
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ALBT(V1,V2,V3,M,AL,BT,IK,EP1,ER,ER2)
      IMPLICIT REAL*16 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      real*8 V1(M),V2(M),V3(M),AL,BT,EP1,ER,ER2,qq,z1,z2
      common/uu/qq,z1,z2
      X1=0.
      X2=0.
      X3=0.
      X12=0.
      X23=0.
      X13=0.
      DO 1 L=M,1,-1
      q1=V1(L)
      q2=V2(L)
      q3=V3(L)
      X1=X1+q1*q1
      X2=X2+q2*q2
      X3=X3+q3*q3
      X12=X12+q1*q2
      X23=X23+q3*q2
 1    X13=X13+q1*q3
      EP=(EP1**2)*X1
      z1=x1
      z2=x2
      D=X2*X1-X12*X12
      IF(D.LT.1q-99)GOTO 2
      xAL=(X13*X12-X23*X1)/D
      xBT=(X23*X12-X13*X2)/D
      AL=xAL
      BT=xBT
      Y1=X13+xAL*X12+xBT*X1
      Y2=X23+xAL*X2+xBT*X12
      IF(qABS(Y1).LT.EP.AND.qABS(Y2).LT.EP)GOTO 3
 2    P1=X12/X2
      P2=X23/X3
      Y12=0.
      Y23=0.
      Y13=0.
      Y11=0.
      Y22=0.
      DO 4 L=M,1,-1
      q1=V1(L)
      q2=V2(L)
      q3=V3(L)
      Y1=q1-P1*q2
      Y2=q2-P2*q3
      Y12=Y12+Y1*Y2
      Y11=Y11+Y1*Y1
      Y22=Y22+Y2*Y2
      Y23=Y23+q3*Y2
 4    Y13=Y13+Y1*q3
      D=Y22*Y11-Y12*Y12
      IF(D.LE.0.)GOTO 5
      xAL=(Y13*Y12-Y23*Y11)/D
      xBT=(Y23*Y12-Y13*Y22)/D
      Q=1q0-P2*xAL
      xAL=(xAL-xBT*P1)/Q
      xBT=xBT/Q
      AL=xAL
      BT=xBT
      Y1=X13+xAL*X12+xBT*X1
      Y2=X23+xAL*X2+xBT*X12
      IF(qABS(Y1).GT.EP.OR.qABS(Y2).GT.EP)then
        IK=2
        Q=X23/X2
        xAL=-2q0*Q
        xBT=Q*Q
        AL=xAL
        BT=xBT
        GOTO 5
      endif
 3    IK=1
      xER=0.
      DO 6 L=M,1,-1
      D=V3(L)+xAL*V2(L)+xBT*V1(L)
 6    xER=xER+D*D
      ER=xER
      IF(ER2.LT.-.5D0)RETURN
      Q=X23/X2
 5    xER2=0
      DO 7 L=M,1,-1
      q2=V2(L)
      q3=V3(L)
      D=q3-Q*q2
 7    xER2=xER2+D*D
      ER2=xER2
      qq=q
      RETURN
      END


      SUBROUTINE PROD2(NIT,VI,VO,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter(ncomplex=20)
c      common/cc/visc,sc,CF(2,ncomplex),ell,gamma,du(3,12),mp,nseq,nf
      common/cc/sc,CF(2,ncomplex),mp,nseq

      DIMENSION V(M)
      DIMENSION VI(M),VO(M)
c      print *,"PROD2:M=",M
      CALL PROD(NIT,VI,VO,x,M)
      DO 3 j=1,nseq-1
      X=DSQRT(1d0/X)
      DO 1 L=1,M
 1    v(L)=VO(L)*X
 3    CALL PROD(NIT,V,VO,X,M)
      X=DSQRT(1d0/X)
      DO 6 L=1,M
 6    VO(L)=VO(L)*X
      RETURN
      END


      SUBROUTINE PROD(NIT,VI,VO,x,M)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      dimension vi(M),vo(M)
c      print *,"PROD:M=",M
      CALL PRODX(NIT,VI,VO,M)
      x=0.
      do 1 i=M,1,-1
 1    x=x+VO(i)*VO(i)
      RETURN
      END


      SUBROUTINE EIGEN(v1,xp,eim,ep,thr,M,name,v2,v3,v4,v5,v6,v7)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      parameter(ncomplex=20,small1=1d-9,thri=.5,
     * nmd=31,ms=3*(2*nmd+1)*(nmd+1)*(nmd+1),
     * small2=1d-12,iter0=1500000,itst=500,init=1,kopt=1)
      DIMENSION V1(M),V2(M),V3(M),V4(M),V5(M),V6(M),V7(M)
      character*100 name
c      common/cc/visc,sc,CF(2,ncomplex),ell,gamma,du(3,12),mp,nseq,nf
      common/cc/sc,CF(2,ncomplex),mp,nseq
      common/uu/qq,z1,z2

c      print *,"M=",M," ",xp,eim,ep,thr,name

      iwri=50
      init0=20
      jj=0
      eim=0.
      EPS=dabs(EP/SC)
      EP1=DSQRT(EPS)
      EP2=EPS*EPS
      EP3=min(1d-21,EP2)
 99   nitextr=itst
      write(6,*)'mp=',mp
      NIT=0

      DO 20 N=1,init0
      CALL PROD(NIT,V1,V2,x1,M)
 20   CALL PROD2(NIT,V2,V1,M)

 1    DO 2 N=1,init
      CALL PROD(NIT,V1,V2,x1,M)
 2    CALL PROD2(NIT,V2,V1,M)
      DO 3 N=1,MP
      CALL PROD(NIT,V1,V2,x1,M)
      CALL PROD(NIT,V2,V3,x1,M)
      X1=0D0
      DO 4 L=M,1,-1
      V1(L)=V3(L)+CF(1,N)*V2(L)+V1(L)*CF(2,N)
 4    X1=X1+V1(L)*V1(L)
      X1=DSQRT(1d0/X1)
      DO 3 L=1,M
 3    V1(L)=V1(L)*X1
c     if(init.eq.0.and.mp.eq.0)then
c     X1=0D0
c     DO 74 L=M,1,-1
c74   X1=X1+V1(L)*V1(L)
c     X1=DSQRT(1d0/X1)
c     DO 73 L=1,M
c73   V1(L)=V1(L)*X1
c     endif
      if(jj+iwri.lt.nit)then
        jj=nit

c Modified by Manuel Baptista (Apr2005)
c
c        print *,"some comp v1: ",(v1(i),i=15,30)

        OPEN(3,STATUS='UNKNOWN',FILE=name,
     *    FORM='UNFORMATTED',ACCESS='DIRECT',RECL=8*m)
        WRITE(3,REC=1)v1
        CLOSE(3)        

c        call vzdeigen_save_ffile(v,m,name)
c
c End of modification

      endif

      DO 7 N=1,kopt
      CALL PROD(NIT,V1,V2,X2,M)
      P=0D0
      DO 93 L=M,1,-1
 93   P=P+V1(L)*V2(L)
      EC=0D0
      DO 94 L=M,1,-1
      ER=V2(L)-P*V1(L)
 94   EC=EC+ER*ER
      IF(EC.LE.small2.or.
     * (EC.LE.small1.and.m.eq.ms))EP3=1d-30
      CALL PROD(NIT,V2,V3,X3,M)
      EA=-1D0
      CALL ALBT(V1,V2,V3,M,AL,BT,IK,EP1,EB,EA)
c     write(6,*)al,bt,ec,p
      D=AL*AL/4D0-BT
      XP=(P-1d0)*sc
      write(6,8)nit,EB,EA,D,EC,XP
 8    FORMAT(I6,' ERRORS:',2E15.7,' D=',3E15.7)
C
C   DOMINATING EIGENVALUES FOUND?
C
      IF(NIT.ge.iter0)then
      write(6,*)NIT,' ITERATIONS: STOP'
c     IF(D.lt.0.)then
c       XP=(-AL/2d0-1d0)*SC
c       eim=DSQRT(-D)*SC
c       EP=DSQRT(EB)*SC
c     else
c       EP=DSQRT(EC)*SC
c     endif
c     return
      stop
      endif

 89   if(ik.eq.2)then
      IF(EA.GE.EP2)GOTO 11
      AA=0D0
      XP=(-AL/2d0-1d0)*sc
      GOTO 12
      endif
      IF(D.LT.0.)GOTO 15
c     IF(D+sc*dsqrt(eb).LT.0.)GOTO 15
c     IF(D+EP1+dsqrt(eb).LT.0.)GOTO 15
      EA=0D0
      P=0D0
      DO 13 L=M,1,-1
 13   P=P+V2(L)*V3(L)
      AA=P/DSQRT(X3*X2)
      P=P/X2
      DO 14 L=M,1,-1
      ER=V3(L)-P*V2(L)
 14   EA=EA+ER*ER
      EA=EA/X2
      IF(EA.GE.EP2)IF(D+EP1+dsqrt(eb))15,11,11
      XP=(P-1d0)*sc
 12   EP=DSQRT(EA)*SC
      X3=DSQRT(1d0/X3)
      DO 16 L=1,M
 16   V1(L)=V3(L)*X3
      write(6,17)XP,EP,AA
 17   FORMAT('EIGENVALUE=',E20.12,', ERROR NORM=',E12.5/
     * 'SCALAR PRODUCT=',E30.20)
      return

c     EB=EPS
c     D=0D0
c18   X4=0D0
c     BB=0D0
c     DO 24 L=M,1,-1
c     IF(DABS(V1(L))-EB)25,25,26
c25   V4(L)=0D0
c     BB=DMAX1(BB,DABS(V1(L)))
c     GOTO 24
c26   V4(L)=V1(L)
c     X4=X4+V4(L)*V4(L)
c24   CONTINUE
c     IF(BB.EQ.0D0)GOTO 21
c     CALL PROD(NIT,V4,V2,X2)
c     P=0D0
c     DO 29 L=M,1,-1
c29   P=P+V2(L)*V4(L)
c     AA=P/DSQRT(X4*X2)
c     P=P/X4
c     ER=0D0
c     DO 27 L=M,1,-1
c     EK=V2(L)-P*V4(L)
c27   ER=ER+EK*EK
c     ER=ER/X4
c23   EB=EB/2D0
c     IF(EB.GT.BB)GOTO 23
c     IF(ER.GT.EA)GOTO 18
c     ER=DSQRT(ER)*SC
c     X4=DSQRT(1d0/X4)
c     DO 28 L=1,M
c28   V1(L)=V4(L)*X4
c     XP=(P-1d0)*sc
c     write(6,22)XP,ER,AA
c22   FORMAT('REDUCED DIMENSION: EIGENVALUE=',E20.12,
c    *', ERROR NORM=',E12.5/'SCALAR PRODUCT=',E30.20)
c21   return
C
 15   IF(EB.GT.EP3) IF(D+DSQRT(EB)+EP1)60,60,11
      write(6,*)'D.LT.-EP1-dsqrt(eb)?',d,ep1,eb
      XP=(-AL/2d0-1d0)*SC
      eim=DSQRT(-D)*SC
      EP=DSQRT(EB)*SC
      write(6,32)XP,eim,EP
 32   FORMAT('EIGENVALUE=',e18.12,'+/-',E18.12,
     * 'I, ERROR NORM=',E13.7)

c      if(mp.le.ncomplex.and.m.eq.ms.and.eim.ge.thri)then
c      CALL GENP(nmd,nmd,nmd,v1,v2)
c      mp=mp+1
c      cf(1,mp)=AL
c      cf(2,mp)=BT
c      io=0
c      name(4:4)=char(ichar(name(4:4))+1)
c      goto 99
c      endif

      return

c     EE=EPS
c33   X5=0D0
c     BB=0D0
c     DO 34 L=M,1,-1
c     IF(DABS(V1(L))-EE)35,35,36
c35   V5(L)=0D0
c     BB=DMAX1(BB,DABS(V1(L)))
c     GOTO 34
c36   V5(L)=V1(L)
c     X5=X5+V5(L)*V5(L)
c34   CONTINUE
c     IF(BB.EQ.0D0)GOTO 37
c     CALL PROD(NIT,V5,V3,X3)
c     CALL PROD(NIT,V3,V4,X4)
c     EE2=-1D0
c     CALL ALBT(V5,V3,V4,M,AL,BT,IE,EP1,X5,X3,X4,ER,EE2)
c     ER=ER/X5
c38   EE=EE/2D0
c     IF(EE.GT.BB)GOTO 38
c     IF(ER.GT.EB.OR.IE.EQ.2)GOTO 33
c     EE=DSQRT(ER)*SC
c     XP=(-AL/2D0-1D0)*SC
c     X5=DSQRT(X5)
c     DO 39 L=1,M
C     V2(L)=V3(L)/X5
c39   V1(L)=V5(L)/X5
c     DD=AL*AL/4D0-BT
c     D=DSQRT(DABS(DD))*SC
c     IF(DD.GT.0D0)D=-D
c     write(6,31)XP,D,EE
c31   FORMAT('REDUCED DIMENSION: EIGENVALUE=',E16.8,'+/-',
c    * E16.8,'I, ERROR NORM=',E16.8)
c37   CONTINUE
C     DO 30 L=1,M
C30   WRITE(2,REC=L+3+M)V2(L)
c     GOTO 21
C
C   EXTRAPOLATING VECTOR COORDINATES
C
 11   x3=1d0/dsqrt(x3)
      do 201 l=1,m
 201  v1(l)=v3(l)*x3
      CALL PROD2(NIT,V1,V2,M)
      CALL PROD2(NIT,V2,V3,M)
      CALL PROD2(NIT,V3,V4,M)
c: x2=x3=x4=1
      D=0D0
      DO 40 L=M,1,-1
 40   D=V3(L)*V4(L)+D
c     D=D/X3
      ER=0D0
      DO 41 L=M,1,-1
      P=V4(L)-D*V3(L)
 41   ER=ER+P*P
c     ER=ER/X3
      CALL PROD2(NIT,V4,V5,M)
      CALL PROD2(NIT,V5,V3,M)
      CALL PROD2(NIT,V3,V2,M)
      CALL PROD2(NIT,V2,V1,M)
c: x1=x2=x3=x4=x5=1
      Y1=0D0
      Y2=0D0
      Y3=0D0
      Y5=0D0
      DO 42 L=M,1,-1
      Y5=V4(L)*V5(L)+Y5
      Y3=V5(L)*V3(L)+Y3
      Y2=V3(L)*V2(L)+Y2
 42   Y1=V2(L)*V1(L)+Y1
c     Y1=Y1/X2
c     Y2=Y2/X3
c     Y3=Y3/X5
c     Y5=Y5/X4
      EA=0D0
      EB=0D0
      EE=0D0
      EK=0D0
      DO 43 L=M,1,-1
      P=V5(L)-Y5*V4(L)
      EA=EA+P*P
      P=V3(L)-Y3*V5(L)
      EB=EB+P*P
      P=V2(L)-Y2*V3(L)
      EE=EE+P*P
      P=V1(L)-Y1*V2(L)
 43   EK=EK+P*P
c     EA=EA/X4
c     EB=EB/X5
c     EE=EE/X3
c     EK=EK/X2
      IF(er.le.ea.or.ea.le.eb.or.eb.le.ee.or.ee.le.ek)
     * nitextr=nitextr+itst
      IF(nitextr.gt.nit)goto 51
      EE2=(EA+EB+ER+EE+EK)/5D0
      ER=DSQRT(ER)
      EA=DSQRT(EA)
      EB=DSQRT(EB)
      EE=DSQRT(EE)
      EK=DSQRT(EK)
      AA=(EA+EB+ER+EE+EK)/5D0
      D=(EA-AA)**2+(EB-AA)**2+(EE-AA)**2+(EK-AA)**2+(ER-AA)**2
      DD=DMAX1(EA,EB,EE,EK,ER)
      P4=(EE2-ER*AA)/D
      P5=(EE2-EA*AA)/D
      P3=(EE2-EB*AA)/D
      P2=(EE2-EE*AA)/D
      P1=(EE2-EK*AA)/D
      ER=(ER-AA)/D
      EA=(EA-AA)/D
      EB=(EB-AA)/D
      EE=(EE-AA)/D
      EK=(EK-AA)/D
      X3=0D0
      DO 44 L=M,1,-1
      P=P5*V5(L)+P4*V4(L)+P3*V3(L)+P2*V2(L)+P1*V1(L)
      V2(L)=EA*V5(L)+ER*V4(L)+EB*V3(L)+EE*V2(L)+EK*V1(L)
      V3(L)=P
 44   X3=X3+P*P
      CALL PROD(NIT,V3,V4,X4,M)
      CALL PROD(NIT,V2,V5,X5,M)
c: x1=x2=1 still
      AL=0D0
      BT=0D0
      P=0D0
      DO 45 L=M,1,-1
      BT=BT+V3(L)*V2(L)
      AL=AL+V3(L)*V4(L)
 45   P=P+V2(L)*V4(L)
      AL=AL/X3
      P=(AL*BT-P)/X3
      EA=0D0
      EB=0D0
      EE=0D0
      DO 46 L=M,1,-1
      P1=V5(L)-AL*V2(L)
      P2=P1+P*V3(L)
      P3=V4(L)-AL*V3(L)
      EA=EA+P1*P3
      EB=EB+P2*P2
 46   EE=EE+V2(L)*P3
      EB=EB+EE*(2D0*P-EE/X3)
      IF(DABS(EA).LT.2D0*DD*DABS(EB))THEN
      EB=-EA/EB
      X3=0D0
      X4=0D0
      DO 47 L=M,1,-1
      V3(L)=EB*V2(L)+V3(L)
      X3=X3+V3(L)*V3(L)
      V4(L)=EB*V5(L)+V4(L)
 47   X4=X4+V4(L)*V4(L)
      ENDIF
      CALL PROD(NIT,V4,V5,X5,M)
      EK2=0D0
      EE2=0D0
      CALL ALBT(V3,V4,V5,M,AK,BK,IK,EP1,EK,EK2)
      EK=EK/X3
      EK2=EK2/X3
      CALL PROD(NIT,V1,V2,X2,M)
      CALL PROD(NIT,V2,V3,X3,M)
c: x1=1 still
      CALL ALBT(V1,V2,V3,M,AL,BT,IE,EP1,EE,EE2)
C     IF(EK2.GE.EE2*THR.OR.(IE.EQ.1.AND.AL*AL.LT.BT*4D0.AND.
C    * (IK.NE.1.OR.EK.GE.EE*THR.OR.AK*AK.GE.BK*4D0)))GOTO 51
      IF(EK2.GE.EE2*THR.OR.(IE.EQ.1.AND.AL*AL.LT.BT*4D0))GOTO 51
      nitextr=nit+itst
      write(6,48)EK2,EE2
 48   FORMAT('NEW VECTOR (EXTRAPOLATED):',E14.8,'/',E14.8)
 49   X5=DSQRT(1d0/X5)
      DO 50 L=1,M
 50   V1(L)=V5(L)*X5
      GOTO 7
C
C    TWO MINOR EIGENVECTORS FOR D UNKNOWN ELIMINATION FOLLOWS
C
 51   CALL PROD2(NIT,V1,V2,M)
      CALL PROD2(NIT,V2,V3,M)
      CALL PROD2(NIT,V3,V4,M)
      CALL PROD2(NIT,V4,V5,M)
c: x1=x2=x3=x4=x5=1
      P=0D0
      DO 52 L=M,1,-1
 52   P=V5(L)*V4(L)+P
c     P=P/X4
      X1=0D0
      X2=0D0
      X3=0D0
      DO 53 L=M,1,-1
      BB=V3(L)-P*V2(L)
      V1(L)=BB-P*(V2(L)-P*V1(L))
      X1=X1+V1(L)*V1(L)
      AA=V4(L)-P*V3(L)
      V2(L)=AA-P*BB
      X2=X2+V2(L)*V2(L)
      V3(L)=(V5(L)-P*V4(L))-P*AA
 53   X3=X3+V3(L)*V3(L)
      EE2=-1D0
      CALL ALBT(V1,V2,V3,M,AA,BB,IE,EP1,EE,EE2)
      CALL PROD2(NIT,V5,V1,M)
      CALL PROD2(NIT,V1,V2,M)
      CALL PROD2(NIT,V2,V3,M)
c: x1=x2=x3=1
      X3=1D0
      X4=0D0
      X5=0D0
      X6=0D0
      DO 54 L=M,1,-1
      V4(L)=V1(L)+AA*V5(L)+BB*V4(L)
      X4=X4+V4(L)*V4(L)
      V6(L)=V2(L)+AA*V1(L)+BB*V5(L)
      X6=X6+V6(L)*V6(L)
      V5(L)=V3(L)+AA*V2(L)+BB*V1(L)
 54   X5=X5+V5(L)*V5(L)
      EK2=0D0
      EE2=0D0
      CALL ALBT(V4,V6,V5,M,AK,BK,IK,EP1,EK,EK2)
      EK=EK/X5
      EK2=EK2/X5
      CALL ALBT(V1,V2,V3,M,AL,BT,IE,EP1,EE,EE2)
C     IF(EK2.GE.EE2*THR.OR.(IE.EQ.1.AND.AL*AL.LT.BT*4D0.AND.
C    * (IK.NE.1.OR.EK.GE.EE*THR.OR.AK*AK.GE.BK*4D0)))GOTO(60,66),IE
      IF(EK2.GE.EE2*THR.OR.(IE.EQ.1.AND.AL*AL.LT.BT*4D0))GOTO(60,66),IE
      IF(nitextr.gt.nit)GOTO(60,66),IE
      nitextr=nit+itst
      write(6,55)EK2,EE2,EK,EE,AA,BB,P
 55   FORMAT('NEW VECTOR (D - UNKNOWN):',
     * E14.8,3('/',E14.8)/'COEFF:',3E15.8)
      GOTO 49
C
C    TWO MINOR EIGENVECTORS FOR D<0 ELIMINATION FOLLOWS
C
 60   x3=dsqrt(1d0/x3)
      do 202 l=1,m
 202  v1(l)=v3(l)*x3
      CALL PROD2(NIT,V1,V2,M)
      CALL PROD2(NIT,V2,V3,M)
      CALL PROD2(NIT,V3,V4,M)
      CALL ALBT(V1,V2,V3,M,AL,BT,IK,EP1,EK,EK2)
      CALL PROD2(NIT,V4,V5,M)
      CALL PROD2(NIT,V5,V6,M)
      CALL PROD2(NIT,V6,V7,M)
c: x4=x5=x6=x7=1
      X1=0D0
      X2=0D0
      X3=0D0
      DO 61 L=M,1,-1
      V1(L)=V3(L)+AL*V2(L)+BT*V1(L)
      V2(L)=V4(L)+AL*V3(L)+BT*V2(L)
      V3(L)=V5(L)+AL*V4(L)+BT*V3(L)
      V4(L)=V6(L)+AL*V5(L)+BT*V4(L)
      V5(L)=V7(L)+AL*V6(L)+BT*V5(L)
      V1(L)=V3(L)+AL*V2(L)+BT*V1(L)
      X1=X1+V1(L)*V1(L)
      V2(L)=V4(L)+AL*V3(L)+BT*V2(L)
      X2=X2+V2(L)*V2(L)
      V3(L)=V5(L)+AL*V4(L)+BT*V3(L)
 61   X3=X3+V3(L)*V3(L)
      EK2=-1D0
      CALL ALBT(V1,V2,V3,M,AL,BT,IK,EP1,EK,EK2)
      CALL PROD2(NIT,V7,V1,M)
      CALL PROD2(NIT,V1,V2,M)
      CALL PROD2(NIT,V2,V3,M)
c: x1=x2=x3=1
      EE2=0D0
      CALL ALBT(V1,V2,V3,M,AA,BB,IE,EP1,EE,EE2)
      X4=0D0
      X5=0D0
      X1=0D0
      DO 62 L=M,1,-1
      V4(L)=V1(L)+AL*V7(L)+BT*V6(L)
      X4=X4+V4(L)*V4(L)
      V5(L)=V2(L)+AL*V1(L)+BT*V7(L)
      X5=X5+V5(L)*V5(L)
      V1(L)=V3(L)+AL*V2(L)+BT*V1(L)
 62   X1=X1+V1(L)*V1(L)
      EK2=0D0
c: x2=x3=x6=x7=1
      CALL ALBT(V4,V5,V1,M,AK,BK,IK,EP1,EK,EK2)
      EK=EK/X1
      EK2=EK2/X1
      IF(IE.EQ.2.OR.AA*AA.GT.BB*4D0)THEN
      IF(EK2.GE.EE2*THR.OR.(IK.EQ.1.AND.AK*AK.LE.BK*4D0))GOTO 66
      ELSE
c     IF(IK.EQ.2.OR.AK*AK.GE.BK*4D0.OR.EK.GE.EE*THR)GOTO 66
      IF(IK.EQ.2.OR.AK*AK.GE.BK*4D0.OR.EK.GE.EE*THR.OR.EK2.GE.EE2*THR)
     * GOTO 66
      ENDIF
      IF(nitextr.gt.nit)GOTO 66
      nitextr=nit+itst
      write(6,63)EK2,EE2,EK,EE,AL,BT
 63   FORMAT('NEW VECTOR (D<0):',E14.8,3('/',E14.8)/'COEFF:',2E15.8)
      X1=DSQRT(1d0/X1)
      DO 65 L=1,M
 65   V1(L)=V1(L)*X1
      GOTO 7
 66   DO 67 L=1,M
 67   V1(L)=V3(L)
 7    CONTINUE
      GOTO 1
      END
