c  Main program and subroutines used to calibrate EOS when adding
c  a new element.
C  ATOMIC UNITS:
C  density:        11.206*A g/cc
C  temperature:    27.21 eV
C  pressure:       2.942E14 erg/cc
c**********************************************************************

      PROGRAM EOSADD
C      SUBROUTINE EOSADD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/PRS/GYC,DLYDLT,G,GP,YC,YW,PE,PN,PS
      COMMON/EOSCA1/BULK0,D00,T00,P00,E00,THU
      DIMENSION TACRIT(7),PACRIT(7),DPACRI(7)
      EXTERNAL FZERGA,FZERHU

      AA=10.82D0
      ZZ=5.D0

c      AA=183.85D0
c      ZZ=74.D0
c      AA=197.D0
c      ZZ=79.D0
      CALL LAUNCH(AA,ZZ)

c      ANU=0.004D0
c      DEL=.07D0
c: solve for the value of GAM:
c  BULK0 = dP/dln(rho) = bulk modulus in a.u.:
c      BULK0=0.0106D0
c       BULK0=0.007D0
c      GAM=ZERAG(FZERGA,.1D0,1.D1,1.D0,1.D-4)

C      GAM=1.526D0

c: bulk modulus:
      D=DZ
      BULK=.5D0*(BARA(1.01D0*D,0.D0,0.D0)-BARA(D/1.01D0,0.D0,0.D0))/
     /LOG(1.01D0)
c: cohesive energy:
      NPOI=200
      D00=DZ
      ECOHE=0.D0
      HD=D00/NPOI
      DO 1100 K=1,NPOI
      D=0.5D0*HD+(K-1)*HD
 1100 ECOHE=ECOHE-HD*BARA(D,0.D0,0.D0)/D**2

      PRINT 9010,GAM,DEL,BULK,ECOHE

c      GOTO 2600

c: construct Hugoniot:
      D00=2.34D0/11.206D0/AA
      T00=300.D0/11605.D0/27.21D0
      P00=BARA(D00,T00,T00)
      E00=ENA(D00,T00,T00)
      UUNIT=SQRT(2.942D4/11.206/AA)
      PRINT 9020,ZZ,AA,D00*11.206D0*AA,T00*27.21D0
      DO 2500 K=1,21
      THU=10.D0**(.1D0*(K-11))/27.21D0
      RAB=D00
      DHU=ZERAG(FZERHU,RAB,10.D0*RAB,2.D0*RAB,1.D-4)
      PHU=BARA(DHU,THU,THU)
      DSHHU=SQRT((PHU-P00)*DHU/(DHU-D00)/D00)
      UHU=(PHU-P00)/DSHHU/D00
      PRINT 9030,THU*27.21D0,DHU*AA*11.206D0,PHU*294.2D0
     *,UHU*UUNIT,DSHHU*UUNIT,4.015D0+1.252D0*UHU*UUNIT
 2500 CONTINUE

c: critical point:
 2600 CONTINUE
      DO 2650 K=1,7
      PACRIT(K)=0.D0
 2650 TACRIT(K)=.3D0+.05D0*(K-1)
      PRINT 9040,TACRIT
      DO 2800 J=1,95
      D=.01D0*DZ*J
      DO 2700 K=1,7
      T=TACRIT(K)/27.21D0
      RAB=BARA(D,T,T)*294.2D0
      DPACRI(K)=RAB-PACRIT(K)
 2700 PACRIT(K)=RAB
      PRINT 9050,DPACRI
 2800 CONTINUE

 9010 FORMAT(' gam,del=',0P2F8.3,' B-mod.(a.u.)=',1PE12.4
     *,' E-coh.(a.u.)=',1PE12.4)
 9020 FORMAT(' HUGONIOT for  Z,A=',0P2F8.3,' ro_0=',1PE12.4,
     *'g/cc, T_0=',1PE12.4,'eV'/
     *'  T(eV)       ro(g/cc)    P(Mbar)     u(km/s)     D(km/s)'
     *,5X,'D_exp(km/s)')
 9030 FORMAT(1P6E12.4)
 9040 FORMAT(' Incr.-of-press. table (Mbar) around critical point'
     *' for temperature columns (eV)'/1P7E11.3/)
 9050 FORMAT(1P7E11.3)
      STOP
      END

C*****************************************************************
      DOUBLE PRECISION FUNCTION FZERGA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/EOSCA1/BULK0,D00,T00,P00,E00,THU
c: bulk modulus:
      GAM=X
      D=DZ
      FZERGA=.5D0*(BARA(1.01D0*D,0.D0,0.D0)-BARA(D/1.01D0,0.D0,0.D0))/
     /LOG(1.01D0)-BULK0
      RETURN
      END

C*****************************************************************
      DOUBLE PRECISION FUNCTION FZERHU(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/EOSCA1/BULK0,D00,T00,P00,E00,THU
c: bulk modulus:
      D=X
      T=THU
      FZERHU=ENA(D,T,T)-E00-.5D0*(BARA(D,T,T)+P00)*
     *(1.D0/D00-1.D0/D)
      RETURN
      END

C*****************************************************************

C  ZERAG = the root X of equation  F(X)=0  in the interval XL<X<XR;
C  error in the value of X is less than EPS*(XR-XL);
      DOUBLE PRECISION FUNCTION ZERAG(F,XL,XR,X0,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F

      NRCMAX=200
      X=X0
      H=.2D0*(XR-XL)
      HMIN=EPS*(XR-XL)
      HMMIN=.1D0*HMIN
      FP=F(X)
      IF(FP.EQ.0.D0) GOTO 77
      NRCAL=1
      KF=0
      KFR=0
      KFL=0
      XP=X

 10   X=X+H
      IF(X.LT.XR) GOTO 14
      KFR=KFR+1
      X=.5D0*(XR+XP)
      GOTO 19
 14   KFR=0
      IF(X.GT.XL) GOTO 19
      KFL=KFL+1
      X=.5D0*(XL+XP)
      GOTO 20
 19   KFL=0

 20   F2=F(X)
      NRCAL=NRCAL+1
      IF(F2.EQ.0.D0) GOTO 77
      F2SGN=F2/DABS(F2)
      IF(FP*F2SGN.GT.0.D0) GOTO 30
      IF(DABS(H).LT.HMIN) GOTO 77
      H=-.3D0*H
      GOTO 40

 30   IF(DABS(FP).GT.DABS(F2)) GOTO 40
      IF(KF.EQ.1) GOTO 33
      KF=1
      X=XP
      H=-H
      GOTO 10

 33   KF=0
      X=XP
      H=.3D0*H
      GOTO 50

 40   XP=X
      KF=0
      FP=F2
C: check for abnormal termination:
 50   IF(DABS(H).LT.HMMIN) GOTO 66
      IF(NRCAL.GT.NRCMAX) GOTO 66
      IF(XR-X.LT.HMMIN.AND.KFR.GT.2) GOTO 66
      IF(X-XL.LT.HMMIN.AND.KFL.GT.2) GOTO 66
      GOTO 10

 66   PRINT 511,NRCAL,XL,X,XR,H,FP,F2,EPS
 77   ZERAG=X
      RETURN

 511  FORMAT(' F(X).NE.0 in ZERAG: NRCAL=',I6/' XL,X,XR=',1P3E15.8/
     *' H=',1PE15.8,' FP,F2=',1P2E15.8,' EPS=',1PE15.8)
      END
