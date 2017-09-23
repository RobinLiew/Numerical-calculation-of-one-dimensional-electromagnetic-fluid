C                   ************************
C                   *                      *
C                   *    T A B I N - 3     *
C                   *                      *
C**********************************************************************
C   Version: 3.1   Date: 22 Nov 06   Author: M.Basko
C----------------------------------------------------------------------
C  TABIN3 is a package for filling EOS and BEAM tables used in DEIRA3.
C  Presently, equation-of-state (EOS) data are available for elements
C  with Z=1,2,3,4,5,6,7,8,11,13,18,22,26,29,42,47,55,74,79,82,83,92;
C  tentative EOS for Z=5,7,8,54; (Z=5 added on 21.01.01);
C  tentative stopping (NESHL and ESHL values) for Z=5,74.
C  Differs from TABIN2: 1) opacity model of 1994; 2) polished BLOCK
C  DATA, and LAEOS, LZEOS added to c/blk /ELEMNT/.
C**********************************************************************
C                  I N S T R U C T I O N S                            *
C**********************************************************************
C  TO FILL EOS AND BEAM TABLES:                                       *
C  1) edit block `Initial Info' to set initial data for each          *
C     of the 5 elements that can be accomodated in the table.         *
C  2) edit lines after   C: fill BEAM tables:                         *
C**********************************************************************

      PROGRAM TABIN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL FARR

C: c/blks for EOS and opacity:
      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/OPAC1/POTOFI(100,5),BOPAC(100,5)
      COMMON/OPAC2/NLAGER,NWWWW,XLAGER(12),WLAGER(12),RWEI(12),PWEI(12)
      COMMON/EPZER/P00(5),E00(5)
      COMMON/URSEDG/NROO(6),ROILL(5),HLROO(5)
     *,NTEMM(3,6),TEMILL(3,5),HLTEMM(3,5)
      COMMON/URSTBL/FARR(30000)
      COMMON/URTBLI/QURS(1000)
C  One should watch that
C  3*(NROO(1)*NTEMM(3,1)+...+NROO(5)*NTEMM(3,5)) < 30000, and
C  4*(NROO(1)+...+NROO(5)) < 1000.

C: c/blks for beam:
      COMMON/BECNST/CBEAM(8)
      COMMON/BEAM1/ABEAM,ZBEAM,CBEAM1(406)
      COMMON/BEAM2/CBEAM2(350),NCBEM2(106),CBEM22(125)

      DIMENSION NRO(5),NTEM(5),ROMIN(5),ROMAX(5),TEMIN(5),TEMAX(5)

C--------------------------------------------------- begin Initial Info
C......................................................................
C Initial information about each element:                             .
C     ZSUB(K)     atomic number of K-th substance;                    .
C     ASUB(K)     atomic weight of K-th substance;                    .
C     NRO(K) >1   number of table grid nodes along the ln(D) axis;    .
C     ROMIN(K)    minimum density in g/cm**3;                         .
C     ROMAX(K)    maximum density in g/cm**3;                         .
C     NTEM(K) >3  number of table grid nodes along the ln(TE) axis;   .
C     TEMIN(K)    minimum electron temperature in keV;                .
C     TEMAX(K)    maximum electron temperature in keV;                .
C                                                                     .
C Comments:                                                           .
C 1. Dependence of PI and EI on ion temperature TI is calculated      .
C    in URSAPB analytically for any TI.GE.0.D0                        .
C                                                                     .
C 2. Set ZSUB(K)=0.D0 for K.GE.I if substances with K.GE.I are not used
C                                                                     .
C 3. For a decent accuracy, table parameters should be chosen such    .
C    that the grid steps in ln(D) and ln(T) do not exceed 0.5.        .
C......................................................................

      INQUIRE(FILE='EOSTA3',EXIST=LEXIST)
      IF(LEXIST) THEN
        WRITE(0,9010)
        PRINT 9010
        GOTO 999
      ENDIF

      INQUIRE(FILE='BEMTA3',EXIST=LEXIST)
      IF(LEXIST) THEN
        WRITE(0,9020)
        PRINT 9020
        GOTO 999
      ENDIF

      ZSUB(1)=1.D0   ***************燃料层的初始化信息
      ASUB(1)=2.D0
      NRO(1)=35
      ROMIN(1)=1.D-4
      ROMAX(1)=1.D+4
      NTEM(1)=41
      TEMIN(1)=3.D-7
      TEMAX(1)=1.D+1

      ZSUB(2)=4.D0
      ASUB(2)=9.D0
      NRO(2)=35
      ROMIN(2)=1.D-4
      ROMAX(2)=1.D+3
      NTEM(2)=45
      TEMIN(2)=3.D-7
      TEMAX(2)=3.D+1

      ZSUB(3)=79.D0
      ASUB(3)=197.D0
      NRO(3)=35
      ROMIN(3)=1.D-4
      ROMAX(3)=1.D+4
      NTEM(3)=45
      TEMIN(3)=3.D-7
      TEMAX(3)=1.D2

      ZSUB(4)=82.D0
      ASUB(4)=207.2D0
      NRO(4)=35
      ROMIN(4)=1.D-4
      ROMAX(4)=1.D+4
      NTEM(4)=45
      TEMIN(4)=3.D-7
      TEMAX(4)=1.D+2

      ZSUB(5)=6.D0
      ASUB(5)=12.D0
      NRO(5)=37
      ROMIN(5)=1.D-5
      ROMAX(5)=3.D+3
      NTEM(5)=45
      TEMIN(5)=3.D-7
      TEMAX(5)=1.D+2       ********************************

C===================================================== end Initial Info

C: parameters of the table grid in atomic units:
C......................................................................
C     DTLL  is the halfwidth of ln(T) interval where half             .
C           of the regular temperature increment H is used;           .
C     TLL0  is the middle of this interval;                           .
C......................................................................
      DTLL=2.D0
      DO 80 K=1,5
      IF(ZSUB(K).LT..5D0) GOTO 80
      IF(NRO(K).LT.2) GOTO 410
      IF(NTEM(K).LT.4) GOTO 420
      NROO(K)=NRO(K)
      ROILL(K)=DLOG(ROMIN(K)/(11.206D0*ASUB(K)))
      HLROO(K)=DLOG(ROMAX(K)/ROMIN(K))/(NRO(K)-1)
      TLL0=1.3333D0*DLOG(1.D0+ZSUB(K))-2.4D0
      TLLL=DLOG(TEMIN(K)/.02721D0)
      TLLR=DLOG(TEMAX(K)/.02721D0)
      TEMILL(1,K)=TLLL
      NTEMM(3,K)=NTEM(K)

      IF(TLLR.LE.TLL0-DTLL.OR.TLLL.GE.TLL0+DTLL.OR.
     *(TLLL.GE.TLL0-DTLL.AND.TLLR.LE.TLL0+DTLL)) GOTO 18
      GOTO 20
 18   HHLT=(TLLR-TLLL)/(NTEMM(3,K)-1)
      NTEMM(1,K)=2
      NTEMM(2,K)=3
      HLTEMM(1,K)=HHLT
      HLTEMM(2,K)=HHLT
      HLTEMM(3,K)=HHLT
      GOTO 80
 20   IF(TLLL.LT.TLL0-DTLL) GOTO 30
      HHLT=(TLLR-2.D0*TLLL+TLL0+DTLL)/(NTEMM(3,K)-1)
      NRAB=(TLLR-TLL0-DTLL)/HHLT
      NTEMM(2,K)=NTEMM(3,K)-NRAB
      IF(NTEMM(2,K).GT.NTEMM(3,K)-1) NTEMM(2,K)=NTEMM(3,K)-1
      IF(NTEMM(2,K).LT.3) NTEMM(2,K)=3
      NTEMM(1,K)=NTEMM(2,K)-1
      HLTEMM(1,K)=(TLLR-TLLL)/(2*NTEMM(3,K)-NTEMM(2,K)-1)
      HLTEMM(2,K)=HLTEMM(1,K)
      HLTEMM(3,K)=2.D0*HLTEMM(1,K)
      GOTO 80
 30   IF(TLLR.GT.TLL0+DTLL) GOTO 40
      HHLT=(2.D0*TLLR-TLLL-TLL0+DTLL)/(NTEMM(3,K)-1)
      NRAB=(TLL0-DTLL-TLLL)/HHLT
      NTEMM(1,K)=NRAB+1
      IF(NTEMM(1,K).LT.2) NTEMM(1,K)=2
      IF(NTEMM(1,K).GT.NTEMM(3,K)-2) NTEMM(1,K)=NTEMM(3,K)-2
      NTEMM(2,K)=NTEMM(1,K)+1
      HLTEMM(3,K)=(TLLR-TLLL)/(NTEMM(3,K)+NTEMM(1,K)-2)
      HLTEMM(2,K)=HLTEMM(3,K)
      HLTEMM(1,K)=2.D0*HLTEMM(3,K)
      GOTO 80
 40   HHLT=(TLLR-TLLL+2.D0*DTLL)/(NTEMM(3,K)-1)
      NRAB=(TLL0-DTLL-TLLL)/HHLT
      NTEMM(1,K)=NRAB+1
      IF(NTEMM(1,K).LT.2) NTEMM(1,K)=2
      IF(NTEMM(1,K).GT.NTEMM(3,K)-2) NTEMM(1,K)=NTEMM(3,K)-2
      NRAB=(TLLR-TLL0-DTLL)/HHLT
      NTEMM(2,K)=NTEMM(3,K)-NRAB
      IF(NTEMM(2,K).GT.NTEMM(3,K)-1) NTEMM(2,K)=NTEMM(3,K)-1
      IF(NTEMM(2,K).LT.NTEMM(1,K)+1) NTEMM(2,K)=NTEMM(1,K)+1
      HLTEMM(2,K)=(TLLR-TLLL)/(NTEMM(1,K)-NTEMM(2,K)+2*NTEMM(3,K)-2)
      HLTEMM(1,K)=2.D0*HLTEMM(2,K)
      HLTEMM(3,K)=HLTEMM(1,K)
 80   CONTINUE

C: fill in EOS tables:
      DO 110 K=1,5
      ZZ=ZSUB(K)
      IF(ZZ.LT..5D0) GOTO 110
      AA=ASUB(K)
      CALL LAUNCH(AA,ZZ)
  100 CALL INUR2(K)
  110 CONTINUE
      CALL INOPAC

      OPEN(11,FILE='EOSTA3',STATUS='NEW',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
      WRITE(11,9999) ASUB,ZSUB,P00,E00
      WRITE(11,9998) NROO
      WRITE(11,9999) ROILL,HLROO
      WRITE(11,9998) NTEMM
      WRITE(11,9999) TEMILL,HLTEMM,FARR,QURS
      WRITE(11,9998) NLAGER,NWWWW
      WRITE(11,9999) POTOFI,BOPAC,XLAGER,WLAGER,RWEI,PWEI
      CLOSE(11)

C: fill in BEAM tables:
      ZB=83.D0
      AB=209.D0
      CALL BELAUN(AB,ZB)

      OPEN(12,FILE='BEMTA3',STATUS='NEW',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
      WRITE(12,9999) CBEAM,ABEAM,ZBEAM,CBEAM1,CBEAM2
      WRITE(12,9998) NCBEM2
      WRITE(12,9999) CBEM22
      CLOSE(12)
 200  CONTINUE

C: tests:
      OPEN(11,FILE='EOSTA3',STATUS='OLD',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
      READ(11,9999) ASUB,ZSUB,P00,E00
      READ(11,9998) NROO
      READ(11,9999) ROILL,HLROO
      READ(11,9998) NTEMM
      READ(11,9999) TEMILL,HLTEMM,FARR,QURS
      READ(11,9998) NLAGER,NWWWW
      READ(11,9999) POTOFI,BOPAC,XLAGER,WLAGER,RWEI,PWEI
      CLOSE(11)
      CALL URTEST

      DO 280 K=1,5
      OPEN(12,FILE='BEMTA3',STATUS='OLD',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
      READ(12,9999) CBEAM,ABEAM,ZBEAM,CBEAM1,CBEAM2
      READ(12,9998) NCBEM2
      READ(12,9999) CBEM22
      CLOSE(12)
 280  CALL BETEST(K)
 300  GOTO 999
 410  PRINT 9110,NRO(K),K
      GOTO 999
 420  PRINT 9120,NTEM(K),K
 
 999  WRITE(0,9090)
      STOP

 9010 FORMAT(/'!!! STOP in TABIN3: file EOSTA3 already exists !!!'/)
 9020 FORMAT(/'!!! STOP in TABIN3: file BEMTA3 already exists !!!'/)
 9090 FORMAT(/'* * *     End of the TABIBN3 run    * * *'
     &//' Press any key to continue')
 9110 FORMAT(/' NRO=',I8,' is less than 2 for substance no.',I2/
     *' NO TABLES FILLED'/)
 9120 FORMAT(/' NTEM=',I8,' is less than 4 for substance no.',I2/
     *' NO TABLES FILLED'/)
 9998 FORMAT(6I10)
 9999 FORMAT(1P5E15.8)
      END
C****************************************************

      SUBROUTINE URTEST
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/UR2OUT/FE(12),FI(6)
      COMMON/PRS/GDIN(6),PE,PN,PS
      COMMON/ENRG/EI,EE(4),EETO,PTO
      DIMENSION RAB(19)

      DO 300 K=1,5
      IF(ZSUB(K).LT..5D0) GOTO 300
      AA=ASUB(K)
      ZZ=ZSUB(K)
      CALL LAUNCH(AA,ZZ)
      PRINT 9601,K

      DO 290 KK=1,5
      PRINT 9602
      DO 290 J=1,20
      RAB(1)=-7.5D0+.5D0*J
      T=10.D0**(RAB(1))/27.21D-3
      DO 250 L=2,19
      D=10.D0**(-5.0D0+.5D0*L)/11.206D0/ASUB(K)
      CALL URSAPB(K,D,T,T,3)
      QA=FE(KK)
      IF(KK.GE.4) QA=FI(KK-3)
      IF(KK.EQ.1) Q=DINA(D,T)
      IF(KK.EQ.2) Q=BARA(D,T,0.D0)
      IF(KK.EQ.3) Q=ENA(D,T,0.D0)
      IF(KK.NE.4) GOTO 200
      Q=BARA(D,T,T)
      Q=PN
      GOTO 230
 200  IF(KK.NE.5) GOTO 230
      Q=ENA(D,T,T)
      Q=EI
 230  RAB(L)=1.D18
      IF(ABS(Q).GT.1.D-18*ABS(QA)) RAB(L)=QA/Q-1.D0
 250  CONTINUE
      PRINT 9604,(RAB(I),I=1,10)
      PRINT 9605,(RAB(I),I=11,19)
 290  CONTINUE
 300  CONTINUE
      RETURN

 9601 FORMAT(///' Tables of relative errors in Y,PE,EE,PI,EI for'
     *,' substance',I3)
 9602 FORMAT(/9X,'Horiz.axis: lg(D) [g/cc] = -4.0,-3.5, ..., +4.5'/
     *9X,'Vertic.axis: lg(T) [keV] = -7.0,-6.5, ..., +2.5'/' lg(T)')
 9604 FORMAT(1X,F6.2,1P9E8.0)
 9605 FORMAT(7X,1P9E8.0)
      END
C****************************************************

      SUBROUTINE BETEST(K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/BECNST/CBEAM(8)
      COMMON/BEAM1/ABEAM,ZBEAM,CBEAM1(406)
      DIMENSION RAB(3,3,3)

      IF(ZSUB(K).LT..5D0) GOTO 500
      PRINT 9611,K

      DO 400 KE=1,3
      E1=5.D0+2.5D0*KE*(KE-1)
      DO 400 KD=1,3
      D=10.D0**(2*KD-3)
      DO 400 KT=1,3
      T=10.D0**(2*KT-5)
 400  RAB(KT,KD,KE)=BESTO2(E1,D,T,T,K)

      AA1=ABEAM
      ZZ1=ZBEAM
      CALL BELAUN(AA1,ZZ1)
      DO 200 KE=1,3
      E1=(5.D0+2.5D0*KE*(KE-1))/(4.96D-5*ABEAM)
      DO 200 KD=1,3
      D=10.D0**(2*KD-3)/(11.206D0*ASUB(K))
      DO 200 KT=1,3
      T=10.D0**(2*KT-2)/27.21D0
 200  RAB(KT,KD,KE)=RAB(KT,KD,KE)/
     *(.04589D0*ABEAM/ASUB(K)*BESTOP(E1,D,T,T,K))-1.D0
      PRINT 9612,RAB
 500  CONTINUE
      RETURN

 9611 FORMAT(//' Error table for the stopping power of substance No.'
     *,I2/)
 9612 FORMAT(1P9G12.3)
      END

C************************************************************
C  BELAUN initializes common blocks /BECNST/,/BEAM1/,/BEAM2/
C  for BESTOP (BESTO2)
      SUBROUTINE BELAUN(ABEAM,ZBEAM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL RA,ESHL,BEPOT,D300
      COMMON/URSAZ/ATARG(5),ZTARG(5)
      COMMON/BECNST/C11,C12,C13,C14,C15,C16,C17,C1R
      COMMON/BEAM1/A1,Z1,B1,BET1,SGM1,AMU1,PTIF1,DINP1,APIN1(4,100)
      COMMON/BEAM2/POT2(350),NSH2(6),NE2(100),ESH2(100),GB0(5),B2(5),
     *BET2(5),SGM2(5),AMU2(5)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/PINT/PTIF,DINP,APIN(4,100)
      COMMON/UR2OUT/YAP,PEAP,EEAP,YLAP,PELAP,EELAP,DYLR,DPELR,DEELR,
     *DYLT,DPELT,DEELT,PIAP,EIAP,DPILR,DEILR,DPILT,DEILT
      COMMON/EOSZ/NA(30)
      COMMON/EOSPAR/RA(331)
      COMMON/NESHL/NESHL(250)
      COMMON/BESHL/ESHL(250)
      COMMON/BEAM3/BEPOT(30),D300(30)
      COMMON/BESOUT/CLBE,SBE,CLFE,SFE,CLFI,SFI,CLNU,SNU,STPE,STPI,
     *Y2,Y300,Z1EF,NIBES

      C11=2.D0/3.D0
      C12=.5D0*(3.D0*3.1416D0**2)**C11
      C13=3.1416D0
      C14=1.781D0**2
      C15=2.71828D0
      C16=1.D0/1.37D2**2

C---------------------------------------------------- begin PROJECTILE
      A1=ABEAM
      Z1=ZBEAM
      R=Z1+.5D0
      IZ=R
      J=0
      DO 10 K=1,30
      IF(NA(K).EQ.IZ) GOTO 20
 10   J=J+11
      PRINT 511,IZ
      STOP
 511  FORMAT(/' BELAUN: no such beam element in c/blk /EOSZ/',I6)

C: assign EOS parameters of the element ZBEAM,ABEAM required for solving
C  the ionization equation for the projectile ion at rest:
 20   BET1=RA(J+2)
      SGM1=RA(J+5)
      AMU1=RA(J+6)
      B1=RA(J+10)
      A=A1
      Z=Z1
      CALL PINTA
      PTIF1=PTIF
      DINP1=DINP
      DO 30 K=1,IZ
      DO 30 J=1,4
 30   APIN1(J,K)=APIN(J,K)
      PRINT 514,A1,Z1,B1,BET1,SGM1,AMU1
 514  FORMAT(//' BEAM  ELEMENT:  A =',F7.2,' Z =',F5.1,' B =',1PE9.3,
     *' BET =',1PE9.3,' SGM =',1PE9.3,' AMU =',1PE9.3)
C====================================================== end PROJECTILE
C-------------------------------------------------------- begin TARGET
C: set atomic data (ioniz. potentials and subshell binding energies),
C  factor G0, and 4 EOS parameters in c/blk /BEAM2/ for each
C  of the target elements:
      IP=0
      ISH=0
      DO 199 KK=1,5
      Z=ZTARG(KK)
      R=Z+.5D0
      IZ=R
      IF(IZ.LT.1) GOTO 199
      A=ATARG(KK)
      CALL PINTA
      DO 110 L=1,IZ
      IP=IP+1
 110  POT2(IP)=POT(L+1)
      IS=0
      IA=0
      DO 120 K=1,30
      IQ=0
      DO 130 J=1,27
      IS=IS+1
      IQ=IQ+NESHL(IS)
      IF(IQ.EQ.NA(K)) GOTO 136
 130  CONTINUE
      PRINT 512,NA(K)
      STOP
 512  FORMAT(/' BELAUN: electron sum over subshells not equal to Z',I6)
 136  IF(NA(K).EQ.IZ) GOTO 140
 120  IA=IA+11
      PRINT 513,IZ
      STOP
 513  FORMAT(/' BELAUN: subshell data not found for Z2=',I8)
 140  IS=IS-J
      NSH2(KK)=J
      DO 150 L=1,J
      IS=IS+1
      ISH=ISH+1
      NE2(ISH)=NESHL(IS)
 150  ESH2(ISH)=ESHL(IS)
      BET2(KK)=RA(IA+2)
      SGM2(KK)=RA(IA+5)
      AMU2(KK)=RA(IA+6)
      B2(KK)=RA(IA+10)

C: calculate parameter G0:
      IF(IZ.NE.1) GOTO 159
      GB0(KK)=1.105D0
      Y300=0.D0
      GOTO 180
 159  D=D300(K)/(11.206D0*A)
      A1=1.D0
      Z1=1.D0
      R=1.D0
      GB0(KK)=1.D0
      C17=1999.D0
      L=0
      R1=DLOG(1.5D0*27.21D0/(C16*BEPOT(K)))
 160  GB0(KK)=GB0(KK)+R
      L=L+1
      F2=BESTOP(1.D0/C16,D,9.5D-4,9.5D-4,KK)
      F2=R1-(1.D0-Y300/Z)*CLBE-Y300*CLFE/Z
      IF(L.EQ.1) F1=1.1D0*F2
      IF(L.GT.30) GOTO 177
      IF(F2*F1) 165,180,170
 165  IF(DABS(R).LT.1.D-4) GOTO 180
 166  R=-.3D0*R
 167  F1=F2
      GOTO 160
 170  IF(DABS(F2)-DABS(F1)) 167,166,171
 171  GB0(KK)=GB0(KK)-R
      R=-R
      GOTO 160
 177  PRINT 516,Z,GB0(KK),F2
 180  C17=0.D0
      A1=ABEAM
      Z1=ZBEAM
      Y20=Y300
 516  FORMAT(' BELAUN: DIVERGENCE  WHEN  CALCULATING  G0;  Z2,G0,F2= '
     *,1P3E12.5)
      PRINT 515,KK,A,Z,B2(KK),BET2(KK),SGM2(KK),AMU2(KK),GB0(KK),Y20
 515  FORMAT(' TARG.ELEM.NO.',I2,': A =',F7.2,' Z =',F5.1,' B ='
     *,1PG9.3,' BET =',1PG9.3,' SGM =',1PG9.3,' AMU =',1PG9.3,' G0 ='
     *,1PG9.3,' Y300 =',1PG9.3)
 199  CONTINUE
C========================================================== end TARGET
      RETURN
      END

C**********************************************************************
C  BESTOP = stopping power in units  458.9*A1/A2 MeV*cm**2/g;
C   E1 = the projectile energy in units  49.6*A1 keV;
C   D  = the target density in a.u. (11.206*A2 gm/cc);
C   TE, TI are the electron and ion temperatures in a.u. (27.21 eV);
C    NSUB = number of element in array ZTARG used as target material.
C---------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION BESTOP(E1,D,TE,TI,NSUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/URSAZ/ATARG(5),ZTARG(5)
      COMMON/BECNST/C11,C12,C13,C14,C15,C16,C17,C1R
      COMMON/BEAM1/A1,Z1,B1,BET1,SGM1,AMU1,PTIF1,DINP1,APIN1(4,100)
      COMMON/BEAM2/POT2(350),NSH2(6),NE2(100),ESH2(100),GB0(5),B2(5),
     *BET2(5),SGM2(5),AMU2(5)
      COMMON/UR2OUT/YAP,PEAP,EEAP,YLAP,PELAP,EELAP,DYLR,DPELR,DEELR,
     *DYLT,DPELT,DEELT,PIAP,EIAP,DPILR,DEILR,DPILT,DEILT
      COMMON/BESOUT/CLBE,SBE,CLFE,SFE,CLFI,SFI,CLNU,SNU,STPE,STPI,
     *Y2,Y300,Z1EF,NIBES
      H(X)=DLOG(1.D0+X/(1.D0 +3.5D0/DSQRT(X)))
      ETA(X)=.353D0+X**2*((X**3+2.34D0)/(X**3+11.D0))

 2    CLBE=0.D0
      SBE=0.D0
      CLFE=0.D0
      SFE=0.D0
      CLFI=0.D0
      SFI=0.D0
      CLNU=0.D0
      SNU=0.D0
      IF(E1-1.D-30) 999,999,10

 10   A2=ATARG(NSUB)
      Z2=ZTARG(NSUB)
      R=Z2+.5D0
      IZ2=R
C......................................................................
C  Here we have to calculate the ionization degree of target element  .
C     # NSUB for density D and two temperatures: the room temperature .
C     9.5D-4 a.u. and given temperature TE;                           .
C  Option #1: use QION when QION.GE.0 for `hot' ionization and        .
C             assume zero `cold' ionization;                          .
C  Option #2: use tabular EOS by invoking URSAPB when QION.LT.0;      .
C  Option #3: use the full scale EOS model of the author when QION.LT.0
C             this last option is rather time consuming;              .
C......................................................................
C-------------------------------------------------------- begin TABURS
      CALL URSAPB(NSUB,D,9.5D-4,9.5D-4,1)
      Y300=YAP
      CALL URSAPB(NSUB,D,TE,TI,1)
      Y2=YAP
C========================================================== end TABURS

 30   CONTINUE
      IF(Y2.LT.1.D-6) Y2=1.D-6
      RY2=(D*Y2)**C11
      VV1=E1*(2.D0+C16*E1)/(1.D0+C16*E1)**2
      DCL=2.D0*DLOG(1.D0+C16*E1)-C16*VV1
C: Z_eff of hydrogen is set to 1:
      Z1EF=1.D0
      YT=1.D0
      IF(Z1.LT.1.5D0) GOTO 67
C: Z_eff of a non-hydrogen projectile ion:
      Z1EF=Z1/(1.D0+(.62D0*Z1**C11/DSQRT(VV1))**1.7D0)**.58824D0
C---------------------------------------------------------- begin IDAR
C: evaluate the `Ionization Degree At Rest' of the beam ion:
      YTN=C1R
      IF(YTN.LT.1.D-8) YTN=.5D0*Z1
      IF(YTN.GT.Z1-1.D-8) YTN=.5D0*Z1
      HM=.2D0*Z1
      NIBES=0
      R=C12*RY2
      R1=AMU1*TE**2
      RMU=R+TE*(TE/(TE+.4D0*R)-1.5D0*DLOG(1.D0+2.5D0*TE/R))
 40   NIBES=NIBES+1
      YT=YTN
      R=Y2*D/YT
      R4=R**BET1
      IF(R*1.D-8-1.D0) 41,43,43
 41   R2=R**SGM1
      IF(R2-1.D-16) 42,44,44
 42   IF(R1-1.D-16) 43,44,44
 43   R3=1.D0
      GOTO 45
 44   R3=R2/(R2+R1)
 45   BB1=B1*R4*R3
      DB1=BB1*(BET1+SGM1*(1.D0-R3))/YT
      IF(YT-Z1+1.D0) 48,50,50
 48   K=YT
      K=K+1
      PO=APIN1(4,K)+YT*(APIN1(3,K)+YT*(APIN1(2,K)+YT*APIN1(1,K)))
      DPO=APIN1(3,K)+2.D0*YT*APIN1(2,K)+3.D0*YT**2*APIN1(1,K)
      GOTO 51
 50   R3=PTIF1/(Z1-YT)**DINP1
      PO=R3*YT
      DPO=R3*(1.D0+DINP1*YT/(Z1-YT))
 51   F2=RMU+PO-BB1
      DF=DPO+DB1
      IF(DABS(F2)-DABS(.01D0*YT*DF)) 66,66,52
 66   IF(DABS(F2)-DABS(1.D-8+.1D0*PO)) 67,67,52
 52   IF(YT-1.D-6) 53,53,54
 53   IF(F2) 54,67,67
 54   IF(Z1-YT-1.D-6) 55,55,56
 55   IF(F2) 67,67,56
 56   IF(NIBES.EQ.1) F1=F2/DABS(F2)
      IF(F1*F2) 57,67,58
 57   HM=.3D0*HM
      IF(HM-.003D0) 67,67,58
 58   IF(DABS(F2)-DABS(DF*HM)) 59,60,60
 59   YTN=YT-F2/DF
      GOTO 61
 60   YTN=YT-F2*HM/DABS(F2)
 61   IF(YTN-1.D-8) 62,62,63
 62   YTN=.5D0*(YT+1.D-8)
      GOTO 65
 63   IF(Z1-YTN-1.D-8) 64,64,65
 64   YTN=.5D0*(YT+Z1-1.D-8)
 65   F1=F2/DABS(F2)
      IF(NIBES.LT.100) GOTO 40
      PRINT 511,D,TE,Y2,YT,YTN,F2,F1,HM
      YT=0.D0
 511  FORMAT(' DIVERG.IN BESTOP,Z1-ION.:',1P8E9.3)
C: `ionization degree at rest' is calculated.
C============================================================ end IDAR
 67   C1R=YT
      IF(Z1EF.LT.YT) Z1EF=YT
C------------------------------------------------------------ begin BE
C: evaluate the stopping power of bound electrons:
      K=1
      KSH=0
      KPT=0
 110  IF(K.EQ.NSUB) GOTO 120
      KSH=KSH+NSH2(K)
      R=ZTARG(K)+.5D0
      N1=R
      KPT=KPT+N1
      K=K+1
      GOTO 110
 120  Q=Z2-Y2
      IF(Q.LT.1.D-8) Q=1.D-8
      IQ=Q
      IF(IQ.GE.IZ2) IQ=IZ2-1
      NOST=0
      KSS=KSH
 130  KSS=KSS+1
      NOST=NOST+NE2(KSS)
      IF(NOST.LT.IQ+1) GOTO 130
      NOST=NOST-NE2(KSS)
      KS=KSS
      IF(NOST.EQ.IQ) KS=KSS-1
      IY=IZ2-IQ+KPT
      R=D**SGM2(NSUB)
      R1=AMU2(NSUB)*TE*TE
      IF(R-1.D-16) 140,150,150
 140  IF(R1-1.D-16) 144,150,150
 144  R=1.D0
      GOTO 160
 150  R=R/(R+R1)
 160  BB2=B2(NSUB)*R*D**BET2(NSUB)
      R=Q-IQ
      IF(IQ.EQ.0) GOTO 170
C: DESH is the increase in the subshell binding energies
C  due to ionization:
      DESH=R*(POT2(IY)-ESH2(KSS)-BB2)+(1.D0-R)*(POT2(IY+1)-ESH2(KS)-BB2)
      IF(IZ2.EQ.1) GOTO 166
      IF(Y2-Y300) 162,162,164
 162  GB=GB0(NSUB)
      GOTO 180
 164  GB=GB0(NSUB)+(1.105D0-GB0(NSUB))*(Y2-Y300)/(Z2-Y300-1.D0)
      GOTO 180
 166  GB=1.105D0
      GOTO 180
 170  DESH=POT2(IY)-ESH2(KSS)-BB2
      GB=1.105D0
 180  IF(DESH.LT.0.D0) DESH=0.D0
      K=KSS
      R=2.D0*VV1/(GB*(ESH2(K)+DESH))
      R=R/DSQRT(1.D0+C14*Z1EF**2/VV1)
      ZLB=(Q-NOST)*H(R)
 186  K=K-1
      IF(K.EQ.KSH) GOTO 190
      R=2.D0*VV1/(GB*(ESH2(K)+DESH))
      R=R/DSQRT(1.D0+C14*Z1EF**2/VV1)
      ZLB=ZLB+NE2(K)*H(R)
      GOTO 186
 190  CLBE=ZLB/Q
      ZLB=ZLB+Q*DCL
      SBE=4.D0*C13*Z1EF**2*ZLB/(A1*VV1)
C: the stopping due to bound electrons is calculated.
C============================================================== end BE

C------------------------------------------------------------ begin FE
C: stopping due to free electrons:
      R=911.5D0*A2*VV1
      IF(R-1.D8*TI) 202,204,204
 202  TI2=TI
      GOTO 210
 204  TI2=1.D-8*R
 210  XI=DSQRT(R/TI2)
      TEF=DSQRT(TE*TE+(1.26D0*C13*RY2)**2)
      R=.5D0*VV1
      IF(R-1.D8*TEF) 212,214,214
 212  TE2=TEF
      GOTO 220
 214  TE2=1.D-8*R
 220  XE=DSQRT(R/TE2)
      VE2=2.D0*TE2*ETA(XE)
      R=4.D0*C13*D*Y2*(1.D0/VE2+Y2/(2.D0*TI2*ETA(XI)))
      R=1.D0/R
      R1=(.75D0*Z1EF/C13)**C11/RY2
      IF(R.LT.R1) R=R1
      R1=4.D0*R*VE2/(1.D0+C14*Z1EF**2/VE2)
      R1=DSQRT(R1)
      CLFE=DLOG(1.D0+R1/(1.D0+.5D0/DSQRT(R1)))
      SFE=4.D0*C13*Z1EF**2*(XE**3/(XE**3+1.33D0))/(A1*VV1)
      SFE=SFE*Y2*(CLFE+DCL)
C============================================================== end FE

      IF(C17.GT.1823.D0) GOTO 999
C------------------------------------------------------------ begin NU
C: stopping due to the nucleus-nucleus Coulomb collisions:
      R2=1.D0/(1823.D0*A2)
      VI2=2.D0*R2*TI2*ETA(XI)
      R3=(2.D0*1823.D0*A2*A1/(A2+A1))**2*VI2
      R1=Z1-Z1EF
      IF(R1.LT.1.D0) R1=1.D0
      R4=(Z1/R1**C11)**2
      R1=Z2-Y2
      IF(R1.LT.1.D0) R1=1.D0
      R4=.26D0/(R4+(Z2/R1**C11)**2)
      R1=DSQRT(R3*R4/(1.D0+C14*(Z1*Z2)**2/VI2))
      CLNU=DLOG(1.D0+R1/(1.D0+.35D0/DSQRT(R1)))
      R2=4.D0*C13*R2*(XI**3/(XI**3+1.33D0))/(A1*VV1)
      R2=R2*DSQRT(1.D0+A2/A1)
      SNU=R2*(Z1*Z2)**2*CLNU
C============================================================== end NU

C------------------------------------------------------------ begin FI
C: stopping due to the ion-ion Coulomb collisions:
      R1=(1.D0+C14*(Z1EF*Y2)**2/VI2)/R3
      IF(R1.LT.R4) R1=R4
      R1=DSQRT(R/R1)
      CLFI=DLOG(1.D0+R1/(1.D0+.5D0/DSQRT(R1)))
      SFI=R2*(Z1EF*Y2)**2*CLFI
C============================================================== end FI
 999  BESTOP=SBE+SFE+SFI+SNU
      RETURN
      END

C**********************************************************************
C                      BLOCK OF DATA
C**********************************************************************
      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL POTIN,P1,P2,P3,P4,P5,P6,P7
      REAL RA
      REAL ESHL
      REAL BEPOT,D300
      COMMON/ZPOTIN/NPOT(36)
      COMMON/POTIN/POTIN(1200)
      COMMON/EOSZ/NA(30)
      COMMON/EOSPAR/RA(331)
      COMMON/NESHL/NESHL(250)
      COMMON/BESHL/ESHL(250)
      COMMON/BEAM3/BEPOT(30),D300(30)

      DIMENSION P1(157),P2(93),P3(129),P4(109),P5(153),P6(165),P7(92)
      EQUIVALENCE (POTIN(1),P1(1)),(POTIN(158),P2(1))
     *,(POTIN(251),P3(1)),(POTIN(380),P4(1)),(POTIN(489),P5(1))
     *,(POTIN(642),P6(1)),(POTIN(807),P7(1))

C: Z of elements for which ionization potentials are available:
      DATA NPOT/1,2,3,4,5,6,7,8,9,10,11,12,13,18,22,26,
     *28,29,36,  40,42,47,  54,55,  74,79,  82,83,  92,7*0/

C: ionization potentials of elements (I1,I2,...,in eV) ordered after
C  c/blk /ZPOTIN/:
      DATA P1/13.5985,   24.5876, 54.418,   5.3918, 75.641, 122.45,
C.. Z=4,5:
     *9.3227, 18.211, 153.90, 217.72,
     *8.2981, 25.155, 37.931, 259.38, 340.23,
C.. Z=6:
     *11.26, 24.384, 47.89, 64.49, 392.09, 490.,
C.. Z=7:
     *14.534, 29.60, 47.45, 77.47, 97.89, 552.07, 667.05,
C.. Z=8:
     *13.618, 35.118, 54.936, 77.414, 113.90, 138.12, 739.34, 871.42,
C.. Z=9:
     *17.423, 34.971, 62.71, 87.14, 114.24, 157.164, 185.19, 953.91
     *,1103.13,
C.. Z=10:
     *21.565, 40.963, 63.46, 97.12, 126.2, 157.9, 207.28, 239.10
     *,1195.8, 1362.21,
C.. Z=11:
     *5.139, 47.287, 71.62, 98.92, 138.40, 172.2, 208.5, 264.2, 299.87
     *,1465.1, 1648.7,
C.. Z=12:
     *7.646, 15.035, 80.144, 109.27, 141.27, 186.5, 224.9, 266.0, 328.2
     *,367.5, 1761.8, 1962.68,
C.. Z=13:
     *5.986, 18.829, 28.448, 119.99, 153.83, 198.48, 241.44, 284.6
     *,330.1, 399.4, 442.0, 2086.0, 2304.2,
C.. Z=18:
     *15.76, 27.63, 40.91, 59.81, 75.0, 91.0, 124.3, 143.5, 422.4
     *,478.7, 539.0, 618.3, 686.1, 755.7, 854.8, 918., 4121., 4426.,
C.. Z=22:
     *6.82, 13.58, 27.49, 43.27, 99.30, 119.53, 140.8, 170.4, 192.1
     *,215.9, 265., 291.5, 787.8, 863.1, 942., 1044., 1131., 1221.
     *,1346., 1426., 6249., 6626.,
C.. Z=26:
     *7.902, 16.19, 30.65, 54.8, 75.5, 100., 128., 151., 235., 262.
     *,290., 331., 361., 392., 457., 489., 1266., 1358., 1456., 1582.
     *,1689., 1799., 1950., 2045., 8828., 9278./
C.. Z=28:
      DATA P2/7.637, 18.17, 35.3, 54.9, 75.5, 108., 134., 164., 193.
     *,225., 321., 352., 384., 430., 464., 499., 571., 607., 1546.
     *,1648., 1756., 1894., 2010., 2131., 2295., 2398., 10290., 10775.,
C.. Z=29:
     *7.726, 20.29, 36.84, 57.4, 79.9, 103., 139., 167., 199., 232.
     *,266., 369., 401., 435., 484., 520., 557., 633., 671., 1698.
     *,1804., 1919., 2060., 2182., 2310., 2478., 2560., 11050., 11200.,
C.. Z=36:
     *14.00, 24.36, 37., 52.5, 64.7, 78.5, 116., 133., 219., 269., 318.
     *,367., 416., 465., 515., 565., 614., 664., 837., 887., 937., 987.
     *,1046., 1097., 1219., 1271., 2728., 2897., 3065., 3234., 3457.
     *,3630., 3872., 4021., 16750., 17400./
C.. Z=40:
      DATA P3/6.634, 13.1, 23., 34.3, 81.5, 98., 118., 137., 161., 181.
     *,226., 248., 402., 462., 522., 582., 643., 703., 766., 826., 887.
     *,948., 1144., 1203., 1262., 1320., 1406., 1467., 1609., 1669.
     *,3547., 3738., 3928., 4118., 4395., 4591., 4869., 5037., 20880.
     *,21610.,
C.. Z=42:
     *7.092, 16.15, 27.2, 46.4, 61.2, 68., 127., 146., 167., 189., 217.
     *,239., 288., 313., 505., 571., 636., 702., 767., 833., 902., 968.
     *,1034., 1101., 1305., 1368., 1431., 1494., 1591., 1655., 1805.
     *,1869., 3990., 4191., 4392., 4593., 4902., 5110., 5407., 5585.
     *,23120., 23890.,
C.. Z=47:
     *7.576, 21.5, 34.8, 60.5, 80., 99.5, 119., 139., 159., 179., 199.
     *,277., 304., 331., 357., 389., 417., 482., 512., 817., 895., 974.
     *,1053., 1131., 1210., 1295., 1374., 1454., 1533., 1761., 1835.
     *,1908., 1982., 2109., 2185., 2357., 2431., 5213., 5441., 5669.
     *,5897., 6302., 6540., 6889., 7092., 29270., 30140./
C.. Z=54:
      DATA P4/12.13, 21.2, 32.1, 46.7, 59.7, 71.8, 98., 112., 171.
     *,202., 233., 264., 294., 325., 358., 390., 421., 452., 573., 608.
     *,643., 678., 726., 762., 853., 891., 1394., 1491., 1587., 1684.
     *,1781., 1877., 1987., 2085., 2183., 2281., 2548., 2637., 2726.
     *,2814., 3001., 3093., 3296., 3386., 7224., 7491., 7758., 8024.
     *,8617., 8899., 9330., 9569., 39250., 40270.,
C.. Z=55:
     *3.894, 25.1, 35.3, 48.1, 60.9, 75.6, 89., 118., 134., 201., 233.
     *,266., 298., 330., 363., 397., 430., 463., 496., 622., 658., 695.
     *,731., 782., 820., 914., 953., 1489., 1588., 1687., 1787., 1886.
     *,1985., 2099., 2199., 2300., 2401., 2673., 2764., 2855., 2946.
     *,3143., 3237., 3446., 3537., 7540., 7812., 8084., 8356., 8982.
     *,9271., 9714., 9958., 40820., 41860./
C.. Z=74:
      DATA P5/7.98, 17., 25.4, 39.3, 53.2, 67., 120., 141., 162., 183.
     *,204., 241., 263., 295., 340., 370., 395., 436., 481., 526., 571.
     *,617., 665., 710., 756., 801., 847., 893., 1154., 1206., 1259.
     *,1312., 1365., 1417., 1483., 1537., 1591., 1645., 1870., 1926.
     *,1981., 2037., 2163., 2223., 2386., 2447., 3734., 3882., 4029.
     *,4177., 4325., 4472., 4684., 4836., 4987., 5139., 5538., 5671.
     *,5803., 5936., 6468., 6611., 6919., 7055., 14760., 15140., 15520.
     *,15900., 17630., 18060., 18800., 19150., 77510., 78990.,
C.. Z=79:
     *9.226, 20.5, 37.4, 54.2, 71., 87.8, 105., 123., 141., 159., 176.
     *,250., 275., 299., 323.5, 365., 392., 433., 487., 517., 546.
     *,600., 654., 709., 763., 818., 872., 931., 986., 1042., 1097.
     *,1152., 1207., 1516., 1575., 1634., 1692., 1751., 1810., 1888.
     *,1948., 2009., 2069., 2325., 2387., 2448., 2510., 2671., 2738.
     *,2924., 2991., 4516., 4676., 4837., 4997., 5158., 5318., 5566.
     *,5731., 5896., 6061., 6500., 6644., 6787., 6931., 7615., 7772.
     *,8111., 8259., 17090., 17500., 17910., 18320., 20570., 21040.
     *,21870., 22260., 89680., 91290./
C.. Z=82:
      DATA P6/7.417, 15.6, 32.3, 43.6, 67.2, 88., 109., 130., 150., 171.
     *,195., 216., 238., 259., 347., 374., 401., 428., 478., 507., 572.
     *,614., 646., 695., 754., 814., 874., 934., 994., 1053., 1118.
     *,1179., 1240., 1300., 1361., 1421., 1758., 1820., 1883., 1946.
     *,2008., 2071., 2157., 2221., 2285., 2350., 2625., 2690., 2755.
     *,2821., 3008., 3079., 3279., 3350., 5026., 5194., 5362., 5530.
     *,5698., 5866., 6139., 6312., 6485., 6658., 7122., 7273., 7423.
     *,7574., 8368., 8534., 8891., 9048., 18580., 19010., 19440., 19860.
     *,22480., 22980., 23880., 24290., 97550., 99250.,
C.. Z=83:
     *7.285, 16.7, 26., 46., 58., 86., 108., 130., 151., 173., 195.
     *,221., 243., 266., 288., 381., 409., 436., 464., 517., 548., 621.
     *,659., 691., 748., 809., 871., 932., 994., 1055., 1117., 1184.
     *,1246., 1308., 1371., 1433., 1496., 1841., 1905., 1969., 2033.
     *,2097., 2161., 2249., 2315., 2380., 2446., 2729., 2795., 2861.
     *,2927., 3125., 3197., 3402., 3474., 5202., 5372., 5543., 5713.
     *,5884., 6054., 6335., 6511., 6687., 6863., 7336., 7489., 7642.
     *,7795., 8629., 8797., 9161., 9321., 19100., 19530., 19960.
     *,20390., 23150., 23650., 24580., 24990., 100300., 102000./
C.. Z=92:
      DATA P7/6.194, 11.6, 18.1, 30.9, 49.9, 68.9, 90.4, 105., 119.
     *,133., 165., 181., 223., 242., 313., 343., 374., 404., 434., 465.
     *,506., 537., 569., 600., 731., 766., 801., 836., 928., 966.
     *,1072., 1112., 1212., 1288., 1365., 1441., 1517., 1593., 1670.
     *,1746., 1834., 1912., 1989., 2066., 2144., 2221., 2630., 2705.
     *,2780., 2856., 2931., 3006., 3146., 3223., 3301., 3379., 3718.
     *,3795., 3872., 3949., 4267., 4352., 4604., 4689., 6899., 7093.
     *,7286., 7480., 7673., 7867., 8240., 8442., 8643., 8844., 9397.
     *,9571., 9744., 9918., 11190., 11390., 11830., 12020., 24020.
     *,24500., 24990., 25470., 29790., 30380., 31550., 32030., 127300.
     *,129300./

C:  Z of elements for which  EOS  is available:
      DATA NA/1,2,3,4,5,6, 7,8, 11,13,18, 22,26,29, 42,47,54,55,
     *74,79,82, 83,92,7*0/

C:  EOS  parameters: [regular: B=5.7*Z**(2/3)]
C:   AL,  BET, GAM, DEL,    SGM,    AMU,   ANU, ALAM,Y1, B,   1/DZ
      DATA RA/
C.. Z=1,2,3,4,6
     *0.,.6667,1.387,2.813, 1.75,  .1554,   0.,  .1,0.,  5.7,  124.5,
     *0.,.65,  1.3, 10.27,  1.8,   .0517,   0.,  .1,0., 10.7,  300.,
     *0.,.64,  1.1754,.061064,1.75,.0265,  .12,  .1,0., 15.306, 97.21,
     *0.,.6,   1.071, .269, 1.8,  5.705E-3,.18,  .1,0., 26.73,  52.8,
     *0.,.6,   1.74 , .07,  1.8,  5.70E-3, .008, .1,0., 30.  ,  51.82,
     *0.,.6,    .65, 5.,    1.6,  2.04E-2,  0.,  .1,0., 34.36,  59.76,
C.. Z=7,8:
     *0.,.67,  1.3,  9.,    1.8,  1.00E-2,  0.,  .1,0., 21.,   152.8,
     *0.,.67,  1.3,  9.,    1.8,  1.00E-2,  0.,  .1,0., 22.8,  122.8,
C.. Z=11,13,18:
     *0.,.65,  1.193, .175, 1.8,  8.754E-3,.006, .1,0., 32.22, 252.7,
     *0.,.59,  1.205, .269, 1.72, 4.97E-3, .035, .1,0., 47.87, 108.8,
     *0.,.68,  1.35, 8.858, 1.8,  7.891E-3,.5,   .1,0., 35.27, 249.,
C.. Z=22,26,29:
     *0.,.67,  1.4,   .2,   1.8,  3.73E-3, .03,  .1,0., 43.63, 114.2,
     *0.,.6667,1.185, .856, 1.82, 3.2E-3,  .025, .1,0., 50.,    78.61,
     *0.,.6667,1.158,1.313, 1.82, 3.32E-3, .04,  .1,0., 53.8,   78.52,
C.. Z=42,47,54,55:
     *0.,.6667,2.,    .138, 1.8,  2.48E-3, 2.E-4,.1,0., 68.87, 104.44,
     *0.,.6667,1.193,1.496, 1.8,  2.787E-3,.05,  .1,0., 74.22, 112.7,
     *0.,.67,  1.,   1.,    1.8,  3.E-3,   .1,   .1,0., 81.4,  544.,
     *0.,.705, 1.,    .17,  1.85, 2.615E-3,.015, .1,0., 63.95, 727.02,
C.. Z=74,79,82,83,92:
     *0.,.6667,1.526, .297, 1.85, 1.373E-3,.004, .1,0.,100.5,  106.2,
     *0.,.6667,1.33, 1.78,  1.85, 1.373E-3,.0125,.1,0.,105.,   113.2,
     *0.,.6667,1.25, 1.24,  1.85, 1.192E-3,.008, .1,0.,107.6,  201.7,
     *0.,.6667,1.24, 1.06,  1.85, 1.136E-3,.005, .1,0.,108.5,  236.,
     *0.,.68,   .93, 4.19,  1.8,  1.39E-3, .15,  .1,0.,109.4,  141.1,
     *78*0./

C:  numbers of electrons in subshells ordered after /EOSZ/:
      DATA NESHL/
C.. Z=1,2,3,4,5,6,7,8:
     *1,  2,  2,1,  2,2, 2,2,1,  2,2,2,  2,2,2,1, 2,2,2,2,
C.. Z=11:
     *2,2,2,4,1,
C.. Z=13:
     *2,2,2,4,2,1,
C.. Z=18:
     *2,2,2,4,2,2,4,
C.. Z=22:
     *2,2,2,4,2,2,4,2,2,
C.. Z=26:
     *2,2,2,4,2,2,4,4,2,2,
C.. Z=29:
     *2,2,2,4,2,2,4,4,6,1,
C.. Z=42:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,1,1,
C.. Z=47:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,1,
C.. Z=54:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,2,4,
C.. Z=55:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,2,4,1,
C.. Z=74 (tentative):
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,6,8,2,4,4,2,
C.. Z=79:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,6,8,2,4,4,6,1,
C.. Z=82:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,6,8,2,4,4,6,2,2,
C.. Z=83:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,2,6,8,2,4,4,6,2,2,1,
C.. Z=92:
     *2,2,2,4,2,2,4,4,6,2,2,4,4,6,6,8,2,2,4,4,6,2,2,4,3,2,1,
     *0/

C:  binding energies of subshells (in a.u.) ordered after /EOSZ/:
      DATA ESHL/
C.. Z=1,2,3,4,5,6:
     *.5, .904, 2.2,.1984, 4.53,.293, 7.3,.44,.3,  10.9,.665,.348,
C.. Z=7,8:
     *14.81, .747, .534, .534,   19.77, 1.047, .50, .50,
C.. Z=11:
     *39.5, 2.64, 1.34, 1.34, .178,
C.. Z=13:
     *57.6, 4.7, 2.98, 2.95, .374, .203,
C.. Z=18:
     *118.,  12., 9.2,  9.1, 1.22, .543, .536,
C.. Z=22:
     *183., 21.1, 17.4, 17.1, 2.75, 1.73, 1.61, .304, .207,
C.. Z=26:
     *262., 31.7, 27., 26.5, 4.03, 2.7, 2.47, .439, .411, .248,
C.. Z=29:
     *330., 40.6, 35.2, 34.5, 4.83, 3.12, 3.03, .264, .255, .235,
C.. Z=42:
     *735.,  106., 96.8, 92.9, 18.9, 15.5, 14.8, 8.75, 8.6, 2.76
     *,1.79, 1.57, .238, .212, .205,
C.. Z=47:
     *938., 140., 130., 123., 26.8,  22.5, 21.4,  14., 13.8, 4.07
     *,2.69, 2.5, .41,  .39,  .234,
C.. Z=54:
     *1270., 199.1, 187.6, 175.5, 41.2, 36.4, 34., 25.4, 24.9, 7.66
     *,5.91, 5.44,  2.57,  2.49,  .867, .456, .404,
C.. Z=55:
     *1323., 210., 197., 184., 45.1, 39.5,  37.,  27.4, 27., 9.08
     *,7.,   6.5,  3.09,  3.,  1.25, .673, .609,  .125,
C.. Z=74 (tentative):
     *2600., 500., 480., 420., 120., 110., 100., 80., 75., 25.
     *,24., 20., 13., 12., 4.4, 3.3, 3.2, 3., 2.2, .4, .3,
C.. Z=79:
     *2966., 528., 505., 438., 127., 116., 101., 84.5, 81.3, 28.5
     *,24.1, 20.5, 13.3, 12.7, 4.49, 3.39, 3.25, 3., 2.4, .395
     *,.342, .276,
C.. Z=82:
     *3234., 583., 559., 479., 142., 131., 113., 95.3, 91.5, 33.3
     *,28.5, 24.1, 16.3, 15.5, 5.91, 5.41, 5.23, 4.21, 3.38, .988
     *,.883, .523, .234,
C.. Z=83:
     *3327., 603., 578., 493., 148., 136., 117.,  99.,  95., 35.1
     *,30.1, 25.4, 17.5, 16.6, 6.45, 6.18, 5.98, 4.68, 3.77, 1.23
     *,1.12, .644,  .3,  .2405,
C.. Z=92:
     *4248., 800., 770., 631., 205., 191., 159., 137., 131., 53.6
     *,47.3, 38.8, 29.,  27.5, 14.6, 14.2, 12.3, 9.8,  7.78, 4.17
     *,3.75, 2.03, 1.27, .886, .233, .19,  .182,
     *0./

C: recommended mean excitation energies (in eV) ordered after /EOSZ/:
      DATA BEPOT/
C.. Z=1,2,3,4,5,6,7,8:
     *15., 42.3, 40., 64., 76., 79., 82., 98.5,
C.. Z=11,13,18:
     *148., 164., 188.,
C.. Z=22,26,29:
     *228., 275., 317.,
C.. Z=42,47,54,55:
     *414., 469.,  489., 490.,
C.. Z=74,79,82,83,92:
     *722., 770., 793., 765., 884.,
     *7*0./

C:  densities (in g/cc) at room temperature and 1 atm pressure
C   ordered after /EOSZ/
      DATA D300/
C.. Z=1,2,3,4,5,6,7,8:
     *9.E-5, 1.8E-4, .534, 1.85, 2.34, 2.25, 1.251E-3, 1.429E-3,
C.. Z=11,13,18:
     *.97, 2.7, 1.66E-3,
C.. Z=22,26,29:
     *4.54, 7.87, 8.93,
C.. Z=42,47,54,55:
     *10.23, 10.5, 5.85E-3, 1.9,
C.. Z=74,79,82,83,92:
     *19.3, 19.3, 11.34, 9.79, 18.6,
     *7*0./
      END

C*********************************************************************
C  PINTA initializes  c/blks /CONST/ and /PINT/, and fills array POT
C  in c/blk /ELEMNT/ for a given element (A,Z)
      SUBROUTINE PINTA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL POTIN
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/CONST/C1,C2,C3,C4,C5
      COMMON/PINT/PTIF,DINP,APIN(4,100)
      COMMON/SGLASP/XI(250),WI(250),YT(250),DV1(250),DV2(250),BD(2)
      COMMON/ZPOTIN/NPOT(36)
      COMMON/POTIN/POTIN(1200)
      EXTERNAL FIZA

C: set a signal that the EOS c/blks are violated and LAUNCH may need
C  to be re-called (when PINTA is called not from LAUNCH):
      LZEOS=0

      C1=2.D0/3.D0
      C2=.5D0*(3.D0*3.1416D0**2)**C1
      C3=(4.D0*3.1416D0/3.D0)**(1.D0/3.D0)
      C5=2.5D0

      R=Z+.5D0
      NZ=R
      NN=0
      POT(1)=0.D0
      XI(1)=0.D0
      IF(NZ.EQ.3) GOTO 121
      IF(NZ.EQ.1) GOTO 111
      IF(NZ.EQ.2) GOTO 131

      NIT=3
      DO 1 K=1,30
      IF(NPOT(K).EQ.NZ) GOTO 11
 1    NN=NN+NPOT(K)
      PRINT 555,Z
      RETURN
 555  FORMAT(/'NO SUCH ELEMENT IN THE TABLE OF IONIZATION POTENTIALS',
     *1PE12.4)
 11   WI(1)=1.D30
      NF=NZ+1
      BD(1)=0.D0
      BD(2)=0.D0
      IT=0
      RI=0.D0
 31   IT=IT+1
      DO 2 I=2,NF
      K=NN+I-1
      POT(I)=(1.D0+RI*((I-1)/Z)**2)*POTIN(K)/27.2D0
      WI(I)=(1.D0+30.D0*((I-1)/Z)**2)/(.2D0+POT(I))
      IF(NZ.LE.10) WI(I)=5.D0/POT(I)
 2    XI(I)=I-1
      CALL SMOOSP(NF,XI,POT,WI,1.D0,BD,APIN,YT,DV1,DV2,3)
      PTIF=YT(NZ)/(Z-1.D0)
      DINP=(DV1(NZ)/PTIF-1.D0)/(Z-1.D0)
      RE=0.D0
      RAA=PTIF*(Z/(1.D0-DINP)-1.D0/(2.D0-DINP))
      DO 17 I=1,NZ
      IF(I.EQ.NZ) GOTO 16
      R=I
      RAA=RAA+GQUAD(FIZA,R-1.D0,R,8)
 16   K=NN+I
      RE=RE+POTIN(K)/27.2D0
 17   CONTINUE
      IF(NZ.GT.10) GOTO 134
      RI=RE/RAA
      APIN(4,1)=0.D0
      PTIF=RI*PTIF
      DO 19 I=1,4
      DO 19 K=1,NF
 19   APIN(I,K)=RI*APIN(I,K)
      RETURN
 134  RI=RI+1.2D0*(RE/RAA-1.D0)
      IF(IT.LE.NIT) GOTO 31
      DO 135 I=2,NF
      K=NN+I-1
 135  POT(I)=POTIN(K)/27.2D0
      APIN(4,1)=0.D0
      RETURN
C: hydrogen:
 111  DINP=1.D-1
      POT(1)=0.D0
      POT(2)=.5D0
      PTIF=(1.D0-DINP)*(1.D0-.5D0*DINP)
      RETURN
C: lithium:
 121  PTIF=1.292857D0
      DINP=.4D0
      APIN(4,1)=0.D0
      APIN(1,1)=0.D0
      APIN(2,1)=.6137271124D0
      APIN(3,1)=.1040373904D0
      APIN(4,2)=.0772671108D0
      APIN(3,2)=-.1277639418D0
      APIN(2,2)=.8455284446D0
      APIN(1,2)=-APIN(4,2)
      POT(1)=0.D0
      POT(2)=5.4D0/27.2D0
      POT(3)=75.6D0/27.2D0
      POT(4)=122.5D0/27.2D0
      RETURN
C: helium:
 131  DINP=.3D0
      PTIF=1.0677D0
      APIN(4,1)=0.D0
      APIN(1,1)=0.D0
      APIN(3,1)=.74739D0
      APIN(2,1)=.32031D0
      POT(1)=0.D0
      POT(2)=24.6D0/27.2D0
      POT(3)=54.4D0/27.2D0
      RETURN
      END

C*********************************************************************
C*********************************************************************
C   The subroutines below belong to the EOS complex; they need
C   the above BLOCK DATA as the input information, plus PINTA.
C*********************************************************************
C  LAUNCH initializes the EOS complex for a given element (A,Z)
      SUBROUTINE LAUNCH(AA,ZZ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL RA
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/CONST/C1,C2,C3,C4,C5
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/PRS/GYC,DLYDLT,G,GP,YC,YW,PE,PN,PS
      COMMON/EOSZ/NA(30)
      COMMON/EOSPAR/RA(331)
      A=AA
      Z=ZZ
      CALL PINTA
      R=Z+.5D0
      NZ=R

      J=0
      DO 1 K=1,30
      IF(NA(K).EQ.NZ) GOTO 7
 1    J=J+11
      PRINT 501,ZZ
      STOP
 501  FORMAT(/10X,'NO  SUCH  ELEMENT  IN  LAUNCH-TABLES, Z=',1PE12.5)
 7    AL=RA(J+1)
      BET=RA(J+2)
      GAM=RA(J+3)
      DEL=RA(J+4)
      SGM=RA(J+5)
      AMU=RA(J+6)
      ANU=RA(J+7)
      ALAM=RA(J+8)
      Y1=RA(J+9)
      B=RA(J+10)
      DZ=1.D0/RA(J+11)
C      DENZ=DZ*A*11.206D0
      EPS=1.D-6
      NS=32
      Y0=DINA(DZ,0.D0)

C: set Z and A indicators of initialized EOS
C  (should be done only after calling PINTA):
      LZEOS=NZ
      LAEOS=A

C      GYC0=GYC
C      G0=.5D0+GYC0/(1.D0-Y1/Y0)
C      V2=1.D0
C      P2=294.2D0*BARA(2.D0*DZ,0.D0,0.D0)
C      DPDN=C1*C2*Y0*(DZ*Y0)**C1*(.2D0+GYC0-.6D0*GAM*DEL/(1.D0+DEL))/V2
C      ESUB=-2626.D0*ENA(DZ,0.D0,0.D0)/A
C      VSND=51.24D0*DSQRT(DPDN/A)
C      ALVT=(3.D0*G0+.4D0*C2*(DZ*Y0)**C1*
C     *Y0*ANU/(1.D0+DEL)/DZ**GAM)/(DPDN*V2)
C      PRINT 531,A,Z
C      PRINT 532,DENZ,VSND
C      PRINT 533,AL,BET,GAM,DEL,SGM,AMU
C      PRINT 534,B,ANU,ALAM,G0,EPS,NS
C      PRINT 535,Y0,Y1,GYC0,DPDN,ESUB,P2
C      PRINT 5351,ALVT
C 531  FORMAT(//40X,'A = ',F8.2,10X,'Z = ',F8.2/)
C 532  FORMAT(40X,'INITIAL  PARAMETERS :'/' NORMAL  DENSITY (G/CM**3,T=0,
C     *P=0) = ',1PG10.3,5X,'SOUND VELOCITY (KM/S,T=0,P=0) =',1PG10.3)
C 533  FORMAT(' ALPHA =',1PG10.4,'  BETA =',1PG10.4,' GAMMA =',1PG10.4
C     *,' DELTA =',1PG10.4,' SIGMA =',1PG10.4,'  MU =',1PG10.4)
C 534  FORMAT('  B =',1PG10.4,'  NU =',1PG10.4,' LAMBDA =',1PG10.4
C     *,' GRUNAIZ =',1PG10.4,'  EPS =',1PG10.4,' N(GAUSS) =',I3/)
C 535  FORMAT(30X,'ENSUING  PARAMETERS:'/' Y0 =',1PG10.4,'  Y1='
C     *,1PG11.4,' DLYCDLN=',1PG11.4,' DPDN=',1PG11.4,' ESUB=',1PG11.4
C     *,'  P(2DZ) =',1PG11.4)
C 5351 FORMAT(/' THERM.EXP.COEFF. (VOL.,IN  A.U.)::',1PG11.4)
      RETURN
      END

C*******************************************************************
C  DINA = degree of ionization as a function of density D1 (in a.u.)
C  and  electron  temperature  T1 (in a.u.)
      DOUBLE PRECISION FUNCTION DINA(D1,T1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/CONST/C1,C2,C3,C4,C5
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/PRS/GYC,DLYDLT,G,GP,YC,YW,PE,PN,PS
      COMMON/PINT/PTIF,DINP,AA(4,100)
      COMMON/FIX/R,M,NP(80)
      D=D1
      T=T1
      IF(D.GT.1.D-30) GOTO 1
      PRINT 547,D,T
      Y=0.D0
      GOTO 133
 1    N=0
      YN=R
      HM=.2D0*Z
      IF(YN-1.D-26) 601,601,602
 601  YN=.5D0*Z
      GOTO 3
 602  IF(Z-YN-1.D-12) 603,603,3
 603  YN=.5D0*Z
 3    M=M+1
      R1=C2*D**C1
      R9=D**SGM
      IF(T) 6,6,7
 6    T=0.D0
      R8=1.D0
      GOTO 8
 7    IF(R9) 7001,7001,7002
 7001 R8=0.D0
      GOTO 8
 7002 R8=R9/(R9+AMU*T*T)
 8    R2=B*R8*D**BET
 11   N=N+1
      Y=YN
      R3=R1*Y**C1
      R5=T/(T+.4D0*R3)
      R6=1.5D0*DLOG(1.D0+2.5D0*T/R3)
      F=R3+T*(R5-R6)-R2
      DF=C1*(R3+T*R5*(R5+.5D0))/Y
      IF(Y-Z+1.D0)22,22,24
 22   K=Y
      K=K+1
      DI=AA(3,K)+Y*(2.D0*AA(2,K)+3.D0*Y*AA(1,K))
      PT=AA(4,K)+Y*(AA(3,K)+Y*(AA(2,K)+Y*AA(1,K)))
      GOTO 26
 24   R4=PTIF/(Z-Y)**DINP
      PT=R4*Y
      DI=R4*(1.D0+DINP*Y/(Z-Y))
 26   F=F+PT
      DF=DF+DI
      R7=DABS(F)
      IF(R7-DABS(EPS*Y*DF)) 610,610,612
 610  IF(R7-EPS*(1.D-26+1.D1*PT)) 100,100,612
 612  IF(Y-2.D-26) 613,613,614
 613  IF(F) 616,100,100
 614  IF(Z-Y-EPS) 615,615,616
 615  IF(F) 100,100,616
 616  IF(N.EQ.1) FB=F/R7
      IF(F*FB) 617,100,619
 617  IF(HM-EPS) 100,100,618
 618  HM=.3D0*HM
 619  IF(R7-DABS(DF*HM)) 620,621,621
 620  YN=Y-F/DF
      GOTO 622
 621  YN=Y-HM*(F/R7)
 622  IF(YN-1.D-26) 623,623,624
 623  YN=.5D0*(Y+1.D-26)
      GOTO 626
 624  IF(Z-YN-1.D-12) 625,625,626
 625  YN=.5D0*(Y+Z-1.D-12)
 626  FB=F/R7
      IF(N-100) 11,200,200
 200  PRINT 511,D,T,HM,Y,F,DF,PT,DI,R
 100  R=Y
 103  IF(T) 105,105,107
 105  DLYDLT=0.D0
      GYC=(BET*R2-C1*R3)/(Y*DI+C1*R3)
      GOTO 133
 107  R11=R5*(R5-.5D0)+R6
      R12=-2.D0*R2*(1.D0-R8)
      R13=R3+T*R5*(R5+.5D0)
      DLYDLT=(T*R11+R12)/(Y*DI+C1*R13)
 133  DINA=Y
      RETURN
 511  FORMAT(' DINA:',1P10E11.4)
 547  FORMAT(/10X,'DENSITY L.E. ZERO IN DINA: D,T=',1P2E15.8)
      END

C*********************************************************************
C  BARA = pressure (in a.u.) as a function of density D (in a.u.),
C  electron temperature TE (in a.u.), and ion temperature TI (in a.u.)
      DOUBLE PRECISION FUNCTION BARA(D,TE,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/CONST/C1,C2,C3,C4,C5
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/PRS/GYC,DLYDLT,G,GP,YC,YW,PE,PN,PS
      COMMON/PNZ/Y,PTI,DPTLY,DDPTLY

      T=TE
      IF(T.LT.0.D0) T=0.D0
      R4=DEL*DZ**GAM+ANU*T
      R6=D**(1.D0/3.D0)
      YC=DINA(D,0.D0)
      IF(T) 3,3,7
 3    YW=YC
      GOTO 8
 7    YW=DINA(D,T)
 8    R1=C2*R6**2*YW**C1
      PE=D*YW*(.4D0*R1+T*(T/(T+.4D0*R1)))
      T=TI
      IF(T.LT.0.D0) T=0.D0
      R2=C3*YC**2*R6
      IF(1.D-30*R2-T) 22,24,24
 22   GP=R2/T
      GOTO 25
 24   GP=1.D30
 25   G=.5D0+GYC
 28   PN=D*T*(1.D0+3.D0*ALAM*G*GP)/(1.D0+ALAM*GP)
      IF(R4) 31,31,32
 31   R3=1.D0
      GOTO 43
 32   IF(D-1.D0) 33,33,34
 33   R3=D**GAM
      R3=R3/(R3+R4)
      GOTO 43
 34   R3=1.D0/(1.D0+R4*D**(-GAM))
 43   PS=.4D0*C2*(DZ/Y0)**(1.D0/3.D0)*Y0**2*D*R6*(1.D0+DEL)*R3
      BARA=PE+PN-PS
      RETURN
      END

C*********************************************************************
C  ENA = internal energy (in a.u. per atom) as a function of
C  density D (in a.u.), electron temperature TE (in a.u.),
C  and ion temperature TI (in a.u.)
      DOUBLE PRECISION FUNCTION ENA(D,TE,TI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/ENRG/EI,EE,EIZ,EPI,ES,EET,PPP
      COMMON/CONST/C1,C2,C3,C4,C5
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/PRS/GYC,DLYDLT,G,GP,YC,YW,PE,PN,PS
      COMMON/FEN/TTT,DDD,R00
      EXTERNAL FIZA,FENA,FURA,FENZA

      TTT=TE
      IF(TTT.LT.0.D0) TTT=0.D0
      DDD=D
      PPP=BARA(DDD,TTT,TI)
      R=DDD**BET
      R00=ANU*TTT+DEL*DZ**GAM
      R7=AMU*TTT**2
      EE=1.5D0*PE/D
      EI=1.5D0*TI*(1.D0+2.D0*ALAM*GP)/(1.D0+ALAM*GP)
      X1=0.D0
      X2=YW
      NGQ=2*NS
      NG=NS
      RZZ=Z+.2D0
      NZZ=RZZ
      IF(NZZ.EQ.1) GOTO 41
      IF(X2-Z+1.D0) 41,41,42
 41   EIZ=GQUAD(FIZA,X1,X2,NGQ)
      GOTO 43
 42   EIZ=GQUAD(FIZA,X1,Z-1.D0,NGQ)+GQUAD(FIZA,Z-1.D0,X2,NG)
 43   CONTINUE
      IF(R7-1.D-18) 2,3,3
 2    X1=0.D0
      X2=R
      EPI=B*(GQUAD(FENZA,X1,X2,NG)-R*YW)
      GOTO 22
 3    IF(DDD-1.D0) 4,4,7
 4    R9=DDD**SGM
      R8=R9/(R9+R7)
      GOTO 5
 7    R8=1.D0/(1.D0+R7*DDD**(-SGM))
 5    X1=-8.D0
      X2=0.D0
      DDD=7.D0*R7**(1.D0/SGM)
      IF(DDD-D) 9,13,13
 9    X3=DLOG(D)-DLOG(DDD)
      S=GQUAD(FENA,X1,X2,NG)+GQUAD(FENA,X2,X3,NG)
      GOTO 15
 13   DDD=D
      S=GQUAD(FENA,X1,X2,NG)
 15   EPI=B*(S-R*YW*R8*(3.D0-2.D0*R8))
 22   X1=0.D0
      X2=D**(1.D0/3.D0)
      ES=1.2D0*C2*DZ**(1.D0/3.D0)*Y0**(5.D0/3.D0)*
     *(1.D0+DEL)*GQUAD(FURA,X1,X2,NG)
      EET=EE+EIZ+EPI-ES
      ENA=EI+EET
      RETURN
      END

C*********************************************************************
C  FIZA = ionization  potential I(Y) (in a.u.)
      DOUBLE PRECISION FUNCTION FIZA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/PINT/PTIF,DINP,APIN(4,100)

      IF(X-Z+1.D0) 1,1,2
 1    K=X
      K=K+1
      FIZA=APIN(4,K)+X*(APIN(3,K)+X*(APIN(2,K)+X*APIN(1,K)))
      RETURN
 2    IF(X-Z) 3,4,4
 3    FIZA=PTIF*X/(Z-X)**DINP
      RETURN
 4    FIZA=1.D30
      RETURN
      END

C*********************************************************************
C  FENZA is an integrand for the cold pressure-ionization energy
      DOUBLE PRECISION FUNCTION FENZA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS

      D=X**(1.D0/BET)
      R=DINA(D,0.D0)
      FENZA=R
      RETURN
      END

C*********************************************************************
C  FENA is an integrand for the pressure-ionization energy
      DOUBLE PRECISION FUNCTION FENA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/PRS/GYC,DLYDLT,G,GP,YC,YW,PE,PN,PS
      COMMON/FEN/TTT,DDD,R00

      D=DDD*DEXP(X)
      R=DINA(D,TTT)
      R9=D/DZ
      R7=AMU*TTT**2
      IF(R7) 2,2,3
 2    R8=1.D0
      GOTO 5
 3    IF(D-1.D0) 4,4,7
 4    R9=D**SGM
      R8=R9/(R9+R7)
      GOTO 5
 7    R8=1.D0/(1.D0+R7*D**(-SGM))
 5    R1=(1.D0-DLYDLT)*(BET+SGM*(1.D0-R8))
      R1=R1+2.D0*(1.D0-R8)*(BET+SGM*(1.D0-2.D0*R8))
      FENA=R*R1*R8*D**BET
      RETURN
      END

C*********************************************************************
C  FURA is an integrand for electron-ion interaction energy
      DOUBLE PRECISION FUNCTION FURA(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/FEN/TTT,DDD,R00

      IF(X-1.D0) 1,1,11
 1    R=X**(3.D0*GAM)
      P=R+R00
      IF(P) 2,2,3
 2    FURA=1.D0
      GOTO 100
 3    FURA=(R/P)*(1.D0+ANU*TTT/P)
      GOTO 100
 11   R=X**(-3.D0*GAM)
      P=1.D0+R*R00
      FURA=(1.D0+ANU*TTT*R/P)/P
 100  RETURN
      END

C*********************************************************************
C  SMOOSP generates a smooth interpolation to a discrete function Y(X);
C  PI are the weights;
      SUBROUTINE SMOOSP(N,X,Y,PI,P,G,A,YT,D1,D2,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),Y(N),PI(N),G(2),A(4,N),YT(N),D1(N),D2(N)

      P6=6.D0*P
      XIP1=X(3)-X(2)
      XI=X(2)-X(1)
      XIM1=1.D0
      DXI=XI+XIP1
      DXIM1=XIM1+XI
      DY=6.D0*(Y(2)-Y(1))/XI
      OT=DXI/XIP1/PI(2)+DXIM1/XIM1/PI(1)
      NM1=N-1
      DO 2 I=2,NM1
      XIM2=XIM1
      XIM1=XI
      XI=XIP1
      IF(I.NE.NM1) XIP1=X(I+2)-X(I+1)
      DYM1=DY
      DY=6.D0*(Y(I+1)-Y(I))/XI
      DXIM1=DXI
      DXI=XI+XIP1
      OTM1=OT
      OT=DXI/XIP1/PI(I+1)+DXIM1/XIM1/PI(I)
      A1=P6/(PI(I-1)*XIM1*XIM2)
      A2=XIM1-P6/XIM1**2*OTM1
      A(1,I)=2.D0*DXIM1+P6*((DXIM1/XIM1/XI)**2/PI(I)+
     *1.D0/PI(I-1)/XIM1**2+1.D0/PI(I+1)/XI**2)
      A(2,I)=XI-P6/XI**2*OT
      A(3,I)=P6/(PI(I+1)*XI*XIP1)
      A(4,I)=DY-DYM1
      IF(I.LE.3) A(4,I)=A(4,I)-A1*G(1)
      IF(I.EQ.2) A(4,I)=A(4,I)-A2*G(1)
      IF(I.LT.4) GOTO 1
      AK=A1/A(1,I-2)
      A2=A2-A(2,I-2)*AK
      A(1,I)=A(1,I)-A(3,I-2)*AK
      A(4,I)=A(4,I)-A(4,I-2)*AK
 1    IF(I.LT.3) GOTO 2
      AK=A2/A(1,I-1)
      A(1,I)=A(1,I)-A(2,I-1)*AK
      A(2,I)=A(2,I)-A(3,I-1)*AK
      A(4,I)=A(4,I)-A(4,I-1)*AK
 2    CONTINUE
      D2(1)=G(1)
      D2(N)=G(2)
      D2(N-1)=(A(4,N-1)-D2(N)*(A(2,N-1)+A(3,N-1)))/A(1,N-1)
      DO 3 I=3,NM1
      J=N-I+1
 3    D2(J)=(A(4,J)-A(2,J)*D2(J+1)-A(3,J)*D2(J+2))/A(1,J)
      IF(K.EQ.0) RETURN
      DD2X=(D2(2)-D2(1))/(X(2)-X(1))
      YT(1)=Y(1)-P/PI(1)*DD2X
      DO 4 I=2,NM1
      DD2XM1=DD2X
      DD2X=(D2(I+1)-D2(I))/(X(I+1)-X(I))
 4    YT(I)=Y(I)-P/PI(I)*(DD2X-DD2XM1)
      YT(N)=Y(N)+P/PI(N)*DD2X
      IF(K.EQ.1) RETURN
      DO 5 I=1,NM1
      A(1,I)=(D2(I+1)-D2(I))/(X(I+1)-X(I))/6.D0
      A(2,I)=D2(I)/2.D0-3.D0*A(1,I)*X(I)
      A(3,I)=-A(1,I)*X(I)**2-(X(I)+X(I+1))*(A(2,I)+A(1,I)*X(I+1))+
     *(YT(I)-YT(I+1))/(X(I)-X(I+1))
      A(4,I)=YT(I)-((A(1,I)*X(I)+A(2,I))*X(I)+A(3,I))*X(I)
 5    CONTINUE
      IF(K.EQ.2) RETURN
      DO 6 I=1,NM1
 6    D1(I)=(A(1,I)*3.D0*X(I)+2.D0*A(2,I))*X(I)+A(3,I)
      D1(N)=(A(1,NM1)*3.D0*X(N)+2.D0*A(2,NM1))*X(N)+A(3,NM1)
      RETURN
      END

C*********************************************************************
C  URSAPB calculates two-temperature EOS from tables (in a.u.);
C  NSUBST is substance number (in array ZSBS(I));
C  output data are located in COMMON/UR2OUT/; depending on JOB value,
C  they are:
C  JOB=1: YI=FAPP(1) (ionization degree) and its derivatives;
C      2: YI,PE=FAPP(2),PI,EI and their derivatives;
C      3: YI,PE,EE=FAPP(3),PI,EI and their derivatives (complete set);
C      4: PE,PI,EI      and their derivatives;
C      5: EE,PI,EI      and their derivatives;
C      6: PE,EE,PI,EI   and their derivatives;
C    > 6: PI,EI         and their derivatives;
C
C  FLAPP(1)=ln(YI); FLAPP(2)=ln(PE+P00); FLAPP(3)=ln(EE+E00);
C  DLFAR(I)=d[FAPP(I)]/d[ln(rho)]; DLFAT(I)=d[FAPP(I)]/d[ln(TE)];
C---------------------------------------------------------------------
      SUBROUTINE URSAPB(NSUBST,D,TE,TI,JOB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL FARR
      COMMON/URTBLI/QURS(1000)
      COMMON/URSTBL/FARR(30000)
      COMMON/URSEDG/NROO(6),ROILL(5),HLROO(5),NTEMM(3,6),
     *TEMILL(3,5),HLTEMM(3,5)
      COMMON/URSAZ/ASBS(5),ZSBS(5)
      COMMON/EPZER/P00(5),E00(5)
      COMMON/UR2OUT/FAPP(3),FLAPP(3),DLFAR(3),DLFAT(3),
     *PIAP,EIAP,DPILR,DEILR,DPILT,DEILT

      J=NSUBST
      ROIL=ROILL(J)
      HLRO=HLROO(J)
      NRO=NROO(J)
      NTE2=NTEMM(1,J)
      NTE3=NTEMM(2,J)
      NTEM=NTEMM(3,J)
      TEMI1=TEMILL(1,J)
      TEMI2=TEMILL(2,J)
      TEMI3=TEMILL(3,J)
      HLT1=HLTEMM(1,J)
      HLT2=HLTEMM(2,J)
      HLT3=HLTEMM(3,J)
      MS=0
      MQ=1
      IF(J.EQ.1) GOTO 303
      J1=J-1
      DO 301 K=1,J1
      MQ=MQ+4*NROO(K)
 301  MS=MS+3*NROO(K)*NTEMM(3,K)
 303  CONTINUE
      ITAB=0
      XLD=DLOG(D)
      XX=(XLD-ROIL)/HLRO
      I=XX
      IF(JOB.GE.7) GOTO 400
      IF(I.GE.1) GOTO 1
      I=0
      GOTO 101
 1    IF(I.LE.NRO-3) GOTO 2
      I=NRO-2
      R=I
      XX=XX-R
      GOTO 101
 2    R=I
      XX=XX-R
 20   XY=DLOG(TE)
      ETAD=1.D0
      ETAU=1.D0
      IF(XY.GT.TEMI2) GOTO 3
      HLTEM=HLT1
      XY=(XY-TEMI1)/HLTEM
      K=XY
      IF(K.LE.0) K=0
      R=K
      XY=XY-R
      IF(K.EQ.NTE2-2) ETAU=HLT2/HLT1
      IF(K.GE.1) GOTO 7
      GOTO 102
 3    IF(XY.GT.TEMI3) GOTO 4
      HLTEM=HLT2
      XY=(XY-TEMI2)/HLTEM
      K=XY
      R=K
      XY=XY-R
      K=K+NTE2-1
      IF(K.EQ.NTE2-1) ETAD=HLT1/HLT2
      IF(K.EQ.NTE3-2) ETAU=HLT3/HLT2
      GOTO 7
 4    HLTEM=HLT3
      XY=(XY-TEMI3)/HLTEM
      K=XY
      IF(K.GE.NTEM-NTE3-1) K=NTEM-NTE3-1
      R=K
      XY=XY-R
      K=K+NTE3-1
      IF(K.EQ.NTE3-1) ETAD=HLT2/HLT3
      IF(K.LE.NTEM-3) GOTO 7
      GOTO 102
 7    CONTINUE
      M1=1
      IF(JOB.GE.4) M1=2
      IF(JOB.EQ.5) M1=3
      M2=M1+JOB-1
      IF(JOB.GE.4) M2=M1
      IF(JOB.GE.6) M2=M1+1

      DO 707 L=M1,M2
      MA=MS+((L-1)*NRO+I)*NTEM+K+1
      MB=MA+NTEM
      F00=FARR(MA)
      F01=FARR(MA+1)
      F10=FARR(MB)
      F11=FARR(MB+1)
      IF(ITAB.EQ.0) GOTO 71
      A20=0.D0
      A21=0.D0
      A02=0.D0
      A12=0.D0
      GOTO 72
 71   CONTINUE
      MA=MA-NTEM
      MB=MB+NTEM
      A20=.25D0*(FARR(MB)+FARR(MA)-F00-F10)
      A21=.25D0*(FARR(MB+1)+FARR(MA+1)-F01-F11)-A20
      ETAU1=1.D0/(1.D0+ETAU)
      ETAD1=1.D0/(1.D0+ETAD)
      MA=MA+NTEM-1
      MB=MB-NTEM-1
      A02=.5D0*(((FARR(MA+3)-F01)/ETAU+ETAU*(F01-F00))*ETAU1-
     *((F00-FARR(MA))/ETAD+ETAD*(F01-F00))*ETAD1)
      A12=.5D0*(((FARR(MB+3)-F11)/ETAU+ETAU*(F11-F10))*ETAU1-
     *((F10-FARR(MB))/ETAD+ETAD*(F11-F10))*ETAD1)-A02
 72   CONTINUE
      A11=F00+F11-F01-F10-A21-A12
 73   A01=F01-F00-A02
      A10=F10-F00-A20
      FLAPP(L)=F00+XX*(A10+XY*(A11+XY*A12)+XX*(A20+A21*XY))+
     +XY*(A01+A02*XY)
      DLFAR(L)=(A10+2.D0*XX*(A20+A21*XY)+XY*(A11+A12*XY))/HLRO
      DLFAT(L)=(A01+2.D0*XY*(A02+A12*XX)+XX*(A11+A21*XX))/HLTEM
      R=0.D0
      IF(L.EQ.2) R=P00(J)
      IF(L.EQ.3) R=E00(J)
      R1=DEXP(FLAPP(L))
      FAPP(L)=R1-R
      DLFAR(L)=DLFAR(L)*R1
      DLFAT(L)=DLFAT(L)*R1
 707  CONTINUE

 400  IF(JOB.EQ.1) GOTO 100
      IF(XLD.GT.ROIL) GOTO 403
      M=MQ
      R1=ROIL
      GOTO 405
 403  R1=ROIL+(NRO-1)*HLRO
      IF(XLD.LT.R1) GOTO 407
      M=MQ+4*(NRO-2)
 405  QL=(QURS(M)*R1*(3.D0*XLD-2.D0*R1)+QURS(M+1)*
     *(2.D0*XLD-R1))*R1+QURS(M+2)*XLD+QURS(M+3)
      DG=0.D0
      DQL=(3.D0*R1*QURS(M)+2.D0*QURS(M+1))*R1+QURS(M+2)
      GOTO 411
 407  M=MQ+4*I
      QL=((QURS(M)*XLD+QURS(M+1))*XLD+QURS(M+2))*XLD+QURS(M+3)
      DQL=(3.D0*QURS(M)*XLD+2.D0*QURS(M+1))*XLD+QURS(M+2)
      DG=3.D0*QURS(M)*XLD+QURS(M+1)
 411  Q=DEXP(QL)
      G1=1.5D0*DQL
      G3=G1+1.D0
      R1=TI+Q
      R2=D*TI/R1
      PIAP=R2*(TI+G3*Q)
      DPILR=R2*(TI*(TI+(1.D0+2.D0*G1*G1/3.D0)*Q)/R1+(G3+3.D0*DG)*Q)
      DPILT=R2*(TI*(TI+(1.D0-G1)*Q)/R1+G3*Q)
 444  EIAP=1.5D0*TI*(R1+Q)/R1
      DEILR=G1*Q*(TI/R1)**2
      DEILT=1.5D0*TI*(1.D0+(Q/R1)**2)
 100  RETURN
 101  ITAB=1
      GOTO 20
 102  ITAB=1
      GOTO 7
      END

C**********************************************************************
C  INUR2 initializes c/blks /URSAZ/,/URSTBL/,/URTBLI/,/EPZER/
C  for URSAPB; should be preceded by CALL LAUNCH(A,Z)
      SUBROUTINE INUR2(NSUBST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL FARR
      COMMON/URTBLI/QURS(1000)
      COMMON/URSTBL/FARR(30000)
      COMMON/URSEDG/NROO(6),ROILL(5),HLROO(5),NTEMM(3,6),
     *TEMILL(3,5),HLTEMM(3,5)
      COMMON/URSAZ/ASBS(5),ZSBS(5)
      COMMON/EPZER/P00(5),E00(5)
      COMMON/CONST/C1,C2,C3,C4,C5
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/BASIC/AL,BET,GAM,DEL,SGM,AMU,ANU,ALAM,Y0,Y1,B,EPS,DZ,NS
      COMMON/ENRG/EI,EE,EIZ,EPI,ES,EET,PPP
      COMMON/SGLASP/XI(250),WI(250),YI(250),DER(250),DER2(250),EDG(2)

      J=NSUBST
      ROIL=ROILL(J)
      HLRO=HLROO(J)
      NRO=NROO(J)
      NTE2=NTEMM(1,J)
      NTE3=NTEMM(2,J)
      NTEM=NTEMM(3,J)
      HLT1=HLTEMM(1,J)
      HLT2=HLTEMM(2,J)
      HLT3=HLTEMM(3,J)
      TEMI1=TEMILL(1,J)
      MT=NRO*NTEM
      TEMI2=TEMI1+HLT1*(NTE2-1)
      TEMI3=TEMI2+HLT2*(NTE3-NTE2)
      TEMILL(2,J)=TEMI2
      TEMILL(3,J)=TEMI3
      ASBS(J)=A
      ZSBS(J)=Z
      MS=0
      MQ=1
      IF(J.EQ.1) GOTO 303
      J1=J-1
      DO 301 K=1,J1
      MQ=MQ+4*NROO(K)
 301  MS=MS+3*NROO(K)*NTEMM(3,K)
 303  CONTINUE
      M1=MS+1
      M2=MS+3*MT
      M3=MQ+4*NRO-1
      P00(J)=-4.D0*BARA(.7D0*DZ,0.D0,0.D0)
      E00(J)=-2.D0*ENA(DZ,0.D0,0.D0)
      PRINT 511,J
      R1=ROIL+HLRO*(NRO-1)
      NR1=1
      PRINT 512,ROIL,R1,HLRO,NR1,NRO
      PRINT 513,TEMI1,TEMI2,HLT1,NR1,NTE2
      PRINT 513,TEMI2,TEMI3,HLT2,NTE2,NTE3
      R1=TEMI3+HLT3*(NTEM-NTE3)
      PRINT 513,TEMI3,R1,HLT3,NTE3,NTEM
      PRINT 514,P00(J),E00(J),M1,M2,MQ,M3
 511  FORMAT(/' Grid parameters (in a.u.) of the EOS table for substance
     * No.',I2)
 512  FORMAT(' LN(D)=',1PG11.4,' ->  ',1PG11.4,'  STEP=',1PG11.4
     *,' NODES:',I3,' -> ',I3)
 513  FORMAT(' LN(T)=',1PG11.4,' ->  ',1PG11.4,'  STEP=',1PG11.4
     *,' NODES:',I3,' -> ',I3)
 514  FORMAT(' P00,E00(IN A.U.)=',1P2E10.3,'; TABLE: FARR',2I6
     *,' QURS',2I5)
      DO 7 K=1,NTEM
      IF(K.GE.NTE2) GOTO 2
      T=DEXP(TEMI1+HLT1*(K-1))
      GOTO 4
 2    IF(K.GT.NTE3) GOTO 3
      T=DEXP(TEMI2+HLT2*(K-NTE2))
      GOTO 4
 3    T=DEXP(TEMI3+HLT3*(K-NTE3))
 4    CONTINUE
      DO 7 I=1,NRO
      D=DEXP(ROIL+HLRO*(I-1))
      MA=MS+(I-1)*NTEM+K
      MB=MA+MT
      MC=MB+MT
      FARR(MA)=DLOG(DINA(D,T))
      FARR(MC)=DLOG(ENA(D,T,0.D0)+E00(J))
      FARR(MB)=DLOG(PPP+P00(J))
 7    CONTINUE
      EDG(1)=0.D0
      EDG(2)=0.D0
      DO 21 K=1,NRO
      XI(K)=ROIL+(K-1)*HLRO
      R1=DINA(DEXP(XI(K)),0.D0)
      WI(K)=100.D0
 21   YI(K)=XI(K)/3.D0+DLOG(C3*ALAM*R1**2)
      CALL SMOOSP(NRO,XI,YI,WI,1.D0,EDG,QURS(MQ),DER,DER,DER2,2)
      RETURN
      END

C*****************************************************************
C  ZERA = the root X of equation  F(X)=0  in the interval XL<X<XR;
C  error in the value of X is less than EPS*(XR-XL);
      DOUBLE PRECISION FUNCTION ZERA(F,XL,XR,X0,EPS)
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
 77   ZERA=X
      RETURN

 511  FORMAT(' F(X).NE.0 in ZERA: NRCAL=',I6/' XL,X,XR=',1P3E15.8/
     *' H=',1PE15.8,' FP,F2=',1P2E15.8,' EPS=',1PE15.8)
      END

C**********************************************************************
C  GQUAD = Gauss quadrature = \int_A^B  F(X) dX
      DOUBLE PRECISION FUNCTION GQUAD(F,A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL F
      DIMENSION XG(222),WG(222),LN(12)
      DIMENSION X1(56),X2(56),X3(56),X4(54),W1(56),W2(56),W3(56),W4(54)
      EQUIVALENCE (XG(1),X1(1)),(XG(57),X2(1)),(XG(113),X3(1)),
     *(XG(169),X4(1))
      EQUIVALENCE (WG(1),W1(1)),(WG(57),W2(1)),(WG(113),W3(1)),
     *(WG(169),W4(1))
C  available N:
      DATA LN/4,8,12,16,20,24,32,40,48,64,80,96/
      DATA X1/3.39981043584856D-1,
     *8.61136311594053D-1,1.8343464249565D-1,
     *5.25532409916329D-1,7.96666477413627D-1,9.60289856497536D-1,
     *1.25233408511469D-1,3.67831498998180D-1,5.87317954286617D-1,
     *7.69902674194305D-1,9.04117256370475D-1,9.81560634246719D-1,
     *9.501250983763744D-2,2.816035507792589D-1,4.580167776572274D-1,
     *6.178762444026437D-1,7.554044083550030D-1,8.656312023878317D-1,
     *9.445750230732326D-1,9.894009349916499D-1,7.652652113349733D-2,
     *2.277858511416451D-1,3.737060887154196D-1,5.108670019508271D-1,
     *6.360536807265150D-1,7.463319064601508D-1,8.391169718222188D-1,
     *9.122344282513259D-1,9.639719272779138D-1,9.931285991850949D-1,
     *6.405689286260563D-2,1.911188674736163D-1,3.150426796961634D-1,
     *4.337935076260451D-1,5.454214713888395D-1,6.480936519369756D-1,
     *7.401241915785544D-1,8.200019859739029D-1,8.864155270044010D-1,
     *9.382745520027328D-1,9.747285559713095D-1,9.951872199970214D-1,
     *4.830766568773832D-2,1.444719615827965D-1,2.392873622521371D-1,
     *3.318686022821276D-1,4.213512761306353D-1,5.068999089322294D-1,
     *5.877157572407623D-1,6.630442669302152D-1,7.321821187402897D-1,
     *7.944837959679424D-1,8.493676137325699D-1,8.963211557660521D-1,
     *9.349060759377397D-1,9.647622555875064D-1/
      DATA W1/6.52145154862546D-1,3.47854845137454D-1,
     *3.62683783378362D-1,3.13706645877887D-1,2.22381034453374D-1,
     *1.01228536290376D-1,2.49147045813403D-1,2.33492536538355D-1,
     *2.03167426723066D-1,1.60078328543346D-1,1.06939325995318D-1,
     *4.717533638651200D-2,1.894506104550685D-1,1.826034150449236D-1,
     *1.691565193950025D-1,1.495959888165767D-1,1.246289712555339D-1,
     *9.515851168249278D-2,6.225352393864789D-2,2.715245941175409D-2,
     *1.527533871307258D-1,1.491729864726037D-1,1.420961093183821D-1,
     *1.316886384491766D-1,1.181945319615184D-1,1.019301198172404D-1,
     *8.327674157670475D-2,6.267204833410906D-2,4.060142980038694D-2,
     *1.761400713915212D-2,1.279381953467522D-1,1.258374563468283D-1,
     *1.216704729278034D-1,1.155056680537256D-1,1.074442701159656D-1,
     *9.761865210411389D-2,8.619016153195328D-2,7.334648141108031D-2,
     *5.929858491543678D-2,4.427743881741981D-2,2.853138862893366D-2,
     *1.234122979998720D-2,9.654008851472780D-2,9.563872007927486D-2,
     *9.384439908080457D-2,9.117387869576388D-2,8.765209300440381D-2,
     *8.331192422694676D-2,7.819389578707031D-2,7.234579410884851D-2,
     *6.582222277636185D-2,5.868409347853555D-2,5.099805926237618D-2,
     *4.283589802222668D-2,3.427386291302143D-2,2.539206530926206D-2/
      DATA X2/9.856115115452683D-1,9.972638618494816D-1,
     *3.877241750605082D-2,1.160840706752552D-1,1.926975807013711D-1,
     *2.681521850072537D-1,3.419940908257585D-1,4.137792043716050D-1,
     *4.830758016861787D-1,5.494671250951282D-1,6.125538896679802D-1,
     *6.719566846141795D-1,7.273182551899271D-1,7.783056514265194D-1,
     *8.246122308333117D-1,8.659595032122595D-1,9.020988069688743D-1,
     *9.328128082786765D-1,9.579168192137917D-1,9.772599499837743D-1,
     *9.907262386994570D-1,9.982377097105592D-1,3.238017096286936D-2,
     *9.700469920946270D-2,1.612223560688917D-1,2.247637903946891D-1,
     *2.873624873554556D-1,3.487558862921607D-1,4.086864819907167D-1,
     *4.669029047509584D-1,5.231609747222330D-1,5.772247260839727D-1,
     *6.288673967765136D-1,6.778723796326639D-1,7.240341309238147D-1,
     *7.671590325157403D-1,8.070662040294426D-1,8.435882616243935D-1,
     *8.765720202742479D-1,9.058791367155697D-1,9.313866907065543D-1,
     *9.529877031604309D-1,9.705915925462473D-1,9.841245837228269D-1,
     *9.935301722663508D-1,9.987710072524261D-1,2.435029266342443D-2,
     *7.299312178779904D-2,1.214628192961206D-1,1.696444204239928D-1,
     *2.174236437400071D-1,2.646871622087674D-1,3.113228719902110D-1,
     *3.572201583376681D-1,4.022701579639916D-1,4.463660172534641D-1/
      DATA W2/1.627439473090567D-2,7.0186100094700097D-3,
     *7.750594797842481D-2,7.703981816424797D-2,7.611036190062624D-2,
     *7.472316905796826D-2,7.288658239580406D-2,7.061164739128678D-2,
     *6.791204581523390D-2,6.480401345660104D-2,6.130624249292894D-2,
     *5.743976909939155D-2,5.322784698393682D-2,4.869580763507223D-2,
     *4.387090818567327D-2,3.878216797447202D-2,3.346019528254785D-2,
     *2.793700698002340D-2,2.224584919416696D-2,1.642105838190789D-2,
     *1.049828453115281D-2,4.521277098533191D-3,6.473769681268392D-2,
     *6.446616443595008D-2,6.392423858464819D-2,6.311419228625403D-2,
     *6.203942315989266D-2,6.070443916589388D-2,5.911483969839564D-2,
     *5.727729210040322D-2,5.519950369998416D-2,5.289018948519367D-2,
     *5.035903555385447D-2,4.761665849249047D-2,4.467456085669428D-2,
     *4.154508294346475D-2,3.824135106583071D-2,3.477722256477044D-2,
     *3.116722783279809D-2,2.742650970835695D-2,2.357076083932438D-2,
     *1.961616045735553D-2,1.557931572294385D-2,1.147723457923454D-2,
     *7.327553901276262D-3,3.153346052305839D-3,4.869095700913972D-2,
     *4.857546744150343D-2,4.834476223480296D-2,4.799938859645831D-2,
     *4.754016571483031D-2,4.696818281621002D-2,4.628479658131442D-2,
     *4.549162792741814D-2,4.459055816375656D-2,4.358372452932345D-2/
      DATA X3/4.894031457070530D-1,5.312794640198945D-1,
     *5.718956462026340D-1,6.111553551723933D-1,6.489654712546573D-1,
     *6.852363130542332D-1,7.198818501716108D-1,7.528199072605319D-1,
     *7.839723589433414D-1,8.132653151227976D-1,8.406292962525804D-1,
     *8.659993981540928D-1,8.893154459951141D-1,9.105221370785028D-1,
     *9.295691721319396D-1,9.464113748584028D-1,9.610087996520537D-1,
     *9.733268277899110D-1,9.833362538846260D-1,9.910133714767443D-1,
     *9.963401167719553D-1,9.993050417357721D-1,1.951138325679400D-2,
     *5.850443715242067D-2,9.740839844158460D-2,1.361640228091439D-1,
     *1.747122918326468D-1,2.129945028576661D-1,2.509523583922721D-1,
     *2.885280548845119D-1,3.256643707477019D-1,3.623047534994873D-1,
     *3.983934058819692D-1,4.338753708317561D-1,4.686966151705445D-1,
     *5.028041118887850D-1,5.361459208971319D-1,5.686712681227098D-1,
     *6.003306228297517D-1,6.310757730468720D-1,6.608598989861198D-1,
     *6.896376443420276D-1,7.173651853620999D-1,7.440002975835973D-1,
     *7.695024201350414D-1,7.938327175046054D-1,8.169541386814635D-1,
     *8.388314735802553D-1,8.594314066631111D-1,8.787225676782138D-1,
     *8.966755794387707D-1,9.132631025717577D-1,9.284598771724458D-1,
     *9.422427613098727D-1,9.545907663436349D-1,9.654850890437993D-1/
      DATA X4/9.749091405857278D-1,9.828485727386291D-1,
     *9.892913024997555D-1,9.942275409656883D-1,9.976498643982377D-1,
     *9.995538226516306D-1,1.627674484960297D-2,4.881298513604973D-2,
     *8.129749546442556D-2,1.136958501106659D-1,1.459737146548969D-1,
     *1.780968823676186D-1,2.100313104605672D-1,2.417431561638400D-1,
     *2.731988125910491D-1,3.043649443544964D-1,3.352085228926254D-1,
     *3.656968614723136D-1,3.957976498289086D-1,4.254789884073005D-1,
     *4.547094221677430D-1,4.834579739205964D-1,5.116941771546677D-1,
     *5.393881083243574D-1,5.665104185613972D-1,5.930323647775721D-1,
     *6.189258401254686D-1,6.441634037849671D-1,6.687183100439162D-1,
     *6.925645366421716D-1,7.156768123489676D-1,7.380306437444001D-1,
     *7.596023411766475D-1,7.803690438674332D-1,8.003087441391408D-1,
     *8.194003107379317D-1,8.376235112281871D-1,8.549590334346015D-1,
     *8.713885059092965D-1,8.868945174024204D-1,9.014606353158523D-1,
     *9.150714231208981D-1,9.277124567223087D-1,9.393703397527552D-1,
     *9.500327177844376D-1,9.596882914487425D-1,9.683268284632642D-1,
     *9.759391745851365D-1,9.825172635630147D-1,9.880541263296238D-1,
     *9.925439003237626D-1,9.959818429872093D-1,9.983643758631817D-1,
     *9.996895038832308D-1/
      DATA W3/4.247351512365359D-2,4.126256324262353D-2,
     *3.995374113272034D-2,3.855015317861563D-2,3.705512854024005D-2,
     *3.547221325688238D-2,3.380516183714161D-2,3.205792835485155D-2,
     *3.023465707240248D-2,2.833967261425948D-2,2.637746971505466D-2,
     *2.435270256871087D-2,2.227017380838325D-2,2.013482315353021D-2,
     *1.795171577569734D-2,1.572603047602472D-2,1.346304789671864D-2,
     *1.116813946013113D-2,8.846759826363948D-3,6.504457968978363D-3,
     *4.147033260562468D-3,1.783280721696433D-3,3.901781365630665D-2,
     *3.895839596276953D-2,3.883965105905197D-2,3.866175977407646D-2,
     *3.842499300695942D-2,3.812971131447764D-2,3.777636436200140D-2,
     *3.736549023873049D-2,3.689771463827601D-2,3.637374990583598D-2,
     *3.579439395341605D-2,3.516052904474759D-2,3.447312045175393D-2,
     *3.373321498461152D-2,3.294193939764540D-2,3.210049867348777D-2,
     *3.121017418811470D-2,3.027232175955798D-2,2.928836958326785D-2,
     *2.825981605727686D-2,2.718822750048638D-2,2.607523576756512D-2,
     *2.492253576411549D-2,2.373188286593010D-2,2.250509024633246D-2,
     *2.124402611578201D-2,1.995061087814200D-2,1.862681420829903D-2,
     *1.727465205626931D-2,1.589618358372569D-2,1.449350804050908D-2,
     *1.306876159240134D-2,1.162411412079783D-2,1.016176604110306D-2/
      DATA W4/8.683945269260858D-3,7.192904768117313D-3,
     *5.690922451403199D-3,4.180313124694895D-3,2.663533589512682D-3,
     *1.144950003186942D-3,3.255061449236317D-2,3.251611871386884D-2,
     *3.244716371406427D-2,3.234382256857593D-2,3.220620479403025D-2,
     *3.203445623199266D-2,3.182875889441101D-2,3.158933077072717D-2,
     *3.131642559686136D-2,3.101033258631384D-2,3.067137612366915D-2,
     *3.029991542082759D-2,2.989634413632839D-2,2.946108995816791D-2,
     *2.899461415055524D-2,2.849741106508539D-2,2.797000761684833D-2,
     *2.741296272602924D-2,2.682686672559176D-2,2.621234073567241D-2,
     *2.557003600534936D-2,2.490063322248361D-2,2.420484179236469D-2,
     *2.348339908592622D-2,2.273706965832937D-2,2.196664443874435D-2,
     *2.117293989219130D-2,2.035679715433332D-2,1.951908114014502D-2,
     *1.866067962741147D-2,1.778250231604526D-2,1.688547986424517D-2,
     *1.597056290256229D-2,1.503872102699494D-2,1.409094177231486D-2,
     *1.312822956696157D-2,1.215160467108832D-2,1.116210209983850D-2,
     *1.016077053500842D-2,9.148671230783387D-3,8.126876925698759D-3,
     *7.096470791153865D-3,6.058545504235962D-3,5.014202742927518D-3,
     *3.964554338444687D-3,2.910731817934946D-3,1.853960788946922D-3,
     *7.967920655520124D-4/
      NI=0
      DO 10 K=1,12
      IF(N.EQ.LN(K)) GOTO 20
 10   NI=NI+LN(K)/2
      PRINT 510,N,LN
      STOP
 510  FORMAT(//' NO  SUCH  N  IN  GQUAD: N=',I6/
     *' WHILE  AVAILABLE  ONLY  N=',12I4)
 20   XA=.5D0*(A+B)
      XB=.5D0*(B-A)
      S=0.D0
      J=N/2
      DO 30 K=1,J
      NI=NI+1
      R=F(XA+XB*XG(NI))+F(XA-XB*XG(NI))
 30   S=S+R*WG(NI)
      GQUAD=S*XB
      RETURN
      END

C*********************************************************************
C  BESTO2 = stopping power in units GeV*mm**2/mg
C  as a function of ion energy E1G (in GeV),
C  density DG (in g/cc), electron TEG and ion TIG temperatures (in keV),
C  and substance number  NSUB
      DOUBLE PRECISION FUNCTION BESTO2(E1G,DG,TEG,TIG,NSUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/URSAZ/ATARG(5),ZTARG(5)
      COMMON/BECNST/C11,C12,C13,C14,C15,C16,C17,C1R
      COMMON/BEAM1/A1,Z1,B1,BET1,SGM1,AMU1,PTIF1,DINP1,APIN1(4,100)
      COMMON/BEAM2/POT2(350),NSH2(6),NE2(100),ESH2(100),GB0(5),B2(5),
     *BET2(5),SGM2(5),AMU2(5)
      COMMON/UR2OUT/YAP,PEAP,EEAP,YLAP,PELAP,EELAP,DYLR,DPELR,DEELR,
     *DYLT,DPELT,DEELT,PIAP,EIAP,DPILR,DEILR,DPILT,DEILT
      COMMON/BESOUT/CLBE,SBE,CLFE,SFE,CLFI,SFI,CLNU,SNU,STPE,STPI,
     *Y2,Y300,Z1EF,NIBES

      H(X)=LOG(1.D0+X/(1.D0 +3.5D0/SQRT(X)))
      ETA(X)=.353D0+X**2*((X**3+2.34D0)/(X**3+11.D0))
      E1=E1G/(4.96D-5*A1)
 2    CLBE=0.D0
      SBE=0.D0
      CLFE=0.D0
      SFE=0.D0
      CLFI=0.D0
      SFI=0.D0
      CLNU=0.D0
      SNU=0.D0
      IF(E1-1.D-8) 999,999,10
 10   A2=ATARG(NSUB)
      Z2=ZTARG(NSUB)
      R=Z2+.5D0
      IZ2=R
      D=DG/(11.206D0*A2)
      TE=TEG/27.21D-3
      TI=TIG/27.21D-3
      CALL URSAPB(NSUB,D,9.5D-4,9.5D-4,1)
      Y300=YAP
      CALL URSAPB(NSUB,D,TE,TI,1)
      Y2=YAP
      IF(Y2.LT.1.D-6) Y2=1.D-6
      RY2=(D*Y2)**C11
      VV1=E1*(2.D0+C16*E1)/(1.D0+C16*E1)**2
      DCL=2.D0*LOG(1.D0+C16*E1)-C16*VV1

C: evaluate Z1EF:
      Z1EF=1.D0
      YT=1.D0
      IF(Z1.LT.1.5D0) GOTO 67
      Z1EF=Z1/(1.D0+(.62D0*Z1**C11/SQRT(VV1))**1.7D0)**.58824D0
      YTN=C1R
      IF(YTN.LT.1.D-8) YTN=.5D0*Z1
      IF(YTN.GT.Z1-1.D-8) YTN=.5D0*Z1
      HM=.2D0*Z1
      NIBES=0
      R=C12*RY2
      R1=AMU1*TE**2
      RMU=R+TE*(TE/(TE+.4D0*R)-1.5D0*LOG(1.D0+2.5D0*TE/R))
 40   NIBES=NIBES+1
      YT=YTN
      R=Y2*D/YT
      R4=R**BET1
      IF(R*1.D-8-1.D0) 41,43,43
 41   R2=R**SGM1
      IF(R2-1.D-16) 42,44,44
 42   IF(R1-1.D-16) 43,44,44
 43   R3=1.D0
      GOTO 45
 44   R3=R2/(R2+R1)
 45   BB1=B1*R4*R3
      DB1=BB1*(BET1+SGM1*(1.D0-R3))/YT
      IF(YT-Z1+1.D0) 48,50,50
 48   K=YT
      K=K+1
      PO=APIN1(4,K)+YT*(APIN1(3,K)+YT*(APIN1(2,K)+YT*APIN1(1,K)))
      DPO=APIN1(3,K)+2.D0*YT*APIN1(2,K)+3.D0*YT**2*APIN1(1,K)
      GOTO 51
 50   R3=PTIF1/(Z1-YT)**DINP1
      PO=R3*YT
      DPO=R3*(1.D0+DINP1*YT/(Z1-YT))
 51   F2=RMU+PO-BB1
      DF=DPO+DB1
      IF(ABS(F2)-ABS(.01D0*YT*DF)) 66,66,52
 66   IF(ABS(F2)-ABS(1.D-8+.1D0*PO)) 67,67,52
 52   IF(YT-1.D-6) 53,53,54
 53   IF(F2) 54,67,67
 54   IF(Z1-YT-1.D-6) 55,55,56
 55   IF(F2) 67,67,56
 56   IF(NIBES.EQ.1) F1=F2/ABS(F2)
      IF(F1*F2) 57,67,58
 57   HM=.3D0*HM
      IF(HM-.003D0) 67,67,58
 58   IF(ABS(F2)-ABS(DF*HM)) 59,60,60
 59   YTN=YT-F2/DF
      GOTO 61
 60   YTN=YT-F2*HM/ABS(F2)
 61   IF(YTN-1.D-8) 62,62,63
 62   YTN=.5D0*(YT+1.D-8)
      GOTO 65
 63   IF(Z1-YTN-1.D-8) 64,64,65
 64   YTN=.5D0*(YT+Z1-1.D-8)
 65   F1=F2/ABS(F2)
      IF(NIBES.LT.100) GOTO 40
      PRINT 511,D,TE,Y2,YT,YTN,F2,F1,HM
      YT=0.D0
 511  FORMAT(' DIVERG.IN BESTOP,Z1-ION.:',1P8G9.3)
 67   C1R=YT
      IF(Z1EF.LT.YT) Z1EF=YT

C: bound-electron stopping power:
      K=1
      KSH=0
      KPT=0
 110  IF(K.EQ.NSUB) GOTO 120
      KSH=KSH+NSH2(K)
      R=ZTARG(K)+.5D0
      N1=R
      KPT=KPT+N1
      K=K+1
      GOTO 110
 120  Q=Z2-Y2
      IF(Q.LT.1.D-8) Q=1.D-8
      IQ=Q
      IF(IQ.GE.IZ2) IQ=IZ2-1
      NOST=0
      KSS=KSH
 130  KSS=KSS+1
      NOST=NOST+NE2(KSS)
      IF(NOST.LT.IQ+1) GOTO 130
      NOST=NOST-NE2(KSS)
      KS=KSS
      IF(NOST.EQ.IQ)KS=KSS-1
      IY=IZ2-IQ+KPT
      R=D**SGM2(NSUB)
      R1=AMU2(NSUB)*TE*TE
      IF(R-1.D-16) 140,150,150
 140  IF(R1-1.D-16) 144,150,150
 144  R=1.D0
      GOTO 160
 150  R=R/(R+R1)
 160  BB2=B2(NSUB)*R*D**BET2(NSUB)
      R=Q-IQ
      IF(IQ.EQ.0) GOTO 170
      DESH=R*(POT2(IY)-ESH2(KSS)-BB2)+(1.D0-R)*
     *(POT2(IY+1)-ESH2(KS)-BB2)
      IF(IZ2.EQ.1) GOTO 166
      IF(Y2-Y300) 162,162,164
 162  GB=GB0(NSUB)
      GOTO 180
 164  GB=GB0(NSUB)+(1.105D0-GB0(NSUB))*(Y2-Y300)/(Z2-Y300-1.D0)
      GOTO 180
 166  GB=1.105D0
      GOTO 180
 170  DESH=POT2(IY)-ESH2(KSS)-BB2
      GB=1.105D0
 180  IF(DESH.LT.0.D0) DESH=0.D0
      K=KSS
      R=2.D0*VV1/(GB*(ESH2(K)+DESH))
      R=R/SQRT(1.D0+C14*Z1EF**2/VV1)
      ZLB=(Q-NOST)*H(R)
 186  K=K-1
      IF(K.EQ.KSH) GOTO 190
      R=2.D0*VV1/(GB*(ESH2(K)+DESH))
      R=R/SQRT(1.D0+C14*Z1EF**2/VV1)
      ZLB=ZLB+NE2(K)*H(R)
      GOTO 186
 190  CLBE=ZLB/Q
      ZLB=ZLB+Q*DCL
      SBE=4.D0*C13*Z1EF**2*ZLB/(A1*VV1)

C: free-electron stoppig power:
      R=911.5D0*A2*VV1
      IF(R-1.D8*TI) 202,204,204
 202  TI2=TI
      GOTO 210
 204  TI2=1.D-8*R
 210  XI=SQRT(R/TI2)
      TEF=SQRT(TE*TE+(1.26D0*C13*RY2)**2)
      R=.5D0*VV1
      IF(R-1.D8*TEF) 212,214,214
 212  TE2=TEF
      GOTO 220
 214  TE2=1.D-8*R
 220  XE=SQRT(R/TE2)
      VE2=2.D0*TE2*ETA(XE)
      R=4.D0*C13*D*Y2*(1.D0/VE2+Y2/(2.D0*TI2*ETA(XI)))
      R=1.D0/R
      R1=(.75D0*Z1EF/C13)**C11/RY2
      IF(R.LT.R1) R=R1
      R1=4.D0*R*VE2/(1.D0+C14*Z1EF**2/VE2)
      R1=SQRT(R1)
      CLFE=LOG(1.D0+R1/(1.D0+.5D0/SQRT(R1)))
      SFE=4.D0*C13*Z1EF**2*(XE**3/(XE**3+1.33D0))/(A1*VV1)
      SFE=SFE*Y2*(CLFE+DCL)
      IF(C17.GT.1823.D0) GOTO 999

C: ion+nuclear stopping power:
      R2=1.D0/(1823.D0*A2)
      VI2=2.D0*R2*TI2*ETA(XI)
      R3=(2.D0*1823.D0*A2*A1/(A2+A1))**2*VI2
      R1=Z1-Z1EF
      IF(R1.LT.1.D0) R1=1.D0
      R4=(Z1/R1**C11)**2
      R1=Z2-Y2
      IF(R1.LT.1.D0) R1=1.D0
      R4=.26D0/(R4+(Z2/R1**C11)**2)
      R1=SQRT(R3*R4/(1.D0+C14*(Z1*Z2)**2/VI2))
      CLNU=LOG(1.D0+R1/(1.D0+.35D0/SQRT(R1)))
      R2=4.D0*C13*R2*(XI**3/(XI**3+1.33D0))/(A1*VV1)
      R2=R2*SQRT(1.D0+A2/A1)
      SNU=R2*(Z1*Z2)**2*CLNU
      R1=(1.D0+C14*(Z1EF*Y2)**2/VI2)/R3
      IF(R1.LT.R4) R1=R4
      R1=SQRT(R/R1)
      CLFI=LOG(1.D0+R1/(1.D0+.5D0/SQRT(R1)))
      SFI=R2*(Z1EF*Y2)**2*CLFI
 999  BESTO2=4.589D-2*A1*(SBE+SFE+SFI+SNU)/A2
      RETURN
      END

C*********************************************************************
C  INOPAC initializes c/blks /OPAC1/ and /OPAC2/ required for
C  calculating Rosseland and Planck mean free paths of photons
      SUBROUTINE INOPAC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ELEMNT/A,Z,POT(100),LAEOS,LZEOS
      COMMON/URSAZ/ATARG(5),ZTARG(5)
      COMMON/OPAC1/POTOFI(100,5),BOPAC(100,5)
      COMMON/OPAC2/NLAGER,NWWWW,XLAGER(12),WLAGER(12),RWEI(12),PWEI(12)

      DIMENSION LL(4),XL(24),WL(24)
      DATA LL/6,12,0,0/
      DATA XL/
     *.22284 66041 79D0, 1.18893 21016 73D0, 2.99273 63260 59D0,
     *5.77514 35691 05D0, 9.83746 74183 83D0, 15.98287 39806 02D0,
     *.11572 21173 58D0, .61175 74845 15D0, 1.51261 02697 76D0,
     *2.83375 13377 44D0, 4.59922 76394 18D0, 6.84452 54531 15D0,
     *9.62131 68424 57D0, 13.00605 49933 06D0, 17.11685 51874 62D0,
     *22.15109 03793 97D0, 28.48796 72509 84D0, 37.09912 10444 67D0,
     *6*0.D0/
      DATA WL/
     *4.58964 673950D-1, 4.17000 830772D-1, 1.13373 382074D-1,
     *1.03991 974531D-2, 2.61017 202815D-4, 8.98547 906430D-7,
     *2.64731 371055D-1, 3.77759 275873D-1, 2.44082 011320D-1,
     *9.04492 222117D-2, 2.01023 811546D-2, 2.66397 354187D-3,
     *2.03231 592663D-4, 8.36505 585682D-6, 1.66849 387654D-7,
     *1.34239 103052D-9, 3.06160 163504D-12, 8.14807 746743D-16,
     *6*0.D0/
      PI=3.1415 92653 58979D0

C:  fill c/blk /OPAC2/:
      NLAGER=12

      NI=0
      DO 10 K=1,4
      IF(NLAGER.EQ.LL(K)) GOTO 20
 10   NI=NI+LL(K)
      PRINT 9010,NLAGER,LL
      STOP

 20   DO 30 J=1,NLAGER
      NI=NI+1
      XLAGER(J)=XL(NI)
      WLAGER(J)=WL(NI)
      RAB=DEXP(-XLAGER(J))
      RWEI(J)=WLAGER(J)*(3.75D0/PI**4)*XLAGER(J)**4/(1.D0-RAB)**2
      PWEI(J)=WLAGER(J)*(15.D0/PI**4)*XLAGER(J)**3/(1.D0-RAB)
 30   CONTINUE

C: fill c/blk /OPAC1/:
      DO 600 KK=1,5
      Z=ZTARG(KK)
      RAB=Z+.5D0
      IZ=RAB
      IF(IZ.LT.1) GOTO 600
      A=ATARG(KK)
      CALL PINTA
      DO 590 I=1,IZ
C: ionization potentials in keV:
      POTOFI(I,KK)=POT(I+1)*2.721D-2
C:  calculate constants b_i:
      IF(I.NE.IZ) GOTO 160
      BOPAC(I,KK)=0.7495D0
      GOTO 590
 160  ALAM=4.8368D0*2.D0*Z**4*1.2504D-5/(POTOFI(I,KK)**2*(IZ+1-I))
C: iterations:
      NITER=0
      X0=ALAM
 200  NITER=NITER+1
      IF(NITER.LT.30) GOTO 220
      PRINT 9040,ALAM,X0,X1
 220  IF(DABS(1.D0-X0).LT.1.D-9) GOTO 240
      RAB=X0*DLOG(X0)/(X0-1.D0)
      GOTO 250
 240  RAB=1.D0
 250  X1=ALAM*(1.D0+1.5D0*RAB)
      IF(DABS(1.D0-X1/X0).LT.1.D-6) GOTO 280
      X0=X1
      IF(NITER.LT.30) GOTO 200
 280  BOPAC(I,KK)=X1
 590  CONTINUE
 600  CONTINUE

      RETURN
 9010 FORMAT(/' NO SUCH NLAGER IN LAGERR TABLES: NLAGER=',I6/
     *' WHILE AVAILABLE ONLY NLAGER=',12I4)
 9020 FORMAT(1P4E19.12)
 9030 FORMAT(/)
 9040 FORMAT(' INOPAC: more than 30 iterations, ALAM,X0,X1=',
     *1P3E16.9)
 9110 FORMAT(1P5E14.7)
      END
