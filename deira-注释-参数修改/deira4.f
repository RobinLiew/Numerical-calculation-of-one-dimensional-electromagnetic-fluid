C                  **************************
C                  *                        *
C                  *     D E I R A - 4      *
C                  *                        *
C**********************************************************************
C   Version: 4.2  Date: 24 Nov 06   Author: M.Basko
C   Last modification: 13 Dec 06
C----------------------------------------------------------------------
C  1-D 3-T hydro code for simulating ion-driven ICF targets
C  Compared to version 4: 
C     1) nuclear scattering has been accounted for in the transport-
C        relaxation of the energy of 14-MeV protons (Nov 03);
C     2) a new option (IFN14=4) for 14-MeV neutron heating has been
C        added (Nov 06): the neutron energy is spread according to the
C        diffusion profile in spherical geometry.
C**********************************************************************
C                  I N S T R U C T I O N S                            *
C**********************************************************************
C >>> Input files: <<<                                                *
C  Minimum configuration: 2 files: 'deira4.for' and 'de4com.fi';      *
C  If tabular equation of state (EOS) is used, or IFOPAC=3 option is  *
C  chosen for opacity, file 'EOSTA3' with the EOS and                 *
C  opacity tables should be prepared by the TABIN3 package.           *
C  If the ion drive is used, file 'BEMTA3' with stopping data should  *
C  be prepared by the TABIN3 package.                                 *
C                                                                     *
C >>> To start a new run: <<<                                         *
C  1. assign appropriate NN and NNZ values (both even numbers)        *
C     in the PARAMETER statement of file "de4com.fi";                 *
C  2. set ISTART=0;                                                   *
C  3. assign run parameters and set the initial target structure by   *
C     editing the INIDAT block; if necessary, edit the initial density*
C     profile in subroutine RO0PR;                                    *
C                                                                     *
C >>> To continue an old run (saved in file 'DEISA3'): <<<            *
C  1. set ISTART=1;                                                   *
C                                                                     *
C >>> To modify EOS: <<<                                              *
C  1. write a new version of subroutine EOS; as an example, see       *
C     EOSIG - equation of state of the ideal Boltzmann gas.           *
C                                                                     *
C >>> To modify target drive: <<<                                     *
C  1. write a new version of subroutine DRIVE.                        *
C                                                                     *
C >>> To modify transport coefficients: <<<                           *
C  1. write a new version of subroutine RELCON.                       *
C                                                                     *
C >>> Boundary conditions: <<<                                        *
C  Following types of boundary conditions are envisaged:              *
C  1) right boundary, r=R(N+1): prescribed boundary pressure PBR(t),  *
C     external rad. temperature TREX(t) and magnetic field HZBR(t);   *
C  2) left boundary, r=R(1):                                          *
C     IFLBND=0  -> fixed center of symmetry, R(1)=U(1)=0;             *
C                  all diffusion fluxes are zero;                     *
C     IFLBND=1  -> a closed void cavity at 0 < r < R(1) with          *
C                  a prescribed boundary pressure PBL(t);             *
C                  all diffusion fluxes are zero;                     *
C     IFLBND=-1 -> open halfspace at r < R(1)[allowed for IGEO=0 only]:
C                  prescribed boundary pressure PBL(t), external rad. *
C                  temperature TRLEX(t) and magnetic field HZBL(t);   *
C                  in this case the artificial t-viscosity is turned  *
C                  off; to compensate -> increase the value of SMU2.  *
C  The time dependence of the boundary values of the relevant physical*
C  quantities must be programmed in the subroutine BNDVAL(t).         *
C.....................................................................*
C  Accuracy and artificial viscosity constants, sensitivity thresholds,
C  and flux limits are controlled via the constants in the BLOCK DATA.*
C                                                                     *
C  Any block in the code between markers      ---------- begin #      *
C  and     =========== end #   can be erased without repercussions.   *

C**********************************************************************
C                             OPTIONS                                 
C**********************************************************************
C >>> Geometry: <<<
C IGEO = 0 -> plane-parallel hydrodynamic flow,
C      = 1 -> cylindrical flow,
C      = 2 -> spherical flow.

C >>> Boundary conditions: <<<
C IFLBND = 0 -> center of symmetry at the left boudary: R(1)=U(1)=0;
C        = 1 -> a closed void cavity at 0 < r < R(1): P_bl=P_bl(t);
C        =-1 -> open halfspace at r < R(1): P_bl=P_bl(t);
C               allowed with IGEO=0 only! artificial t-viscosity off!

C >>> Magnetic field: <<<
C IFHZ = 0 -> no equation for the magnetic field is solved;
C      = 1 -> a diffusion equation for z-component of the magnetic
C             field is solved; meaningful for IGEO=0,1 only.

C >>> Burn: <<<
C IFBURN =-1 -> burn bypassed, no check for ignition;
C IFBURN = 0 -> burn bypassed, check for ignition;
C        = 1 -> intermediate value when turning the burn on;
C        = 2 -> burn calculated, no check for being stopped.

C >>> Charged fusion products: <<<
C IFAL(IFP3,IFP14) = 0 -> ignored;
C                  = 1 -> diffusion;
C                  = 2 (or any other) -> local;
C NALPRI = 1 -> print EAL, HIAL;
C        = 2 -> print EP3, HIP3;
C        = 3 -> print EP14, HIP14;

C >>> Neutron heating of the central fuel "sphere": <<<
C IFN14(IFN2) = 0 -> ignored;
C             = 1 -> 1-st scattering;
C             = 2 -> local;
C             = 3 -> spread uniformly;
C IFN14       = 4 -> sphere, spread according to the diffusion profile;
C IFN2        = 4 -> reserved for the future;
C IFN14(IFN2) > 4 -> not allowed.

C >>> Neutron heating of possible outer fuel layers: <<<
C IFN14(IFN2) = 0 -> ignored;
C             > 0 -> local;
C Neutron heating of the non-fuel layers is not accounted for.

C >>> Opacity: <<<
C IFOPAC < 2 - default (zero) values of CAPROU and HIEROU;
C        = 2 - calculate CAPROU and HIEROU from fast formulas of DEIRA-2
C        = 3 - calculate CAPROU and HIEROU by integrating absorption
C              cross-section (model DEIRA-3);
C        > 3 - default (zero) values of CAPROU and HIEROU.

C**********************************************************************
C**********************************************************************
      PROGRAM DEIRA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL FARR
      LOGICAL LEXIST,LOPEN
      CHARACTER FNAME*12,FSTAT*3
      INCLUDE 'de4com.fi'

C: c/blks for the equation of state, filled from file 'EOSTA3':
      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/OPAC1/POTOFI(100,5),BOPAC(100,5)
      COMMON/OPAC2/NLAGER,NWWWW,XLAGER(12),WLAGER(12),RWEI(12),PWEI(12)
      COMMON/EPZER/P00(5),E00(5)
      COMMON/URSEDG/NROO(6),ROILL(5),HLROO(5)
     *,NTEMM(3,6),TEMILL(3,5),HLTEMM(3,5)
      COMMON/URSTBL/FARR(30000)
      COMMON/URTBLI/QURS(1000)
C  Here one should watch that
C  3*(NROO(1)*NTEMM(3,1)+...+NROO(5)*NTEMM(3,5)) < 30000, and
C  4*(NROO(1)+...+NROO(5)) < 1000

C: c/blks for the beam, filled from file 'BEMTA3':
      COMMON/BECNST/CBEAM(8)
      COMMON/BEAM1/ABEAM,ZBEAM,CBEAM1(406)
      COMMON/BEAM2/CBEAM2(350),NCBEM2(106),CBEM22(125)

C: c/blks for profile and history plots:
C: rates of energy transfer (powers) across layer interfaces:
C      COMMON/LEAKPW/PLKW(1:NNZ+1),PLKE(1:NNZ+1),PLKR(1:NNZ+1)
C     *,PLKAL(1:NNZ+1),PLKP3(1:NNZ+1),PLKP14(1:NNZ+1)
C      REAL APROF(10,1:NN)
c      PARAMETER(NNHIS=300,JRHIS=40)
C      REAL TIMHIS(1:NNHIS),AHIS(10,1:NNHIS),RHIS(JRHIS,1:NNHIS)

C: Array of prescribed times for printout (profile plots):
      DIMENSION TIMPRI(15)
      DATA TIMPRI/1.0D0,2.5D0,13*1.D2/

C: assign default values:
      CALL DEFVAL

C: set control parameters for profile plots (at TIME=TPRIN);
C: the K order in APROF(K,J) is TI,TE,TR,RO,P,U,...
C      NPROF=15
C      IPROF=0
C      FNAME='p.dat'
C      FSTAT='OLD'
C      INQUIRE(FILE=FNAME,EXIST=LEXIST)
C      IF(.NOT.LEXIST) FSTAT='NEW'
C      OPEN(21,FILE=FNAME,STATUS=FSTAT,ACCESS='SEQUENTIAL'
C     *,FORM='FORMATTED')
C      WRITE (21,9120) N
C: set control parameters for history plots:
C      THIPLO=0.D0
C      DTHIPL=1.D-2
CC      NDTPLO=40
C      IHIS=0
C      FNAME='ahis.dat'
C      FSTAT='OLD'
C      INQUIRE(FILE=FNAME,EXIST=LEXIST)
C      IF(.NOT.LEXIST) FSTAT='NEW'
C      OPEN(22,FILE=FNAME,STATUS=FSTAT,ACCESS='SEQUENTIAL'
C     *,FORM='FORMATTED')
C      FNAME='rhis.dat'
C      FSTAT='OLD'
C      INQUIRE(FILE=FNAME,EXIST=LEXIST)
C      IF(.NOT.LEXIST) FSTAT='NEW'
C      OPEN(23,FILE=FNAME,STATUS=FSTAT,ACCESS='SEQUENTIAL'
C     *,FORM='FORMATTED')

C: print DEIRA units and explanations to the output:
      CALL PRINUN
C------------------------------------------------------------ begin RUN
 1000 CONTINUE
      IFOPAC=2
C-------------------------------------------------------- begin READEOS
      INQUIRE(FILE='EOSTA3',EXIST=LEXIST)
      IF(.NOT.LEXIST) GOTO 1020
      IFOPAC=3
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
      CALL PRINUS
 1020 CONTINUE
C========================================================== end READEOS
C......................................................................
C  If a non-tabular EOS is used and file 'EOSTA3' does not exist, the .
C  values of ASUB(I), ZSUB(I) should be explicitly assigned here - to .
C  be used when calculating transport coefficients in KINBUR, RELCON; .
C......................................................................

C: Assign the values of the fit parameter FITCAP for the heat and
C  electrical conductivity at room temperature:
      DO 1030 K=1,5
      RAB=ZSUB(K)+.5D0
      IZSUBK=RAB
      IF(IZSUBK.EQ.29) FITCAP(K)=3.D0
      IF(IZSUBK.EQ.79) FITCAP(K)=6.D0
      IF(IZSUBK.EQ.74) FITCAP(K)=10.D0
 1030 CONTINUE

C------------------------------------------------------- begin READBEAM
      INQUIRE(FILE='BEMTA3',EXIST=LEXIST)
      IF(.NOT.LEXIST) GOTO 1040
      OPEN(12,FILE='BEMTA3',STATUS='OLD',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
      READ(12,9999) CBEAM,ABEAM,ZBEAM,CBEAM1,CBEAM2
      READ(12,9998) NCBEM2
      READ(12,9999) CBEM22
      CLOSE(12)
      CALL PRINBE
 1040 CONTINUE
C========================================================= end READBEAM

C--------------------------------------------------------- begin INIDAT
C: variables not present in the "hard-core" c/blks
C  and not saved in 'DEISA3' (for any run):
      ISTART=0
      TFIN=5.D0
      NTFIN=20 000
      TPRIN=TIMPRI(1)
C     DTPRIN=1.D18

      IF(ISTART.GE.1) GOTO 1180
C: variables from the "hard-core" c/blks (for a new run):
      NRUN=2163
      NTPRI=10000
      NJPRI=6
      NALPRI=1
      NTSAVE=2000
C: define geometry, type of the left boundary, and target structure:******几何结构、靶的结构等设定
      IGEO=2                *******=0平板流，=1柱形流，=2球形流
      IFLBND=0              *******右边界，r=R(N+1)：压强PBR(t)，外部辐射温度TREX(t)，磁场HZBR(t)2、左边界，r=R(1)：
                          *******IFLBND=0 固定中心对称，R(1)=U(1)=0,扩散flux为零
                          *******IFLBND=1 封闭空腔，0 < r < R(1), 边界压强PBR(t)扩散flux为零
                          *******IFLBND=-1 开放半空间， r < R(1),只用于平板结构，压强PBR(t)，外部辐射温度TREX(t)，磁场HZBR(t)，人工粘性关闭（增加SMU2进行补偿）
      NZ=5                  ********靶层数
      NZEVEN=1              ********网格结构开始的层
      QZEVEN=1.D0             *********网格结构开始层的cell质量几何级数因子0.93 < QZEVEN < 1.07

C:: set parameters of target layers:  ******靶层的参数设定
      RZ0(1)=0.D0

      I=1
      NMESH(I)=10
      NSUB(I)=1            *******第I层物质的量
      RZ0(I+1)=4.280D0    *****第I+1层左边界的位置
      ROZ0(I)=3.2D-4     ******靶层I的左边初始密度
      IFMFU0(I)=1         *******分配燃料，等于0该层没有燃料，等于1该层有燃料
      XD0(I)=1.D0          *******第I层初始D燃料的比例
      XT0(I)=1.D0           *******第I层初始T燃料的比例

      I=2
      NMESH(I)=65
      NSUB(I)=1
      RZ0(I+1)=4.500D0
      ROZ0(I)=.225D0
      IFMFU0(I)=1
      XD0(I)=1.D0
      XT0(I)=1.D0

      I=3
      NMESH(I)=32
      NSUB(I)=3
      RZ0(I+1)=4.512D0
      ROZ0(I)=19.5D0

      I=4
      NMESH(I)=24
      NSUB(I)=2
      RZ0(I+1)=4.822D0
      ROZ0(I)=1.90D0

      I=5
      NMESH(I)=25
      NSUB(I)=3
      RZ0(I+1)=4.912D0
      ROZ0(I)=19.5D0

C: define whether the target is magnetized (targets with magnetic
C  field along z-axis are physically meaningful for IGEO=0,1 only!):
      IFHZ=0            ********等于0不加磁场，等于1加磁场
      HZ0=0.D0          ********初始磁场
C: define how to treat fusion products (for details see the OPTIONS):
      IFAL=1            ********3.5-MeV alpha输运选择标志=0忽略；=1扩散；=2(any other)local
      IFP3=1            *********3-MeV proton输运选择标志=0忽略；=1扩散；=2(any other)local
      IFP14=1            ********14-MeV proton输运选择标志=0忽略；=1扩散；=2(any other)local
      IFN14=4           ******对于14-MeV和2.45-MeV中子加热的选择标志=0忽略；=1第一次散射；=2local；=3均匀                         *******传播=4存储；>4不允许（中心）
      IFN2=3
C: set beam-drive parameters, initial time step DT0:
      EIONB=10.D0
      WDRIV=720.D0        ******驱动源的功率power of the driving beam [in units TW*cm**(IGEO-2)]
      TDRFIN=.85D0         ********停止靶驱动的时间time drive finish
      TDRCAL=1.D0/FLOOR     ******time to call DRIVE; incremented in DRIVE
      DT0=1.D-8              ******初始时间步长
C initial target state is defined.
      GOTO 1190

 1180 CONTINUE
C-------------------------------------------------------- begin READSAV
C: take initial data from file 'DEISA3':
      INQUIRE(FILE='DEISA3',EXIST=LEXIST)
      IF(.NOT.LEXIST) GOTO 2010
      OPEN(13,FILE='DEISA3',STATUS='OLD',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
      CALL RECORD('R',13)
      REWIND 13
      PRINT 9030,NRUN,NTSLOI
C========================================================== end READSAV
 1190 CONTINUE
C......................................................................
C Here non-standard values of the constants assigned in DEFVAL can be .
C set up; example:
C      TMIN=3.D-5
C......................................................................

      CALL START(ISTART)

C......................................................................
C Here non-standard initial values of variables can be assigned;      .
C values assigned in START can be overwritten; example:               .
C      IFBURN=-1, DT=other than 1.D-12, ...                           .
C Also, certain values taken from 'DEISA3' can be overwritten!        .
C......................................................................

C: set initial values of the print, terminate, and save flags:
      IPRIN=0
      JPRIN=0
      IFIN=0
      ISAVE=0
C=========================================================== end INIDAT

C: print the initial target state:
      WRITE(0,9004) NRUN
      IF(ISTART.LE.0) PRINT 9005,NRUN
      CALL PRINST
      CALL PRINOO
      CALL PRINIJ

C------------------------------------------------------- begin MAINLOOP
C>>! Never change the order of the CALL st-ts in the HARDCORE block !<<
C----------------------------------------------------------------------
 1200 CONTINUE
C------------------------------------------------------- begin HARDCORE
      CALL URSOS
      TVACCD=TVACCD+VINCD+TEINCD
      IF(TVACCD.GE.CZDRIV) CALL DRIVE
      CALL KINBUR
c:  the new time step DT is determined between here and label 1340:
      CALL STEP
      DTPHYS=DT

C: ensure CALL DRIVE when TIME=TDRFIN:
 1240 IF(TIME+DT-TDRFIN) 1280,1270,1250
 1250 IF(TIME-TDRFIN) 1260,1270,1280
 1260 DT=TDRFIN-TIME
      IF(DT.LT.1.D-3*DTPHYS) DT=1.D-3*DTPHYS
 1270 TVACCD=1.D0+1.1D0*DABS(CZDRIV)
C: ensure CALL DRIVE when TIME=TDRCAL:
 1280 IF(TIME+DT-TDRCAL) 1320,1310,1290
 1290 IF(TIME-TDRCAL) 1300,1310,1320
 1300 DT=TDRCAL-TIME
      IF(DT.LT.1.D-3*DTPHYS) DT=1.D-3*DTPHYS
 1310 TVACCD=1.D0+1.1D0*DABS(CZDRIV)
C: ensure printout when TIME=TPRIN:
 1320 IF(TIME.GE.TPRIN) GOTO 1330
      IF(TIME+DT.LT.TPRIN) GOTO 1340
      DT=TPRIN-TIME
      IF(DT.LT.1.D-3*DTPHYS) DT=1.D-3*DTPHYS
      IPRIN=1
C 1330 TPRIN=TPRIN+DTPRIN
 1330 JPRIN=JPRIN+1
      TPRIN=TIMPRI(JPRIN+1)
 1340 CONTINUE

      CALL EBALNC
      CALL UPSLOI
      NTSLOI=NTSLOI+1
      TIME=TIME+DT
C: hard-core variables are recalculated!
C========================================================= end HARDCORE

C: memorize data for history plots:
C      IF(TIME.LT.THIPLO) GOTO 1390
CC     IF(NDTPLO*(NTSLOI/NDTPLO).NE.NTSLOI) GOTO 1390
C      IHIS=IHIS+1
C      IF(IHIS.GT.NNHIS) GOTO 1390
C      TIMHIS(IHIS)=TIME*10.d0
C      DO J=1,JRHIS
C       JR=1+J*10
C       IF(JR.LE.N2) THEN
C         RHIS(J,IHIS)=R(JR)*.1D0
C       ELSE
C         RHIS(J,IHIS)=R(N2)*.1D0
C       ENDIF
C      ENDDO
C      AHIS(1,IHIS)=-U(90)
C      AHIS(2,IHIS)=-U(70)
C      THIPLO=THIPLO+DTHIPL
C 1390 CONTINUE

C---------------------------------------------------------- begin TCOND
C  to terminate the run: GOTO 1400; to continue: GOTO 1410
      IF(NTSLOI.GE.NTFIN) GOTO 1400
      IF(TIME.GE.TFIN) GOTO 1400
      IF(TI(1).GT.1.D0) GOTO 1410
      IF(U(NFU).LT.1.D-3) GOTO 1410
      IF (RORFU.GT.5.D0) GOTO 1410

 1400 IFIN=1
      ISAVE=1
      IPRIN=1
 1410 CONTINUE
C============================================================ end TCOND

C--------------------------------------------------------- begin PSCOND
C: basic save and print conditions:
      IF(NTSAVE*(NTSLOI/NTSAVE).EQ.NTSLOI) ISAVE=1
      IF(NTPRI*(NTSLOI/NTPRI).EQ.NTSLOI) IPRIN=1
C.........................................................
C Here additinal printout conditions can be introduced   .
C.........................................................

C: condition of extremum of the fuel <RO*R>:
C RORFU=<RO*R> of the central fuel sphere, calculated in 'KINBUR';
      IF(RORFU.GE.1.01D0*RORP) GOTO 1540
      IF(RORFU.GT..99D0*RORP) GOTO 1580
      RORP=RORFU
      IF(JRORP.EQ.-1) GOTO 1580
C: local maximum of RORFU:
      JRORP=-1
      IPRIN=1
      GOTO 1580
 1540 RORP=RORFU
      IF(JRORP.EQ.1) GOTO 1580
C: local minimum of RORFU:
      JRORP=1
C     IPRIN=1
 1580 CONTINUE
C=========================================================== end PSCOND
      IF(IPRIN.EQ.0) GOTO 1700
C------------------------------------------------------- begin PRINTOUT
 1600 CALL PRINOO
      CALL PRINIJ
      IPRIN=0

C: write profile plots:
C      IPROF=IPROF+1
C      IF(IPROF.GT.NPROF) GOTO 1650
C      CALL URSOS
C      V(N1)=1.D18
C      DO 1630 JJ=1,200
C      J=1+2*(JJ-1)
C      IF(J.GT.N) J=N
C      APROF(1,JJ)=.05D0*(R(J)+R(J+1))
C      APROF(2,JJ)=.5D0*(U(J)+U(J+1))
C      APROF(3,JJ)=1.D0/V(J)
C      APROF(4,JJ)=TI(J)
C      APROF(5,JJ)=TE(J)
C      APROF(6,JJ)=TR(J)
C      APROF(7,JJ)=PE(J)+PI(J)+ASBOL*TR(J)**4/3.D0
C 1630 CONTINUE
C      WRITE (21,9129) NRUN,TIME
C      WRITE (21,9120) 
C      WRITE (21,9130) APROF
C      PRINT 9110, IPROF
C 1650 CONTINUE

C: write history plots:
C      WRITE(22,9131) (TIMHIS(I),(AHIS(J,I),J=1,10),I=1,IHIS)
C      REWIND 22
C      WRITE(23,9121) (TIMHIS(I),(RHIS(J,I),J=1,JRHIS),I=1,IHIS)
C      REWIND 23
C========================================================= end PRINTOUT
 1700 CONTINUE
C----------------------------------------------------------- begin SAVE
      IF(ISAVE.EQ.0) GOTO 1790
      INQUIRE(FILE='DEISA3',EXIST=LEXIST,OPENED=LOPEN)
      IF(LOPEN) GOTO 1720
      IF(LEXIST) GOTO 1710
      OPEN(13,FILE='DEISA3',STATUS='NEW',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
 1710 OPEN(13,FILE='DEISA3',STATUS='OLD',ACCESS='SEQUENTIAL'
     *,FORM='FORMATTED')
 1720 CALL RECORD('W',13)
      REWIND 13
      NSAVE=NSAVE+1
      PRINT 9045,NSAVE,NTSLOI,TIME
      ISAVE=0
 1790 CONTINUE
C============================================================= end SAVE

 1800 IF(IFIN.EQ.0) GOTO 1200
C========================================================= end MAINLOOP
      CALL PRINST
      CALL PRINOO
      PRINT 9010,NRUN
      WRITE(0,9085) NRUN
C============================================================== end RUN
      STOP
 2010 PRINT 9070
      STOP
 9004 FORMAT('* * *     DEIRA run ',I8,'     * * *'/)
 9005 FORMAT(1X,52('*'),' DEIRA4: begin run',I8)
 9010 FORMAT(1X,51('*'),' DEIRA4: end of run',I8/)
 9020 FORMAT(8I10)
 9030 FORMAT(1X,27('*'),'DEIRA4: Continue run=',I8
     *,'  from time step',I6)
 9045 FORMAT(' ________SAVE record',I3,', time step',I8
     *,', TIME=',1PG15.8)
 9070 FORMAT(' Could not open file DEISA3 for ISTART.GE.1')
 9085 FORMAT(/'* * *     End of run NRUN=',I8,'     * * *'
     &//' Press any key to continue')
 9110 FORMAT(' ________profile PLOT record',I5)
 9121 FORMAT(1P41E10.3)
 9120 FORMAT('  R,j+1/2,cm U,j+1/2    rho',8X,'T_i',8X,'T_e',8X,'T_r'
     &,8X,'P_e+P_i+P_r')
 9129 FORMAT(/'RUN=',I6,'  TIME = ',1PG12.4)
 9130 FORMAT (1P10E10.3)
 9131 FORMAT (1P11E10.3)
 9998 FORMAT(6I10)
 9999 FORMAT(1P5E15.8)
      END

C**********************************************************************
      SUBROUTINE START(ISTART)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine constructs the mesh and assigns initial values of all
c     the scalar and field variables.

c     Called by:  DEIRA
c     Calls    :  DM0PR,ZERA0
c =====================================================================
      INCLUDE 'de4com.fi'
      EXTERNAL FZER0PR
      COMMON/RO0FZER/RRAB0,HHRAB,IIRAB,JJRAB
      COMMON/URSAZ/ASUB(5),ZSUB(5)

      DIMENSION EE(1:NN+1),EI(1:NN+1),EZALL(1:20*NNZ+22)
      EQUIVALENCE (QE(1),EE(1)),(QI(1),EI(1)),(EZDR(1),EZALL(1))

      IF(ISTART.GT.0) GOTO 2000

c: check for admissable values of IGEO:
      IF(IGEO.EQ.0.OR.IGEO.EQ.1.OR.IGEO.EQ.2) GOTO 120
      PRINT 9010,IGEO
      STOP 'Inadmissible IGEO'
c: check for admissable combination of IFLBND,IGEO and RZ0(1) values:
 120  IF(IFLBND.EQ.0.AND.ABS(RZ0(1)).LT.FLOOR) GOTO 130
      IF(IFLBND.EQ.1.AND.RZ0(1).GT.FLOOR) GOTO 130
      IF(IFLBND.EQ.-1.AND.IGEO.EQ.0) GOTO 130
      PRINT 9020,IGEO,IFLBND,RZ0(1)
      STOP 'Inadmissible IFLBND'
 130  CONTINUE
C: initial values of scalar variables:
      N1MAX=NN+1
      NZMAX=NNZ
      IFBURN=0
      NTSLOI=0
      NTBAD=0
      NTVBAD=0
      NDRIV=0
      NSAVE=0
      NTFLLE=0
      NTFLLI=0
      NTFLLR=0
      NFU=0
      NJZ(1)=1
      JRORP=1
      TIME=0.D0
      RORFUM=0.D0
      VINCR=0.D0
      TEINCR=0.D0
      TIINCR=0.D0
      TRINCR=0.D0
      HZINCR=0.D0
      VINCD=0.D0
      TEINCD=0.D0
      AMSFU=0.D0
      JF=20*NZMAX+22
      DO 200 J=1,JF
  200 EZALL(J)=0.D0
      DO I=1,NZ
        IF(ROZ0(I).LE.0.D0) THEN
          PRINT 9040,NZ
          PRINT 9050,(ROZ0(J),J=1,NZ)
          STOP
        ENDIF
        IF(ROZ01(I).LE.0.D0) ROZ01(I)=ROZ0(I)
      ENDDO

C: geometry:
      I1GEO=IGEO+1
      S1GEO=I1GEO
      CSURF=1.D0
      IF(IGEO.GE.1) CSURF=2*IGEO*PINUM

C: target structure:
      DO 400 I=1,NZ          *******NZ为靶层总数
      IF(IFMFU0(I).EQ.1) GOTO 360
      XD0(I)=0.D0
      XT0(I)=0.D0
      XHE0(I)=0.D0
      XH0(I)=0.D0
      XB0(I)=0.D0
c: "reduced" mass (=\int\rho r^IGEO dr) of layer I:
 360  RAB0=RZ0(I)
      RAB1=RZ0(I+1)
      AMZ(I)=DM0PR(RAB0,RAB1,I,96)
cc      AMZ(I)=ROZ0(I)*(RZ0(I+1)**I1GEO-
cc     *RZ0(I)**I1GEO)/S1GEO
      J=NMESH(I)+NJZ(I)
C: parameters of the central fuel sphere:
      IF(IFMFU0(I).NE.1.OR.NJZ(I).NE.NFU+1) GOTO 400
      NFU=J-1
      AMSFU=AMSFU+AMZ(I)
  400 NJZ(I+1)=J

      N=J-1
      N1=J
      N2=J+1
      IF(N1.LE.N1MAX) GOTO 450
      IF(NZ.LE.NZMAX) GOTO 450
      PRINT 9130,N,NZ
      STOP
C: make neutron heating "fool-proof":
 450  IF(IFN2.LE.0) IFN2=0
      IF(IFN14.LE.0) IFN14=0
      IF(IFLBND.EQ.-1) GOTO 454
      IF(AMSFU.GT.FLOOR) GOTO 460
 454  IF(IFN2.GT.0) IFN2=2
      IF(IFN14.GT.0) IFN14=2

C: construct the mesh in layer I=NZEVEN:
  460 JI=NJZ(NZEVEN)
      JF=NJZ(NZEVEN+1)-1
      R(JI)=RZ0(NZEVEN)
C:: set the denominator of the geometric progression for cell masses:
      QZ=QZEVEN
      IF(QZ.EQ.1.D0) H=AMZ(NZEVEN)/(JF-JI+1)
      IF(QZ.NE.1.D0) H=AMZ(NZEVEN)*(QZ-1.D0)/(QZ**(JF-JI+1)-1.D0)
      IIRAB=NZEVEN
      JJRAB=1
      DO 500 J=JI,JF
      RAB0=R(J)
      RAB1=RZ0(NZEVEN+1)
      RAB=0.5D0*(RAB0+RAB1)
      RRAB0=RAB0
      HHRAB=H
      R(J+1)=ZERA0(FZER0PR,RAB0,RAB1,RAB)
CC      DR3=H*S1GEO/ROZ0(NZEVEN)
CC      IF(IGEO.NE.0) R(J+1)=(R(J)**I1GEO+DR3)**(1.D0/S1GEO)
CC      IF(IGEO.EQ.0) R(J+1)=R(J)+DR3
      DM(J)=H
 500  H=H*QZ

C: construct the mesh in remaining layers:
      DO 800 IDIR=1,2
C  IDIR=1 - construct outward, IDIR=2 - inward;
      ID=3-2*IDIR
      JJRAB=ID
      J=NJZ(NZEVEN+1)-1
      IF(IDIR.EQ.2) J=NJZ(NZEVEN)
      H=DM(J)

      DO 800 IP=1,NZ
      I=NZEVEN+IP*ID
      IF(I.LE.0.OR.I.GT.NZ) GOTO 800
      IIRAB=I
      JI=NJZ(I)
      JF=NJZ(I+1)-1
      QZ=ZN(JF-JI+1,AMZ(I),H)

      DO 750 K=JI,JF
      DR3=H*S1GEO/ROZ0(I)
      J=K
      IF(IDIR.EQ.2) J=JI+JF-K
      DM(J)=H
      J=J+IDIR-1
      J1=J+ID
      IF(J1.EQ.1) THEN
        R(J1)=RZ0(1)
      ELSE
CC        IF(IGEO.EQ.0) R(J1)=R(J)+ID*DR3
CC        IF(IGEO.NE.0) R(J1)=(R(J)**I1GEO+ID*DR3)**(1.D0/S1GEO)
        IF(IDIR.EQ.1) THEN
          RAB0=R(J)
          RAB1=RZ0(I+1)
        ELSE
          RAB0=RZ0(I)
          RAB1=R(J)
        ENDIF
        RAB=0.5D0*(RAB0+RAB1)
        RRAB0=R(J)
        HHRAB=H
        R(J1)=ZERA0(FZER0PR,RAB0,RAB1,RAB)
      ENDIF
 750  H=H*QZ
      H=H/QZ
 800  CONTINUE
      R(N2)=R(N1)

C: initial values of field variables:
      DO 1000 I=1,NZ
      IF(ZSUB(NSUB(I)).GT..5D0.AND.ASUB(NSUB(I)).GT..5D0) GOTO 802
      PRINT 9030,NSUB(I),I,ASUB(NSUB(I)),ZSUB(NSUB(I))
      STOP 'STOP in START: zero ZSUB or ASUB'  
 802  IF(IFMFU0(I).EQ.1) GOTO 820
C: non-fuel layer:
      XMOL(I)=1.D0
      AMOL(I)=ASUB(NSUB(I))
      ZMOL(I)=ZSUB(NSUB(I))
      Z2MOL(I)=ZSUB(NSUB(I))**2
      SMOL(I)=1.D0/(SQRT(ASUB(NSUB(I)))*ZSUB(NSUB(I))**2)
      GOTO 850
C: fuel layer:
 820  XMOL(I)=XD0(I)+XT0(I)+XHE0(I)+XH0(I)+XB0(I)
      AMOL(I)=2.014102D0*XD0(I)+3.016049D0*XT0(I)+3.016029D0*XHE0(I)+
     &1.007825D0*XH0(I)+11.009305D0*XB0(I)
      ZMOL(I)=XD0(I)+XT0(I)+XH0(I)+2.D0*XHE0(I)+5.D0*XB0(I)
      Z2MOL(I)=XD0(I)+XT0(I)+XH0(I)+4.D0*XHE0(I)+25.D0*XB0(I)
      SMOL(I)=XD0(I)/(SQRT(2.014102D0))+XT0(I)/(SQRT(3.016049D0))+
     &XH0(I)/SQRT(1.007825D0)+XHE0(I)/(4.D0*SQRT(3.016029D0))+
     &XB0(I)/(25.D0*SQRT(11.009305D0))
c: convert "reduced" layer masses into real masses:
 850  AMZ(I)=CSURF*AMZ(I)
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 1000 J=JI,JF
      U(J+1)=0.D0
      RS1=R(J)
      IF(IGEO.GE.1) RS1=R(J)**I1GEO
      RS1P=R(J+1)
      IF(IGEO.GE.1) RS1P=R(J+1)**I1GEO
      V(J)=(RS1P-RS1)/(S1GEO*DM(J))
      TE(J)=TE0(I)
      IF(TE(J).LT.TMIN) TE(J)=TMIN
      TI(J)=TI0(I)
      IF(TI(J).LT.TMIN) TI(J)=TMIN
      TR(J)=TR0(I)
      IF(TR(J).LT.TMIN) TR(J)=TMIN
      HZ(J)=HZ0
      EAL(J)=0.D0
      EP3(J)=0.D0
      EP14(J)=0.D0
      IFDMFU(J)=IFMFU0(I)
      XD(J)=XD0(I)
      XT(J)=XT0(I)
      XHE(J)=XHE0(I)
      XH(J)=XH0(I)
      XB(J)=XB0(I)
      ISHLJ(J)=I
 1000 CONTINUE
      ISHLJ(N1)=ISHLJ(N)
      U(1)=0.D0

      DT=DT0
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(N1)**IGEO
      RAB=(.25D0*UCBET*DT)*(RSP/(DM(N)*V(N)))
      RAB=TREX*DSQRT(DSQRT(RAB/(CZDTVT+RAB)))
      IF(TR(N).LT.RAB) TR(N)=RAB
      IF(IFLBND.NE.-1) GOTO 1020
      RAB=.25D0*UCBET*DT/(DM(N)*V(N))
      RAB=TRLEX*DSQRT(DSQRT(RAB/(CZDTVT+RAB)))
      IF(TR(1).LT.RAB) TR(1)=RAB
 1020 CONTINUE

C: initial energy in layers - for balance; only when  ISTART.LE.0:
      CALL URSOS
      DM(N1)=DM(N)
      DO 1200 I=1,NZ
      EZ0(I)=0.D0
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 1200 J=JI,JF
      JM1=J-1
      IF(JM1.LT.1) JM1=1
1200  EZ0(I)=EZ0(I)+CSURF*(DM(J)*(EE(J)+EI(J)+
     +(ASBOL*TR(J)**4+HZ(J)**2/(8.D0*PINUM)+
     +EAL(J)+EP3(J)+EP14(J))*V(J))+
     +.125D0*((DM(JM1)+DM(J))*U(J)**2+
     +(DM(J)+DM(J+1))*U(J+1)**2))

      CALL DRIVE
      DO 1230 J=1,N
      HZOLD(J)=0.D0
      TEOLD(J)=0.D0
1230  VOLD(J)=0.D0
      CALL BNDVAL(TIME)
      PBLOLD=PBL
      PBROLD=PBR
      HBLOLD=HZBL
      HBROLD=HZBR
      PBLSOL=PBLSUM
      PBRSOL=PBRSUM
      CALL KINBUR
      RORP=RORFU
      GOTO 3000
C======================================================================

C: initialization of the run continued from 'DEISA3':
 2000 CALL URSOS
      CALL DRIVE
      DO 2030 J=1,N
      HZOLD(J)=0.D0
      TEOLD(J)=0.D0
 2030 VOLD(J)=0.D0
      CALL KINBUR
C======================================================================

 3000 RETURN

 9010 FORMAT(' Inadmissable value of IGEO =',I8)
 9020 FORMAT(' Inadmissable combination of IGEO=',I8
     *,' IFLBND=',I8,' and RZ0(1)=',1PE16.8)
 9030 FORMAT(/' Substance #',I2,' in shell #',I2,' has A,Z=',1P2E12.4)
 9040 FORMAT(' Non-positive nitial density: NZ=',I4,'  ROZ0(I)=')
 9050 FORMAT(1P5E14.6)
 9110 FORMAT(8I10)
 9120 FORMAT(2I10,5G10.4)
 9130 FORMAT(' N=',I4,' NZ=',I3,' BEYOND THE LIMIT')
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION RO0PR(RRAB,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine defines the initial density profile as a function of
c     radius RRAB in target layer I.

c     Called by:  DM0PR
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'

c: initial density profile in layer I:
      RO0PR=ROZ0(I)+(ROZ01(I)-ROZ0(I))*((RRAB-RZ0(I))/
     &(RZ0(I+1)-RZ0(I)))**BRO0PR(I)

      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION DM0PR(R0,R1,I,NGAUSS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the "reduced" mass in the interval R0<r<R1
c     according to the initial density profile in target layer I.

c     Called by:  START,FZER0PR
c     Calls    :  GQUAD
c =====================================================================
      EXTERNAL FDM0PR
      COMMON/RO0FDM/IIRABDM,JJRABDM

      IIRABDM=I
      DM0PR=GQUAD(FDM0PR,R0,R1,NGAUSS)
      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION FDM0PR(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine defines the function to be integrated when
c     the "reduced" mass interval DM is calculated in DM0PR.

c     Called by:  GQUAD
c     Calls    :  RO0PR
c =====================================================================
      COMMON/RO0FDM/IIRABDM,JJRABDM
      COMMON/PAGEO/IGEO,I1GEO,S1GEO,CSURF

      I=IIRABDM
      RAB=1.D0
      IF(IGEO.GE.1) RAB=RR**IGEO
      FDM0PR=RAB*RO0PR(RR,I)
      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION FZER0PR(RR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine defines the function to be zeroed when costructing
c     the m-mesh in target layer I.

c     Called by:  ZERA0
c     Calls    :  DM0PR
c =====================================================================
      COMMON/RO0FZER/RRAB0,HHRAB,IIRAB,JJRAB

      I=IIRAB
      RRRAB0=RRAB0
      IF(JJRAB.EQ.1) THEN
        FZER0PR=DM0PR(RRRAB0,RR,I,16)/HHRAB-1.D0
      ELSE
        FZER0PR=DM0PR(RR,RRRAB0,I,16)/HHRAB-1.D0
      ENDIF
      RETURN
      END

C**********************************************************************
      SUBROUTINE PRINUN
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine performs all the regular printouts.

c     Called by:  DEIRA
c     Calls    :  URSOS,KINBUR
c =====================================================================
      REAL FARR
      INCLUDE 'de4com.fi'

C: c/blks for the equation of state, filled from file 'EOSTA3':
      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/EPZER/P00(5),E00(5)
      COMMON/URSEDG/NROO(6),ROILL(5),HLROO(5)
     *,NTEMM(3,6),TEMILL(3,5),HLTEMM(3,5)
      COMMON/URSTBL/FARR(30000)
      COMMON/URTBLI/QURS(1000)

C: c/blks for the beam, filled from file 'BEMTA3':
      COMMON/BECNST/CBEAM(8)
      COMMON/BEAM1/ABEAM,ZBEAM,CBEAM1(406)
      COMMON/BEAM2/CBEAM2(350),NCBEM2(106),CBEM22(125)

      DIMENSION YI(1:NN+1),P(1:NN+1),EE(1:NN+1),EI(1:NN+1)
     *,EALL(1:3*NN+3),HIALL(1:3*NN+3)
     *,EMT(31),EMZ(1:NNZ,14),EZALL(1:19*NNZ+19),EMTL(8)
      CHARACTER*6 LGEO(3),LHZ(2),LNEU(5),LAL(3),LAL1(3),LAL2(3),LAL3(3)
      CHARACTER*17 LLBND(3)
      EQUIVALENCE (US(1),YI(1)),(PE(1),P(1))
     *,(QE(1),EE(1)),(QI(1),EI(1))
     *,(EAL(1),EALL(1)),(EP3(1),EALL(NN+2)),(EP14(1),EALL(2*NN+3))
     *,(HIAL(1),HIALL(1)),(HIP3(1),HIALL(NN+2)),(HIP14(1),HIALL(2*NN+3))
     *,(EZDR(1),EZALL(1))
      DATA LGEO/'PLANAR','CYLIND','SPHERE'/
     *,LHZ/'NO    ','YES   '/
     *,LNEU/'IGNORE','1 SCAT','LOCAL ','UNIFOR','D_PROF'/
     *,LAL/'IGNORE','DIFFUS','LOCAL '/
     *,LAL1/' E ALF','  E P3',' E P14'/
     *,LAL2/'  HIAL','  HIP3',' HIP14'/
     *,LAL3/'  HIER','  HIEI','  HZ  '/
     *,LLBND/'"fixed center",  ','"void cavity",   '
     *,'"open halfspace",'/

      PRINT 9210
      PRINT 9212
      PRINT 9213
      RETURN

C**********************************************************************
      ENTRY PRINUS
      PRINT 9219
      JF=0
      K=0

      DO 100 I=1,5
      IF(ZSUB(I).LT..5D0) GOTO 100
      PRINT 9220,I,ZSUB(I),ASUB(I)
      RAB=11.206D0*ASUB(I)*DEXP(ROILL(I))
      RAB1=RAB*DEXP(HLROO(I)*(NROO(I)-1))
      JI=JF+1
      JF=JF+3*NROO(I)*NTEMM(3,I)
      J=K+1
      K=K+4*NROO(I)
      PRINT 9230,RAB,RAB1,JI,JF
      RAB=.02721D0*DEXP(TEMILL(1,I))
      RAB1=RAB*DEXP(HLTEMM(1,I)*(NTEMM(1,I)-1)+
     +HLTEMM(2,I)*(NTEMM(2,I)-NTEMM(1,I))+
     *HLTEMM(3,I)*(NTEMM(3,I)-NTEMM(2,I)))
      PRINT 9240,RAB,RAB1,J,K
      RAB=2.942D0*P00(I)
      RAB1=2.942D0*E00(I)/(11.206D0*ASUB(I))
      PRINT 9250,RAB,RAB1
 100  CONTINUE
      RETURN

C**********************************************************************
      ENTRY PRINBE
      PRINT 9260,ZBEAM,ABEAM
      RETURN

C**********************************************************************
      ENTRY PRINST
C: print basic run parameters:
      J=1
      IF(IFLBND.EQ.1) J=2
      IF(IFLBND.EQ.-1) J=3
      PRINT 9270,NRUN,LGEO(I1GEO),LLBND(J),RZ0(1),NTPRI,NJPRI,NTSAVE
      PRINT 9280,LNEU(IFN14+1),LNEU(IFN2+1)
     *,LAL(IFAL+1),LAL(IFP3+1),LAL(IFP14+1)
      PRINT 9300,EIONB,WDRIV,TDRFIN
      PRINT 9305,LHZ(IFHZ+1),HZ0,IFOPAC
      PRINT 9310

      DO 200 I=1,NZ
      J=NJZ(I+1)-1
      RAB=CSURF*DM(J)
      RAB1=DM(J)/DM(J-1)
      J=J+1-NJZ(I)
      PRINT 9320,I,NSUB(I),J,AMZ(I),RAB,RAB1
     *,ROZ0(I),RZ0(I+1),XD0(I),XT0(I),XHE0(I),XH0(I),XB0(I)
 200  CONTINUE
      RETURN

C**********************************************************************
      ENTRY PRINOO
C: print general information on the target:
      CALL URSOS
C: evaluate mass and energy characteristics:
C  save  YI=QSDT -> US, P -> PE
C  EMT(K) is the total mass (energy) of sort 'K';
C  EMZ(I,K) is mass (energy) of sort 'K' in target layer 'I'
      DO 250 K=1,31
 250  EMT(K)=0.D0
      EBLNC=0.D0
      TOTN2=0.D0
      TOTN14=0.D0
      TOTDHE3=0.D0
      DM(N1)=DM(N)

      DO 700 I=1,NZ
      DO 300 K=1,14
 300  EMZ(I,K)=0.D0
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 500 J=JI,JF
      YI(J)=QSDT(J)
      P(J)=PE(J)+PI(J)+ASBOL*TR(J)**4/3.D0
      RAB=CSURF*DM(J)
      RAB1=0.D0
      IF(IFDMFU(J).EQ.1) RAB1=RAB/AMOL(I)
      EMZ(I,1)=EMZ(I,1)+RAB
      EMZ(I,2)=EMZ(I,2)+2.D0*XD(J)*RAB1
      EMZ(I,3)=EMZ(I,3)+3.D0*XT(J)*RAB1
      EMZ(I,4)=EMZ(I,4)+3.D0*XHE(J)*RAB1
      EMZ(I,5)=EMZ(I,5)+(R(J+1)-R(J))/V(J)
      JM1=J-1
      IF(JM1.LT.1) JM1=1
      EMZ(I,7)=EMZ(I,7)+CSURF*.125D0*
     *((DM(JM1)+DM(J))*U(J)**2+
     +(DM(J)+DM(J+1))*U(J+1)**2)
      EMZ(I,8)=EMZ(I,8)+RAB*EE(J)
      EMZ(I,9)=EMZ(I,9)+RAB*EI(J)
      RAB=RAB*V(J)
      EMZ(I,10)=EMZ(I,10)+RAB*ASBOL*TR(J)**4
      EMZ(I,11)=EMZ(I,11)+RAB*HZ(J)**2/(8.D0*PINUM)
      EMZ(I,12)=EMZ(I,12)+RAB*EAL(J)
      EMZ(I,13)=EMZ(I,13)+RAB*EP3(J)
 500  EMZ(I,14)=EMZ(I,14)+RAB*EP14(J)
      DO 550 K=7,14
 550  EMZ(I,6)=EMZ(I,6)+EMZ(I,K)
      DO 600 K=1,14
 600  EMT(K)=EMT(K)+EMZ(I,K)
      DO 650 K=15,23
      JF=I+(NZMAX+1)*(K-15)
 650  EMT(K)=EMT(K)+EZALL(JF)
      TOTN2=TOTN2+TNUN2(I)
      TOTN14=TOTN14+TNUN14(I)
      TOTDHE3=TOTDHE3+TNUDHE3(I)
 700  EBLNC=EBLNC+EZ0(I)

c: energies that flowed out through the right boundary:
      DO 750 K=24,31
      JF=NZ+1+(NZMAX+1)*(K-15)
 750  EMT(K)=EZALL(JF)
c: energies that flowed in through the left boundary:
      DO 760 K=1,8
      JF=1+(NZMAX+1)*(K+8)
 760  EMTL(K)=EZALL(JF)

      CALL KINBUR
      PRINT 9400,NTSLOI,TIME,DT,NDRIV,NTBAD
     *,NTVBAD,NTFLLE,NTFLLI,NTFLLR
C: thermonuclear energy gain:
      RAB=EMT(15)+DABS(EMT(15))+EBLNC+DABS(EBLNC)
     *+ERIN+DABS(ERIN)
      IF(RAB.GT.0.D0) RAB=2.D0*EMT(17)/RAB
C: energy balance:
      EBLNC=EMT(6)-EBLNC-EMT(15)-
     *EMT(18)-EMT(19)-EMT(20)-EMT(21)-EMT(22)-EMT(23)+
     +EMT(24)+EMT(25)+EMT(26)+EMT(27)+EMT(28)+EMT(29)+EMT(30)+EMT(31)-
     *EMTL(1)-EMTL(2)-EMTL(3)-EMTL(4)-EMTL(5)-EMTL(6)-EMTL(7)-EMTL(8)
      PRINT 9410,EMT
      PRINT 9420,PBLOLD,PBROLD,HBLOLD,HBROLD,TRLEX,TREX,EREX,ERIN
     *,TR(N),RORFU,RORFUM,TOTN2,TOTN14,EBLNC,RAB,TOTDHE3
      IF(IFN14.EQ.4) PRINT 9404,JTAU140,TAU140
      WRITE(0,9402) NTSLOI,TIME,DT,RAB
      RETURN

C**********************************************************************
      ENTRY PRINIJ
      PRINT 9412, EMTL
C: print information on each layer:
      DO 1000 I=1,NZ
      RAB=TIME*1.D1
      AEOS=ASUB(NSUB(I))
      IF(IFMFU0(I).NE.1) GOTO 810
      IF(XHE0(I).GT.XBFLOO) GOTO 810
      IF(XB0(I).GT.XBFLOO) GOTO 810
C: adjust the value of A for a mixture of hydrogen isotopes:
      AEOS=AMOL(I)/XMOL(I)
 810  CONTINUE
      PRINT 9430,I,ZSUB(NSUB(I)),AEOS,RAB
      EBLNC=EMZ(I,6)-EZ0(I)-EZDR(I)-EZCL(I)-EZN14(I)-EZN2(I)-
     *EZAL(I)-EZP3(I)-EZP14(I)+
     +EOUTW(I+1)+EOUTE(I+1)+EOUTI(I+1)+EOUTR(I+1)+EOUTHZ(I+1)+
     +EOUTAL(I+1)+EOUTP3(I+1)+EOUP14(I+1)-
     *EOUTW(I)-EOUTE(I)-EOUTI(I)-EOUTR(I)-EOUTHZ(I)-
     *EOUTAL(I)-EOUTP3(I)-EOUP14(I)

      PRINT 9410,(EMZ(I,K),K=1,14),EZDR(I),EZJL(I),EZFUS(I)
     *,EZCL(I),EZN14(I),EZN2(I),EZAL(I),EZP3(I)
     *,EZP14(I),EOUTW(I+1),EOUTE(I+1),EOUTI(I+1),EOUTR(I+1)
     *,EOUTHZ(I+1),EOUTAL(I+1),EOUTP3(I+1),EOUP14(I+1)
      PRINT 9435,TNUN2(I),TNUN14(I),TNUDHE3(I),EBLNC
      JHI1=1
      IF(IFHZ.EQ.1) JHI1=3
      PRINT 9440,LAL3(JHI1)
      IF(I.NE.1) GOTO 820
      JRAB=0
      PRINT 9451,JRAB,R(1),U(1)
 820  JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 950 J=JI,JF
      IF(NJPRI*(J/NJPRI).NE.J.AND.J.NE.JI.AND.J.NE.JF) GOTO 950
      RAB=1.D0/V(J)
      IF(IFHZ.EQ.1) THEN
        PRINT 9450,J,R(J+1),U(J+1),RAB,TE(J),TI(J),TR(J),HZ(J),P(J)
      ELSEIF(JHI1.EQ.1) THEN
        PRINT 9450,J,R(J+1),U(J+1),RAB,TE(J),TI(J),TR(J),HIER(J),P(J)
      ELSEIF(JHI1.EQ.2) THEN
        PRINT 9450,J,R(J+1),U(J+1),RAB,TE(J),TI(J),TR(J),HIEI(J),P(J)
      ENDIF
 950  CONTINUE
      IF(XB0(I).LT.XBFLOO) GOTO 960
      PRINT 9444,LAL1(NALPRI),LAL2(NALPRI)
      GOTO 962
 960  PRINT 9442,LAL1(NALPRI),LAL2(NALPRI)
 962  CONTINUE
      DO 980 J=JI,JF
      IF(NJPRI*(J/NJPRI).NE.J.AND.J.NE.JI.AND.J.NE.JF) GOTO 980
      K=J+(NALPRI-1)*N1MAX
      RAB1=QE(J)+QI(J)
      IF(XB0(I).LT.XBFLOO) GOTO 976
      PRINT 9452,J,YI(J),CAPE(J),CAPI(J),CAPR(J)
     *,EALL(K),HIALL(K),XH(J),XB(J),RAB1
      GOTO 980
 976  PRINT 9452,J,YI(J),CAPE(J),CAPI(J),CAPR(J)
     *,EALL(K),HIALL(K),XD(J),XT(J),RAB1
 980  CONTINUE
 1000 CONTINUE
C  Electron QE(J) and ion QI(J) heating rates
C  do not include the energy deposition by al,p3,p14;
C  QE=QDRIV + QCHARGEDLOCAL/E + QNEUTRON/E
C  QI=QCHARGEDLOCAL/I + QNEUTRON/I
      RETURN

 9210 FORMAT(' DEIRA UNITS: ','[time]',13X,'= 1.E-8 sec'/
     +14X,'[length]           = 0.1 cm = 1 mm'/
     +14X,'[velocity]         = 1.E+7 cm/s'/
     +14X,'[mass]             = 1 mg/mm**(2-IGEO)'/
     +14X,'[energy]           = 1.E+11 ergs/mm**(2-IGEO)'
     +' = 10 kJ/mm**(2-IGEO)'/
     +14X,'[density]          = 1 g/cm**3'/
     +14X,'[pressure]         = 1.E+14 ergs/cm**3'/
     +14X,'[temperature]      = 1 keV = 1.16E+7 K'/
     +14X,'[magnetic field]   = 1.E+7 Gauss'/
     +14X,'[power]            = 1 TW/mm**(2-IGEO)'/
     +14X,'[spec. power]      = 1.E+22 ergs/(g*sec) = 1 TW/mg'/
     +14X,'[heat cond.coeff.] = 1.E+20 ergs/(sec*cm*keV)'/
     +14X,'[elec.resistivity] = 1.E-8 sec'/
     +14X,'[<RO*R>]           = 1 mg/mm**2 = 0.1 g/cm**2')
 9212 FORMAT(' EXPLANATIONS TO THE OUTPUT:'/
     +' MASS: TOT=total target (layer) mass; D=mass of deuterium; '/
     +'       T=mass of tritium; HEL3=mass of helium-3.'/
     +' ENERGY: TOT=total target (layer) energy; KIN=kinetic energy;'/
     +'         EL=electron, IONS=ion, R=radiation, H=magnetic field'/
     +'         components of the internal energy;'/
     +'         AL=amount of energy in fast (3.5-MeV) alpha particles;'/
     +'         P3=energy in 3-MeV, P14=in 14-MeV protons;'/
     +' HEATING: DR=ext.drive energy deposited in the target (layer);'/
     +'          JL=Joule heating; TN=total amount of thermonuclear'/
     +'          energy released; CH=energy deposited locally by'/
     +'          charged fusion products; N14(N2)=energy deposited by'/
     +'          14(2)-MeV neutrons;'/
     +' DEP TO: AL(P3)(P14)=energy released in the form of 3.5-MeV'/
     +'         alphas (3-MeV protons) (14-MeV protons)')
 9213 FORMAT(
     +' LEAK: various components of energy transferred outward across'/
     +'       the outer (right-hand) boundary of the target (layer);'/
     +' L.INFLUX: the same components of energy transferred inward'/
     +'           across the left boundary of the target;'/
     +' P_bl,r[TR(N)][H_bl,r]=boundary press.[rad.temp.][mag.field];'/
     +' TR_lex,ex=temperature of the external radiation field, driving'/
     +'           the target from the left (right) boundary'/
     +' EREX=time integral of (sigma*TRLEX**4*area+sigma*TREX**4*area)'/
     +'     =energy of external radiation field that enters the target;'
     +/' ERIN=EOUTR(1)-EOUTR(NZ+1)=net radiative energy input into the'/
     +'                            target;'/
     +' RORFU(RORFU_max)=<RO*R>(maximum <RO*R>) of the central fuel'/
     +' NN2[NN14]=total number of 2.5-MeV [14-MeV] neutrons generated'/
     +' (E)BLNC=absolute energy balance in the target (layer);'/
     +' GAIN=thermonuclear energy gain of the target;'/
     +' column Q gives QE+QI, where QE=QDRIV+QCHARGEDLOCAL_e+'/
     +'               +QNEUTRON_e, QI=QCHARGEDLOCAL_i + QNEUTRON_i')
 9219 FORMAT(/' PARAMETERS OF THE EOS TABLES:')
 9220 FORMAT(' Substance',I2,': Z=',F5.1,5X,'A=',F6.2)
 9230 FORMAT(14X,'RO_min=',1PE9.2,5X,'RO_max=',1PE9.2
     *,8X,'FARR:',I6,' -',I6)
 9240 FORMAT(14X,'T_min=',1PE10.2,5X,'T_max=',1PE10.2
     *,8X,'QURS:',I6,' -',I6)
 9250 FORMAT(14X,'P00=',1PE12.2,5X,'E00=',1PE12.2)
 9260 FORMAT(' Ion beam: ZBEAM =',F5.1,5X,'ABEAM =',F6.2)
 9270 FORMAT(1X,'Run',I8,2X,A6,': left boundary = ',A17,
     *' R(1)=',1PG12.5/
     *' Print every',I7,'-th time step, every',I4,'-th cell,'
     *1X,'save every',I7,'-th step.')
 9280 FORMAT(1X,'Fusion products: n14-',A6,'  n2-',A6,
     *'  alpha-',A6,'  p3-',A6,'  p14-',A6,'.')
 9300 FORMAT(1X,'Driver: EIONB=',1PG11.4,'GeV,  WDRIV='
     *,1PG11.4,'TW,  TDRFIN=',1PG13.6)
 9305 FORMAT(1X,'Magnetic z-field: ',A6,' HZ0=',1PE12.5,'; IFOPAC=',I2)
 9310 FORMAT(/'  I SUB  NJ    MASS  right DM_j  PROGRESS ',2X
     *,'DENSITY   R(I+1)     DEUT TRIT HEL3 H    B11')
 9320 FORMAT(I3,2I4,1P2E10.3,1PE11.4,1PE10.3,1PG12.5,0P5F5.2)
 9322 FORMAT(I3,2I4,1P2E10.3,1PE11.4,1PE10.3,1PG12.5,0P2F5.2)
 9400 FORMAT(/1X,78('=')/1X,'STEP',I6,2X,'TIME='
     *,1PG15.8,2X,'DT=',1PG10.3,4X,'DRIVE CALLS =',I5/
     *3X,'bad stps=',I6,' v.bad stps=',I6,
     *2X,'fl.lim: e-',I7,' i-',I7,' r-',I7)
 9402 FORMAT('STEP',I6,2X,'TIME=',1PG15.8,2X,'DT=',1PG10.3
     &,' GAIN=',1PE10.3)
 9404 FORMAT('JTAU140=',I4,'  TAU140=',1PE10.3)
 9410 FORMAT(' MASS: TOT=',1PE11.4,' D=',1PE10.3
     *,' T=',1PE10.3,' HEL3=',1PE10.3,' RO*R=',1PE10.3
     */' ENRG: TOT=',1PE11.4,' KIN=',1PE11.4,' EL=',1PE11.4
     *,' IONS=',1PE11.4/7X,'R=',1PE11.4,' H=',1PE11.4
     *,' AL=',1PE11.4,' P3=',1PE11.4,' P14=',1PE11.4
     */' HEATING: DR=',1PE10.3,' JL=',1PE10.3,' TN=',1PE10.3
     *,' CH=',1PE10.3,' N14=',1PE10.3
     */9X,' N2=',1PE10.3,3X,'DEP TO: AL=',1PE10.3
     *,' P3=',1PE10.3,' P14=',1PE10.3
     */' LEAK: WRK=',1PE11.4,' E-CON=',1PE11.4,' I-CON=',1PE11.4
     *,' R-CON=',1PE11.4/7X,'H=',1PE11.4,' AL=',1PE11.4,' P3=',1PE11.4
     *,' P14=',1PE11.4)
 9412 FORMAT(/1X,71('-')/' L.INFLUX: W=',1PE11.4,' E-CON=',1PE11.4
     *,' I-CON=',1PE11.4,' R-CON=',1PE11.4/11X,'H=',1PE11.4
     *,' AL=',1PE11.4,' P3=',1PE11.4,' P14=',1PE11.4)
 9420 FORMAT(' BNDRS: P_bl=',1PE10.3,' P_br=',1PE10.3,' H_bl=',1PE10.3
     *,' H_br=',1PE10.3/8X,'TR_lex=',1PE10.3,' TR_ex=',1PE10.3
     *,' EREX=',1PE10.3,' ERIN=',1PE10.3/' TR(N)=',1PE10.3
     *,' RORF=',1PE9.2,' RORF_m=',1PE9.2,' NN2=',1PE10.3
     *,' NN14=',1PE10.3/' BLNC=',1PE11.4,' GAIN=',1PE10.3,' NUDHE3='
     &,1PE10.3)
 9430 FORMAT(/1X,16('*'),' Target layer',I3,' ****** (Z,A=',F5.2,','
     *,F6.2,'), time=',F7.2,' ns')
 9435 FORMAT(7X,'NN2=',1PE11.4,' NN14=',1PE11.4,' NUDHE3=',1PE11.4
     &,' EBLNC=',1PE11.4)
 9440 FORMAT(/3X,'J',4X,'R(J+1)   U(J+1)   DENSITY',3X
     *,'T ELEC   T ION    T RAD   ',A6,'   P(e+i+r)')
 9442 FORMAT(/3X,'J','  IONIZ   CAPE     CAPI     CAPR',4X
     *,A6,1X,A6,6X,'XD     XT  Q(dr+cl+n)')
 9444 FORMAT(/3X,'J','  IONIZ   CAPE     CAPI     CAPR',4X
     *,A6,1X,A6,6X,'XH     XB  Q(dr+cl+n)')
 9450 FORMAT(I4,1PE12.5,1P7E9.2)
 9451 FORMAT(I4,1PE12.5,1PE9.2)
 9452 FORMAT(I4,0PF7.3,1P5E9.2,0P2F7.3,1PE9.2)
      END

C**********************************************************************
      SUBROUTINE STEP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine evaluates the next value of the time step DT.

c     Called by:  DEIRA
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'

      DTM=DTMAX
C: constraints due to the boundary pressure:
      IF(IFLBND.EQ.0) GOTO 60
      US2=(PBLOLD-PE(1)-PI(1))*V(1)
      IF(US2.LE.0.D0) GOTO 60
      RAB=.5D0*CZDT*(R(2)-R(1))/DSQRT(US2)
      IF(RAB.LT.DTM) DTM=RAB
 60   US2=(PBROLD-PE(N)-PI(N))*V(N)
      IF(US2.LE.0.D0) GOTO 100
      RAB=.5D0*CZDT*(R(N1)-R(N))/DSQRT(US2)
      IF(RAB.LT.DTM) DTM=RAB
 100  CONTINUE

C: constraints due to the Courant cond., artif.visc., heating:
      DO 300 J=1,N
      DEL=R(J+1)-R(J)
      IF(IFLBND.NE.0) GOTO 140
      IF(J.EQ.1) DEL=R(2)/S1GEO
C: COURANT condition and linear viscosity:
 140  RAB=CZDT*DEL/(1.D0+SMU1+TMU1)
      IF(DTM*US(J).GT.RAB) DTM=RAB/US(J)

C: quadratic viscosity:
      DV=U(J)-U(J+1)
      IF(DV.LE.0.D0) GOTO 240
      RAB=.25D0*CZDT*DEL
      IF(RAB.LT.DTM*DV*(SMU2+TMU2)) DTM=RAB/(DV*(SMU2+TMU2))

C: heating:
 240  RAB=CZDTQ*(EEMIN+EET(J)*TE(J)+EIT(J)*TI(J))
      RAB1=DABS(QE(J))+DABS(QI(J))
      IF(IFBURN.EQ.2) RAB1=RAB1+EAL(J)*HIAL(J)+
     +EP3(J)*HIP3(J)+EP14(J)*HIP14(J)
      IF(IFHZ.EQ.1) RAB1=RAB1+QJL(J)
      IF(DTM*RAB1.GT.RAB) DTM=RAB/RAB1
 300  CONTINUE

      RAB1=1.D0+CZDDT
C: if too small a time step is forced by the coditions above, try
C      DTMIN=1.D-18
C      IF(TE(1).LT.3.D-3) DTMIN=.3D0*DTM
C      IF(DTMIN.GT.DT*RAB1) DTMIN=DT*RAB1

C: correct the previous time-step value:
      RAB=VINCR
      IF(RAB.LT.TEINCR) RAB=TEINCR
      IF(RAB.LT.TIINCR) RAB=TIINCR
      IF(RAB.LT.TRINCR) RAB=TRINCR
      IF(IFHZ.EQ.1.AND.RAB.LT.HZINCR) RAB=HZINCR
C: the number of 'bad' and 'very bad' steps:
      IF(RAB.GT.2.D0*RAB1*CZDTVT) NTBAD=NTBAD+1
      IF(RAB.GT.1.D1*RAB1*CZDTVT) NTVBAD=NTVBAD+1
      IF(RAB.LT.CZDTVT/RAB1) DT=DT*RAB1
      IF(RAB.GT.CZDTVT*RAB1) DT=DT*CZDTVT/RAB
      IF(DT.GT.DTM) DT=DTM
 400  IF(DT.LT.DTMIN) DT=DTMIN
      RETURN
      END

C**********************************************************************
      SUBROUTINE KINBUR
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates all the transport coefficients,
c     the thermonuclear burn and energy deposition rates.

c     Called by:  DEIRA,PRINUN,
c     Calls    :  RELCON,ELLICK
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/CAHIOU/YIONOU,CAPEOU,CAPIOU,HIEIOU,ETHZOU
     *,CAPROU,HIEROU,ETV0OU,ETV1OU

      DIMENSION YI(1:NN+1),HIALL(1:3*NN+3),CPRO(3),IFPRO(3)
      EQUIVALENCE (QSDT(1),YI(1)),(HIAL(1),HIALL(1))
     *,(HIP3(1),HIALL(NN+2)),(HIP14(1),HIALL(2*NN+3))
     *,(IFAL,IFPRO(1)),(IFP3,IFPRO(2)),(IFP14,IFPRO(3))
      DIMENSION RORAFJ(1:NN+1)
      DATA CPRO/.077D0,.12D0,.13D0/

      RORFU=0.D0
      RORAFU=0.D0
      RORAFJ(1)=0.D0
C: "zero option" for the heating rates:
      DO 100 J=1,N
      QE(J)=QDRIV(J)
 100  QI(J)=0.D0
C-------------------------------------------------------------- begin 1
C: check the beginning of thermonuclear burn:
C IFBURN < 0 - no check, no burn anticipated !
C IFBURN=0 - not burning yet;  IFBURN=2 - burning;
      IF(IFBURN.EQ.2.OR.IFBURN.LE.-1) GOTO 200
      DO 190 J=1,N
      IF(IFDMFU(J).NE.1) GOTO 190
      IF(TI(J).GT.TBURN0) IFBURN=1
 190  CONTINUE
 200  CONTINUE
C================================================================ end 1

C: main loop:
      DO 2000 I=1,NZ
      Z=ZSUB(NSUB(I))
      AMEAN=AMOL(I)/XMOL(I)
      AMOL2=AMOL(I)**2
      WN2(I)=0.D0
      WN14(I)=0.D0
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 2000 J=JI,JF
      VJ=V(J)
      TEJ=TE(J)
      TIJ=TI(J)
      TRJ=TR(J)
      HZJ=HZ(J)
      YIJ=YI(J)
      IF(YIJ.LT.1.D-9) YIJ=1.D-9
      IF(YIJ.GT.Z) YIJ=Z
      YSTAR=YIJ
      IF(YIJ.LT.1.D0) YSTAR=1.D0
      FYSTAR=1.D0+YIJ-YSTAR

      IF(J.GT.NFU) GOTO 350
      RAB=(R(J+1)-R(J))/VJ
      RORFU=RORFU+RAB
      RORAFU=RORAFU+RAB/AMEAN
      RORAFJ(J+1)=RORAFJ(J)+RAB/AMEAN
 350  CONTINUE

      YAE=YIJ*ZMOL(I)/(Z*AMOL(I))
      RYA=YAE/VJ
C: prepare for the flux limiting:
      CAPEB(J)=FLLE*1280.D0*RYA*TEJ*DSQRT(TEJ)
      CAPIB(J)=FLLI*29.97D0*TIJ*DSQRT(TIJ/AMEAN)/(VJ*AMEAN)
      CAPRB(J)=FLLR*(ASBOL*UCBET*TRJ**2)*TRJ**2

C-------------------------------------------------------------- begin 2
C: check the condition for recalculating the transport coefficients
C  in mesh cell J; if 'TRUE', goto 400:
      IF(VOLD(J).LT.FLOOR) GOTO 400
      IF(TEOLD(J).LT.FLOOR) GOTO 400
      RAB=VJ/VOLD(J)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      IF(RAB.GT.1.D0+CZVKIN) GOTO 400
      RAB=(TEJ+TRJ)/TEOLD(J)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      IF(RAB.GT.1.D0+CZTKIN) GOTO 400
      IF(IFHZ.NE.1) GOTO 390
      IF(HZOLD(J).LT.FLOOR) GOTO 400
      RAB=HZJ/HZOLD(J)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      IF(RAB.GT.1.D0+CZVKIN) GOTO 400
 390  IF(IFBURN.LE.0) GOTO 2000
      IF(IFBURN.EQ.2) GOTO 1700
 400  VOLD(J)=VJ
      HZOLD(J)=HZJ
      TEOLD(J)=TEJ+TRJ
C================================================================ end 2
C: "zero option" for temperature relaxation, heat conduction, viscosity
C  and resistivity coefficients:
      CAPE(J)=0.D0
      CAPI(J)=0.D0
      HIEI(J)=1.D0/FLOOR
      ETHZ(J)=0.D0
      CAPR(J)=0.D0
      HIER(J)=0.D0
      ETVI0(J)=0.D0
      ETVI1(J)=0.D0
C: "zero option" for energy dissipation by charged products:
      HIAL(J)=1.D0/FLOOR
      HIP3(J)=1.D0/FLOOR
      HIP14(J)=1.D0/FLOOR

C-------------------------------------------------------------- begin 3
C:  electron and radiative temperature relaxation and heat conduction
C   coefficients, electrical resistivity, ion viscosity:
      CALL RELCON(I,VJ,TEJ,TIJ,TRJ,HZJ,YIJ)
      CAPE(J)=CAPEOU
      CAPI(J)=CAPIOU
      HIEI(J)=HIEIOU
      ETHZ(J)=ETHZOU
      CAPR(J)=CAPROU
      HIER(J)=HIEROU
      ETVI0(J)=ETV0OU
      ETVI1(J)=ETV1OU
C================================================================ end 3
      IF(IFBURN.LE.0) GOTO 2000
C--------------------------------------------------------------- begin 5
C: coefficients of energy dissipation by charged fusion products
C  HIAL, HIP3, HIP14:
      EFE=.026D0*RYA**.6666666667D0
      TFE=DSQRT(TEJ**2+(.6666666667D0*EFE)**2)
      SQTFE=DSQRT(TFE)

      DO 1630 K=1,3
      IF(IFPRO(K).NE.1) GOTO 1630
      U0=UPRO(K)
      XE=U0/(187.5D0*SQTFE)
      XE3=XE*XE*XE
      BY=.353D0+XE*XE*(2.34D0+XE3)/(11.D0+XE3)
      AAL=APRO(K)
      ZAL=ZPRO(K)
      CL=7.207D0+.5D0*DLOG((AMEAN/(1.D0+AMEAN/AAL))**2*TFE*
     *BY*U0*U0/(RYA*(1.D0+(39.D0*ZAL*YSTAR/U0)**2)))
      UVK=(YSTAR*ZMOL(I)/(AMOL(I)*Z))**2*(143.D0/U0)**3*
     *DSQRT(1.D0+AMEAN/AAL)*CL*FYSTAR*ZAL**2/AAL
      RAB=139.D0*TFE*BY/DSQRT(RYA)
      CL=DLOG(1.D0+RAB/(1.D0+.5D0/DSQRT(RAB)))
      UVK=UVK+YAE*(1746.D0/U0)**3*(ZAL**2/AAL)*CL/(1.D0+1.33D0/XE3)
      JK=J+(K-1)*N1MAX
      HIALL(JK)=UVK*(2.5D0-1.5D0/(1.D0+2.4D-3/(5.D-4+XE3)+
     +CPRO(K)*XE3))
      IF(K.EQ.3) THEN
c  nuclear-scattering contribution to HIP14:
        HIP14N(I)=7.95D0/AMEAN
        HIALL(JK)=HIALL(JK)+HIP14N(I)
      ENDIF
 1630 CONTINUE
C================================================================ end 5

C  !!! write over YI: QSDT(J) -> YI(J) !!!
C: "zero option" for thermonuclear burn rates:
 1700 QSDT(J)=0.D0
      QSDD(J)=0.D0
      QSDHE(J)=0.D0
      QSBH(J)=0.D0
C-------------------------------------------------------------- begin 6
C: thermonuclear burn rates:
      IF(IFDMFU(J).NE.1) GOTO 2000
      IF(TIJ.LT.TBURN0) GOTO 2000
      T13=TIJ**.333333333D0
      TM23=1.D0/T13**2
      QSDT(J)=1.58D4*TM23*((1.D0+.16D0*TIJ)*
     *DEXP(-19.98D0/T13-(TIJ/10.34D0)**2)+
     +.0108D0*DEXP(-45.07D0/TIJ))
      QSDD(J)=81.4D0*TM23*(1.D0+.01D0*TIJ)*DEXP(-18.81D0/T13)
      QSDHE(J)=1.3D4*TM23*(1.D0+5.D-4*TIJ**2)*
     *DEXP(-31.72D0/T13-(TIJ/27.14D0)**2)+
     +40.5D0*DEXP(-148.2D0/TIJ)/DSQRT(TIJ)
      QSBH(J)=5.05D4*TM23*(1.D0+.008D0*T13+.063D0*T13*T13+
     +.0034D0*TIJ+.0056D0*TIJ*T13+7.9D-4*TIJ*T13*T13)*
     *DEXP(-53.423D0/T13-(TIJ/174.06D0)**2)+
     +(59.2D0*DEXP(-149.33D0/TIJ)+6.51D4*DEXP(-618.45D0/TIJ))/
     /(TIJ*DSQRT(TIJ))+3.336D2*TM23*DEXP(-1094.D0/TIJ)
C: local component of the energy deposition by charged fusion products:
      FET=7.D0/(7.D0+TEJ)
      FEHE3=5.6D0/(5.6D0+TEJ)
      RAB=.5D0*XD(J)**2*QSDD(J)/(VJ*AMOL2)
      QE(J)=QE(J)+RAB*(.97D4*FEHE3+.79D4*FET)
      QI(J)=QI(J)+RAB*(.97D4*(1.D0-FEHE3)+
     +.79D4*(1.D0-FET))
C================================================================ end 6
 2000 CONTINUE
C: end of main loop.

C: calculate the Joule heating QJL(J) (to be used in STEP only):
      IF(IFHZ.NE.1) GOTO 20100
      HZASTP=HZ(1)
      IF(IFLBND.EQ.-1) HZASTP=HBLOLD
      RAB=(UCBET/(4.D0*PINUM))**2
      DO 20020 J=1,N
      HZAST=HZASTP
      IF(J.NE.N) HZASTP=HZ(J)+(HZ(J+1)-HZ(J))*
     *((R(J+1)-R(J))/(R(J+2)-R(J)))
      IF(J.EQ.N) HZASTP=HBROLD
20020 QJL(J)=RAB*ETHZ(J)*V(J)*((HZASTP-HZAST)/(R(J+1)-R(J)))**2
20100 CONTINUE

C: impose flux limits;  CAPE, CAPR  are non-limited conduction
C  coefficients at the centers of mesh cells;   CAPEB, CAPRB  are
C  flux limited coefficients at the interfaces of mesh cells:
      RLE=CAPEB(1)
      RLI=CAPIB(1)
      RLR=CAPRB(1)
      DO 2100 J=2,N
      RLEM=RLE
      RLIM=RLI
      RLRM=RLR
      RLE=CAPEB(J)
      RLI=CAPIB(J)
      RLR=CAPRB(J)

c: calculate RAB = elec.cond.coefficient at the interface r=R(J):
      RAB=.5D0*(CAPE(J-1)+CAPE(J))
c: impose flux limit:
      RAB1=TE(J)-TE(J-1)
      IF(RAB1.GT.0.D0) RLEM=RLE
      RLEM=RLEM*.5D0*(R(J+1)-R(J-1))
      IF(RAB*DABS(RAB1).LE.RLEM) GOTO 2050
      RAB=RLEM/DABS(RAB1)
      NTFLLE=NTFLLE+1
 2050 CAPEB(J)=RAB

c: calculate RAB = ion cond.coefficient at the interface r=R(J):
      RAB=.5D0*(CAPI(J-1)+CAPI(J))
c: impose flux limit:
      RAB1=TI(J)-TI(J-1)
      IF(RAB1.GT.0.D0) RLIM=RLI
      RLIM=RLIM*.5D0*(R(J+1)-R(J-1))
      IF(RAB*DABS(RAB1).LE.RLIM) GOTO 2060
      RAB=RLIM/DABS(RAB1)
      NTFLLI=NTFLLI+1
 2060 CAPIB(J)=RAB

c: calculate RAB = rad.cond.coefficient at the interface r=R(J):
      RAB=.5D0*(CAPR(J-1)+CAPR(J))
c: impose flux limit:
      RAB1=TR(J)-TR(J-1)
      IF(RAB1.GT.0.D0) RLRM=RLR
      RLRM=RLRM*.5D0*(R(J+1)-R(J-1))
      IF(RAB*DABS(RAB1).LE.RLRM) GOTO 2070
      RAB=RLRM/DABS(RAB1)
      NTFLLR=NTFLLR+1
 2070 CAPRB(J)=RAB
 2100 CONTINUE
      CAPEB(1)=CAPE(1)
      CAPIB(1)=CAPI(1)
      CAPRB(1)=CAPR(1)
      CAPEB(N1)=0.D0
      CAPIB(N1)=0.D0
      CAPRB(N1)=0.D0

      IF(IFBURN.EQ.1) IFBURN=2
      IF(RORFU.GT.RORFUM) RORFUM=RORFU
      IF(IFBURN.LE.0) GOTO 3100
C......................................................................

C             - NEUTRONS -

      IF(IFN2.EQ.0.AND.IFN14.EQ.0) GOTO 3100
C------------------------------------------------------------- begin 61
C: prepare for the diffusion-like spread of 14-MeV neutrons:
      IF(IFN14.NE.4) GOTO 2190
      IF(IGEO.NE.2) THEN
       WRITE(0,9010) IFN14,IGEO
       PRINT 9010,IFN14,IGEO
       STOP
      ENDIF

C: calculate TAU140 and the total N14 power (=Q14):
      JRAB0=1
      IRAB0=1
      RAB0= -1.D120
      Q14=0.D0
      I=1
      DO J=1,NFU
       IF(J.GE.NJZ(I+1)) I=I+1
       RAB=QSDT(J)*XD(J)*XT(J)*DM(J)/(V(J)*AMOL(I)**2)
       RAB00=RAB/(R(J+1)-R(J))
       IF(RAB00.GT.RAB0) THEN
        JRAB0=J
        IRAB0=I
        RAB0=RAB00
       ENDIF
       Q14=Q14+RAB
      ENDDO

      I=IRAB0
      DO J=JRAB0,NFU
       IF(J.GE.NJZ(I+1)) I=I+1
       RAB=QSDT(J)*XD(J)*XT(J)*DM(J)/(V(J)*AMOL(I)**2*(R(J+1)-R(J)))
       IF(RAB.LT..5D0*RAB0) THEN
        JTAU140=J-1
        GOTO 2120
       ENDIF
      ENDDO
      JTAU140=NFU
 2120 TAU140=RORAFJ(JTAU140+1)/20.D0

      RABC2=0.5D0*(TAU140-1.D0+(TAU140+1.D0)*EXP(-2.D0*TAU140))
      
 2190 CONTINUE
C--------------------------------------------------------------- end 61

C-------------------------------------------------------------- begin 7
C: lay away for the uniformly distributed neutron heating to be
C  calculated below:
      IF(IFN2.NE.3.AND.IFN14.NE.3) GOTO 2300
      Q2=0.D0
      Q14=0.D0
      I=1
      DO 2200 J=1,NFU
      IF(J.GE.NJZ(I+1)) I=I+1
      IF(IFN14.EQ.3) Q14=Q14+QSDT(J)*XD(J)*XT(J)*
     *DM(J)/(V(J)*AMOL(I)**2)
 2200 IF(IFN2.EQ.3) Q2=Q2+QSDD(J)*(XD(J)/AMOL(I))**2*
     *DM(J)/V(J)
      IF(IFN14.EQ.3) Q14=Q14*13.6D4*RORAFU/
     *((RORAFU+80.D0)*AMSFU)
      IF(IFN2.EQ.3) Q2=Q2*2.37D4*.5D0*RORAFU/
     *((RORAFU+26.D0)*AMSFU)
C================================================================ end 7

C: main cycle for neutron heating:
 2300 CONTINUE
      DO 3000 I=1,NZ
      AMEAN=AMOL(I)/XMOL(I)
      AMOL2=AMOL(I)**2
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 3000 J=JI,JF
      IF(IFDMFU(J).NE.1) GOTO 3000
      QN2=0.D0
      QN14=0.D0
      IF(J.GT.NFU) GOTO 2700

C: neutron heating of the central fuel sphere:
C-------------------------------------------------------------- begin 8
C: approximation of the  1-st scattering:
      IF(IFN2.NE.1.AND.IFN14.NE.1) GOTO 2500
      IF(IGEO.EQ.0) GOTO 2500
      RJ12=.5D0*(R(J)+R(J+1))
      IK=1

      DO 2490 K=1,NFU
      IF(K.GE.NJZ(IK+1)) IK=IK+1
      RK12=.5D0*(R(K)+R(K+1))
      IF(K.EQ.J) RK12=(R(K)+2.D0*R(K+1))/3.D0
      IF(IGEO.EQ.2) GJK=DLOG((RK12+RJ12)/DABS(RK12-
     *RJ12))/(RK12*RJ12)
      IF(IGEO.EQ.1) GJK=(2.D0/(RJ12+RK12))*
     *ELLICK(2.D0*DSQRT(RK12*RJ12)/(RK12+RJ12))
      IF(IFN14.EQ.1) QN14=QN14+(1100.D0/AMEAN)*QSDT(K)*
     *XD(K)*XT(K)*GJK*DM(K)/(V(K)*AMOL(IK)**2)
 2490 IF(IFN2.EQ.1) QN2=QN2+(325.D0/AMEAN)*QSDD(K)*
     *(XD(K)/AMOL(IK))**2*GJK*DM(K)/V(K)
C================================================================ end 8

C: branching for 14-MeV neutrons:
 2500 IF(IFN14.EQ.0) GOTO 2600
      GOTO (2600,2520,2540,2560),IFN14
C..local:
 2520 QN14=13.6D4*QSDT(J)*XD(J)*XT(J)/(V(J)*AMOL2)
      GOTO 2600
C..uniform distribution:
 2540 QN14=Q14
      GOTO 2600
C..diffusion-profile spread:
 2560 CONTINUE
      RABTAU=(RORAFJ(J)+RORAFJ(J+1))/40.D0
      IF(J.EQ.1) THEN
       RABTAUR=RORAFJ(J+1)/(20.D0*R(J+1))
      ELSE
       RABTAUR=(RORAFJ(J)/R(J)+RORAFJ(J+1)/R(J+1))/40.D0
      ENDIF

      IF(RABTAU.LT.TAU140) THEN
       RABPHI=1.D0-.5D0*(1.D0+TAU140)*(EXP(RABTAU-TAU140)-
     & EXP(-RABTAU-TAU140))/RABTAU
      ELSE
       RABPHI=RABC2*EXP(TAU140-RABTAU)/RABTAU
      ENDIF
      RABPHI=RABPHI*3.D0/TAU140**3
      QN14=13.6D4*Q14*RABPHI*RABTAUR**2/(20.D0*AMEAN)
 
C: branching for 2-Mev neutrons:
 2600 IF(IFN2.EQ.0) GOTO 2900
      GOTO (2900,2620,2640,2660),IFN2
C..local:
 2620 QN2=2.37D4*.5D0*QSDD(J)*XD(J)**2/(V(J)*AMOL2)
      GOTO 2900
C..uniform distribution:
 2640 QN2=Q2
      GOTO 2900
C..reserve:
 2660 CONTINUE
      GOTO 2900

C: neutron heating in shells detached from the central fuel sphere:
 2700 IF(IFN14.GT.0) QN14=13.6D4*QSDT(J)*XD(J)*XT(J)/
     *(V(J)*AMOL2)
      IF(IFN2.GT.0) QN2=2.37D4*.5D0*QSDD(J)*
     *XD(J)**2/(V(J)*AMOL2)

C: distribute neutron heating over mesh cells:
 2900 RAB=(1.D0+AMEAN)**2*TE(J)
      FEN14=750.D0/(750.D0+RAB)
      FEN2=130.D0/(130.D0+RAB)
      QE(J)=QE(J)+QN14*FEN14+QN2*FEN2
      QI(J)=QI(J)+QN14*(1.D0-FEN14)+QN2*(1.D0-FEN2)
C..for energy balance:
      WN14(I)=WN14(I)+QN14*DM(J)
      WN2(I)=WN2(I)+QN2*DM(J)
 3000 CONTINUE
 3100 CONTINUE
      RETURN
 9010 FORMAT(/'STOP in KINBUR: IFN14=',I6,'  IGEO=',I6
     &/'whereas (IFN14=4 and IGEO.NE.2) is unacceptable!')
      END

C**********************************************************************
      SUBROUTINE BNDVAL(TYME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the values of the boundary pressure, 
c     external radiation temperature, and magnetic field at TIME=TYME.

c     Called by:  
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'

C: boundary pressure:      ***********添加磁压的区间************************************************************************************************************
      PBL=0.D0
      PBR=0.D0

C: external X-ray drive temperature:
      TRLEX=0.D0
      TREX=0.D0*TYME

C: boundary magnetic z-field:
      HZBL=HZ0
      HZBR=HZ0

C: total external pressure:
      PBLSUM=PBL+ASBOL*TRLEX**4/3.D0
      PBRSUM=PBR+ASBOL*TREX**4/3.D0
      IF(IFHZ.NE.1) GOTO 290
      PBLSUM=PBLSUM+HZBL**2/(8.D0*PINUM)
      PBRSUM=PBRSUM+HZBR**2/(8.D0*PINUM)
 290  CONTINUE
      RETURN
      END

C**********************************************************************
      SUBROUTINE UPSLOI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine recalculates basic dynamic variables from
c     the old time TIME to the new time TIME+DT;

c     Called by:  DEIRA
c     Calls    :  
c =====================================================================
      INCLUDE 'de4com.fi'

      DIMENSION VN(1:NN+1),DVDT(1:NN+1),SCAVIS(1:NN+1),TENVIS(1:NN+1)
     *,IFPRO(3),QAL(1:NN+1),QP3(1:NN+1),QP14(1:NN+1)
     *,A0P(1:NN+1),B0P(1:NN+1),A1P(1:NN+1),B1P(1:NN+1),A3P(1:NN+1)
     *,B3P(1:NN+1),AA11(1:NN+1),AA12(1:NN+1),AA13(1:NN+1),BB1(1:NN+1)
     *,AA21(1:NN+1),AA22(1:NN+1),AA23(1:NN+1),BB2(1:NN+1)
     *,AA31(1:NN+1),AA32(1:NN+1),AA33(1:NN+1),BB3(1:NN+1)
     *,EALL(1:3*NN+3),HIALL(1:3*NN+3),QALL(1:3*NN+3),EOUALL(1:3*NNZ+3)
      EQUIVALENCE (PE(1),A0P(1),A1P(1),A3P(1),AA11(1))
     *,(EET(1),AA12(1))
     *,(EEV(1),AA13(1))
     *,(QE(1),BB1(1))
     *,(US(1),SCAVIS(1),VN(1))
     *,(PI(1),B0P(1),DVDT(1),AA21(1))
     *,(EIT(1),AA22(1))
     *,(EIV(1),AA23(1))
     *,(QI(1),BB2(1))
     *,(QJL(1),TENVIS(1),B3P(1),AA31(1))
     *,(QSDT(1),QAL(1),QALL(1),B1P(1),AA32(1))
     *,(QSDD(1),QP3(1),QALL(NN+2),AA33(1))
     *,(QSDHE(1),QP14(1),QALL(2*NN+3),BB3(1))
     *,(EAL(1),EALL(1)),(EP3(1),EALL(NN+2)),(EP14(1),EALL(2*NN+3))
     *,(HIAL(1),HIALL(1)),(HIP3(1),HIALL(NN+2)),(HIP14(1),HIALL(2*NN+3))
     *,(IFAL,IFPRO(1)),(IFP3,IFPRO(2)),(IFP14,IFPRO(3))
     *,(EOUTAL(1),EOUALL(1)),(EOUTP3(1),EOUALL(NNZ+2))
     *,(EOUP14(1),EOUALL(2*NNZ+3))
      COMMON/LEAKPW/PLKW(1:NNZ+1),PLKE(1:NNZ+1),PLKR(1:NNZ+1)
     *,PLKAL(1:NNZ+1),PLKP3(1:NNZ+1),PLKP14(1:NNZ+1)
      DIMENSION PLKALL(1:3*NNZ+3)
      EQUIVALENCE (PLKAL(1),PLKALL(1)),(PLKP3(1),PLKALL(NNZ+2))
     *,(PLKP14(1),PLKALL(2*NNZ+3))

C--------------------------------------------------------- begin NEWXDT
C: solve the thermonuclear burn equations:
      IF(IFBURN.LE.0) GOTO 110
      DO 100 I=1,NZ
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 100 J=JI,JF
      IF(IFDMFU(J).NE.1) GOTO 100
      XDJ=XD(J)
      QDTJ=QSDT(J)
      QDDJ=QSDD(J)
      QDHEJ=QSDHE(J)
C: injection of energy into fast products:
C !!! write over QSDT,QSDD,QSDHE:  QAL(J) -> QSDT(J)
C     QP3(J) -> QSDD(J),  QP14(J) -> QSDHE(J)   !!!
      RAB=1.D0/(V(J)*AMOL(I)**2)
      QAL(J)=RAB*(3.4D4*XDJ*XT(J)*QDTJ+
     +3.54D4*XDJ*XHE(J)*QDHEJ+8.38D4*XH(J)*XB(J)*QSBH(J))
      QP3(J)=2.91D4*RAB*XDJ*XDJ*QDDJ/2.D0
      QP14(J)=14.16D4*RAB*XDJ*XHE(J)*QDHEJ
C: new abundances for DD group:
      RAB=DT/(V(J)*AMOL(I))
      RAB1=RAB*XDJ*QDDJ
      XT(J)=(XT(J)+.5D0*RAB1*XDJ)/(1.D0+RAB*XDJ*QDTJ)
      XHE(J)=(XHE(J)+.5D0*RAB1*XDJ)/(1.D0+RAB*XDJ*QDHEJ)
      XD(J)=XDJ*(1.D0-2.D0*RAB1-RAB*(XT(J)*QDTJ+
     +XHE(J)*QDHEJ))
C: new abundances for BH group:
      RAB=(DT*QSBH(J))/(V(J)*AMOL(I))
      IF(XB(J).GT.XH(J)) GOTO 80
      RAB1=XB(J)
      XB(J)=RAB1/(1.D0+RAB*XH(J))
      XH(J)=XH(J)+XB(J)-RAB1
      GOTO 100
 80   RAB1=XH(J)
      XH(J)=RAB1/(1.D0+RAB*XB(J))
      XB(J)=XB(J)+XH(J)-RAB1
 100  CONTINUE
 110  CONTINUE
C  burn equations are solved.
C=========================================================== end NEWXDT

C---------------------------------------------------------- begin NEWRU
C: solve the Euler equation:
C: boundary conditions, initial assignments:
      DT5=.5D0*DT
      CALL BNDVAL(TIME+DT5)
      DM(N1)=0.D0
      IFVOCL=0
      R1OLD=R(1)
      U1OLD=U(1)

C: lay away for the work at boundaries:
      WN14(1)=0.D0
      WN2(1)=0.D0
      IF(IFLBND.EQ.0) GOTO 260
      WN14(1)=U(1)*DT5
      IF(IFLBND.EQ.1) GOTO 240
      WN2(1)=PBLSUM
      GOTO 260
 240  WN2(1)=PBL+ASBOL*TR(1)**4/3.D0
      IF(IFHZ.EQ.1) WN2(1)=WN2(1)+HZ(1)**2/(8.D0*PINUM)
      IF(IFBURN.EQ.2) WN2(1)=WN2(1)+2.D0*(EAL(1)+EP3(1)+EP14(1))/3.D0
 260  CONTINUE

      DO 280 I=2,NZ
      J=NJZ(I)
      K=J-1
      WN14(I)=U(J)*DT5
      PO=.5D0*(PE(K)+PE(J)+PI(K)+PI(J)+
     +ASBOL*(TR(K)**4+TR(J)**4)/3.D0)
      IF(IFHZ.EQ.1) PO=PO+(HZ(K)**2+HZ(J)**2)/(16.D0*PINUM)
      IF(IFBURN.EQ.2) PO=PO+(EAL(K)+EAL(J)+
     +EP3(K)+EP3(J)+EP14(K)+EP14(J))/3.D0
 280  WN2(I)=PO
      WN14(NZ+1)=U(N1)*DT5
      WN2(NZ+1)=PBRSUM

C: check for central void closure (only when IFLBND=1):
      IF(IFLBND.EQ.0) GOTO 310
      RAB=R(1)+U(1)*DT5
      IF(IFLBND.NE.1) GOTO 300
      IF(RAB.GT..5D0*R(1)) GOTO 300
C - central void is closed at this time step !!!
      IFVOCL=1
      RAB=.5D0*R(1)
C: radii at half-step:
 300  R(1)=RAB
 310  CONTINUE
      DO 320 J=2,N1
 320  R(J)=R(J)+U(J)*DT5

C: progon coefficients for velocity: C  !!!  write over PE and PI  !!!
C  !!!  A0P(J) -> PE(J), B0P(J) -> PI(J)  !!!
      RS=0.D0
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(1)**IGEO
      DIV=0.D0
      RABDM=0.D0
      RABSV=0.D0
      RABTV=0.D0
      PO=0.D0
      IF(IFLBND.EQ.0) GOTO 330
c: left boundary pressure for "open halfspace":
      PO=PBLSUM
      IF(IFLBND.EQ.-1) GOTO 330
c: left boundary pressure for "closed cavity":
      PO=PBL+ASBOL*TR(1)**4/3.D0
      IF(IFHZ.EQ.1) PO=PO+HZ(1)**2/(8.D0*PINUM)
      IF(IFBURN.EQ.2) PO=PO+2.D0*(EAL(1)+EP3(1)+EP14(1))/3.D0
 330  CONTINUE

      DO 600 J=1,N1
      JM1=J-1
      IF(JM1.LT.1) JM1=1
      RSM=RS
      RS=RSP
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(J+1)**IGEO
      DIVM=DIV
      DIV=U(J+1)*RSP-U(J)*RS
      RABDMM=RABDM
      RABDM=DM(J)
      RABSVM=RABSV
      RABTVM=RABTV
      RBYRP=0.D0
      IF(IFLBND.NE.-1) RBYRP=R(J)/R(J+1)
C: calculate scalar viscosity SCAVIS(J)=(\eta_sca,j)/(\Delta m_j) and
C  tensor viscosity TENVIS(J)=(\eta_ten,j)/(\Delta m_j r_j+1^s+2):
      RABSV=0.D0
      RABTV=0.D0
      IF(J.EQ.N1) GOTO 490
      RABS=0.D0
      RABT=0.D0
      RABDU=U(J)-U(J+1)
      IF(RABDU.LE.0.D0) GOTO 450
      RABS=(SMU1*US(J)+SMU2*RABDU)/((1.D0-SIG)*RS+SIG*RSP)
      IF(IFLBND.NE.-1) RABT=(TMU1*US(J)+TMU2*RABDU)*
     *(SIGT+(1.D0-SIGT)*(RS/RSP)*RBYRP*RBYRP)
C: US(J) is free; add physical (ion) viscosity:
 450  IF(IGEO.GE.1) GOTO 460
      RABSF=ETVI0(J)/3.D0+ETVI1(J)
      RABTF=0.D0
      GOTO 480
 460  IF(IGEO.GE.2) GOTO 470
      RABSF=ETVI0(J)/3.D0
      RABTF=ETVI1(J)
      GOTO 480
 470  RABSF=0.D0
      RABTF=1.3333333333D0*ETVI0(J)
 480  RABSV=(RABS+RABSF/DM(J))/V(J)
      IF(IFLBND.NE.-1) RABTV=(RABT+RABTF*(SIGT*RSP+
     +(1.D0-SIGT)*RS*(RS/RSP)*RBYRP*RBYRP)/DM(J))/V(J)
C  !!! write over US(J): SCAVIS(J) -> US(J)  !!!
C  !!! write over QJL(J): TENVIS(J) -> QJL(J)  !!!
 490  SCAVIS(J)=RABSV
      TENVIS(J)=RABTV

C: clear out PE(J) and PI(J) and calculate pressure PO:
      EEV(J)=EEV(J)+PE(J)
      EIV(J)=EIV(J)+PI(J)
      POM=PO
      IF(J.LE.N) GOTO 520
      PO=PBRSUM
      GOTO 540
 520  RAB=TR(J)**2
      PO=PE(J)+PI(J)+(ASBOL*RAB)*RAB/3.D0
      IF(IFHZ.EQ.1) PO=PO+HZ(J)**2/(8.D0*PINUM)
      IF(IFBURN.EQ.2) PO=PO+2.D0*(EAL(J)+EP3(J)+EP14(J))/3.D0

 540  RABTAU=RS*(DT/(RABDM+RABDMM))
      IF(J.EQ.1.AND.IFLBND.EQ.0) RABTAU=0.D0
      RPBRS2=0.D0
      IF(IFLBND.EQ.-1) GOTO 570
      IF(J.EQ.1.AND.IFLBND.EQ.0) GOTO 570
      RPBRS2=(RSP/RS)/RBYRP**2
 570  RBYRM=0.D0
      IF(IFLBND.EQ.-1) GOTO 590
      IF(J.GE.3) RBYRM=R(J)/R(JM1)
      IF(J.EQ.2.AND.IFLBND.EQ.1) RBYRM=R(J)/R(JM1)
 590  RABOMG=1.D0+RABTAU*(RABSV*RS+RABSVM*(RS-RSM*A0P(JM1))+
     +RABTV*RPBRS2+RABTVM*(1.D0-RBYRM*A0P(JM1)))
      A0P(J)=RABTAU*(RABSV*RSP+RABTV*RPBRS2*RBYRP)/RABOMG
      B0P(J)=(U(J)+RABTAU*(2.D0*(POM-PO)+RABSV*DIV-RABSVM*
     *(DIVM-RSM*B0P(JM1))+RABTV*RPBRS2*(U(J+1)*RBYRP-U(J))-
     *RABTVM*(U(J)-RBYRM*(U(JM1)+B0P(JM1)))))/RABOMG
 600  CONTINUE

C: reverse progon, new velocities, dV/dt, viscous dissipation:
      RS=1.D0
      IF(IGEO.GE.1) RS=R(N1)**IGEO
      UO=U(N1)
      U(N1)=B0P(N1)
      UU=.5D0*(UO+U(N1))
      U(N2)=U(N1)

      DO 800 K=1,N
      J=N1-K
      RSP=RS
      RS=1.D0
      IF(IGEO.GE.1) RS=R(J)**IGEO
      UO=U(J)
      U(J)=A0P(J)*U(J+1)+B0P(J)
c: check for void closure:
      IF(J.NE.1) GOTO 740
      IF(IFLBND.EQ.0.OR.IFLBND.EQ.-1) GOTO 740
      IF(IFVOCL.EQ.1) GOTO 730
      RABR1=R(1)+U(1)*DT5
      IF(RABR1.LT.1.D-2*ABS(R(1))) GOTO 730
      GOTO 740
c: correct the velocity U(1) for exact void closure:
 730  U(1)=-R(1)/DT5
      IFVOCL=1
 740  UUP=UU
      UU=.5D0*(UO+U(J))
C: calculate dV/dt:
C  !!!  write over B0P(J)=PI(J): DVDT(J) -> B0P(J)=PI(J)  !!!
      DVDT(J)=(UUP*RSP-UU*RS)/DM(J)
C: redefine SCAVIS(J) and TENVIS(J):
      SCAVIS(J)=(DVDT(J)*SCAVIS(J))*DM(J)
      RABSI=0.D0
      IF(IFLBND.EQ.-1) GOTO 790
      IF(IFLBND.EQ.1.OR.J.NE.1) GOTO 770
      RABSI=UUP/R(J+1)
      GOTO 780
 770  RABSI=UUP/R(J+1)-UU/R(J)
 780  TENVIS(J)=(RABSI*TENVIS(J))*RSP*R(J+1)**2
C: viscous dissipation - all into ions(!):
 790  QI(J)=QI(J)+DVDT(J)*SCAVIS(J)+(RABSI/DM(J))*TENVIS(J)
 800  CONTINUE
C: add kinetic energy dissipation due to void closure:
      IF(IFLBND.EQ.1.AND.IFVOCL.EQ.1) QI(1)=QI(1)+
     +.125D0*U1OLD**2/DT

C: work at boundaries = outward leak of internal energy:
      PLKW(1)=0.D0
      IF(IFLBND.EQ.0) GOTO 860
      RS=1.D0
      IF(IGEO.GE.1) RS=R(1)**IGEO
      PLKW(1)=CSURF*(WN14(1)/DT+.5D0*U(1))*(WN2(1)*RS)
      EOUTW(1)=EOUTW(1)+CSURF*(WN14(1)+DT5*U(1))*(WN2(1)*RS)
 860  CONTINUE

      DO 900 I=2,NZ+1
      J=NJZ(I)
      K=J-1
      RS=1.D0
      IF(IGEO.GE.1) RS=R(J)**IGEO
      RAB=0.D0
      IF(I.EQ.NZ+1) GOTO 880
      RAB=RS*(SCAVIS(J)+SCAVIS(K))
      IF(IFLBND.NE.-1) RAB=RAB+(TENVIS(J)+TENVIS(K))/R(J)
 880  PLKW(I)=CSURF*(WN14(I)/DT+.5D0*U(J))*
     *(WN2(I)*RS-.5D0*RAB)
 900  EOUTW(I)=EOUTW(I)+CSURF*(WN14(I)+DT5*U(J))*
     *(WN2(I)*RS-.5D0*RAB)

C: calculate new radii and specific volumes VN(J):
      VINCR=0.D0
      VINCD=0.D0
      IF(IFLBND.EQ.0) GOTO 930
      IF(IFVOCL.NE.1) GOTO 940
c: change the type of boundary condition in the case of void closure:
      PRINT 9110,NTSLOI,TIME,R1OLD,U1OLD
 9110 FORMAT(/1X,'!!!  Central void closure at NTSLOI=',I8
     *,'  TIME=',1PE16.8,'  !!!'/1X,'void closure: R_1,old=',1PE16.8
     *,',  U_1,old=',1PE16.8)
      IFLBND=0
 930  R(1)=0.D0
      RSP=0.D0
      U(1)=0.D0
      GOTO 960
 940  RAB=U(1)*DT5
      R(1)=R(1)+RAB
      RAB1=R(1)+RAB
      RSP=RAB1
      IF(IGEO.GE.1) RSP=RAB1**I1GEO
 960  CONTINUE

      DO 1000 J=1,N
      RS=RSP
      RAB=U(J+1)*DT5
      R(J+1)=R(J+1)+RAB
      RAB1=R(J+1)+RAB
      RSP=RAB1
      IF(IGEO.GE.1) RSP=RAB1**I1GEO
C: check new radii for monotonousness:
      IF(RSP.GT.RS) GOTO 980
      PRINT 9001,NTSLOI,J,DT,R,U
      STOP
C !!!  write over US: VN(J) -> US(J)  !!!
 980  VN(J)=(RSP-RS)/(S1GEO*DM(J))
C: maximum change in specific volume:
      RAB=VN(J)/V(J)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      RAB=RAB-1.D0
      IF(RAB.GT.VINCR) VINCR=RAB
      IF(RAB.GT.VINCD.AND.QDRIV(J).GT.0.D0) VINCD=RAB
 1000 CONTINUE
      R(N2)=R(N1)
C  the Euler equation is solved.
C  in c/blk /XXX/ vacant are: PE,WN14,WN2
C============================================================ end NEWRU

C--------------------------------------------------------- begin NEWEAL
C: solve diffusion equations for the energy of charged fusion products;
C  add their energy deposition to QE(J), QI(J):
      IF(IFBURN.LE.0) GOTO 2000
      DO 1600 K=1,3
      IF(IFPRO(K).EQ.0) GOTO 1600
      J0=N1MAX*(K-1)
C------------------------------------------------------------- begin 11
      IF(IFPRO(K).EQ.1) GOTO 1200
C: local energy deposition:
      DO 1100 J=1,N
      JAL=J0+J
      EALL(JAL)=QALL(JAL)/HIALL(JAL)
      RAB=TE(J)/TPRO(K)
      QE(J)=QE(J)+QALL(JAL)/(1.D0+RAB)
 1100 QI(J)=QI(J)+RAB*QALL(JAL)/(1.D0+RAB)
      GOTO 1600
C=============================================================== end 11
 1200 CONTINUE
C------------------------------------------------------------- begin 12
C: diffusion of the energy of charged products:
C: boundary conditions:
      J=J0+N1
      EALL(J)=0.D0
      HIALL(J)=HIALL(J-1)
      V(N1)=V(N)
C: progon coefficients:
      FP=0.D0
      PO=.125D0*UPRO(K)**2
      ZBYAAL=964.9D0*ZPRO(K)/APRO(K)
      RAB=HIALL(J0+1)/V(1)
      IF(IFHZ.EQ.1) RAB=RAB+.5D0*(ZBYAAL*HZ(1))**2/RAB
      IF(K.EQ.3) RAB=RAB+.25D0*HIP14N(1)/V(1)
      DALP=PO/RAB

      DO 1400 J=1,N
      JM1=J-1
      IF(JM1.LT.1) JM1=1
      JAL=J0+J
      F=FP
      DAL=DALP
      RAB=HIALL(JAL+1)/V(J+1)
      IF(IFHZ.EQ.1) RAB=RAB+.5D0*(ZBYAAL*HZ(J+1))**2/RAB
      IF(K.EQ.3) RAB=RAB+.25D0*HIP14N(ISHLJ(J+1))/V(J+1)
      DALP=PO/RAB
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(J+1)**IGEO
      FP=RSP*(DT/(R(J+2)-R(J)))*(DAL+DALP)
      POM=VN(J)+DT*HIALL(JAL)+(FP+F*(1.D0-A1P(JM1)))/DM(J)
      IF(J.EQ.1.AND.IFLBND.EQ.-1)
     *POM=POM+2.D0*DAL*(DT/(R(2)-R(1)))/DM(1)
C  !!!  write over A0P=PE: A1P(J) -> A0P(J)=PE(J)  !!!
      A1P(J)=FP/(POM*DM(J))
C  !!!  write over QAL=QSDT: B1P(J) -> QAL(J)=QSDT(J)  !!!
 1400 B1P(J)=(V(J)*EALL(JAL)+F*B1P(JM1)/DM(J)+
     +DT*(QALL(JAL)-2.D0*EALL(JAL)*DVDT(J)/3.D0))/POM

C: reverse progon, contribution to the heating:
      DO 1500 JF=1,N
      J=N1-JF
      JAL=J0+J
      EALL(JAL)=A1P(J)*EALL(JAL+1)+B1P(J)
      RAB=TE(J)/TPRO(K)
      IF(K.EQ.3) THEN
        RAB1=EALL(JAL)*HIP14N(ISHLJ(J))
        RAB2=HIALL(JAL)-HIP14N(ISHLJ(J))
        IF(RAB2.LT.0.D0) THEN
          WRITE(*,9040) J,HIALL(JAL),HIP14N(ISHLJ(J))
          STOP
        ENDIF
        RAB2=EALL(JAL)*RAB2
        RAB3=(1.D0+AMOL(ISHLJ(J))/XMOL(ISHLJ(J)))**2*TE(J)/750.D0
        QE(J)=QE(J)+RAB1/(1.D0+RAB3)+RAB2/(1.D0+RAB)
        QI(J)=QI(J)+RAB1*RAB3/(1.D0+RAB3)+RAB2*RAB/(1.D0+RAB)
      ELSE
        RAB1=EALL(JAL)*HIALL(JAL)
        QE(J)=QE(J)+RAB1/(1.D0+RAB)
        QI(J)=QI(J)+RAB1*RAB/(1.D0+RAB)
      ENDIF
 1500 CONTINUE

C: outward leak of energy from layers - for balance:
      JF=(NZMAX+1)*(K-1)+1
      PLKALL(JF)=0.D0
      IF(IFLBND.NE.-1) GOTO 1520
      JAL=J0+1
      RAB=HIALL(JAL)/V(1)
      IF(IFHZ.EQ.1) RAB=RAB+.5D0*(ZBYAAL*HZ(1))**2/RAB
      DAL=PO/RAB
      PLKALL(JF)=-2.D0*CSURF*DAL*EALL(JAL)/(R(2)-R(1))
      EOUALL(JF)=EOUALL(JF)-2.D0*CSURF*DAL*EALL(JAL)*(DT/(R(2)-R(1)))
 1520 CONTINUE

      DO 1590 I=2,NZ+1
      J=NJZ(I)
      JAL=J0+J
      JF=(NZMAX+1)*(K-1)+I
      RAB=HIALL(JAL)/V(J)
      IF(IFHZ.EQ.1) RAB=RAB+.5D0*(ZBYAAL*HZ(J))**2/RAB
      DALP=PO/RAB
      RAB=HIALL(JAL-1)/V(J-1)
      IF(IFHZ.EQ.1) RAB=RAB+.5D0*(ZBYAAL*HZ(J-1))**2/RAB
      DAL=PO/RAB
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(J)**IGEO
      PLKALL(JF)=CSURF*RSP*(DAL+DALP)*
     *(EALL(JAL-1)-EALL(JAL))/(R(J+1)-R(J-1))
 1590 EOUALL(JF)=EOUALL(JF)+CSURF*RSP*(DAL+DALP)*
     *(EALL(JAL-1)-EALL(JAL))*(DT/(R(J+1)-R(J-1)))
C=============================================================== end 12
 1600 CONTINUE
 2000 CONTINUE
C  diffusion equations for charged products are solved;
C  in c/blk /XXX/ vacant are: PE,QSDT,QSDD,QSDHE,QJL,WN14,WN2.
C=========================================================== end NEWEAL

C---------------------------------------------------------- begin NEWHZ
C: solve diffusion equation for z-component of magnetic field HZ(J):
      IF(IFHZ.NE.1) GOTO 4490
C: magnetic field at the outer boundary (HZ_tilde, for Joule heating:
      QJL(N1)=.5D0*(HZBR+HBROLD)
      HTILBL=.5D0*(HZBL+HBLOLD)
C: boundary conditions:
      ETHZ(N1)=ETHZ(N)
      HZ(N1)=HZBR
C: progon coefficients:
      FP=0.D0
      RAB=UCBET**2/(4.D0*PINUM)
      DO 4200 J=1,N
      JM1=J-1
      IF(JM1.LT.1) JM1=1
      F=FP
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(J+1)**IGEO
      FP=RSP*(DT/(R(J+2)-R(J)))*(RAB*(ETHZ(J)+ETHZ(J+1)))
      POM=VN(J)+(FP+F-F*A3P(JM1))/DM(J)
      IF(J.EQ.1.AND.IFLBND.EQ.-1)
     *POM=POM+2.D0*RAB*ETHZ(1)*(DT/(R(2)-R(1)))/DM(1)
C !!!  write over PE=A0P=A1P: A3P(J) -> PE(J)=A0P(J)=A1P(J)  !!!
      A3P(J)=FP/(POM*DM(J))
C !!!  write over QJL: B3P(J) -> QJL(J)  !!!
      RAB1=V(J)*HZ(J)+F*B3P(JM1)/DM(J)
      IF(J.EQ.1.AND.IFLBND.EQ.-1)
     *RAB1=RAB1+2.D0*RAB*ETHZ(1)*HZBL*(DT/(R(2)-R(1)))/DM(1)
 4200 B3P(J)=RAB1/POM

      HZINCR=0.D0
C: reverse progon:
      DO 4250 JF=1,N
      J=N1-JF
      HZO=HZ(J)
      HZ(J)=A3P(J)*HZ(J+1)+B3P(J)
      IF(HZ(J).GE.0.D0) GOTO 4230
      PRINT 9020,TIME,NTSLOI,J,HZ(J)
      HZ(J)=0.D0
 4230 RAB=HZO+HZMIN
C: memorize H_tilde(J) -> QJL(J): C !!!  write over B3P(J)  !!!
      QJL(J)=.5D0*(HZO+HZ(J))
C: maximum change in magnetic field:
      RAB=(HZ(J)+HZMIN)/RAB
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      RAB=RAB-1.D0
      IF(RAB.GT.HZINCR) HZINCR=RAB
 4250 CONTINUE

      RAB9=.5D0*(UCBET/(4.D0*PINUM))**2
C: outward leak of magnetic energy - for balance:
      IF(IFLBND.NE.-1) GOTO 4280
      EOUTHZ(1)=EOUTHZ(1)+DT*CSURF*
     *(RAB9*2.D0*ETHZ(1))*(HTILBL+QJL(1))*
     *((HZBL-HZ(1))/(R(2)-R(1)))

 4280 CONTINUE
      DO 4300 I=2,NZ+1
      J=NJZ(I)
      RS=1.D0
      IF(IGEO.GE.1) RS=R(J)**IGEO
 4300 EOUTHZ(I)=EOUTHZ(I)+DT*CSURF*RS*
     *(RAB9*(ETHZ(J-1)+ETHZ(J)))*(QJL(J-1)+QJL(J))*
     *((HZ(J-1)-HZ(J))/(R(J+1)-R(J-1)))

C: calculate the Joule heating QJLJ and add to QE(J) and EZJL(I):
      FP=0.D0
      DHTILP=0.D0
      DO 4350 I=1,NZ
      JI=NJZ(I)
      JF=NJZ(I+1)-1
      DO 4350 J=JI,JF
      F=FP
      DHTIL=DHTILP
      RSP=1.D0
      IF(IGEO.GE.1) RSP=R(J+1)**IGEO
      FP=RSP*(RAB9*(ETHZ(J)+ETHZ(J+1)))*
     *((HZ(J)-HZ(J+1))/(R(J+2)-R(J)))
      DHTILP=QJL(J)-QJL(J+1)
      QJLJ=F*(DHTIL/DM(J))+FP*(DHTILP/DM(J))
      IF(J.EQ.1.AND.IFLBND.EQ.-1)
     *QJLJ=QJLJ+2.D0*RAB9*ETHZ(1)*((HZBL-HZ(1))/(R(2)-R(1)))*
     *((HTILBL-QJL(1))/DM(1))
      IF(QJLJ.LT.0.D0) QJLJ=0.D0
      QE(J)=QE(J)+QJLJ
 4350 EZJL(I)=EZJL(I)+DT*CSURF*(QJLJ*DM(J))
 4490 CONTINUE
C  diffusion equation for z-component of magnetic field is solved;
C: in c/blk /XXX/ vacant are: PE,QSDT,QSDD,QSDHE,QJL,WN14,WN2.
C============================================================ end NEWHZ

C--------------------------------------------------------- begin NEWTTT
C: solve the energy equations for electrons, ions and radiation:
      DO 2100 J=1,N
      QE(J)=QE(J)-EEV(J)*DVDT(J)
 2100 QI(J)=QI(J)-EIV(J)*DVDT(J)
C  vacant are: PE,EEV,EIV,QSDT,QSDD,QSDHE,QJL

C: boundary coditions for matrix progon:
      TRGR=TR(N)
      TRGR1=TR(1)
      RAB=1.D0
      IF(IGEO.GE.1) RAB=R(N1)**IGEO
      QRGR=.25D0*UCBET*ASBOL*RAB*(DT/DM(N))
      QRGR1=.25D0*UCBET*ASBOL*(DT/DM(1))
      CAPEB(N1)=0.D0
      CAPIB(N1)=0.D0
      CAPRB(N1)=0.D0
C: matrix progon:
      FEP=0.D0
      FRP=0.D0
      FIP=0.D0
      DO 2700 J=1,N
      JM1=J-1
      IF(JM1.LT.1) JM1=1
      RAB=DM(JM1)/DM(J)
      FE=FEP*RAB
      FR=FRP*RAB
      FI=FIP*RAB
      RAB1=1.D0
      IF(IGEO.GE.1) RAB1=R(J+1)**IGEO
      RAB=2.D0*RAB1*(DT/(R(J+2)-R(J)))/DM(J)
      FEP=RAB*CAPEB(J+1)
      FRP=RAB*CAPRB(J+1)
      FIP=RAB*CAPIB(J+1)
      HIEIDT=DT*HIEI(J)
      HIERDT=DT*HIER(J)
      G11X=EET(J)+FEP+FE-FE*AA11(JM1)
      G11=G11X+HIERDT+HIEIDT
      G12X=-FE*AA12(JM1)
      G12=G12X-HIERDT
      G13X=-FE*AA13(JM1)
      G13=G13X-HIEIDT
      Q1=EET(J)*TE(J)+QE(J)*DT+FE*BB1(JM1)
      G21X=-FR*AA21(JM1)
      G21=G21X-HIERDT
      G22X=4.D0*ASBOL*TR(J)**3*VN(J)+FRP+FR-FR*AA22(JM1)
      IF(J.EQ.1.AND.IFLBND.EQ.-1) G22X=G22X+QRGR1*TR(1)**3
      IF(J.EQ.N) G22X=G22X+QRGR*TR(N)**3
      G22=G22X+HIERDT
      G23=-FR*AA23(JM1)
      Q2=ASBOL*TR(J)**4*(3.D0*VN(J)+V(J)-DT*DVDT(J)/3.D0)+FR*BB2(JM1)
      IF(J.EQ.N) Q2=Q2+TREX**2*(QRGR*TREX**2)
      IF(J.EQ.1.AND.IFLBND.EQ.-1) Q2=Q2+TRLEX**2*(QRGR1*TRLEX**2)
      G31X=-FI*AA31(JM1)
      G31=G31X-HIEIDT
      G32=-FI*AA32(JM1)
      G33X=EIT(J)+FIP+FI-FI*AA33(JM1)
      G33=G33X+HIEIDT
      Q3=EIT(J)*TI(J)+QI(J)*DT+FI*BB3(JM1)
C: calculate DET=|G|; modify row 3 by adding row 2:
      G3111=G31X+G11X+HIERDT
      G3212=G32+G12
      G3313=G33X+G13X
      G2111=G21X+G11X+HIEIDT
      G2212=G22X+G12X
      DET=G3111*(G12*G23-G13*G22)+G3212*(G13*G21-G11*G23)+
     +G3313*(G11*G2212-G12*G2111)
      IF(DABS(DET).GT.0.D0) GOTO 2510
      PRINT 9030,TIME,NTSLOI,J
      STOP
C: invert matrix G: R=G**-1:
 2510 R11=(G22*G33-G23*G32)/DET
      R12=(G13*G32-G12*G33)/DET
      R13=(G12*G23-G13*G22)/DET
      R21=(G23*G31-G21*G33)/DET
      R22=(G11*G3313-G13*G3111)/DET
      R23=(G13*G21-G11*G23)/DET
      R31=(G21*G32-G22*G31)/DET
      R32=(G12*G31-G11*G32)/DET
      R33=(G11*G2212-G12*G2111)/DET
C  !!!  AA11(J) -> PE(J);   AA12(J) -> EET(J);  AA13(J) -> EEV(J);  !!!
C  !!!  AA21(J) -> PI(J);   AA22(J) -> EIT(J);  AA23(J) -> EIV(J);  !!!
C  !!!  AA31(J) -> QJL(J);  AA32(J) -> QSDT(J); AA33(J) -> QSDD(J); !!!
C  !!!  BB1(J) ->  QE(J);   BB2(J) ->  QI(J);   BB3(J) -> QSDHE(J); !!!
      AA11(J)=R11*FEP
      AA12(J)=R12*FRP
      AA13(J)=R13*FIP
      AA21(J)=R21*FEP
      AA22(J)=R22*FRP
      AA23(J)=R23*FIP
      AA31(J)=R31*FEP
      AA32(J)=R32*FRP
      AA33(J)=R33*FIP
      BB1(J)=R11*Q1+R12*Q2+R13*Q3
      BB2(J)=R21*Q1+R22*Q2+R23*Q3
 2700 BB3(J)=R31*Q1+R32*Q2+R33*Q3

C: reverse progon, reassign densities:
      TEINCR=0.D0
      TIINCR=0.D0
      TRINCR=0.D0
      TEINCD=0.D0
      DO 2800 JF=1,N
      J=N1-JF
      RAB=TE(J)+TMIN
      TE(J)=AA11(J)*TE(J+1)+AA12(J)*TR(J+1)+AA13(J)*TI(J+1)+BB1(J)
      RAB1=TE(J)
      IF(RAB1.LT.TMIN) RAB1=TMIN
      RAB=RAB/(RAB1+TMIN)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      RAB=RAB-1.D0
      IF(RAB.GT.TEINCR) TEINCR=RAB
      IF(RAB.GT.TEINCD.AND.QDRIV(J).GT.0.D0) TEINCD=RAB
      RAB=TR(J)+TMIN
      TR(J)=AA21(J)*TE(J+1)+AA22(J)*TR(J+1)+AA23(J)*TI(J+1)+BB2(J)
      RAB1=TR(J)
      IF(RAB1.LT.TMIN) RAB1=TMIN
      RAB=RAB/(RAB1+TMIN)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      RAB=RAB-1.D0
      IF(RAB.GT.TRINCR) TRINCR=RAB
      RAB=TI(J)+TMIN
      TI(J)=AA31(J)*TE(J+1)+AA32(J)*TR(J+1)+AA33(J)*TI(J+1)+BB3(J)
      RAB1=TI(J)
      IF(RAB1.LT.TMIN) RAB1=TMIN
      RAB=RAB/(RAB1+TMIN)
      IF(RAB.LT.1.D0) RAB=1.D0/RAB
      RAB=RAB-1.D0
      IF(RAB.GT.TIINCR) TIINCR=RAB
 2800 V(J)=VN(J)

C: outward leak of energy from layers - for balance:
      PLKE(1)=0.D0
      PLKR(1)=0.D0
      IF(IFLBND.NE.-1) GOTO 2890
      PLKR(1)=(CSURF*.25D0*UCBET*ASBOL)*(TRLEX**4-TRGR1**3*TR(1))
      EOUTR(1)=EOUTR(1)+(CSURF*QRGR1*DM(1))*(TRLEX**4-TRGR1**3*TR(1))
 2890 CONTINUE
      DO 2900 I=2,NZ
      J=NJZ(I)
      RAB=1.D0
      IF(IGEO.GE.1) RAB=R(J)**IGEO
      RSLK=(CSURF/(R(J+1)-R(J-1)))*RAB
      RS=RSLK*DT
      PLKE(I)=RSLK*2.D0*CAPEB(J)*(TE(J-1)-TE(J))
      PLKR(I)=RSLK*2.D0*CAPRB(J)*(TR(J-1)-TR(J))
      EOUTE(I)=EOUTE(I)+RS*2.D0*CAPEB(J)*(TE(J-1)-TE(J))
      EOUTI(I)=EOUTI(I)+RS*2.D0*CAPIB(J)*(TI(J-1)-TI(J))
 2900 EOUTR(I)=EOUTR(I)+RS*2.D0*CAPRB(J)*(TR(J-1)-TR(J))
      EOUTR(NZ+1)=EOUTR(NZ+1)+(CSURF*QRGR*DM(N))*
     *(TRGR**3*TR(N)-TREX**4)
      PLKR(NZ+1)=(CSURF*(QRGR/DT)*DM(N))*(TRGR**3*TR(N)-TREX**4)
      EREX=EREX+CSURF*QRGR*DM(N)*TREX**4
      IF(IFLBND.EQ.-1) EREX=EREX+CSURF*QRGR1*DM(1)*TRLEX**4
      ERIN=EOUTR(1)-EOUTR(NZ+1)

C: lower bound on temperatures:
      DO 2950 J=1,N
      IF(TE(J).LT.TMIN) TE(J)=TMIN
      IF(TR(J).LT.TMIN) TR(J)=TMIN
 2950 IF(TI(J).LT.TMIN) TI(J)=TMIN
C: the energy equations are solved.
C=========================================================== end NEWTTT
c: memorize the boundary values used:
      PBLOLD=PBL
      PBROLD=PBR
      PBLSOL=PBLSUM
      PBRSOL=PBRSUM
      HBLOLD=HZBL
      HBROLD=HZBR

      RETURN
 9001 FORMAT(/' UPSLOI: radii are not monotonous, NTSLOI=',
     *I6,' J=',I3,' DT=',G17.11/5G20.11)
 9020 FORMAT(/' UPSLOI: HZ_new < 0,  TIME,NTSLOI,J,HZ_new = ',
     *1PE12.5,I7,I5,1PE12.5)
 9030 FORMAT(/' UPSLOI, new TTT: DET=0 IN matrix progon; TIME,NTSLOI,J='
     *,1PE12.5,I7,I4)
 9040 FORMAT(/'STOP in NEWEAL: J=',I6,' HIP14=',1PE12.4,
     &'  <  HIP14N=',1PE12.4)  
      END

C**********************************************************************
      SUBROUTINE EBALNC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates various quantities needed for the energy 
c     balance.

c     Called by:  DEIRA
c     Calls    :  
c =====================================================================
      INCLUDE 'de4com.fi'

      CDT=DT*CSURF
      DO 300 I=1,NZ
c: neutron energy deposition in shell I:
      EZN2(I)=EZN2(I)+CDT*WN2(I)
      EZN14(I)=EZN14(I)+CDT*WN14(I)
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 300 J=JI,JF
      RAB=DM(J)
c: driver energy deposition in shell I:
      EZDR(I)=EZDR(I)+(RAB*QDRIV(J))*CDT
      IF(IFBURN.LE.0) GO TO 300
      IF(IFDMFU(J).NE.1) GOTO 300
      RAB1=RAB/(V(J)*AMOL(I)**2)
      RAB2=RAB1*XD(J)
      XQD=XD(J)*QSDD(J)/2.D0
      XQT=XT(J)*QSDT(J)
      XQHE=XHE(J)*QSDHE(J)
      XQBH=XH(J)*XB(J)*QSBH(J)
c: local energy deposition by slow charged fusion products in shell I:
      EZCL(I)=EZCL(I)+(1.76D4*RAB2*XQD)*CDT
c: total fusion energy generated in shell I:
      EZFUS(I)=EZFUS(I)+((1.7D5*XQT+7.04D4*XQD+
     +1.77D5*XQHE)*RAB2+8.38D4*XQBH*RAB1)*CDT
c: energy of fast charged fusion products generated in shell I:
      EZAL(I)=EZAL(I)+((3.4D4*XQT+3.54D4*XQHE)*RAB2+
     +8.38D4*XQBH*RAB1)*CDT
      EZP3(I)=EZP3(I)+(2.91D4*XQD*RAB2)*CDT
      EZP14(I)=EZP14(I)+(14.16D4*XQHE*RAB2)*CDT
c: total number of fusion neutrons generated in shell I:
      TNUN2(I)=TNUN2(I)+(6.022D20*CDT)*RAB2*XQD
      TNUN14(I)=TNUN14(I)+(6.022D20*CDT)*RAB2*XQT
      TNUDHE3(I)=TNUDHE3(I)+(6.022D20*CDT)*RAB2*XQHE
 300  CONTINUE
      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION ELLICK(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the full first elliptic integral.

c     Called by:  KINBUR
c     Calls    :  none
c =====================================================================

      Z=1.D0-X
      IF(Z.GT.1.D0) Z=1.D0
      IF(Z.LT.1.D-12) Z=1.D-12
      ELLICK=1.3862944D0+Z*(.1119723D0+7.25296D-2*Z)+
     +(.5D0+Z*(.1213478D0+2.88729D-2*Z))*DLOG(1.D0/Z)
      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION ZN(N,HS,H1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the denominator of a geometric progression
c     given the number of terms N,the sum HS and the 1-st term H1;
c     sets Q1=1.D0  for the increasing, and Q1=0.D0 for the decreasing
c     progressions.

c     Called by:  START
c     Calls    :  none
c =====================================================================

      RAB=HS/H1
      RABN=N
      IF (RAB.NE.RABN) GOTO 5
      ZN=1.D0
      RETURN
 5    Q1=1.D0
      IF (RAB.LT.RABN) Q1=0.D0
      DQ=0.5D0
      DDQ=0.5D0
 10   IF(DDQ.LT.1.D-11) GO TO 3
      X=Q1+DQ
      IF(N*DLOG(X).GT.4.D1) GO TO 1
      AF=(X**N-1.D0)/(X-1.D0)-RAB
      IF(AF) 2,3,1

 1    DDQ=0.5D0*DDQ
      DQ=DQ-DDQ
      GO TO 10
 2    DDQ=0.5D0*DDQ
      DQ=DQ+DDQ
      GO TO 10
 3    ZN=Q1+DQ
      RETURN
      END

C**********************************************************************
      SUBROUTINE URSOS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine fills in the following EOS arrays:
c     YI(J)=ionization degree, PE(J)=electron pressure
c     (in units 10**14 ergs/cm**3), PI(J)=ion pressure, EE(J)=electron
c     specific energy (in units 10**14 ergs/g), EI(J)=ion specific
c     energy, EEV(J)=d(EE(J))/d(V(J)), EET(J)=d(EE(J))/d(TE(J)),
c     EIV(J)=d(EI(J))/d(V(J)), EIT(J)=d(EI(J))/d(TI(J)), US(J)=sound
c     speed (in units 10**7 cm/s) for given values of TE(J),TI(J),TR(J)
c     (in keV), and V(J)=1/RHO(J) (in units cm**3/g).

c     Called by:  DEIRA,START,PRINUN
c     Calls    :  EOS
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/EOSOUT/YOU,PEOU,PIOU,EEOU,EEVOU,EETOU,EIOU,EIVOU,EITOU,USOU

      DIMENSION YI(1:NN+1),EE(1:NN+1),EI(1:NN+1)
      EQUIVALENCE (QSDT(1),YI(1)),(QE(1),EE(1)),(QI(1),EI(1))

      DO 1200 I=1,NZ
      JI=NJZ(I)
      JF=NJZ(I+1)-1

      DO 1200 J=JI,JF
      VJ=V(J)
      TEJ=TE(J)
      TIJ=TI(J)
      CALL EOS(I,VJ,TEJ,TIJ)
      YI(J)=YOU
      PE(J)=PEOU
      PI(J)=PIOU
      EE(J)=EEOU
      EI(J)=EIOU
      EEV(J)=EEVOU
      EET(J)=EETOU
      EIV(J)=EIVOU
      EIT(J)=EITOU
      RAB=USOU**2+.444444444D0*ASBOL*TR(J)**4*VJ
      IF(IFHZ.EQ.1) RAB=RAB+VJ*HZ(J)**2/(4.D0*PINUM)
      US(J)=DSQRT(RAB)
 1200 CONTINUE
      RETURN
      END

C**********************************************************************
      SUBROUTINE DRIVE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the specific external heating (in addition
c     to a possible radiation drive, set by the boundary radiation 
C     temperatures TREX and TRLEX) of target shells in units 
c     [TW/mg] - array QDRIV(J);

c     Called by:  DEIRA
c     Calls    :  BESTO2
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/BECNST/CBEAM(8)
      COMMON/BEAM1/ABEAM,ZBEAM,CBEAM1(406)
      COMMON/BEAM2/CBEAM2(350),NCBEM2(106),CBEM22(125)

      TVACCD=0.D0
      NDRIV=NDRIV+1
      DO 100 J=1,N
  100 QDRIV(J)=0.D0
      IF(TIME.GE.TDRFIN) RETURN

C: here - alternative drive energy input;
C  the value of TDRCAL may be changed here;

      IF(EIONB.LE.0.D0) GOTO 410
      EPS=1.D-5*ABEAM
      EBK=EIONB
C: forward-backward loop:
      DO 400 ID=1,2
      IF(ID.EQ.2.AND.IFLBND.EQ.-1) GOTO 410
C: loop over target layers:
      DO 400 IZ=1,NZ
      I=(NZ-IZ+1)*(2-ID)+IZ*(ID-1)
C: loop over mesh cells:
      JI=NJZ(I)
      JF=NJZ(I+1)-1
      DO 400 JJ=JI,JF
      J=(JI+JF-JJ)*(2-ID)+JJ*(ID-1)
C: integrate over the cell, RUNGE-KUTTA:
      EBN=EBK
      RAB=0.D0
      ROJ=1.D0/V(J)
      DO 300 K=1,4
      EBR=EBN-RAB
      IF(K.EQ.4) EBR=EBR-RAB
      RAB=ROJ*BESTO2(EBR,ROJ,TE(J),TI(J),NSUB(I))*
     *(R(J+1)-R(J))/2.D0
 300  EBK=EBK-(K*(5-K)-2)*RAB/6.D0
      IF(EBK.LE.EPS) EBK=0.D0
      QDRIV(J)=QDRIV(J)+WDRIV*(EBN-EBK)/(EIONB*CSURF*DM(J))
      IF(EBK.LE.EPS) GOTO 410
 400  CONTINUE
 410  CONTINUE
      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION BESTO2(E1G,DG,TEG,TIG,NSUB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the stopping power in units GeV*mm**2/mg
c     as a function of the ion energy E1G (in GeV), 
c     the density DG (in g/cc), the electron, TEG, and the ion, TIG, 
c     temperatures (in keV),and the substance number  NSUB.

c     Called by:  DRIVE
c     Calls    :  none
c =====================================================================

      COMMON/URSAZ/ATARG(5),ZTARG(5)
      COMMON/BECNST/C11,C12,C13,C14,C15,C16,C17,C1R
      COMMON/BEAM1/A1,Z1,B1,BET1,SGM1,AMU1,PTIF1,DINP1,APIN1(4,100)
      COMMON/BEAM2/POT2(350),NSH2(6),NE2(100),ESH2(100),GB0(5),B2(5),
     *BET2(5),SGM2(5),AMU2(5)
      COMMON/UR2OUT/YAP,PEAP,EEAP,YLAP,PELAP,EELAP,DYLR,DPELR,DEELR
     *,DYLT,DPELT,DEELT,PIAP,EIAP,DPILR,DEILR,DPILT,DEILT
      COMMON/BESOUT/CLBE,SBE,CLFE,SFE,CLFI,SFI,CLNU,SNU,STPE,STPI,
     *Y2,Y300,Z1EF,NIBES

      H(X)=DLOG(1.D0+X/(1.D0 +3.5D0/DSQRT(X)))
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
      A2=ATARG(NSUB)
      Z2=ZTARG(NSUB)
      IF(E1-1.D-8) 999,999,10
 10   R=Z2+.5D0
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
      DCL=2.D0*DLOG(1.D0+C16*E1)-C16*VV1

C: evaluate Z1EF:
      Z1EF=1.D0
      YT=1.D0
      IF(Z1.LT.1.5D0) GOTO 67
      Z1EF=Z1/(1.D0+(.62D0*Z1**C11/DSQRT(VV1))**1.7D0)**.58824D0
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

C: free-electron stoppig power:
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
      R1=DSQRT(R3*R4/(1.D0+C14*(Z1*Z2)**2/VI2))
      CLNU=DLOG(1.D0+R1/(1.D0+.35D0/DSQRT(R1)))
      R2=4.D0*C13*R2*(XI**3/(XI**3+1.33D0))/(A1*VV1)
      R2=R2*DSQRT(1.D0+A2/A1)
      SNU=R2*(Z1*Z2)**2*CLNU
      R1=(1.D0+C14*(Z1EF*Y2)**2/VI2)/R3
      IF(R1.LT.R4) R1=R4
      R1=DSQRT(R/R1)
      CLFI=DLOG(1.D0+R1/(1.D0+.5D0/DSQRT(R1)))
      SFI=R2*(Z1EF*Y2)**2*CLFI
 999  BESTO2=4.589D-2*A1*(SBE+SFE+SFI+SNU)/A2
      RETURN
      END

C**********************************************************************
      SUBROUTINE RELCON(ISHELL,VJJ,TEJJ,TIJJ,TRJJ,HZJJ,YIJJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the values of the transversal (with 
c     respect to the magnetic field) electron heat conduction 
c     coefficient CAPEOU, of the ion, CAPIOU, and the radiation, CAPROU, 
c     heat conduction coefficients [in 10**20 ergs/(cm s keV)], 
c     of the electron-ion, HIEIOU, and the electron-radiation, HIEROU, 
c     temperature relaxation coefficients [in 10**22 ergs/(g s keV)], 
c     and of the transversal electrical resistivity ETHZOU [10**-8 s]; 
c     all the output information is transferred via the c/blk /CAHIOU/; 
c     YIONOU is the actual ionization degree used to calculate 
c     the above mentioned coefficients;
c     the input parameters are: the target layer number ISHELL, the
c     specific volume VJJ [cm**3/g], the electron, TEJJ, the ion, TIJJ, 
c     and the radiation, TRJJ, temperatures [in keV], the magnetic 
c     z-field HZJJ [in 10**7 Gauss], the ionization degree YIJJ.

c     Called by:  KINBUR
c     Calls    :  
c =====================================================================
      INCLUDE 'de4com.fi'
C: only IFHZ,AMOL(I),XMOL(I),ZMOL(I),Z2MOL(I),SMOL(I),FITCAP(I)
C  are used from the above c/blks.

      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/OPAC1/POTOFI(100,5),BOPAC(100,5)
      COMMON/OPAC2/NLAGER,NWWWW,XLAGER(12),WLAGER(12),RWEI(12),PWEI(12)
      COMMON/CAHIOU/YIONOU,CAPEOU,CAPIOU,HIEIOU,ETHZOU
     *,CAPROU,HIEROU,ETV0OU,ETV1OU
      FI1(X,B)=X/(210.D0*X+5.D0*B+.25D0*(B+3.D0)/X)
      FI2(X,B)=X/(X+.5D0*(B+2.5D0)*(1.D0+1.D0/(12.D0*X)))
      GA11(W)=3.25D0+2.D0*W
      GA12(W)=(58.725625D0+W*(375.74D0+W*(427.68D0+129.6D0*W)))/49.D0
      GA71(W)=(366.625625D0+W*(412.69D0+132.52D0*W))/49.D0
      GA72(W)=0.31D0+W*(12.08D0+5.76D0*W)/7.D0
      GA51(W)=(226.895625D0+123.705D0*W)/49.D0
      GA52(W)=(3.3201D0+W*(34.1064D0+W*(95.7888D0+W*41.472D0)))/49.D0

      NSUBI=NSUB(ISHELL)
      Z2=ZSUB(NSUBI)
      A2M=AMOL(ISHELL)/XMOL(ISHELL)
      RAB=Z2+.5
      IZ2=RAB

C: common quantities:
      YION=YIJJ
      IY=YION+1.D0
      IF(IY.LT.1) IY=1
      IF(IY.GT.IZ2) IY=IZ2
      YIONOU=YION
      VJ=VJJ
      TEJ=TEJJ
      TIJ=TIJJ
      TRJ=TRJJ
      HZJ=HZJJ
      TE32=TEJ*DSQRT(TEJ)
      YAE=(YION/Z2)*ZMOL(ISHELL)/AMOL(ISHELL)
      Y2AE=(YION/Z2)**2*Z2MOL(ISHELL)/AMOL(ISHELL)
      RYA=YAE/VJ
      RY2A=Y2AE/VJ
      EFE=.026D0*RYA**.666666667D0
      TFE=DSQRT(TEJ**2+(.666666667D0*EFE)**2)
      TET=TEJ/EFE
      TET32=TET*DSQRT(TET)

C-------------------------------------------------------------- begin 2
C: coefficient of radiative heat conduction CAPROU,
C  coefficient of electron-radiation temperature relaxation HIEROU:
C.......................................................................
C  IFOPAC<2 - default (zero) values of CAPROU and HIEROU;              .
C  IFOPAC=2 - calculate CAPROU and HIEROU from fast formulas of DEIRA-2.
C  IFOPAC=3 - calculate CAPROU and HIEROU by integrating absorption    .
C             cross-section (model DEIRA-3);                           .
C  IFOPAC>3 - default (zero) values of CAPROU and HIEROU;              .
C.......................................................................
      IF(IFOPAC.EQ.2) GOTO 300
      IF(IFOPAC.EQ.3) GOTO 1000
 200  CONTINUE
      CAPROU=0.D0
      HIEROU=0.D0
      GOTO 3000
C------------------------------------------------------------- begin 21
C: fast formulas from DEIRA-2:
 300  CONTINUE
C..Compton:
      CCS=.04D0*TRJ*RYA
      HICS=.94D0*YAE*ASBOL*TRJ**4
C..free-free:
      CFF=RY2A/(TRJ+9.3D0*TET32*TRJ*(TRJ/TEJ))
      HIFF=1.12D5*Y2AE*(TEJ+TRJ)/(1.D0+2.19D0*TET32)
C..free-bound:
      TRT=TEJ/TRJ
      IF(DABS(TRT-1.D0).GT.1.D-4) GOTO 500
      TRJ=TEJ*.9999D0
      TRT=TEJ/TRJ
 500  DY=Z2-YION
      IF(DY.GT.1.D0) GOTO 700
C_hydrogen-like ions:
      PTY=.0136D0*Z2*Z2
      RAB=(EFE-PTY)/TEJ
      DY=0.D0
      IF(RAB.LT.-42.D0) DY=1.D0
      IF(DABS(RAB).LT.42.D0)
     *DY=1.D0/(1.D0+.665D0*TET32*DEXP(RAB))
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      CPH=DY*FI1(OMR,-1.397D0)
      HIPH=DY*((TRT*FI2(OME,-1.397D0)-
     *FI2(OMR,-1.397D0))/(TRT-1.D0))
      GOTO 990
 700  IF(DY.GT.2.D0) GOTO 800
C_helium-like ions:
      PTY=.0136D0*Z2*Z2
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      RAB=2.D0-DY
      CPH=RAB*FI1(OMR,-1.397D0)
      HIPH=RAB*(TRT*FI2(OME,-1.397D0)-FI2(OMR,-1.397D0))
      PTY=.0136D0*(Z2-.65D0)**2
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      BY=-1.397D0+13.25D0/(Z2-.4D0)**2
      RAB=2.D0*(DY-1.D0)
      CPH=CPH+RAB*FI1(OMR,BY)
      HIPH=(HIPH+RAB*(TRT*FI2(OME,BY)-FI2(OMR,BY)))
     */(TRT-1.D0)
      GOTO 990
 800  IF(DY.GT.3.D0) GOTO 900
C_lithium-like ions:
      PTY=.0136D0*(Z2-.65D0)**2
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      BY=-1.397D0+13.25D0/(Z2-.4D0)**2
      RAB=2.D0*(3.D0-DY)
      CPH=RAB*FI1(OMR,BY)
      HIPH=RAB*(TRT*FI2(OME,BY)-FI2(OMR,BY))
      PTY=.0036D0*(Z2-1.77D0)**2
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      RAB=(1.D0-1.77D0/Z2)**4
      BY=5.53D0*DLOG(11.D0/RAB)/RAB
      RAB=2.D0*(DY-2.D0)
      CPH=CPH+RAB*FI1(OMR,BY)
      HIPH=(HIPH+RAB*(TRT*FI2(OME,BY)-FI2(OMR,BY)))
     */(TRT-1.D0)
      GOTO 990
C_highly ionized ions:
 900  LY=YION
      YRJ=LY
      PTY=.0036D0*(YRJ+1.3D0)**2
      IF(Z2-YRJ.GT.10.001D0) PTY=.0063D0*(1.D0+YRJ)*DSQRT(1.D0+YRJ)+
     +1.5D-7*(1.D0+YRJ)**4
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      RAB=(.0272D0*Z2*Z2/PTY)**2/(Z2-YRJ)
      BY=.245D0*RAB*DLOG(RAB)
      RAB1=YRJ+1.D0-YION
      CPH=2.D0*FI1(OMR,BY)*RAB1
      HIPH=2.D0*(TRT*FI2(OME,BY)-FI2(OMR,BY))*RAB1
      YRJ=LY+1
      PTY=.0036D0*(YRJ+1.3D0)**2
      IF(Z2-YRJ.GT.10.001D0) PTY=.0063D0*(1.D0+YRJ)*DSQRT(1.D0+YRJ)+
     +1.5D-7*(1.D0+YRJ)**4
      OME=(TEJ/PTY)**2
      OMR=(TRJ/PTY)**2
      RAB=(.0272D0*Z2*Z2/PTY)**2/(Z2-YRJ)
      BY=.245D0*RAB*DLOG(RAB)
      RAB1=1.D0-RAB1
      CPH=2.D0*FI1(OMR,BY)*RAB1+CPH
      HIPH=(2.D0*(TRT*FI2(OME,BY)-FI2(OMR,BY))*RAB1+HIPH)/(TRT-1.D0)
 990  RAB=Z2**4
      CPH=CPH*1.2D0*RAB/(VJ*A2M*TRJ**2)
      HIPH=HIPH*759.D0*RAB/A2M
      CAPROU=5484.D0*TRJ**2*(TRJ**2/(CCS+CFF+CPH))
      HIEROU=HICS+HIFF+HIPH
      GOTO 3000
C=============================================================== end 21

C------------------------------------------------------------- begin 22
C: model DEIRA-3: absorption cross-section is numerically integrated;
 1000  CONTINUE
      EFETE=EFE/TEJ
      EFTE32=EFETE*DSQRT(EFETE)
      IF(EFETE.LT.2.D1) GOTO 1020
      EXMUTE=1.D30
      AMUTE=EFETE-2.6587D0/EFTE32
      AMU=AMUTE*TEJ
      GOTO 1060
 1020 RAB=1.D0+2.6587D0/EFTE32
      IF(EFETE.GT.1.D-4) GOTO 1040
      EXMUTE=2.D0/RAB
      GOTO 1050
 1040 EXMUTE=(1.D0+DEXP(EFETE))/RAB
 1050 AMUTE=DLOG(EXMUTE)
      AMU=AMUTE*TEJ
 1060 CONTINUE

 1070 DIY=IY-YION
      IF(IY.NE.IZ2) GOTO 1090
C: a special formula for DIY for hydrogen-like ions:
      RAB=(EFE-0.0136D0*Z2**2)/TEJ
      IF(RAB.LT.5.D1) GOTO 1080
      DIY=0.D0
      GOTO 1085
 1080 RAB=158.54D0*TE32*A2M*VJ*DEXP(RAB)
      DIY=2.D0*Z2/(1.D0+Z2+RAB+DSQRT((Z2-1.D0)**2+
     +RAB*(RAB+2.D0*(1.D0+Z2))))
 1085 YIONOU=Z2-DIY
 1090 AVY2=DIY*(IY-1)**2+(1.D0-DIY)*IY**2

      IF(DABS(1.D0-TRJ/TEJ).LT.1.D-2) TRJ=.99D0*TEJ
C: calculate necessary integrals:
      CPHOT0=4.1174D0*4.8368D0*Z2**4
      SUMR=0.D0
      SUMPE=0.D0
      SUMPR=0.D0

      DO 1970 L=1,NLAGER
      EFO=XLAGER(L)*TRJ
      IFLAG=0
C-------------------------------------------------------- begin SIGMAA
 1100 CONTINUE
C: calculate SIGMAA(EFO):
C: photo-absorption:
      SIGPH=0.D0
      IION=IY
      DFRAC=DIY
 1110 RAB=EFO/POTOFI(IION,NSUBI)
      CPHOT=CPHOT0
      IF(IION.NE.IZ2) CPHOT=CPHOT0*2.D0
      IF(RAB.GT.1.D0) GOTO 1150
      SIGPH=SIGPH+DFRAC*(CPHOT/(POTOFI(IION,NSUBI)**3*
     *BOPAC(IION,NSUBI)))*RAB**2
      GOTO 1200
 1150 SIGPH=SIGPH+DFRAC*(CPHOT/POTOFI(IION,NSUBI)**3)/
     *(RAB*(RAB**2+BOPAC(IION,NSUBI)-1.D0))
 1200 IF(IION.EQ.IY+1) GOTO 1300
      IF(IY.EQ.IZ2) GOTO 1300
      IION=IY+1
      DFRAC=1.D0-DIY
      GOTO 1110
 1300 CONTINUE
C%      SIGPH=0.D0

C: free-free absorption:
      EFOTE=EFO/TEJ
C.. the Gaunt-factor:
C%      GFF=1.D0
C%      GOTO 1580
      IF(EFO.GE.2.D0*(AMU-.71D0*TEJ)) GOTO 1550
      RAB=2.D0*AMU/EFO
      GFF=.55133D0*DLOG(RAB*(1.D0+DSQRT(1.D0-1.D0/RAB**2)))
      GOTO 1580
 1550 RAB=.71D0*TEJ/EFO
      GFF=.55133D0*DLOG(1.D0+2.D0*RAB*(1.D0+DSQRT(1.D0+1.D0/RAB)))
 1580 CONTINUE

C: the logarithm term in cross-section:
      AMUEFO=AMUTE-EFOTE
      IF(AMUTE.LT.4.6D0) GOTO 1660
C..AMUTE > 4.6:
      IF(AMUEFO.LT.4.6D0) GOTO 1620
      RAB=EFOTE
      GOTO 1800
 1620 IF(AMUEFO.LT.-4.6D0) GOTO 1640
      RAB=AMUTE-DLOG(1.D0+DEXP(AMUEFO))
      GOTO 1800
 1640 RAB=AMUTE
      GOTO 1800

 1660 IF(AMUTE.LT.-4.6D0) GOTO 1720
      RAB=DLOG((1.D0+EXMUTE)/(1.D0+DEXP(AMUEFO)))
      GOTO 1800

C..AMUTE < -4.6:
 1720 IF(EFOTE.LT.4.6D0) GOTO 1740
      RAB=EXMUTE
      GOTO 1800
 1740 IF(EFOTE.LT.-4.6D0) GOTO 1760
      RAB=EXMUTE*(1.D0-DEXP(-EFOTE))
      GOTO 1800
 1760 RAB=EXMUTE*EFOTE
 1800 SIGFF=1463.7D0*(RAB*TEJ/EFO**3)*AVY2*GFF
C%      SIGFF=0.D0
      SIGMAA=SIGPH+SIGFF*Z2MOL(ISHELL)/(Z2*Z2*XMOL(ISHELL))
C: SIGMAA(EFO) calculated.
C========================================================== end SIGMAA
      IF(IFLAG.NE.0) GOTO 1900
C: add Thomson scattering: ?
C%      SIGSCA=0.D0
      SIGSCA=0.6652D0*YION
      SUMR=SUMR+RWEI(L)/(SIGMAA+SIGSCA*ZMOL(ISHELL)/(Z2*XMOL(ISHELL)))
      SUMPR=SUMPR+PWEI(L)*SIGMAA
      EFO=XLAGER(L)*TEJ
      IFLAG=1
      GOTO 1100
 1900 SUMPE=SUMPE+PWEI(L)*SIGMAA
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ begin TEST
C%      PRINT 9201,XLAGER(L),RWEI(L),PWEI(L),SIGPH,SIGFF
C% 9201 FORMAT(1PE12.5,' RW,PW=',1P2E12.5,' PH,FF=',1P2E12.5)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ end TEST
 1970 CONTINUE

      RAB=16.6D0*A2M
      RAB1=RAB*VJ
      ROSSL=RAB1*SUMR
C      PLANL=RAB1/SUMPE
      CAPROU=5484.D0*(ROSSL*TRJJ)*TRJJ**2
      HIEROU=1.29D0*TRJJ**4*YION*ZMOL(ISHELL)/(Z2*AMOL(ISHELL))+
     +(4113.D0/RAB)*(TEJ**4*SUMPE-TRJ**4*SUMPR)/(TEJ-TRJ)
      GOTO 3000

C......................................................................
C  If the Rosseland, ROSSL(V,T), and the Planck, PLANL(V,T), mean     .
C  free paths (in units of mm) are provided from some other source,   .
C  CAPROU and HIEROU can be evaluated as:                             .
C  CAPROU=5484.D0*ROSSL(VJ,TRJ)*TRJ**3                                .
C  HIEROU=4113.D0*VJ*TEJ**3/PLANL(VJ,TEJ)                             .
C......................................................................
C=============================================================== end 22
 3000 CONTINUE
C================================================================ end 2

C-------------------------------------------------------------- begin 3
C: coefficients of the electron CAPEOU and ion CAPIOU heat conduction,
C  coefficient of the electron-ion temperature relaxation HIEIOU,
C  electrical resistivity ETHZOU,
C  coefficients of the ion viscosity ETV0OU and ETV1OU;
C: calculate Coulomb logarithms:
      YSTAR=YION
      IF(YION.LT.1.D0) YSTAR=1.D0
      FYSTAR=1.D0+YION-YSTAR
      DEE2=RYA/TFE
      RAB=315.5D0*TEJ/DSQRT(DEE2*(1.D0+27.56D0*TFE))
      RAB=RAB*DSQRT(RAB)
      CLEE=RAB/(1.D0+1.D0/(8.5D0*RAB))
      IF(CLEE.GT.1.D-4) CLEE=DLOG(1.D0+CLEE)
      CLEE=.6666666667D0*CLEE
      DEI2=DEE2+(YSTAR/Z2)**2*FYSTAR*Z2MOL(ISHELL)/(AMOL(ISHELL)*VJ*TIJ)
      RAB=631.D0*TFE/DSQRT(DEI2*(YSTAR**2+27.56D0*TFE))
      CLEI=FITCAP(NSUBI)*RAB/(1.D0+1.D0/(FITCAP(NSUBI)*6.5D0*RAB))
      IF(CLEI.GT.1.D-4) CLEI=DLOG(1.D0+CLEI)
      RAB=315.5D0*TIJ/DSQRT(DEI2*(YSTAR**4+1.512D-2*TIJ/A2M))
      CLII=RAB*RAB
      IF(CLII.GT.1.D-4) CLII=DLOG(1.D0+CLII)
      CLII=.5D0*CLII
C: e-a collision frequency:
      COLFEA=0.D0
      IF(YION.LT.1.D0) COLFEA=7.987D9*DSQRT(TFE)*
     *DEXP(-(EFE/1.3D-2)**2)*(1.D0-YION)/(VJ*A2M)
C: electron heat conduction:
      BCOLF=3.394D0
      COLFEE=3.914D5*RYA*CLEE/(TEJ*DSQRT(DSQRT(TEJ**2+
     +(BCOLF*EFE)**2)))
      BCOLF=0.3402D0
      RAB=DSQRT(TEJ**2+(BCOLF*EFE)**2)
      COLFEN=COLFEA+5.536D5*CLEI*YSTAR**2*FYSTAR*Z2MOL(ISHELL)/
     *(RAB*DSQRT(RAB)*Z2*Z2*AMOL(ISHELL)*VJ)
      WW=COLFEE/COLFEN
      XHZ2=0.D0
      IF(IFHZ.EQ.1) XHZ2=(1.759D6*HZJ/COLFEN)**2
      IF(XHZ2.GT.1.D4) GOTO 3240
      RAB=(XHZ2*GA11(WW)+GA12(WW))/
     *(XHZ2*(XHZ2+GA71(WW))+(GA72(WW))**2)
      GOTO 3250
 3240 RAB=(GA11(WW)+GA12(WW)/XHZ2)/
     *(XHZ2+GA71(WW)+(GA72(WW))**2/XHZ2)
 3250 CAPEOU=1.697D5*RYA*TEJ*RAB/COLFEN
C: ion heat conduction and viscosity:
      COLFIA=0.D0
      IF(YION.LT.1.D0) COLFIA=1.871D7*DSQRT(TIJ/A2M)*
     *(1.D0-YION)/(VJ*A2M)
      RAB=Z2MOL(ISHELL)/(VJ*Z2*Z2*AMOL(ISHELL)*DSQRT(A2M))
      COLFIN=9166.D0*CLII*YSTAR**4*FYSTAR*RAB/(TIJ*DSQRT(TIJ))+COLFIA
      XHZ2=0.D0
      IF(IFHZ.EQ.1) XHZ2=(964.9D0*HZJ*YSTAR/(A2M*COLFIN))**2
      IF(XHZ2.GT.1.D4) GOTO 3340
      RAB=(2.D0*XHZ2+2.645D0)/(XHZ2*(XHZ2+2.7D0)+.677D0)
      RAB1=(4.8D0*XHZ2+2.23D0)/(XHZ2*(16.D0*XHZ2+16.12D0)+2.33D0)
      GOTO 3350
 3340 RAB=(2.D0+2.645D0/XHZ2)/(XHZ2+2.7D0+.677D0/XHZ2)
      RAB1=(4.8D0+2.23D0/XHZ2)/(16.D0*XHZ2+16.12D0+2.33D0/XHZ2)
 3350 CAPIOU=93.89D0*TIJ*(RAB/COLFIN)*
     *Z2*Z2*SMOL(ISHELL)/(VJ*AMOL(ISHELL)*DSQRT(A2M))
      ETV0OU=9.263D0*TIJ/(VJ*A2M*COLFIN)
      ETV1OU=9.649D0*TIJ*RAB1/(VJ*A2M*COLFIN)
C: electron-ion temperature relaxation:
      BCOLF=0.8271D0
      RAB=DSQRT(TEJ**2+(BCOLF*EFE)**2)
      COLFEN=COLFEA+5.536D5*CLEI*YSTAR**2*FYSTAR*Z2MOL(ISHELL)/
     *(RAB*DSQRT(RAB)*Z2*Z2*AMOL(ISHELL)*VJ)
      HIEIOU=1.588D-2*YAE*COLFEN/A2M
C: electrical resistivity:
      BCOLF=0.3665D0
      RAB=DSQRT(TEJ**2+(BCOLF*EFE)**2)
      COLFEN=COLFEA+5.536D5*CLEI*YSTAR**2*FYSTAR*Z2MOL(ISHELL)/
     *(RAB*DSQRT(RAB)*Z2*Z2*AMOL(ISHELL)*VJ)
      WW=COLFEE/COLFEN
      XHZ2=0.D0
      IF(IFHZ.EQ.1) XHZ2=(1.759D6*HZJ/COLFEN)**2
      IF(XHZ2.GT.1.D4) GOTO 3440
      RAB=(XHZ2*GA51(WW)+GA52(WW))/
     *(XHZ2*(XHZ2+GA71(WW))+(GA72(WW))**2)
      GOTO 3450
 3440 RAB=(GA51(WW)+GA52(WW)/XHZ2)/
     *(XHZ2+GA71(WW)+(GA72(WW))**2/XHZ2)
 3450 ETHZOU=(6.556D-17*COLFEN)*(1.D0-RAB)/RYA
C================================================================ end 3
      RETURN
      END

c**********************************************************************
      SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ)
c      SUBROUTINE EOSTAB(ISHELL,VJ,TEJ,TIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c   >>> A version of the SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ) <<<
c ---------------------------------------------------------------------
c     This routine calculates the equation of state (with possible 
c     contribution from radiation excluded) in the target layer ISHELL
c     as a function of VJ=specific volume (in units cm**3/g), and of the
c     electron, TEJ, and ion, TIJ, temperatures (in keV).
c
c     The output results are in c/blk /EOSOUT/: YOU=ionization degree,
c     PEOU= (PIOU=) electron (ion) pressure (in 10**14 erg/cm**3);
c     EEOU= (EIOU=) electron (ion) specific energy (in 10**14 erg/g);
c     EEVOU=d(EE)/dV, EIVOU=d(EI)/dV, EETOU=d(EE)/dTE, EITOU=d(EI)/dTI;
c     USOU=sound speed (in 10**7 cm/s) due to matter pressure only.

c     Called by:  URSOS
c     Calls    :  URSAPB
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/UR2OUT/YAP,PEAP,EEAP,FLAP(3),DYLR,DPELR,DEELR
     *,DYLT,DPELT,DEELT,PIAP,EIAP,DPILR,DEILR,DPILT,DEILT
      COMMON/EOSOUT/YOU,PEOU,PIOU,EEOU,EEVOU,EETOU,EIOU,EIVOU,EITOU,USOU

      NSUBI=NSUB(ISHELL)
C: atomic units for specific volume and energy:
      AEOS=ASUB(NSUBI)
      IF(IFMFU0(ISHELL).NE.1) GOTO 200
      IF(XHE0(ISHELL).GT.XBFLOO) GOTO 200
      IF(XB0(ISHELL).GT.XBFLOO) GOTO 200
C: adjust the value of A for a mixture of hydrogen isotopes:
      AEOS=AMOL(ISHELL)/XMOL(ISHELL)
 200  VU=1.D0/(11.206D0*AEOS)
      EU=2.942D0*VU

      ROA=VU/VJ
      TEA=TEJ/.02721D0
      TIA=TIJ/.02721D0
      CALL URSAPB(NSUBI,ROA,TEA,TIA,3)
      IF(YAP.GT.ZSUB(NSUBI)) YAP=ZSUB(NSUBI)
      YOU=YAP
      PEOU=PEAP*2.942D0
      PIOU=PIAP*2.942D0
      EEOU=EEAP*EU
      EIOU=EIAP*EU
      EEVOU=-DEELR*2.942D0*ROA
      EIVOU=-DEILR*2.942D0*ROA
      EITOU=DEILT*EU/TIJ
      RAB=DEELT*EU/TEJ
C: a lower limit on the electron heat capacity:
      RAB1=EEMIN/TEJ
      IF(RAB1.GT.VU) RAB1=VU
      IF(RAB.LT.RAB1) RAB=RAB1
      EETOU=RAB
C: sound speed due to matter pressure only:
      RAB=VJ*(DPELR+DPILR+(EU*DPELT**2/(TEJ*EETOU)+
     +DPILT**2/DEILT)/ROA)*2.942D0
      IF(RAB.LT.0.D0) RAB=0.D0
      USOU=DSQRT(RAB)

      RETURN
      END

C**********************************************************************
      SUBROUTINE URSAPB(NSUBST,D,TE,TI,JOB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the two-temperature equation of state 
c     from tables (in a.u.);
c     NSUBST is the substance number (according to c/blk /URSAZ/);
c     output data are stored in COMMON/UR2OUT/; 
c     depending on the JOB value, they are:
c  JOB=1: YI=FAPP(1) (ionization degree) and its derivatives;
c      2: YI,PE=FAPP(2),PI,EI and their derivatives;
c      3: YI,PE,EE=FAPP(3),PI,EI and their derivatives (complete set);
c      4: PE,PI,EI      and their derivatives;
c      5: EE,PI,EI      and their derivatives;
c      6: PE,EE,PI,EI   and their derivatives;
c    > 6: PI,EI         and their derivatives;
c
c     FLAPP(1)=ln(YI); FLAPP(2)=ln(PE+P00); FLAPP(3)=ln(EE+E00);
c     DLFAR(I)=d[FAPP(I)]/d[ln(rho)]; DLFAT(I)=d[FAPP(I)]/d[ln(TE)];

c     Called by:  EOS
c     Calls    :  none
c =====================================================================
      REAL FARR
      COMMON/URTBLI/QURS(1000)
      COMMON/URSTBL/FARR(30000)
      COMMON/URSEDG/NROO(6),ROILL(5),HLROO(5),NTEMM(3,6)
     *,TEMILL(3,5),HLTEMM(3,5)
      COMMON/URSAZ/ASBS(5),ZSBS(5)
      COMMON/EPZER/P00(5),E00(5)
      COMMON/UR2OUT/FAPP(3),FLAPP(3),DLFAR(3),DLFAT(3)
     *,PIAP,EIAP,DPILR,DEILR,DPILT,DEILT

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
      SUBROUTINE EOSIG(ISHELL,VJ,TEJ,TIJ)
c      SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c   >>> A version of the SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ) <<<
c ---------------------------------------------------------------------
c     This routine calculates the equation of state of the ideal 
c     Boltzmann gas P=CEOS*RO*T, E=CEOS*T/(GEOS-1) with the adiabatic
c     index GEOS; pressure and energy are divided equally between 
c     electrons and ions; 
c     ISHELL is the sequential number of the target layer;
c     VJ is the specific volume (in units cm**3/g);
c     TEJ (TIJ) is the electron (ion) temperature (in keV).
c
c     The output results are in c/blk /EOSOUT/: YOU=ionization degree,
c     PEOU= (PIOU=) electron (ion) pressure (in 10**14 erg/cm**3);
c     EEOU= (EIOU=) electron (ion) specific energy (in 10**14 erg/g);
c     EEVOU=d(EE)/dV, EIVOU=d(EI)/dV, EETOU=d(EE)/dTE, EITOU=d(EI)/dTI;
c     USOU=sound speed (in 10**7 cm/s) due to matter pressure only.

c     Called by:  URSOS
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/EOSOUT/YOU,PEOU,PIOU,EEOU,EEVOU,EETOU,EIOU,EIVOU,EITOU,USOU

      DIMENSION CEOS(10),GEOS(10)
      DATA CEOS/10*1.D0/,GEOS/10*1.6666666666667D0/

      NSUBI=NSUB(ISHELL)

      YOU=1.D0
      PEOU=.5D0*CEOS(NSUBI)*TEJ/VJ
      PIOU=.5D0*CEOS(NSUBI)*TIJ/VJ
      EEOU=.5D0*CEOS(NSUBI)*TEJ/(GEOS(NSUBI)-1.D0)
      EIOU=.5D0*CEOS(NSUBI)*TIJ/(GEOS(NSUBI)-1.D0)
      EEVOU=0.D0
      EIVOU=0.D0
      EETOU=EEOU/TEJ
      EITOU=EIOU/TIJ
C: sound speed due to matter pressure only:
      RAB=VJ*(PEOU+PIOU)*GEOS(NSUBI)
      IF(RAB.LT.0.D0) RAB=0.D0
      USOU=DSQRT(RAB)

      RETURN
      END

C**********************************************************************
      SUBROUTINE EOSLIN(ISHELL,VJ,TEJ,TIJ)
c      SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c   >>> A version of the SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ) <<<
c ---------------------------------------------------------------------
c     This routine calculates the so called "linear equation of state":
c           P=P_c(RO)+GGEOS*RO*CVEOS*T;
c           E=E_c(RO)+CVEOS*T;
c           P_c=RO0EOS*C0EOS**2*[(RO/RO0EOS)**(GGEOS+1)-1]/(GGEOS+1);
c           E_c=C0EOS**2*{(RO/RO0EOS)**GGEOS/[GGEOS*(GGEOS+1)]+
c              +RO0EOS/[(GGEOS+1)*RO]-1/GGEOS};
c           C**2=[P(GGEOS+1)+RO0EOS*C0EOS**2]/RO;

c     Here parameters are:
c     RO0EOS   - normal density; must be >0;
c     C0EOS    - sound speed at normal density and T=0; if zero, one has
c                an ideal Boltzmann gas with \gamma= 1+GGEOS;    
c     GGEOS    - Grunaizen coefficient;
c     CVEOS    - heat capacity;
c     pressure and energy are equally divided between electrons and ions;
c
c     The output results are in c/blk /EOSOUT/: YOU=ionization degree,
c     PEOU= (PIOU=) electron (ion) pressure (in 10**14 erg/cm**3);
c     EEOU= (EIOU=) electron (ion) specific energy (in 10**14 erg/g);
c     EEVOU=d(EE)/dV, EIVOU=d(EI)/dV, EETOU=d(EE)/dTE, EITOU=d(EI)/dTI;
c     USOU=sound speed (in 10**7 cm/s) due to matter pressure only.

c     Called by:  URSOS
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/EOSOUT/YOU,PEOU,PIOU,EEOU,EEVOU,EETOU,EIOU,EIVOU,EITOU,USOU

      DIMENSION RO0EOS(10),C0EOS(10),GGEOS(10),CVEOS(10)
      DATA RO0EOS/1.D0, 2*10.D0, 7*0.D0/
     &,C0EOS/0.D0, 2*2.5D-2, 7*0.D0/
     &,GGEOS/.6666666666667D0, 2*2.5D0, 7*.6666666666667D0/
     &,CVEOS/14.5D0, 2*.1D0, 7*1.D0/
c  Substance #1 -> D_2;  #2,3 -> Pb;

      NSUBI=NSUB(ISHELL)

      YOU=1.D0
      RAB=1.D0/(VJ*RO0EOS(NSUBI))
      RAB1=GGEOS(NSUBI)+1.D0
      RAB2=C0EOS(NSUBI)**2
      EETOU=.5D0*CVEOS(NSUBI)
      EITOU=.5D0*CVEOS(NSUBI)

      PCRAB=RO0EOS(NSUBI)*RAB2*(RAB**RAB1-1.D0)/RAB1     
      PEOU=.5D0*PCRAB+GGEOS(NSUBI)*EETOU*TEJ/VJ
      PIOU=.5D0*PCRAB+GGEOS(NSUBI)*EITOU*TIJ/VJ

      EEVOU=-.5D0*PCRAB
      EIVOU=-.5D0*PCRAB
      ECRAB=RAB2*(1.D0/(RAB1*RAB)+(RAB**GGEOS(NSUBI)/RAB1-1.D0)
     &/GGEOS(NSUBI))
      EEOU=.5D0*ECRAB+EETOU*TEJ
      EIOU=.5D0*ECRAB+EITOU*TIJ
C: sound speed due to matter pressure only:
      RAB=VJ*((PEOU+PIOU)*RAB1+RO0EOS(NSUBI)*RAB2)
      IF(RAB.LT.0.D0) RAB=0.D0
      USOU=DSQRT(RAB)

      RETURN
      END

c**********************************************************************
      SUBROUTINE EOSFER(ISHELL,VJ,TEJ,TIJ)
c      SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c   >>> A version of the SUBROUTINE EOS(ISHELL,VJ,TEJ,TIJ) <<<
c ---------------------------------------------------------------------
c     This routine calculates the equation of state of a mixture of
c     the Boltzmann gas of ions and the Fermi gas of electrons (with 
c     possible contribution from radiation excluded) in the target 
c     layer ISHELL as a function of VJ=specific volume (in cm**3/g), 
c     and of the electron, TEJ, and ion, TIJ, temperatures (in keV).
c
c     The output results are in c/blk /EOSOUT/: YOU=ionization degree,
c     PEOU= (PIOU=) electron (ion) pressure (in 10**14 erg/cm**3);
c     EEOU= (EIOU=) electron (ion) specific energy (in 10**14 erg/g);
c     EEVOU=d(EE)/dV, EIVOU=d(EI)/dV, EETOU=d(EE)/dTE, EITOU=d(EI)/dTI;
c     USOU=sound speed (in 10**7 cm/s) due to matter pressure only.

c     Called by:  URSOS
c     Calls    :  FEOS
c =====================================================================
      INCLUDE 'de4com.fi'

      COMMON/URSAZ/ASUB(5),ZSUB(5)
      COMMON/FEOSOU/FEFOU,AMUFOU,PEFOU,EEFOU,DPELTF,DEELTF
     *,DPELVF,DEELVF,PIFOU,EIFOU,DPILTF,DEILTF,DPILVF,DEILVF
      COMMON/EOSOUT/YOU,PEOU,PIOU,EEOU,EEVOU,EETOU,EIOU,EIVOU,EITOU,USOU

      NSUBI=NSUB(ISHELL)
c   here we assume a fully ionized plasma:
      YOU=ZSUB(NSUBI)

c   AEOS is the mean atomic mass (per one nucleus of matter):
      AEOS=AMOL(ISHELL)/XMOL(ISHELL)
c   YEOS is the mean number of electrons per one nucleus of matter:
      YEOS=(YOU/ZSUB(NSUBI))*(ZMOL(ISHELL)/XMOL(ISHELL))
c   atomic units for specific volume and energy:
      VU=1.D0/(11.206D0*AEOS)
      EU=2.942D0*VU

      VA=VJ/VU
      TEA=TEJ/.02721D0
      TIA=TIJ/.02721D0
      
      CALL FEOS(YEOS,VA,TEA,TIA)
      PEOU=PEFOU*2.942D0
      PIOU=PIFOU*2.942D0
      EEOU=EEFOU*EU
      EIOU=EIFOU*EU
      EEVOU=DEELVF*2.942D0/VA
      EIVOU=DEILVF*2.942D0/VA
      EITOU=DEILTF*EU/TIJ
      RAB=DEELTF*EU/TEJ
C: a lower limit on the electron heat capacity:
      RAB1=EEMIN/TEJ
      IF(RAB1.GT.VU) RAB1=VU
      IF(RAB.LT.RAB1) RAB=RAB1
      EETOU=RAB
C: sound speed due to matter pressure only:
      RAB=VJ*(-DPELVF-DPILVF+(EU*DPELTF**2/(TEJ*EETOU)+
     +DPILTF**2/DEILTF)*VA)*2.942D0
      IF(RAB.LT.0.D0) RAB=0.D0
      USOU=DSQRT(RAB)

      RETURN
      END

C**********************************************************************
      SUBROUTINE FEOS(YA,VA,TEA,TIA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the equation of state of a mixture of
c     the Boltzmann gas of ions and the Fermi gas of electrons (with 
c     possible contribution from radiation excluded) according to the
c     interpolation of Fermi integrals by Basko; all the input and 
c     output quantities are in atomic units.
c
c     The input parameters are:
c     YA - ionization degree = number of free electrons per atom,
c     VA - specific volume,  TEA and TIA - electron and ion temperatures.
c     The output parameters are in the c/blk /EOSFOU/:
c     FEFOU  - electron free energy;
c     AMUFOU - electron chemical potential;
c     PEFOU  - electron pressure;
c     EEFOU  - electron energy;
c     DPELTF - dP_e/dlnT_e;   DPELVF - dP_e/dlnV;
c     DEELTF - dE_e/dlnT_e;   DEELVF - dE_e/dlnV;
c     PIFOU  - electron pressure;
c     EIFOU  - electron energy;
c     DPILTF - dP_i/dlnT_i;   DPILVF - dP_i/dlnV;
c     DEILTF - dE_i/dlnT_e;   DEILVF - dE_i/dlnV;
c: Atomic units:
c     density: 11.206*A g/cc;       temperature: e^2/a_0 = 27.21 eV;
c     pressure:        e^2/a_0^4 = 294.2 Mbar;
c     internal energy: 27.21 eV/atom = (26.26/A)*10^12 erg/g.

c     Called by:  EOSFER
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4par.fi'

      COMMON/FEOSOU/FEFOU,AMUFOU,PEFOU,EEFOU,DPELTF,DEELTF
     *,DPELVF,DEELVF,PIFOU,EIFOU,DPILTF,DEILTF,DPILVF,DEILVF

      DATA TFLOOR/1.D-14/
      SAVE TFLOOR
      
c: initial assignments:
      Y=YA
      V=VA
      TE=TEA
      TI=TIA
      IF(Y.GT.0.D0.AND.V.GT.0.D0.AND.TE.GE.0.D0.AND.TI.GE.0.D0)GOTO 100
      PRINT 9010,Y,V,TE,TI
      STOP
 9010 FORMAT(' Exception in FEOS: Y,V,TE,TI (in a.u.) = ',1P4E15.7)
c----------------------------------------------------- begin ion part:
 100  PIFOU=TI/V
      EIFOU=1.5D0*TI
      DPILTF=PIFOU
      DPILVF=-PIFOU
      DEILTF=EIFOU
      DEILVF=0.D0
c======================================================= end ion part:

c------------------------------------------------ begin electron part:
      EFE=.5D0*(3.D0*PINUM**2*Y/V)**.666666666667D0
      IF(TE.LT.TFLOOR*EFE) GOTO 1000
      PSI=EFE/TE
      APSI=AFERFIT*PSI
      APSI1=1.D0/(1.D0+APSI)
      P0=Y*TE/V

c  if chemical potential AMUFOU and free energy FEFOU are not needed,
c  the three lines below may be commented out:
c      RL15=1.5D0*LOG(1.D0+1.D0/APSI) 
c      FEFOU=Y*TE*(0.6D0*PSI-RL15)
c      AMUFOU=TE*(PSI+APSI1-RL15)

      PEFOU=P0*(0.4D0*PSI+APSI1)
      EEFOU=1.5D0*V*PEFOU
      DPELTF=P0*(1.D0+2.D0*APSI)*APSI1**2
      DPELVF=-P0*(.666666666667D0*PSI+(1.D0+APSI/3.D0)*APSI1**2)
      DEELTF=1.5D0*V*DPELTF
      DEELVF=EEFOU+1.5D0*V*DPELVF
      GOTO 8888

c: cold asymptotics:
 1000 CONTINUE
      FEFOU=0.6D0*Y*EFE
      AMUFOU=EFE
      PEFOU=0.4D0*Y*EFE/V
      EEFOU=1.5D0*V*PEFOU
      DPELTF=(2.D0*Y/(AFERFIT*V))*TE*(TE/EFE)
      DPELVF=-0.666666666667D0*Y*EFE/V
      DEELTF=1.5D0*V*DPELTF
      DEELVF=EEFOU+1.5D0*V*DPELVF
      GOTO 8888
c================================================== end electron part:

 8888 RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION TFEOS(YA,VA,PA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the temperature given the value of the 
c     total pressure (Boltzmann ions + Fermi electrons) in a.u.; 
c     works in combination with FEOS.

c:    Atomic units:
c     density: 11.206*A g/cc;       temperature: e^2/a_0 = 27.21 eV;
c     pressure:        e^2/a_0^4 = 294.2 Mbar;
c     internal energy: 27.21 eV/atom = (26.26/A)*10^12 erg/g.

c     Called by:  ...
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4par.fi'

c: initial assignments:
      Y=YA
      V=VA
      IF(Y.GT.0.D0.AND.V.GT.0.D0.AND.PA.GT.0.D0) GOTO 100
      PRINT 9010,Y,V,PA
      STOP
 9010 FORMAT(' Exception in TFEOS: Y,V,P (in a.u.) = ',1P3E15.7)

 100  EFE=.5D0*(3.D0*PINUM**2*Y/V)**.666666666667D0
      PT=PA-0.4D0*Y*EFE/V
      IF(PT.GT.0.D0) GOTO 300
      IF(PT.LT.0.D0) GOTO 200
      TFEOS=0.D0
      GOTO 8888
 200  PRINT 9020,YA,VA,PA,PT
 9020 FORMAT('  P < P_fermi for Y,V,P,PT=',1P4E15.8)
      STOP
      
 300  RAB=V*PT-AFERFIT*EFE
      RAB1=SQRT(RAB**2+4.D0*AFERFIT*(Y+1)*EFE*V*PT)
      IF(RAB.LT.0.D0) GOTO 400
      TFEOS=0.5D0*(RAB1+RAB)/(Y+1.D0)
      GOTO 8888
 400  TFEOS=2.D0*AFERFIT*EFE*V*PT/(RAB1-RAB)

 8888 RETURN
      END

C**********************************************************************
      SUBROUTINE RECORD(JOB,JCHAN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine writes or reads on disk hard-core variables. 

c     Called by:  ...
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'
      CHARACTER JOB*1

      IF(JOB.EQ.'W') GOTO 200
      READ(JCHAN,*) DUMP,DM,R,U,V,TE,TI,TR,HZ
     *,EAL,EP3,EP14,XD,XT,XHE,XH,XB,IFDMFU
     *,IFLBND,IFHZ,IFBURN,NTSLOI,NTBAD,NTVBAD
     *,NTFLLE,NTFLLI,NTFLLR,NT0000,TIME,DT,RORFU,RORAFU,RORFUM
     *,VINCR,TEINCR,VINCD,TEINCD,TVACCD,TIINCR,TRINCR,HZINCR
     *,RORP,JRORP,IFOPAC
     *,PBL,PBR,TRLEX,TREX,HZBL,HZBR,PBLSUM,PBRSUM
     *,PBLOLD,PBROLD,HBLOLD,HBROLD,PBLSOL,PBRSOL
     *,NRUN,NTPRI,NJPRI,NALPRI,NTSAVE,NSAVE,N0000
     *,IFN14,IFN2,IFAL,IFP3,IFP14
     *,EIONB,WDRIV,TDRFIN,TDRCAL,NDRIV,NDRAA
     *,IGEO,I1GEO,S1GEO,CSURF
     *,N,N1,N2,NZ,NJZ,NZEVEN,NFU,AMSFU,QZEVEN
     *,UCBET,ASBOL,SMU1,SMU2,TMU1,TMU2,SIG,SIGT
     *,APRO,ZPRO,UPRO,TPRO,FITCAP,NSUB,AMOL,XMOL,ZMOL,Z2MOL,SMOL
     *,RZ0,ROZ0,TE0,TI0,TR0,XD0,XT0
     *,XHE0,XH0,XB0,AMZ,EZ0,NMESH,IFMFU0,HZ0,DT0
     *,EZDR,EZJL,EZFUS,EZCL,EZN14,EZN2,EZAL,EZP3,EZP14
     *,EOUTW,EOUTE,EOUTI,EOUTR,EOUTHZ,EOUTAL,EOUTP3,EOUP14
     *,TNUN2,TNUN14,EREX,ERIN
     *,N1MAX,NZMAX,DTMIN,DTMAX,EEMIN,TMIN,HZMIN
     *,FLLE,FLLI,FLLR,TBURN0,XBFLOO
     *,CZDT,CZDTQ,CZDTVT,CZDDT,CZDRIV,CZVKIN,CZTKIN
       RETURN

 200  WRITE(JCHAN,*) DUMP,DM,R,U,V,TE,TI,TR,HZ
     *,EAL,EP3,EP14,XD,XT,XHE,XH,XB,IFDMFU
     *,IFLBND,IFHZ,IFBURN,NTSLOI,NTBAD,NTVBAD
     *,NTFLLE,NTFLLI,NTFLLR,NT0000,TIME,DT,RORFU,RORAFU,RORFUM
     *,VINCR,TEINCR,VINCD,TEINCD,TVACCD,TIINCR,TRINCR,HZINCR
     *,RORP,JRORP,IFOPAC
     *,PBL,PBR,TRLEX,TREX,HZBL,HZBR,PBLSUM,PBRSUM
     *,PBLOLD,PBROLD,HBLOLD,HBROLD,PBLSOL,PBRSOL
     *,NRUN,NTPRI,NJPRI,NALPRI,NTSAVE,NSAVE,N0000
     *,IFN14,IFN2,IFAL,IFP3,IFP14
     *,EIONB,WDRIV,TDRFIN,TDRCAL,NDRIV,NDRAA
     *,IGEO,I1GEO,S1GEO,CSURF
     *,N,N1,N2,NZ,NJZ,NZEVEN,NFU,AMSFU,QZEVEN
     *,UCBET,ASBOL,SMU1,SMU2,TMU1,TMU2,SIG,SIGT
     *,APRO,ZPRO,UPRO,TPRO,FITCAP,NSUB,AMOL,XMOL,ZMOL,Z2MOL,SMOL
     *,RZ0,ROZ0,TE0,TI0,TR0,XD0,XT0
     *,XHE0,XH0,XB0,AMZ,EZ0,NMESH,IFMFU0,HZ0,DT0
     *,EZDR,EZJL,EZFUS,EZCL,EZN14,EZN2,EZAL,EZP3,EZP14
     *,EOUTW,EOUTE,EOUTI,EOUTR,EOUTHZ,EOUTAL,EOUTP3,EOUP14
     *,TNUN2,TNUN14,EREX,ERIN
     *,N1MAX,NZMAX,DTMIN,DTMAX,EEMIN,TMIN,HZMIN
     *,FLLE,FLLI,FLLR,TBURN0,XBFLOO
     *,CZDT,CZDTQ,CZDTVT,CZDDT,CZDRIV,CZVKIN,CZTKIN
      RETURN
      END

C**********************************************************************
      DOUBLE PRECISION FUNCTION GQUAD(F,A,B,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates Gaussian quadrature = \int_A^B  F(x) dx
c     Called by:  ...
c     Calls    :  none
c =====================================================================
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

C*****************************************************************
      DOUBLE PRECISION FUNCTION ZERA0(F,XL,XR,X0)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine calculates the root X of equation  F(X)=0 on the 
c     interval XL<X<XR; absolute accuracy of X is better than 1.D-14*XL

c     Called by:  ...
c     Calls    :  none
c =====================================================================

      NRCMAX=300
      X=X0
      H=.2D0*(XR-XL)
      HMIN=1.D-15*MAX(ABS(XL),ABS(XR))
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
 50   IF(DABS(H).LT.HMMIN) GOTO 77
      IF(NRCAL.GT.NRCMAX) GOTO 66
c      IF(XR-X.LT.HMMIN.AND.KFR.GT.2) GOTO 66
c      IF(X-XL.LT.HMMIN.AND.KFL.GT.2) GOTO 66
      GOTO 10

 66   PRINT 511,NRCAL,XL,X,XR,H,FP,F2
 77   ZERA0=X
      RETURN

 511  FORMAT(' F(X).NE.0 in ZERA0: NRCAL=',I6/' XL,X,XR=',1P3E15.8/
     *' H=',1PE15.8,' FP,F2=',1P2E15.8)
      END

C**********************************************************************
      SUBROUTINE DEFVAL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c ---------------------------------------------------------------------
c     This routine assigns default values of various parameters,
c     variables and constants; they may be overwritten in DEIRA, START;
c     this routine replaces the original BLOCK DATA.

c     Called by:  DEIRA
c     Calls    :  none
c =====================================================================
      INCLUDE 'de4com.fi'

C: physical constants:
      UCBET=2.998D3
      ASBOL=1.372D0
      APRO(1)=4.D0
      APRO(2)=1.D0
      APRO(3)=1.D0
      ZPRO(1)=2.D0
      ZPRO(2)=1.D0
      ZPRO(3)=1.D0
      UPRO(1)=130.3D0
      UPRO(2)=240.5D0
      UPRO(3)=530.2D0
      TPRO(1)=20.D0
      TPRO(2)=60.D0
      TPRO(3)=300.D0
      DO K=1,5
       FITCAP(K)=1.D0
      ENDDO
c: artificial viscosity:
      SMU1=0.1D0
      SMU2=0.D0
      SIG=1.D0
      TMU1=0.D0
      TMU2=2.D0
      SIGT=0.1D0

c: limiting and accuracy constants:
      DTMIN=1.D-18
      DTMAX=1.D0
      HZMIN=1.D-28
      TMIN=3.D-7
      TBURN0=0.1D0
      EEMIN=1.D-4
      XBFLOO=1.D-9
      FLLE=0.5D0
      FLLI=0.5D0
      FLLR=0.5D0
      CZDT=0.5D0
      CZDTQ=0.2D0
      CZDTVT=0.2D0
      CZDDT=0.2D0
      CZDRIV=0.2D0
      CZVKIN=0.1D0
      CZTKIN=0.03D0

c: control parameters:
      IFOPAC=1

      IGEO=0
      NRUN=0
      NTPRI=10 000 000
      NJPRI=1
      NALPRI=1
      NTSAVE=10 000 000
      IFN14=0
      IFN2=0
      IFAL=0
      IFP3=0
      IFP14=0

      NZ=1
      NZEVEN=1
      QZEVEN=1.D0

c: initial target state:
      HZ0=0.D0
      DT0=1.D-8
      DO I=1,NNZ
        ROZ01(I)=-1.D0
        BRO0PR(I)=1.D0
        TE0(I)=0.D0
        TI0(I)=0.D0
        TR0(I)=0.D0
        XD0(I)=0.D0
        XT0(I)=0.D0
        XHE0(I)=0.D0
        XH0(I)=0.D0
        XB0(I)=0.D0
        IFMFU0(I)=0
C mychanging
      TE0(1)=0.14D0    !changing
      TI0(1)=0.14D0
      !TR0(1)=0.D0
      !TE0(2)=0.D0    !changing
      !TI0(2)=0.D0
      !TR0(2)=0.D0
      ENDDO
      RETURN
      END
