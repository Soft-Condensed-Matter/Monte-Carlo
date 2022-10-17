!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%***************************************************************************%%
!%%***************************************************************************%%
!%%**  PROGRAM      MONTE CARLO MOLECULAR SIMULATION                        **%%
!%%**  AUTHOR       ALEXIS TORRES CARBAJAL @Alpixels                        **%%
!%%**  LICENSE      LGPL-V3                                                 **%%
!%%**                                                                       **%%
!%%**  ENSEMBLE     NVT                                                     **%%
!%%**  ALGORITHM    METROPOLIS                                              **%%
!%%**  DATE         SEPTEMBER 27, 2022                                      **%%
!%%**                                                                       **%%
!%%**  OBS          HARD-SPHERE FLUID IN THE CANONICAL ENSEMBLE             **%%
!%%**               THREE DIMENSIONAL MONO-COMPONENT                        **%%
!%%**                                                                       **%%
!%%***************************************************************************%%
!%%***************************************************************************%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MCVAR
 IMPLICIT NONE
 INTEGER, PARAMETER:: D       = KIND(1.0D0) !PRESICION
 INTEGER, PARAMETER:: NPARTX  = 2000        !MAXIMUM NUMBER OF PARTICLES
 INTEGER, PARAMETER:: OUTDISP = 1000        !SHOW DATA ON DISPLAY
 INTEGER, PARAMETER:: OUTSAMP = 10          !TAKE A SAMPLE
 INTEGER, PARAMETER:: OUTFILM = 100         !TAKE A PICTURE
 INTEGER, PARAMETER:: NBINX   = 10000       !MAXIMUM NUMBER OF BINS

 REAL(D), PARAMETER:: SIGMA   = 1.0D0       !PARTICLE DIAMETER
 REAL(D), PARAMETER:: ESTAR   = 1.0D0       !PARTICLE INTERACTION WELL
 REAL(D), PARAMETER:: DB      = 0.01D0      !HISTOGRAM BIN WIDTH
 REAL(D), PARAMETER:: PI      = DACOS(-1.0D0)!PI NUMBER 
 REAL(D), PARAMETER:: CUT     = 1.1D0       !CUT-OFF RADIUS
 REAL(D), PARAMETER:: DSIGMA  = 0.02D0      !HEAVISIDE 

 INTEGER:: NC,NPART,MCCYC,NACC,MCTRY 
 INTEGER:: IMC,JMC,MCMV,IDUMM,NB,ISAM
 INTEGER:: IFRAME,SSTAT
 INTEGER:: HGZ(NBINX,2)

 REAL(D):: RHOSTAR,TEMP,BETA,VOL,PRESS,FRACC
 REAL(D):: RCUT,ENER,VIR,DR,PHI,LAMBDA
 REAL(D):: BOXX,BOXY,BOXZ,IBOXX,IBOXY,IBOXZ
 REAL(D):: RX(NPARTX),RY(NPARTX),RZ(NPARTX)
 REAL(D):: SUMEN,SUMVR,SUMPR,SUM2EN,SUM2VR,SUM2PR
 REAL(D):: SIGMA12,WX,WY,WZ,PXX,PYY,PZZ,EPOT,PVNKT
 REAL(D):: SUMPRESS,SUMPVNKT,SM2PRESS,SM2PVNKT
END MODULE MCVAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   MAIN PROGRAM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PROGRAM MCNVT
 USE MCVAR
 IMPLICIT NONE

 CALL BEGINMC                               !BEGIN SIMULATION

 DO IMC=1,MCCYC                             !MONTE CARLO CYCLES

    DO JMC=1,MCMV                           !MONTE CARLO STEPS
       CALL MCMOVE                          !MONTE CARLO MOVEMENT
    ENDDO

    CALL COMPUTE                            !COMPUTE OBSERVABLES
 ENDDO

 CALL FINISHMC                              !END SIMULATION

 STOP
END PROGRAM MCNVT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE BEGINMC
 USE MCVAR
 IMPLICIT NONE

! PRINT*,'************************************************'
! PRINT*,'**             MONTE CARLO: NVT               **'
! PRINT*,'************************************************'

 CALL INPMC                                 !SIMULATION INPUT DATA
 CALL CNFCC                                 !SET UP INITIAL CONFIGURATION
 CALL SETMC                                 !SET UP SIMULATION PARAMETERS
 CALL ETOT                                  !SYSTEM TOTAL ENERGY
 CALL DSPLY                                 !SHOW DATA ON DISPLAY

 OPEN(UNIT=13,FILE="MCStatus.dat")
 OPEN(UNIT=17,FILE="MCMovie.xyz")

 RETURN
END SUBROUTINE BEGINMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE INPMC
 USE MCVAR
 IMPLICIT NONE
 
 OPEN(UNIT=10,FILE="MC.inp",STATUS='OLD')

 READ(10,*)NC                               !NUMBER OF CELLS
 READ(10,*)PHI                              !PACKING
 READ(10,*)TEMP                             !SYSTEM TEMPERATURE
 READ(10,*)MCCYC                            !MONTE CARLO CYCLES
 READ(10,*)SSTAT                            !EQUILIBRIUM CYCLES
 READ(10,*)DR                               !MAXIMUM DISPLACEMENT
 READ(10,*)LAMBDA                           !ATTRACTIVE INTERACTION WELL

 CLOSE(UNIT=10)

 RETURN
END SUBROUTINE INPMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE CNFCC
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: M,I,IX,IY,IZ,LAT
 REAL(D):: CELLX,CELLY,CELLZ
 REAL(D):: CELLUX,CELLUY,CELLUZ

 NPART=4*NC**3                              !NUMBER OF PARTICLES
 RHOSTAR=6.0D0*PHI/PI                       !SYSTEM DENSITY
 BOXX=(DBLE(NPART)*SIGMA**3/RHOSTAR)**(1.0D0/3.0D0)  !BOX LENGHT

 M=0
 CELLX=BOXX/DBLE(NC)
 CELLY=BOXX/DBLE(NC)
 CELLZ=BOXX/DBLE(NC)
 CELLUX=0.5D0*CELLX
 CELLUY=0.5D0*CELLY
 CELLUZ=0.5D0*CELLZ

 RX(1)=0.0D0
 RY(1)=0.0D0
 RZ(1)=0.0D0

 RX(2)=CELLUX
 RY(2)=CELLUY
 RZ(2)=0.0D0

 RX(3)=0.0D0
 RY(3)=CELLUY
 RZ(3)=CELLUZ

 RX(4)=CELLUX
 RY(4)=0.0D0
 RZ(4)=CELLUZ


 DO IX=1,NC
    DO IY=1,NC
       DO IZ=1,NC
          DO LAT=1,4
             RX(LAT + M)=RX(LAT) + CELLX*DBLE(IX - 1)
             RY(LAT + M)=RY(LAT) + CELLY*DBLE(IY - 1)
             RZ(LAT + M)=RZ(LAT) + CELLZ*DBLE(IZ - 1)
          ENDDO
          M=M+4
       ENDDO
    ENDDO
 ENDDO

 DO I=1,NPART
    RX(I)=RX(I) - 0.5D0*BOXX
    RY(I)=RY(I) - 0.5D0*BOXX
    RZ(I)=RZ(I) - 0.5D0*BOXX
 ENDDO

 OPEN(UNIT=11,FILE="MCPic.xyz")             !SAVE INITIAL CONFIGURATION
 WRITE(11,*)NPART
 WRITE(11,*)'FRAME',1

 DO I=1,NPART
    WRITE(11,*)'C',RX(I),RY(I),RZ(I)
 ENDDO

 CLOSE(UNIT=11)

 RETURN
END SUBROUTINE CNFCC
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SETMC
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: I

 BOXY=BOXX                                  !BOX Y SIZE LENGHT 
 BOXZ=4.0D0*BOXX                            !BOX Z SIZE LENGHT 
 IBOXX=1.0D0/BOXX                           !INVERSE BOX X SIZE LENGHT 
 IBOXY=1.0D0/BOXY                           !INVERSE BOX Y SIZE LENGHT 
 IBOXZ=1.0D0/BOXZ                           !INVERSE BOX Z SIZE LENGHT 
 BETA=1.0D0/TEMP                            !INVERSE OF SYSTEM TEMPERATURE
 VOL=BOXX*BOXY*BOXZ                         !SIMULATION BOX VOLUME

 RCUT=CUT*LAMBDA*SIGMA                      !CUT-OFF RADIUS

 SUMPRESS=0.0D0                             !PRESS ACCUMULATOR
 SUMPVNKT=0.0D0
 SM2PRESS=0.0D0
 SM2PVNKT=0.0D0

 SUMEN=0.0D0                                !ENERGY AVERAGE
 SUMVR=0.0D0                                !VIRIAL AVERAGE
 SUMPR=0.0D0                                !PRESSION AVERAGE

 SUM2EN=0.0D0                               !ENERGY SQUARE AVERAGE
 SUM2VR=0.0D0                               !VIRIAL SQUARE AVERAGE
 SUM2PR=0.0D0                               !PRESSION SQUARE AVERAGE 

 NACC=0                                     !ACCPETED TRIES
 ISAM=0                                     !NUMBER OF SAMPLE
 MCTRY=0                                    !NUMBER OF TRIES
 IFRAME=0                                   !NUMBER OF FRAME

 MCMV=NPART                                 !NUMBER OF MONTE CARLO MOVEMENTS

 CALL SYSTEM_CLOCK(IDUMM)                   !SEED ACCORDING MACHINE TIME
 IDUMM=-IDUMM                               !SEED FOR RANDOM NUMBER GENERATOR

 NB=INT(BOXZ/DB) + 1                        !NUMBER OF BINS IN HISTROGRAM

 DO I=1,NB                                  !SET TO ZERO g(r) HISTOGRAM
    HGZ(I,1)=0
    HGZ(I,2)=0
 ENDDO

 RETURN
END SUBROUTINE SETMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE DSPLY
 USE MCVAR
 IMPLICIT NONE 

 OPEN(UNIT=12,FILE="WazapMC.txt")

 ENER=ENER/DBLE(NPART)
 VIR=VIR/DBLE(NPART)

 WRITE(6,01)  RHOSTAR,NPART,SIGMA,ESTAR,TEMP,LAMBDA,BOXX,DR,RCUT,ENER,VIR,MCCYC
 WRITE(12,01) RHOSTAR,NPART,SIGMA,ESTAR,TEMP,LAMBDA,BOXX,DR,RCUT,ENER,VIR,MCCYC

 CLOSE(UNIT=12)

 01 FORMAT(1X,//,'***  SIMULATION PARAMETERS  *** ',// &
          'REDUCED DENSITY               ',F15.8,/ &
          'NUMBER OF PARTICLES           ',I10  ,/ &
          'SOLVENT DIAMETER              ',F15.8,/ &
          'SOLVENT INTERACTION           ',F15.8,/ &
          'TEMPERATURE                   ',F15.8,/ &
          'INTERACTION RANGE             ',F15.8,/ &
          'SIMULATION X BOX SIZE         ',F15.8,/ &
          'MAXIMUM DISPLACEMENT          ',F15.8,/ &
          'CUT OFF RADIUS                ',F15.8,/ &
          'INITIAL ENERGY                ',F15.8,/ &
          'INITIAL VIRIAL                ',F15.8,/ &
          'SIMULATION CYCLES             ',I10  ,//)

 RETURN
END SUBROUTINE DSPLY
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ETOT
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: I,JB
 REAL(D):: ENI
 REAL(D):: RXI,RYI,RZI

 VIR=0.0D0                                  !VIRIAL
 ENER=0.0D0                                 !TOTAL ENERGY

 DO I=1,NPART-1                             !LOOP OVER PAIR PARTICLES
    RXI=RX(I)                               !I PARTICLE X POSITION 
    RYI=RY(I)                               !I PARTICLE Y POSITION 
    RZI=RZ(I)                               !I PARTICLE Z POSITION 

    JB=I+1                                  !NEIGBOUR PARTICLE ID

    CALL ENERI(RXI,RYI,RZI,I,JB,ENI)   

    ENER=ENER + ENI
 ENDDO

 RETURN
END SUBROUTINE ETOT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ENERI(RXI,RYI,RZI,I,JB,EN)
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: I,J,JB
 REAL(D):: EN,BUR,DBUR
 REAL(D):: RXI,RYI,RZI,DX,DY,DZ,RIJ

 EN=0.0D0

 DO J=JB,NPART
    IF(J .NE. I)THEN
      DX=RXI - RX(J)                        !X SEPARATION BETWEEN I & J PARTICLE
      DY=RYI - RY(J)                        !Y SEPARATION BETWEEN I & J PARTICLE
      DZ=RZI - RZ(J)                        !Z SEPARATION BETWEEN I & J PARTICLE

      DX=DX - ANINT(DX*IBOXX)*BOXX          !X MINIMUM IMAGE CONDITION
      DY=DY - ANINT(DY*IBOXY)*BOXY          !Y MINIMUM IMAGE CONDITION
      DZ=DZ - ANINT(DZ*IBOXZ)*BOXZ          !Z MINIMUM IMAGE CONDITION

      RIJ=DSQRT(DX*DX + DY*DY + DZ*DZ)      !SEPATATION DISTANCE
      IF(RIJ .LT. RCUT)THEN
        IF(RIJ .LT. SIGMA)THEN              !HARD-SPHERE CONTRIBUTION
          BUR=10E5
          DBUR=10E5
          EN=EN + BUR
        ELSEIF(RIJ .LE. LAMBDA)THEN         !SQUARE WELL CONTRIBUTION
          BUR=-ESTAR                        !SQUARE-WELL DEPTH
          DBUR=0.0D0
          EN=EN + BUR
        ENDIF
      ENDIF
    ENDIF
 ENDDO
 
 RETURN
END SUBROUTINE ENERI
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MCMOVE 
 USE MCVAR
 IMPLICIT NONE
 LOGICAL:: ACCEPT
 INTEGER:: JB,O
 REAL(D):: RAN1,ARG
 REAL(D):: RXO,RYO,RZO,ENO
 REAL(D):: RXN,RYN,RZN,ENN

 MCTRY=MCTRY+1                              !TRY NUMBER OF MC MOVEMENTS
 JB=1                                       !PARTICLE ONE ID
 O=INT(DBLE(NPART)*RAN1(IDUMM)) + 1         !CHOOSE A RANDOM ID PARTICLE
 RXO=RX(O)                                  !X OLD POSITION
 RYO=RY(O)                                  !Y OLD POSITION
 RZO=RZ(O)                                  !Z OLD POSITION

 CALL ENERI(RXO,RYO,RZO,O,JB,ENO)           !OLD ENERGY

 RXN=RXO + (RAN1(IDUMM) - 0.5D0)*DR         !X NEW POSITION
 RYN=RYO + (RAN1(IDUMM) - 0.5D0)*DR         !Y NEW POSITION
 RZN=RZO + (RAN1(IDUMM) - 0.5D0)*DR         !Z NEW POSITION

 CALL ENERI(RXN,RYN,RZN,O,JB,ENN)           !NEW ENERGY

 ARG=ENN - ENO                              !ENERGY DIFFERENCE

 IF(ARG .LE. 0.0D0)THEN                     !METROPOLIS CRITERIA
   ACCEPT= .TRUE.
 ELSEIF(RAN1(IDUMM) .LT. EXP(-BETA*ARG))THEN
   ACCEPT= .TRUE.
 ELSE
   ACCEPT= .FALSE.   
 ENDIF

 IF(ACCEPT)THEN
   NACC=NACC + 1                            !MC MOVEMENT ACCEPTED

   ENER=ENER + ARG                          !NEW ENERGY

   RXN=RXN - ANINT(RXN*IBOXX)*BOXX          !X PERIODIC BOUNDARY CONDITIONS
   RYN=RYN - ANINT(RYN*IBOXY)*BOXY          !Y PERIODIC BOUNDARY CONDITIONS
   RZN=RZN - ANINT(RZN*IBOXZ)*BOXZ          !Z PERIODIC BOUNDARY CONDITIONS

   RX(O)=RXN                                !UPDATE NEW X POSITION 
   RY(O)=RYN                                !UPDATE NEW Y POSITION
   RZ(O)=RZN                                !UPDATE NEW Z POSITION
 ENDIF

 RETURN 
END SUBROUTINE MCMOVE 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COMPUTE
 USE MCVAR
 IMPLICIT NONE

 EPOT=ENER/DBLE(NPART)
 CALL ADJUST

 IF(MOD(IMC,OUTDISP) .EQ. 0)THEN
   WRITE(6,02)IMC,EPOT,DR,FRACC
   WRITE(13,02)IMC,EPOT,DR,FRACC
 ENDIF

 IF(IMC .GT. SSTAT)THEN
   IF(MOD(IMC,OUTSAMP) .EQ. 0)CALL SAMPLE
!   IF(MOD(IMC,OUTFILM) .EQ. 0)CALL MCFILM
 ENDIF

 02 FORMAT(T10,I7,3(F12.5))
 RETURN 
END SUBROUTINE COMPUTE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ADJUST
 USE MCVAR
 IMPLICIT NONE

 IF(MOD(IMC,10) .EQ. 0)THEN
   FRACC=100.0D0*DBLE(NACC)/DBLE(MCTRY)
   IF(FRACC .GT. 40.0D0)THEN
     DR=DR*1.050D0
   ELSE
     DR=DR*0.950D0
   ENDIF
   IF(DR .GT. 0.5D0*BOXX)THEN
      DR=0.5D0*BOXX
   ENDIF
   NACC=0
   MCTRY=0
 ENDIF

 RETURN
END SUBROUTINE ADJUST
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE SAMPLE
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: I,KK
 REAL(D):: Z,RCM

 ISAM=ISAM + 1                              !NUMBER OF SAMPLE
 RCM=0.0D0

 DO I=1,NPART
    RCM=RCM + RZ(I)
 ENDDO
 
 RCM=RCM/DBLE(NPART)                        !Z MASS CENTER POSITION

 DO I=1,NPART
    Z=RZ(I) + 0.5D0*BOXZ
    KK=INT(Z/DB) + 1
    IF(KK .LE. NBINX)HGZ(KK,1)=HGZ(KK,1) + 1

    !PROFILE ACCORDING TO THE CENTER OF MASS
    Z=Z - RCM
    KK=INT(Z/DB) + 1
    IF(KK .LE. NBINX)HGZ(KK,2)=HGZ(KK,2) + 1
 ENDDO
 
 RETURN
END SUBROUTINE SAMPLE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE FINISHMC
 USE MCVAR
 IMPLICIT NONE

 CALL MCGZ                                  !COMPUTE DENSITY PROFILE
 CALL MCAV                                  !COMPUTE AVERAGES

 CLOSE(UNIT=13)

 RETURN
END SUBROUTINE FINISHMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MCGZ
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: I
 REAL(D):: Z,Z1,DVOL,GZ1,GZ2

 OPEN(UNIT=17,FILE="MCGz.dat")
 
 DO I=1,NB
    Z=DBLE(I)*DB
    Z1=DBLE(I-1)*DB
    DVOL=(Z - Z1)*BOXX*BOXY
    GZ1=DBLE(HGZ(I,1))/(DVOL*DBLE(ISAM))
    GZ2=DBLE(HGZ(I,2))/(DVOL*DBLE(ISAM))
    Z=(DBLE(I)-0.5D0)*DB - 0.5*BOXZ
    WRITE(17,107)Z,GZ1,GZ2
 ENDDO

 CLOSE(UNIT=17)
 107 FORMAT(3(F15.8))

 RETURN
END SUBROUTINE MCGZ
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MCAV
 USE MCVAR 
 IMPLICIT NONE
 INTEGER:: I

 OPEN(UNIT=15,FILE="MCCnf.dat")
 OPEN(UNIT=16,FILE="MCAvr.dat")

 DO I=1,NPART
    WRITE(15,105)RX(I),RY(I),RZ(I)
 ENDDO

 SUMPRESS=SUMPRESS/DBLE(ISAM)
 SUMPVNKT=SUMPVNKT/DBLE(ISAM)
 SM2PRESS=SM2PRESS/DBLE(ISAM)
 SM2PVNKT=SM2PVNKT/DBLE(ISAM)
 SM2PRESS=SM2PRESS - SUMPRESS*SUMPRESS
 SM2PVNKT=SM2PVNKT - SUMPVNKT*SUMPVNKT

 WRITE(6,02) SUMPRESS,SM2PRESS,SUMPVNKT,SM2PVNKT
 WRITE(16,02)SUMPRESS,SM2PRESS,SUMPVNKT,SM2PVNKT


 CLOSE(UNIT=15)
 CLOSE(UNIT=16)

 02 FORMAT(1X,//,'***  SIMULATION RESULTS  *** ',// &
          'P=    ',F15.8,'+/-',F15.8,/ &
          'Z=    ',F15.8,'+/-',F15.8,//)

 105 FORMAT(3(F16.8))

 RETURN
END SUBROUTINE MCAV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE MCFILM
 USE MCVAR
 IMPLICIT NONE
 INTEGER:: I
                                            !VMD AND XMAKEMOL FORMAT
 IFRAME=IFRAME+1
 WRITE(17,*)NPART
 WRITE(17,*)'FRAME',IFRAME

 DO I=1,NPART                               !NANOPARTICLES POSITIONS
    WRITE(17,*)'C',RX(I),RY(I),RZ(I)
 ENDDO

 RETURN
END SUBROUTINE MCFILM
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%   FUNCTIONS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION RAN1(IDUMM)
 INTEGER, PARAMETER:: D = KIND(1.0D0)
 INTEGER:: IDUMM,IA,IM,IQ,IR,NTAB,NDIV
 REAL(D):: RAN1,AM,EPS,RNMX
 PARAMETER( IA=16807, IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
            NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.0-EPS )
 INTEGER:: J,K,IV(NTAB),IY
 SAVE IV,IY
 DATA IV /NTAB*0/,IY /0/

 IF(IDUMM .LE. 0 .OR.  IY .EQ. 0)THEN
   IDUMM=MAX(-IDUMM,1)
   DO J=NTAB+8,1,-1
      K=IDUMM/IQ
      IDUMM=IA*(IDUMM-K*IQ)-IR*K
      IF(IDUMM .LT. 0)IDUMM=IDUMM + IM
      IF(J .LE. NTAB)IV(J)=IDUMM
   ENDDO
   IY=IV(1)
 ENDIF

 K=IDUMM/IQ
 IDUMM=IA*(IDUMM-K*IQ) - IR*K
 IF(IDUMM .LT. 0)IDUMM=IDUMM + IM
 J=1 + IY/NDIV
 IY=IV(J)
 IV(J)=IDUMM
 RAN1=MIN(AM*IY,RNMX)

 RETURN
END FUNCTION RAN1
