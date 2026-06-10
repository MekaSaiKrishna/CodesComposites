C     ! Incremental Orthotropic Elasticity
C     ! -----------------------------------------------------
C     ! 
C     ! -----------------------------------------------------
C     ! Last Debugged: 11/17/2025 
C     ! meka1@purdue.edu
C     !! Note: Double Exlcamation == Notes to myself

!     ordering scheme [eps11 eps22 eps33 gam12 gam13 gam23]
!                     [sig11 sig22 sig33 sig12 sig13 sig23]

!     remaining issues:  
C     
C-----------------------------------------------------------------------
C       SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
C      1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
C      2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
C      3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
C      4 KSPT, KSTEP, KINC)

C       IMPLICIT NONE
C C     
C       INCLUDE 'ABA_PARAM.INC'
C C     
C       CHARACTER*8 CMNAME
C C     
C       DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
C      1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
C      2 PREDEF(1), DPRED(2), PROPS(NPROPS), COORDS(3), DROT(3, 3),
C      3 DFGRD0(3, 3), DFGRD1(3, 3), TIME(2)

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)


      real*8::eps(6)

      eps=STRAN+DSTRAN


C-----------------------------------------------------------------------
C     INCREMENTAL FORMULATION
      call  elasortho(NSTATV, NPROPS,NOEL,
     1       TIME(1), TIME(2), PREDEF(1), DTIME, CELENT,
     2       PROPS, eps, STRESS, STATEV, DDSDDE, DSTRAN, 1)
      RETURN
      END

C       include 'standard_support.for'


C-----------------------------------------------------------------------
      SUBROUTINE elasortho(nstatv, nprops, noel, stepTime, totalTime,
     1 predef, dt, L, props, eps, sig, statev, ddsdde, deps, isImplicit)

      IMPLICIT NONE
      
C     STANDARD VARIABLES
      real*8::sig(6)
      real*8::eps(6), deps(6)
      real*8::ddsdde(6,6)
      real*8::statev
      real*8::props 
      real*8::stepTime,totalTime,dt,L,predef
      integer::nstatv,nprops,noel,isImplicit     

      DIMENSION STATEV(NSTATV),PROPS(NPROPS)

      INTEGER::i,j
      REAL*8::Dco(6,6)

      REAL*8::EMOD1,EMOD2,EMOD3,ENU12,ENU13,ENU23,EG12,EG13,EG23
      REAL*8::ZERO,ONE,TWO
C     
      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0)
C     
C-----------------------------------------------------------------------
C     UMAT FOR ORTHOTROPIC MATERIALS 
C----------------------------------------------------------------------- 
C     PROPS(1) - EMOD1
C     PROPS(2) - EMOD2
C     PROPS(3) - EMOD3
C     PROPS(4) - ENU12
C     PROPS(5) - ENU13
C     PROPS(6) - ENU23
C     PROPS(7) - EG12
C     PROPS(8) - EG13
C     PROPS(9) - EG23
C-----------------------------------------------------------------------
C     ELASTIC PROPERTIES
      EMOD1=PROPS(1)
      EMOD2=PROPS(2)
      EMOD3=PROPS(3)
      ENU12=PROPS(4)
      ENU13=PROPS(5)
      ENU23=PROPS(6)
      EG12=PROPS(7)
      EG13=PROPS(8)
      EG23=PROPS(9)
    
      write(*,*) 'predef: ',predef
     
      CALL STIFFNESS_MATRIX(EMOD1,EMOD2,EMOD3,ENU12, ENU13, ENU23, 
     1      EG12, EG13, EG23, Dco)

C     STRESS UPDATE
      DO i=1,6
            STATEV(i) = sig(i)
      END DO

      DO i=1,6
            DO j=1,6
                  sig(i) = STATEV(i) + Dco(i,j)*deps(j)
                  STATEV(i) = sig(i)
            END DO
      END DO
      

C     JACOBIAN UPDATE
      IF(isIMPLICIT .EQ. 1) ddsdde = Dco

      RETURN
      END

C=======================================================================     
C     DEFINE FUNCTIONS
      SUBROUTINE STIFFNESS_MATRIX(EMOD1,EMOD2,EMOD3,ENU12,ENU13,ENU23, 
     1      EG12,EG13,EG23,Dco)
      IMPLICIT NONE
      REAL*8::EMOD1,EMOD2,EMOD3,ENU12,ENU13,ENU23,EG12,EG13,EG23
      REAL*8::ENU21,ENU32,ENU31,DELTA
      REAL*8::D1111,D2222,D3333,D1122,D1133,D2233,D1212,D1313,D2323
      REAL*8::Dco(6,6)
      REAL*8::ZERO,ONE,TWO

      PARAMETER (ZERO=0.D0, ONE=1.D0, TWO=2.D0)

       ENU21=(EMOD2/EMOD1)*ENU12
       ENU32=(EMOD3/EMOD2)*ENU23
       ENU31=(EMOD3/EMOD1)*ENU13
C     
       DELTA=ONE-(ENU12*ENU21)-(ENU23*ENU32)-(ENU13*ENU31)-(TWO*ENU12*ENU32*ENU13)
C
       D1111=EMOD1*(ONE-(ENU23*ENU32))/DELTA
       D2222=EMOD2*(ONE-(ENU13*ENU31))/DELTA
       D3333=EMOD3*(ONE-(ENU12*ENU21))/DELTA
       D1122=EMOD1*(ENU12+(ENU32*ENU13))/DELTA
       D1133=EMOD3*(ENU13+(ENU12*ENU23))/DELTA
       D2233=EMOD3*(ENU23+(ENU21*ENU13))/DELTA
       D1212=EG12
       D1313=EG13
       D2323=EG23
C
C     ELASTIC STIFFNESS
C
      Dco(1, 1)=D1111
      Dco(1, 2)=D1122
      Dco(2, 1)=D1122
      Dco(1, 3)=D1133
      Dco(3, 1)=D1133
      Dco(2, 2)=D2222
      Dco(2, 3)=D2233
      Dco(3, 2)=D2233
      Dco(3, 3)=D3333
      Dco(4, 4)=D1212
      Dco(5, 5)=D1313
      Dco(6, 6)=D2323  

      RETURN
      END