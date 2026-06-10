C***********************************************************************
C     ! Isotropic TOTAL SCA 3D 
C     ! ----------------------------------------------------------------
C     ! Contains the code concerned with the formulation for SCA in 
C     ! orthotropic composite material, this is used in the 
C     ! 'main_UD_tot.for' file
C     ! ----------------------------------------------------------------
C     ! Last Debugged: 03/09/2025 
C     ! meka1@purdue.edu
C***********************************************************************
      SUBROUTINE elas3dortho(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, DSTRAN, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)
          
      IMPLICIT NONE
      
C     STANDARD VARIABLES
      INTEGER::NSTATV,NPROPS,NOEL
      REAL*8::stepTIME,totalTIME,DT,L
      REAL*8::PROPS, SIG(6)
      REAL*8::EPS(6), DSTRAN(6)
      REAL*8::STATEV
      REAL*8::DDSDDE(6,6)      
      REAL*8::DTEMP
      INTEGER::isIMPLICIT
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)

C     DEFINED VARIABLES
      REAL*8::Dco(6,6),Sco(6,6)
      REAL*8::G1C_F
      REAL*8::E,nu,XT
      REAL*8::isCRACKED, MODE, NONDIM
      INTEGER::counter
      REAL*8::EPS_MAX,EPS_OLD
      REAL*8::FI_FT
      REAL*8::SLOPE
      REAL*8::temp_sig
      REAL*8::PEAK, TEMP
      REAL*8::SECANT
      INTEGER::ZONE

C     DEFINED FUNCTIONS


C***********************************************************************
C       !State Variables List
C-----------------------------------------------------------------------
C       (1) Crack Flag
C       (2) MODE
C-----------------------------------------------------------------------
C       !Property Variables List
C-----------------------------------------------------------------------
C       (1) Elastic Modulus
C       (2) Poisson's ratio
C       (3) Fracture Toughness-Fibre  (Tension)
C       (4) Tensile Strength      (Longitudinal)

C***********************************************************************
C     MODE = 1: Longitudinal Tension
C-----------------------------------------------------------------------
C     !Property List
      E      = props(1)
      nu     = props(2)

C     !Fracture Energy Values
      G1C_F  = props(3)

C     !Strength Values
      XT     = props(4)



      RETURN
      END
     
C=======================================================================     
C     DEFINE FUNCTIONS
      SUBROUTINE STIFFNESS_MATRIX(E, nu, Dco)
      IMPLICIT NONE
      REAL*8::E,nu
      REAL*8::Dco(6,6),Sco(6,6)

C     Isotropic Material 
      Sco = reshape(
     1   (1.0d0/E)*(/ 1.0d0,-1.0d0*nu,-1.0d0*nu,   0.0d0,   0.0d0,   0.0d0, 
     2            -1.0d0*nu,    1.0d0,-1.0d0*nu,   0.0d0,   0.0d0,   0.0d0, 
     3            -1.0d0*nu,-1.0d0*nu,    1.0d0,   0.0d0,   0.0d0,   0.0d0, 
     4                0.0d0,    0.0d0,    0.0d0,1.0d0+nu,   0.0d0,   0.0d0,
     5                0.0d0,    0.0d0,    0.0d0,   0.0d0,1.0d0+nu,   0.0d0,
     6                0.0d0,    0.0d0,    0.0d0,   0.0d0,   0.0d0,1.0d0+nu/),(/6,6/))
      
C     Stiffness Matric of the Continuum (Dco)
      call matrixInverse(Sco,Dco,6,6)

      RETURN
      END

      
      