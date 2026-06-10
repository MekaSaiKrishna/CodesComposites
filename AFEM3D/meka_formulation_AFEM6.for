C***********************************************************************
C     ! Isotropic Incremental SCA 3D - GLOBAL
C     ! ----------------------------------------------------------------
C     ! Contains the code concerned with the formulation for SCA in 
C     ! orthotropic composite material, this is used in the 
C     ! 'main_UD_tot.for' file
C     ! ----------------------------------------------------------------
C     ! Last Debugged: 11/12/2025 
C     ! meka1@purdue.edu
C***********************************************************************
      SUBROUTINE sca3dortho(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, DSTRAN, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)
          
      IMPLICIT NONE
      
C     STANDARD VARIABLES
      INTEGER::NSTATV,NPROPS,NOEL
      REAL*8::stepTIME,totalTIME,DT,L
      REAL*8::PROPS,SIG(6)
      REAL*8::EPS(6),DSTRAN(6)
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

      REAL*8::EPSN, EPSF, ALF1, ALF2
      REAL*8::SIG_CORR,SIG_COMP
      REAL*8::EPS_STAR, DEPS1, DEPS2
      INTEGER::NTENS, K1

C     DEFINED FUNCTIONS


C***********************************************************************
C       !State Variables List
C-----------------------------------------------------------------------
C       (1) Crack Flag
C       (2) MODE
C       (3) EPS_MAX: stores the max strain
C       (4) ZONE: store the loading zone in the traction separation law
C       (5) : store the stress at the maximum strain encountered (EPS_MAX)
C       (6) SLOPE: store the slope value based on the zone 
C       (7) : stores the max stress - 
C       (8) SECANT: store instantaneous secant modulus
C       (9) NONDIM: Strain energy ratio parameter
C       (10)CELENT: Length of element
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

      ! PARAMETERS

      NTENS = 6
      
      EPSN = XT/E                    !Global Strain @ Peak
      EPSF = (2.0d0*G1C_F)/(XT*L)    !Global Strain @ Failure

      ALF1 = E                       ! Pre-Peak Global Slope
      ALF2 = (-1.0d0*XT)/(EPSF-EPSN) ! Post-Peak Global Softening Slope

C-----------------------------------------------------------------------
      isCRACKED = STATEV(1)  ! CRACK FLAG
      MODE      = STATEV(2)
      EPS_MAX   = STATEV(3)  ! Stores the max strain
      ZONE      = STATEV(4)  ! Allot zone in Traction-Separation Law
C       SLOPE     = STATEV(6)  ! SLOPE OF THE SEGMENT
      NONDIM    = STATEV(9)
      STATEV(10)= L
C-----------------------------------------------------------------------
C     Calculate Stiffness Matrix of Continuum (Sco)      
      CALL STIFFNESS_MATRIX(E, nu, Dco)

C-----------------------------------------------------------------------
C     INITIALIZE STRAIN
      DO K1=1,NTENS
            EPS(K1) = EPS(K1) + DSTRAN(K1)
      END DO

C-----------------------------------------------------------------------      
C     CHECK FOR INITIATION OF FAILURE
      IF (isCRACKED .LT. 0.5) THEN
            SIG = SIG + MATMUL(Dco,DSTRAN)
            IF(isIMPLICIT .EQ.1) DDSDDE = Dco

            ZONE = 1
            STATEV(4) = ZONE
            STATEV(8) = ALF1

C-----------------------------------------------------------------------
C           First time when strain overshoots, CUTBACK to the Peak Point

C           LONGITUDINAL TENSION - FAILURE INDEX
            FI_FT = (SIG(1)/XT)**2

            IF ((SIG(1) .GT. 0.0D0) .AND. (FI_FT .GE. 1.0D0)) THEN 
                  MODE      = 1
                  isCRACKED = 1.0D0
                  PEAK      = SIG(1)

                  ! Computed Stress 'just' before overshoot
                  SIG_COMP   = PEAK - ALF1*DSTRAN(1)
                  STATEV(11) = SIG_COMP
                  
                  ! STRAIN in REGION-2 that corresponds to SIG before overshoot
                  EPS_STAR = EPSN - (XT-SIG_COMP)/ALF2

                  DEPS1       = EPS_STAR-EPS(1)
                  DEPS2       = DSTRAN(1) - DEPS1

                  ! Corrected Stress @ overshoot with ALF2
C                   SIG_CORR    = SIG_COMP - ALF2*(DEPS1)
                  SIG_CORR    = SIG_COMP !+ ALF2*(DEPS2)

                  STATEV(12) = SIG_CORR

            END IF


C-----------------------------------------------------------------------
      END IF

      STATEV(1) = isCRACKED 
      STATEV(2) = MODE

C     STORE THE PEAK STRESS VALUE OF TRACTION SEPARATION LAW
      STATEV(7) = DMAX1(SIG_COMP,STATEV(7))

C     STORE THE MAX STRAIN OBSERVED 
      STATEV(3) = DMAX1(EPS(1),EPS_MAX)

C     BEFORE FAILURE -- LINEAR ELASTIC RESPONSE
      IF (isCRACKED .LT. 0.5) THEN
            CALL STIFFNESS_MATRIX(E,nu, Dco)
            IF(isIMPLICIT .EQ.1) DDSDDE = Dco


      ELSE IF (isCRACKED .GE. 0.5) THEN
            EPS_OLD = EPS_MAX
C             write(*,*) 'STATEV(3): ',STATEV(3) ! STORE MAX STRAIN
C             write(*,*) 'EPS(1): ',EPS(1)       ! CURRENT   STRAIN
C             write(*,*) 'EPS_OLD: ',EPS_OLD     ! PREVIOUS  STRAIN

            IF (MODE .EQ. 1) THEN
                  !BEFORE COMPLETE FAILURE
                  IF (EPS(1) .LE. EPSF) THEN
                        ! POST-PEAK
                        IF (STATEV(3) .GE. XT/E) THEN
                              ! LOADING in POST-PEAK
                              IF (EPS(1) .EQ. STATEV(3)) THEN 
                                    ! POST-PEAK
                                    ZONE = 2
                                    STATEV(4) = ZONE

                                    ! Stress in Zone-2
                                    SIG(1) = STATEV(12) + ALF2*DSTRAN(1)
                                    STATEV(12) = SIG(1)

                                    ! STORE THE STRESS AT THE MAX DISPLACEMENT ENCOUNTERED
                                    temp_sig = SIG(1)
                                    STATEV(5) = temp_sig

                                    ! SECANT STIFFNESS
                                    SECANT = SIG(1)/EPS(1)

                                    STATEV(8) = ALF2

                              END IF
      
                              ! UNLOADING/RELOADING in POST-PEAK
                               IF (EPS(1) .LT. STATEV(3)) THEN
                                     ZONE = 3
                                     STATEV(4) = ZONE
      
                                     ! Slope in Zone-3
                                     SLOPE = STATEV(5)/STATEV(3)                         
                                     STATEV(6) = SLOPE
                  
                                     ! Stress in Zone-3
                                     SIG(1) = STATEV(5) + SLOPE*(EPS(1)-STATEV(3))

                                     ! SECANT STIFFNESS
                                     SECANT = SIG(1)/EPS(1)

                                     STATEV(8) = SECANT

                              END IF
                        END IF
                  ELSE !AFTER COMPLETE FAILURE
                        SIG(1) = 0.0d0
                  END IF

            END IF

            IF(isIMPLICIT .EQ.1) THEN
                  DDSDDE      = Dco 
                  DDSDDE(1,1) = STATEV(8)
C                   DDSDDE(1,2) = 0.0d0
C                   DDSDDE(1,3) = 0.0d0
C                   DDSDDE(2,1) = 0.0d0
C                   DDSDDE(3,1) = 0.0d0
            END IF
      END IF

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

      
      
