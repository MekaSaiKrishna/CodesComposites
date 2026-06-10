C***********************************************************************
C     ! Isotropic TOTAL SCA 3D 
C     ! ----------------------------------------------------------------
C     ! Contains the code concerned with the formulation for SCA in 
C     ! orthotropic composite material, this is used in the 
C     ! 'main_UD_tot.for' file
C     ! ----------------------------------------------------------------
C     ! Last Debugged: 02/26/2025 
C     ! meka1@purdue.edu
C***********************************************************************
      SUBROUTINE sca3dortho(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)
          
      IMPLICIT NONE
      
C     STANDARD VARIABLES
      INTEGER::NSTATV,NPROPS,NOEL
      REAL*8::stepTIME,totalTIME,DT,L
      REAL*8::PROPS,EPS(6),SIG(6)
      REAL*8::STATEV
      REAL*8::DDSDDE(6,6)      
      REAL*8::DTEMP
      INTEGER::isIMPLICIT
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)

C     DEFINED VARIABLES
      REAL*8::Dco(6,6),Sco(6,6)
      REAL*8::G1C_F
      REAL*8::E,nu,XT
      REAL*8::isCRACKED, MODE
      INTEGER::counter
      REAL*8::EPS_MAX,EPS_OLD
      REAL*8::FI_FT
      REAL*8::SLOPE
      REAL*8::temp_sig
      REAL*8::PEAK, TEMP
      INTEGER::ZONE

C     DEFINED FUNCTIONS


C***********************************************************************
C       !State Variables List
C-----------------------------------------------------------------------
C       (1) Crack Flag
C       (2) Increment Counter
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

C-----------------------------------------------------------------------
      isCRACKED = STATEV(1)  ! CRACK FLAG
      MODE      = STATEV(2)
      EPS_MAX   = STATEV(3)  ! Stores the max strain
      ZONE      = STATEV(4)  ! Allot zone in Traction-Separation Law
C       SLOPE     = STATEV(6)  ! SLOPE OF THE SEGMENT
C-----------------------------------------------------------------------
C     Calculate Stiffness Matrix of Continuum (Sco)      
      CALL STIFFNESS_MATRIX(E, nu, Dco)
      
C     CHECK FOR INITIATION OF FAILURE
      IF (isCRACKED .LT. 0.5) THEN
            SIG = MATMUL(Dco,EPS)
            IF(isIMPLICIT .EQ.1) DDSDDE = Dco

            ZONE = 1
            STATEV(4) = ZONE

C           LONGITUDINAL TENSION - FAILURE INDEX
            FI_FT = (SIG(1)/XT)**2

            IF ((SIG(1) .GT. 0.0D0) .AND. (FI_FT .GE. 1.0D0)) THEN 
                  MODE      = 1
                  isCRACKED = 1.0D0
C                   TEMP      = SIG(1)
                  PEAK      = SIG(1)
            END IF

            ! STORE THE PEAK STRESS VALUE
            IF (TEMP .GT. 0.0d0) THEN
                  PEAK = TEMP
            END IF
      END IF

      STATEV(1) = isCRACKED 
      STATEV(2) = MODE

C     STORE THE PEAK STRESS VALUE OF TRACTION SEPARATION LAW
      STATEV(7) = DMAX1(PEAK,STATEV(7))

C     STORE THE MAX STRAIN OBSERVED 
      STATEV(3) = DMAX1(EPS(1),EPS_MAX)

C     BEFORE FAILURE -- LINEAR ELASTIC RESPONSE
      IF (isCRACKED .LT. 0.5) THEN
            CALL STIFFNESS_MATRIX(E,nu, Dco)
            SIG = MATMUL(Dco,EPS)

      ELSE IF (isCRACKED .GE. 0.5) THEN
            EPS_OLD = EPS_MAX
C             write(*,*) 'STATEV(3): ',STATEV(3) ! STORE MAX STRAIN
C             write(*,*) 'EPS(1): ',EPS(1)       ! CURRENT   STRAIN
C             write(*,*) 'EPS_OLD: ',EPS_OLD     ! PREVIOUS  STRAIN
            IF (MODE .EQ. 1) THEN
                  !BEFORE COMPLETE FAILURE
                  IF (EPS(1) .LE. (2.0d0*(G1C_F/L)/XT)) THEN
                        ! POST-PEAK
                        IF (STATEV(3) .GE. XT/E) THEN
                              ! LOADING in POST-PEAK
                              IF (EPS(1) .EQ. STATEV(3)) THEN 
                                    ! POST-PEAK
                                    ZONE = 2
                                    STATEV(4) = ZONE

                                    PEAK = STATEV(7)
                                    ! Slope in Zone-2
                                    SLOPE = -1.0d0*PEAK/((2.0d0*(G1C_F/L)/PEAK)-(PEAK/E))

                                    STATEV(6) = SLOPE

                                    ! Stress in Zone-2
                                    SIG(1) = PEAK + SLOPE*(EPS(1)-(PEAK/E))

                                    ! STORE THE STRESS AT THE MAX DISPLACEMENT ENCOUNTERED
                                    temp_sig = SIG(1)
                                    STATEV(5) = temp_sig

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
                              END IF
                        END IF
                  ELSE !AFTER COMPLETE FAILURE
                        SIG(1) = 0.0d0
                  END IF
            END IF

            IF(isIMPLICIT .EQ.1) THEN
                  DDSDDE      = Dco 
                  DDSDDE(1,2) = 0.0d0
                  DDSDDE(1,3) = 0.0d0
                  DDSDDE(2,1) = 0.0d0
                  DDSDDE(3,1) = 0.0d0
            END IF
      END IF

C       IF(isIMPLICIT .EQ.1) DDSDDE = Dco

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

      
      
