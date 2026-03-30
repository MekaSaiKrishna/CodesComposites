C***********************************************************************
C  claude_AFEM2D_PEv3.for
C  SCA Isotropic Incremental UMAT  --  2D Plane Strain
C  Restructured from v2: calcEcr_PE eliminated (inlined); Dcr called
C  exactly twice per cracked increment (pre-solve / post-update).
C  STATEV(9:10) and STATEV(14:15) are written in the same pass from
C  the same Ecr value, so the one-increment lag cannot occur.
C  Last Modified: 03/30/2026  |  meka1@purdue.edu
C***********************************************************************
      SUBROUTINE sca2diso_PE(NTENS, NDI, NSHR, DFGRD1, NSTATV, NPROPS,
     1      NOEL, stepTIME, totalTIME, DT, PNEWDT, L, PROPS, EPS,
     2      DSTRAN, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)

      IMPLICIT NONE

      INTEGER  :: NTENS, NDI, NSHR, NSTATV, NPROPS, NOEL
      REAL*8   :: DFGRD1(3,3)
      REAL*8   :: stepTIME, totalTIME, DT, PNEWDT, L, Leff
      REAL*8   :: PROPS(NPROPS)
      REAL*8   :: SIG(NTENS), EPS(NTENS), DSTRAN(NTENS)
      REAL*8   :: STATEV(NSTATV)
      REAL*8   :: DDSDDE(NTENS,NTENS)
      REAL*8   :: DTEMP
      INTEGER  :: isIMPLICIT

      INTEGER  :: INDEX, I, K1
      REAL*8   :: E, nu, GIc, XT, G, eta
      REAL*8   :: sigcr0, taucr0, sigcr0_r, BL, MODE
      REAL*8   :: FI_T, FI_C, FI_MAX, J_FAIL, J_MIN, PMIN, J_MAX, PMAX
      REAL*8   :: PRINS(3), T(3,3), PRINS_prev(3), T_prev(3,3)
      REAL*8   :: SIG_prev(NTENS)
      REAL*8   :: FI_prev, alpha_init, alpha
      REAL*8   :: Dco(NTENS,NTENS), Dcocr(NTENS,NTENS), Dcr(2,2)
      REAL*8   :: N(4,2)
      REAL*8   :: Ecr_MAX(2), Ecr(2), Ecr_save(2), dEcr(2)
      REAL*8   :: scr_prev(2), scr(2)
      REAL*8   :: NT_Dco(2,4), NT_Dco_N(2,2)
      REAL*8   :: Y(2,2), invY(2,2), detY
      REAL*8   :: Ecr_n_ult

      REAL*8, PARAMETER :: ZERO     = 0.0D0
      REAL*8, PARAMETER :: ONE      = 1.0D0
      REAL*8, PARAMETER :: TWO      = 2.0D0
      REAL*8, PARAMETER :: FI_CRIT  = 1.0D0
      REAL*8, PARAMETER :: R_RES    = 1.0D-8
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20

C  ===================================================================
C  PROPS layout  (NPROPS >= 5, PROPS(6) optional)
C    PROPS(1) = E      Young's modulus                        [MPa]
C    PROPS(2) = nu     Poisson's ratio                        [-]
C    PROPS(3) = GIc    Mode-I fracture energy                 [N/mm]
C    PROPS(4) = XT     Tensile initiation strength            [MPa]
C    PROPS(5) = INDEX  Failure criterion index (integer)      [-]
C    PROPS(6) = eta    Viscous regularization coefficient     [MPa.s]
C               (optional; default 0; mirrors ABAQUS cohesive viscosity)
C  ===================================================================
C  STATEV layout  (NSTATV >= 17)
C    STATEV(1)  = isCRACKED   crack flag: 0 = intact, 1 = active
C    STATEV(2)  = MODE        crack mode: +1 = tensile, -1 = compressive
C    STATEV(3)  = Leff        effective element size stored for output [mm]
C    STATEV(4)  = FI_MAX      peak failure index reached (pre-peak)    [-]
C    STATEV(5)  = cos(theta)  crack-normal cosine (global x-axis)      [-]
C    STATEV(6)  = sin(theta)  crack-normal sine   (global x-axis)      [-]
C    STATEV(7)  = Ecr_MAX(1)  max normal  crack strain reached         [m/m]
C    STATEV(8)  = Ecr_MAX(2)  max shear   crack strain reached         [m/m]
C    STATEV(9)  = Ecr(1)      current normal  crack strain             [m/m]
C    STATEV(10) = Ecr(2)      current shear   crack strain             [m/m]
C    STATEV(11) = dEcr(1)     normal  crack strain increment this step [m/m]
C    STATEV(12) = dEcr(2)     shear   crack strain increment this step [m/m]
C    STATEV(13) = sig_init    max principal stress at initiation       [MPa]
C    STATEV(14) = scr(1)      current normal  crack traction           [MPa]
C    STATEV(15) = scr(2)      current shear   crack traction           [MPa]
C    STATEV(16) = sigcr0      normal  traction at crack initiation     [MPa]
C    STATEV(17) = taucr0      shear   traction at crack initiation     [MPa]
C  ===================================================================

      E     = PROPS(1)
      nu    = PROPS(2)
      GIc   = PROPS(3)
      XT    = PROPS(4)
      INDEX = INT(PROPS(5))
      eta   = ZERO
      IF (NPROPS .GE. 6) eta = PROPS(6)

      G    = E / (TWO*(ONE + nu))
      Leff = ONE*L
C     Leff = 0.02D0   ! uncomment to fix element size for DCB models

      MODE      = STATEV(2)
      sigcr0    = STATEV(16)
      taucr0    = STATEV(17)
      STATEV(3) = Leff

      CALL STIFFNESS_MATRIX_2D_PE(E, nu, Dco)

      DO K1 = 1, NTENS
         EPS(K1) = EPS(K1) + DSTRAN(K1)
      END DO

      IF (sigcr0 .GT. 1.0D-12) THEN
         BL = TWO*GIc*E / (sigcr0**2)
         IF (Leff .GT. BL) THEN
            WRITE(*,'(A)')          ' *** BAZANT LIMIT EXCEEDED ***'
            WRITE(*,'(A,1PG12.5)')  '   2GIcE/sigcr0^2 = ', BL
            WRITE(*,'(A,1PG12.5)')  '   Leff           = ', Leff
            CALL myExit_PE()
         END IF
      END IF

C  ------------------------------------------------------------------
C  PRE-PEAK BLOCK
C  ------------------------------------------------------------------
      IF (STATEV(1) .LT. 0.5D0) THEN

         STATEV(4)  = ZERO
         STATEV(16) = ZERO
         STATEV(17) = ZERO
         SIG_prev   = SIG

         SIG    = SIG + MATMUL(Dco, DSTRAN)
         DDSDDE = Dco

         CALL mysprind_PE(SIG, PRINS, T, 1, 3, 1)

         FI_T   = ZERO
         FI_C   = ZERO
         FI_MAX = ZERO

         PMAX  = MAX(PRINS(1), PRINS(2), PRINS(3))
         J_MAX = 1
         DO I = 1, 3
            IF (PRINS(I) .EQ. PMAX) J_MAX = I
         END DO
         IF (PMAX .GT. ZERO) FI_T = PMAX / XT

         PMIN  = MIN(PRINS(1), PRINS(2), PRINS(3))
         J_MIN = 1
         DO I = 1, 3
            IF (PRINS(I) .EQ. PMIN) J_MIN = I
         END DO
         IF (PMIN .LT. ZERO) FI_C = ABS(PMIN) / XT

         IF (FI_T .GE. FI_C) THEN
            MODE   =  ONE
            FI_MAX =  FI_T
            J_FAIL =  DBLE(J_MAX)
         ELSE
            MODE   = -ONE
            FI_MAX =  FI_C
            J_FAIL =  DBLE(J_MIN)
         END IF

         STATEV(2) = MODE
         STATEV(4) = FI_MAX
         STATEV(5) = T(1, INT(J_FAIL))
         STATEV(6) = T(2, INT(J_FAIL))

         IF (FI_MAX .GE. FI_CRIT) THEN

            IF (FI_MAX .GT. 1.05D0*FI_CRIT) THEN
               CALL mysprind_PE(SIG_prev, PRINS_prev, T_prev, 1, 3, 1)
               IF (MODE .GT. ZERO) THEN
                  FI_prev = MAX(PRINS_prev(1),PRINS_prev(2),
     1                         PRINS_prev(3)) / XT
               ELSE
                  FI_prev = ABS(MIN(PRINS_prev(1),PRINS_prev(2),
     1                             PRINS_prev(3))) / XT
               END IF
               FI_prev = MAX(FI_prev, ZERO)
               IF ((FI_MAX - FI_prev) .GT. 1.0D-12) THEN
                  alpha_init = (FI_CRIT - FI_prev) / (FI_MAX - FI_prev)
                  alpha_init = MAX(MIN(alpha_init, 0.95D0), 0.10D0)
                  PNEWDT = MIN(PNEWDT, alpha_init)
               END IF
               RETURN
            END IF

            STATEV(1)     = ONE
            CALL getN_PE(STATEV(5), STATEV(6), N)
            STATEV(16:17) = MATMUL(TRANSPOSE(N), SIG)
            STATEV(13)    = STATEV(4) * XT
            STATEV(14:15) = STATEV(16:17)
            STATEV(7)     = ZERO
            STATEV(8)     = ZERO
            STATEV(9)     = ZERO
            STATEV(10)    = ZERO
            STATEV(11)    = ZERO
            STATEV(12)    = ZERO

         END IF

      END IF

      sigcr0 = STATEV(16)
      taucr0 = STATEV(17)

C  ------------------------------------------------------------------
C  POST-PEAK BLOCK
C  ------------------------------------------------------------------
      IF (STATEV(1) .GE. ONE) THEN

         Ecr_MAX  = STATEV(7:8)
         Ecr      = STATEV(9:10)
         dEcr     = STATEV(11:12)
         scr_prev = STATEV(14:15)

         CALL getN_PE(STATEV(5), STATEV(6), N)

         sigcr0_r  = sigcr0 * R_RES
         Ecr_n_ult = TWO*GIc / (Leff*(sigcr0 - sigcr0_r))

C        Step A: Dcr at current (old) Ecr — needed to form system matrix Y.
C        scr output from this call is intentionally discarded; only Dcr used.
         CALL calcDcr_PE(Ecr, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                   nu, eta, DT, Dcr, scr_prev, scr, dEcr)

C        Step B: Solve  (Dcr + N^T Dco N) dEcr = N^T Dco DSTRAN
         NT_Dco   = MATMUL(TRANSPOSE(N), Dco)
         NT_Dco_N = MATMUL(NT_Dco, N)
         Y   = Dcr + NT_Dco_N
         detY = Y(1,1)*Y(2,2) - Y(1,2)*Y(2,1)
         IF (ABS(detY) .LT.
     1       SING_TOL*(ABS(Y(1,1))*ABS(Y(2,2)) + ONE)) THEN
            WRITE(*,'(A,1PG12.4)') ' WARN: Y near-singular, det=', detY
            PNEWDT = MIN(PNEWDT, 0.5D0)
            RETURN
         END IF
         CALL matrixInverse(Y, invY, 2, 2)
         dEcr = MATMUL(invY, MATMUL(NT_Dco, DSTRAN))

C        Step C: Failure overshoot — cut step if Ecr would exceed
C        1.05*Ecr_n_ult. Alpha floor of 0.10 prevents infinite loops.
         IF (Ecr_MAX(1) .LT. Ecr_n_ult .AND. dEcr(1) .GT. ZERO) THEN
            IF ((Ecr(1) + dEcr(1)) .GT. 1.05D0*Ecr_n_ult) THEN
               alpha = MAX((Ecr_n_ult - Ecr(1)) / dEcr(1), 0.10D0)
               PNEWDT = MIN(PNEWDT, alpha)
            END IF
         END IF

C        Step D: Update Ecr; save pre-update value for clamp below.
         Ecr_save   = Ecr
         Ecr        = Ecr + dEcr
         Ecr_MAX(1) = MAX(Ecr(1), Ecr_MAX(1))
         Ecr_MAX(2) = MAX(Ecr(2), Ecr_MAX(2))

C        Step E: Interpenetration clamp (tensile mode only).
C        dEcr(1) set to the exact delta that closes the crack to zero.
         IF (MODE .GT. ZERO .AND. Ecr(1) .LT. ZERO) THEN
            dEcr(1) = -Ecr_save(1)
            Ecr(1)  =  ZERO
         END IF

C        Step F: Dcr and scr at NEW Ecr.
C        STATEV(9:10) and STATEV(14:15) will both be written from this
C        same Ecr in Step H — the one-increment lag is impossible.
         CALL calcDcr_PE(Ecr, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                   nu, eta, DT, Dcr, scr_prev, scr, dEcr)

C        Step G: Cracked consistent tangent D_cocr.
         CALL calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)

C        Step H: Write STATEV — all crack variables from the same state.
         STATEV(7:8)   = Ecr_MAX
         STATEV(9:10)  = Ecr
         STATEV(11:12) = dEcr
         STATEV(14:15) = scr

C        Step I: Stress update.
         SIG = SIG + MATMUL(Dcocr, DSTRAN)
         IF (isIMPLICIT .EQ. 1) DDSDDE = Dcocr

      END IF

      RETURN
      END SUBROUTINE sca2diso_PE

C***********************************************************************
C  STIFFNESS_MATRIX_2D_PE  —  plane-strain isotropic D_co  (4x4)
C***********************************************************************
      SUBROUTINE STIFFNESS_MATRIX_2D_PE(E, nu, Dco)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: E, nu
      REAL*8, INTENT(OUT) :: Dco(4,4)
      REAL*8 :: c, G
      c = E / ((1.0D0 - 2.0D0*nu)*(1.0D0 + nu))
      G = E / (2.0D0*(1.0D0 + nu))
      Dco      = 0.0D0
      Dco(1,1) = c*(1.0D0 - nu)
      Dco(1,2) = c*nu
      Dco(1,3) = c*nu
      Dco(2,1) = c*nu
      Dco(2,2) = c*(1.0D0 - nu)
      Dco(2,3) = c*nu
      Dco(3,1) = c*nu
      Dco(3,2) = c*nu
      Dco(3,3) = c*(1.0D0 - nu)
      Dco(4,4) = G
      RETURN
      END SUBROUTINE STIFFNESS_MATRIX_2D_PE

C***********************************************************************
C  getN_PE  —  crack transformation matrix N  (4x2)
C  Maps local crack strains [eps_cr_n, gam_cr_t] to global Voigt space.
C***********************************************************************
      SUBROUTINE getN_PE(c, s, N)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: c, s
      REAL*8, INTENT(OUT) :: N(4,2)
      N(1,1) =  c*c
      N(1,2) = -c*s
      N(2,1) =  s*s
      N(2,2) =  c*s
      N(3,1) =  0.0D0
      N(3,2) =  0.0D0
      N(4,1) =  2.0D0*s*c
      N(4,2) =  c*c - s*s
      RETURN
      END SUBROUTINE getN_PE

C***********************************************************************
C  mysprind_PE  —  thin wrapper around ABAQUS SPRIND
C***********************************************************************
      SUBROUTINE mysprind_PE(S, PS, AN, LSTR, NDI, NSHR)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NDI, NSHR, LSTR
      REAL*8,  INTENT(IN)  :: S(NDI+NSHR)
      REAL*8,  INTENT(OUT) :: PS(3), AN(3,3)
      CALL SPRIND(S, PS, AN, LSTR, NDI, NSHR)
      RETURN
      END SUBROUTINE mysprind_PE

C***********************************************************************
C  myExit_PE  —  thin wrapper around ABAQUS XIT
C***********************************************************************
      SUBROUTINE myExit_PE()
      CALL XIT
      RETURN
      END SUBROUTINE myExit_PE

C***********************************************************************
C  calcDcr_PE  —  crack constitutive matrix D_cr (2x2) and traction scr
C
C  Regime detection (based on Ecr_eff = max(|Ecr|, |Ecr_MAX|)):
C    Ecr_eff < Ecr_n_ult :
C      Ecr(1) < Ecr_MAX(1)  ->  SECANT   : total-form, D_sec = scr_max/Ecr_MAX
C      Ecr(1) = Ecr_MAX(1)  ->  SOFTENING: total-form + viscous correction
C    Ecr_eff >= Ecr_n_ult   ->  RESIDUAL : D = RES_K, scr = sigcr0_r
C
C  Shear: secant total-form when in secant branch; incremental otherwise.
C  Viscous: eta/DT added to softening slope (mirrors ABAQUS cohesive).
C
C  Arguments:
C    Ecr     (IN)    crack strain at which D_cr is evaluated
C    Ecr_MAX (IN)    maximum crack strain history
C    scr_old (IN)    crack traction from previous increment (shear only)
C    dEcr    (IN)    crack strain increment (for viscous term)
C    Dcr     (OUT)   2x2 crack stiffness matrix
C    scr     (OUT)   2-component crack traction vector
C***********************************************************************
      SUBROUTINE calcDcr_PE(Ecr, Ecr_MAX, GIc, sigcr0, taucr0,
     1           Leff, nu, eta, DT, Dcr, scr_old, scr, dEcr)
      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: Ecr(2), Ecr_MAX(2)
      REAL*8, INTENT(IN)  :: GIc, sigcr0, taucr0, Leff, nu, eta, DT
      REAL*8, INTENT(IN)  :: scr_old(2), dEcr(2)
      REAL*8, INTENT(OUT) :: Dcr(2,2), scr(2)

      REAL*8 :: sigcr0_r, term, slope_soft, visc_term
      REAL*8 :: Ecr_eff(2), scr_max, Ecr_n_ult
      REAL*8, PARAMETER :: R_RES   = 1.0D-8
      REAL*8, PARAMETER :: RES_K   = 1.0D-8
      REAL*8, PARAMETER :: ECR_MIN = 1.0D-14

      IF (sigcr0 .LT. 1.0D-12) THEN
         Dcr      = 0.0D0
         Dcr(1,1) = RES_K
         Dcr(2,2) = RES_K
         scr      = scr_old
         RETURN
      END IF

      sigcr0_r   = sigcr0 * R_RES
      term       = 1.0D0 / (2.0D0*(1.0D0 + nu))
      slope_soft = -(sigcr0**2) / (2.0D0*GIc/Leff)
      Ecr_n_ult  = 2.0D0*GIc / (Leff*(sigcr0 - sigcr0_r))
      visc_term  = 0.0D0
      IF (eta .GT. 0.0D0 .AND. DT .GT. 0.0D0) visc_term = eta / DT

      scr_max    = MAX(sigcr0 + slope_soft*Ecr_MAX(1), 0.0D0)
      Ecr_eff(1) = MAX(ABS(Ecr(1)), ABS(Ecr_MAX(1)))
      Ecr_eff(2) = MAX(ABS(Ecr(2)), ABS(Ecr_MAX(2)))

      Dcr = 0.0D0

      IF (Ecr_eff(1) .LT. Ecr_n_ult) THEN
         IF (Ecr(1) .LT. Ecr_MAX(1)) THEN
            IF (ABS(Ecr_MAX(1)) .GT. ECR_MIN) THEN
               Dcr(1,1) = scr_max / Ecr_MAX(1)
            ELSE
               Dcr(1,1) = slope_soft + visc_term
            END IF
            scr(1) = Dcr(1,1) * Ecr(1)
         ELSE
            Dcr(1,1) = slope_soft + visc_term
            scr(1)   = sigcr0 + slope_soft*Ecr(1) + visc_term*dEcr(1)
            scr(1)   = MAX(scr(1), sigcr0_r)
         END IF
      ELSE
         Dcr(1,1) = RES_K
         scr(1)   = sigcr0_r
      END IF

      IF (Ecr_eff(1) .LT. Ecr_n_ult) THEN
         Dcr(2,2) = Dcr(1,1) * term * 1.0D-3
      ELSE
         Dcr(2,2) = RES_K * term
      END IF

      IF (Ecr(1) .LT. Ecr_MAX(1) .AND.
     1    Ecr_eff(1) .LT. Ecr_n_ult) THEN
         scr(2) = Dcr(2,2) * Ecr(2)
      ELSE
         scr(2) = scr_old(2) + Dcr(2,2)*dEcr(2)
      END IF

      RETURN
      END SUBROUTINE calcDcr_PE

C***********************************************************************
C  calcDcocr_PE  —  cracked-continuum tangent D_cocr (4x4)
C
C  D_cocr = D_co - D_co N (D_cr + N^T D_co N)^{-1} N^T D_co
C
C  Near-singular fallback: returns D_co and requests PNEWDT = 0.5.
C***********************************************************************
      SUBROUTINE calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)
      IMPLICIT NONE
      REAL*8, INTENT(IN)    :: Dcr(2,2), Dco(4,4), N(4,2)
      REAL*8, INTENT(OUT)   :: Dcocr(4,4)
      REAL*8, INTENT(INOUT) :: PNEWDT
      REAL*8 :: Y22(2,2), invY22(2,2), A42(4,2), A44(4,4), detY
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20

      Y22  = Dcr + MATMUL(TRANSPOSE(N), MATMUL(Dco, N))
      detY = Y22(1,1)*Y22(2,2) - Y22(1,2)*Y22(2,1)

      IF (ABS(detY) .LT.
     1    SING_TOL*(ABS(Y22(1,1))*ABS(Y22(2,2)) + 1.0D0)) THEN
         WRITE(*,'(A,1PG12.4)')
     1       ' WARN: Dcocr near-singular, det=', detY
         Dcocr  = Dco
         PNEWDT = MIN(PNEWDT, 0.5D0)
         RETURN
      END IF

      CALL matrixInverse(Y22, invY22, 2, 2)
      A42   = MATMUL(Dco, MATMUL(N, invY22))
      A44   = MATMUL(A42, MATMUL(TRANSPOSE(N), Dco))
      Dcocr = Dco - A44

      RETURN
      END SUBROUTINE calcDcocr_PE
