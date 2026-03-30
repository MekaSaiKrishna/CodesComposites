C***********************************************************************
C     ! Isotropic Incremental SCA 2D Plane Strain  [Enhanced v2]
C     ! ----------------------------------------------------------------
C     ! Based on: claude_AFEM2D_PEv1.for (03/11/2026)
C     !
C     ! Changes from v1 (all targeted bug fixes, no behavioural rewrites):
C     !
C     !  FIX-A [sca2diso_PE, local declarations]:
C     !       Removed unused local variables J (INTEGER) and
C     !       Sco(NTENS,NTENS) (REAL*8). Neither was ever referenced.
C     !       Reason: dead declarations can confuse compilers and readers;
C     !       with IMPLICIT NONE they serve no purpose.
C     !
C     !  FIX-B [sca2diso_PE, pre-peak init block]:
C     !       Added  STATEV(4) = ZERO  to the pre-peak initialisation.
C     !       Was: STATEV(16) and STATEV(17) were zeroed but STATEV(4)
C     !       (FI_MAX history) was not, leaving it with the value from
C     !       the previous converged increment when ABAQUS re-enters the
C     !       pre-peak block after a step cut.
C     !       Reason: although STATEV(4) is always overwritten before use
C     !       in the normal path, the initiation-overshoot RETURN exits
C     !       before the  STATEV(4) = FI_MAX  line at crack activation,
C     !       so stale FI_MAX can persist across retried increments.
C     !
C     !  FIX-C [calcEcr_PE, OS-2new block]:
C     !       Removed the  IF (alpha .LE. 0.10D0) PNEWDT = 1.0d0
C     !       override that followed the MIN(PNEWDT, alpha) assignment.
C     !       Was: when alpha_raw < 0.10, alpha was clamped to 0.10 by
C     !       MAX(alpha, 0.10D0), then PNEWDT = MIN(PNEWDT, 0.10) was
C     !       immediately overridden by the hard assignment PNEWDT = 1.0.
C     !       This completely disabled step control precisely when the
C     !       step was overshooting Ecr_n_ult by the largest margin.
C     !       Reason: the intent (avoid infinite cut loops) is already
C     !       served by the 0.10 lower bound on alpha; the override is
C     !       therefore both redundant and dangerous. It also violates
C     !       the "always MIN, never assign" PNEWDT rule and can undo
C     !       a valid PNEWDT reduction made earlier in the same call.
C     !
C     !  FIX-D [calcDcocr_PE]:
C     !       Restored PNEWDT as INTENT(INOUT) argument. Near-singular
C     !       tangent now sets PNEWDT = MIN(PNEWDT, 0.5D0) and prints a
C     !       warning before falling back to Dco.
C     !       Was: PNEWDT was removed from v1 with the comment that
C     !       "cutting the step here is expensive and rarely helps."
C     !       Reason: while the elastic-tangent fallback is correct, a
C     !       silent recovery gives ABAQUS no signal that something went
C     !       wrong. The step cut nudges Newton-Raphson toward a state
C     !       where the tangent is better conditioned and leaves a trace
C     !       in the ABAQUS message file. The call site in sca2diso_PE
C     !       is updated to pass PNEWDT.
C     !
C     ! All other code is identical to claude_AFEM2D_PEv1.for.
C     ! Changed lines are marked  !<Fv2-X>  for traceability.
C     ! ----------------------------------------------------------------
C     ! Last Modified: 03/30/2026
C     ! meka1@purdue.edu
C***********************************************************************
      SUBROUTINE sca2diso_PE(NTENS, NDI, NSHR, DFGRD1, NSTATV, NPROPS,
     1      NOEL, stepTIME, totalTIME, DT, PNEWDT, L, PROPS, EPS,
     2      DSTRAN, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)

      IMPLICIT NONE

C     --- ABAQUS standard arguments ---
      INTEGER     :: NTENS, NDI, NSHR, NSTATV, NPROPS, NOEL
      REAL*8      :: DFGRD1(3,3)
      REAL*8      :: stepTIME, totalTIME, DT, PNEWDT, L, Leff
      REAL*8      :: PROPS(NPROPS)
      REAL*8      :: SIG(NTENS), EPS(NTENS), DSTRAN(NTENS)
      REAL*8      :: STATEV(NSTATV)
      REAL*8      :: DDSDDE(NTENS,NTENS)
      REAL*8      :: DTEMP
      INTEGER     :: isIMPLICIT

C     --- Local variables ---
      INTEGER     :: INDEX, I, K1                               !<Fv2-A> removed unused J
      REAL*8      :: E, nu, GIc, XT, G, eta
      REAL*8      :: sigcr0, taucr0, sigcr0_r, BL
      REAL*8      :: isCRACKED, MODE
      REAL*8      :: FI_T, FI_C, FI_MAX, J_FAIL, J_MIN, PMIN, J_MAX, PMAX
      REAL*8      :: PRINS(3), T(3,3)
      REAL*8      :: Dco(NTENS,NTENS)                           !<Fv2-A> removed unused Sco

      REAL*8      :: N(4,2)
      REAL*8      :: Ecr_MAX(2), Ecr_OLD(2)
      REAL*8      :: Ecr(2), dEcr(2)
      REAL*8      :: scr_old(2), scr(2)
      REAL*8      :: Dcocr(NTENS,NTENS)
      REAL*8      :: Dcr(2,2)
      REAL*8      :: scr_start(2)  ! snapshot of scr_old for Fix-GAP second Dcr call
      REAL*8      :: SIG_prev(NTENS)     ! stress at start of increment
      REAL*8      :: PRINS_prev(3), T_prev(3,3)
      REAL*8      :: FI_prev, alpha_init

C     --- Numerical constants ---
      REAL*8, PARAMETER :: ZERO    = 0.0D0
      REAL*8, PARAMETER :: ONE     = 1.0D0
      REAL*8, PARAMETER :: TWO     = 2.0D0
      REAL*8, PARAMETER :: TINY_   = 1.0D-8 ! seed value for crack strains
      REAL*8, PARAMETER :: RES_    = 1.0D-5 ! residual stiffness after full failure
      REAL*8, PARAMETER :: FI_CRIT = 1.0D0  ! failure index threshold

C***********************************************************************
C     READ MATERIAL PROPERTIES
C***********************************************************************
      E      = PROPS(1)      ! Young's modulus
      nu     = PROPS(2)      ! Poisson's ratio
      GIc    = PROPS(3)      ! Mode I fracture toughness
      XT     = PROPS(4)      ! Tensile (and compressive) strength
      INDEX  = INT(PROPS(5)) ! Failure criteria index
C     PROPS(6) = viscous regularization parameter (optional, default 0)
C     Typical value: 1e-5 to 1e-3 * characteristic_time
C     Mirrors ABAQUS cohesive element viscosity parameter.
      IF (NPROPS .GE. 6) THEN
         eta = PROPS(6)
      ELSE
         eta = 0.0D0
      END IF

      G      = E/(TWO*(ONE + nu))
      Leff   = ONE*L   ! characteristic element size (CELENT for CPE4)
C       Leff   = 0.02  ! override for fixed-mesh DCB model (toggle if needed)

      MODE      = STATEV(2)
      sigcr0    = STATEV(16)
      taucr0    = STATEV(17)
      STATEV(3) = Leff

C***********************************************************************
C     STEP 1: Elastic stiffness matrix D_co
C***********************************************************************
      CALL STIFFNESS_MATRIX_2D_PE(E, nu, Dco)

C***********************************************************************
C     STEP 2: Update total strain tensor
C***********************************************************************
      DO K1 = 1, NTENS
         EPS(K1) = EPS(K1) + DSTRAN(K1)
      END DO

C***********************************************************************
C     STEP 3: Bazant size limit check (guard against sigcr0 = 0)
C***********************************************************************
      IF (sigcr0 .GT. 1.0D-12) THEN
         BL = TWO*GIc*E / (sigcr0**2)
         IF (Leff .GT. BL) THEN
            WRITE(*,'(A)')         ' '
            WRITE(*,'(A)')         ' *** BAZANT LIMIT EXCEEDED ***'
            WRITE(*,'(A,1PG12.5)') '  2*GIc*E/sigcr0^2 = ', BL
            WRITE(*,'(A,1PG12.5)') '  Current Leff     = ', Leff
            WRITE(*,'(A,1PG12.5)') '  Recommended Leff <= ', BL
            CALL myExit_PE()
         END IF
      END IF

C=======================================================================
C     PRE-PEAK BLOCK  (no crack yet)
C=======================================================================
      IF (STATEV(1) .LT. 0.5D0) THEN

C        Initialise local crack stresses (not yet initiated)
         STATEV(16) = ZERO
         STATEV(17) = ZERO
         STATEV(4)  = ZERO   !<Fv2-B> clear FI_MAX history each pre-peak call

C        Save stress at start of increment for overshoot detection
         SIG_prev = SIG

C        --- Linear elastic stress update ---
         SIG    = SIG + MATMUL(Dco, DSTRAN)
         DDSDDE = Dco

C        --- Principal stress computation ---
         CALL mysprind_PE(SIG, PRINS, T, 1, 3, 1)

C        Failure indices
         FI_T   = ZERO
         FI_C   = ZERO
         FI_MAX = ZERO

C        Tensile check: most positive principal stress
         PMAX = MAX(PRINS(1), PRINS(2), PRINS(3))
         IF (PMAX .GT. ZERO) FI_T = PMAX / XT
         J_MAX = 1
         DO I = 1, 3
            IF (PRINS(I) .EQ. PMAX) J_MAX = I
         END DO

C        Compressive check: most negative (magnitude) principal stress
         PMIN = MIN(PRINS(1), PRINS(2), PRINS(3))
         J_MIN = 1
         DO I = 1, 3
            IF (PRINS(I) .EQ. PMIN) J_MIN = I
         END DO
         IF (PMIN .LT. ZERO) FI_C = ABS(PMIN) / XT

C        Dominant mode
         IF (FI_T .GE. FI_C) THEN
            MODE   =  ONE
            FI_MAX = FI_T
            J_FAIL = DBLE(J_MAX)
         ELSE
            MODE   = -ONE
            FI_MAX = FI_C
            J_FAIL = DBLE(J_MIN)
         END IF

C        Store failure mode and crack orientation
         STATEV(2) = MODE
         STATEV(4) = FI_MAX
         STATEV(5) = T(1, INT(J_FAIL))
         STATEV(6) = T(2, INT(J_FAIL))

C        ---------------------------------------------------------------
C        INITIATION OVERSHOOT PREVENTION
C        If the elastic predictor overshoots the failure surface by more
C        than 5%, compute the fractional step alpha at which FI = 1.0
C        and request ABAQUS to retry with PNEWDT = alpha. This ensures
C        the stress at crack onset is within 5% of XT, preventing a
C        large jump in crack strain at the first post-peak increment.
C        ---------------------------------------------------------------
         IF (FI_MAX .GE. FI_CRIT) THEN
            IF (FI_MAX .GT. 1.05D0 * FI_CRIT) THEN
C              Compute FI at start of increment (before DSTRAN)
               CALL mysprind_PE(SIG_prev, PRINS_prev, T_prev, 1, 3, 1)
               IF (MODE .GT. ZERO) THEN
                  FI_prev = MAX(PRINS_prev(1), PRINS_prev(2),
     1                         PRINS_prev(3)) / XT
                  FI_prev = MAX(FI_prev, ZERO)
               ELSE
                  FI_prev = ABS(MIN(PRINS_prev(1), PRINS_prev(2),
     1                             PRINS_prev(3))) / XT
                  FI_prev = MAX(FI_prev, ZERO)
               END IF

C              Linear interpolation: alpha at which FI crosses FI_CRIT
               IF ((FI_MAX - FI_prev) .GT. 1.0D-12) THEN
                  alpha_init = (FI_CRIT - FI_prev)
     1                       / (FI_MAX - FI_prev)
                  alpha_init = MAX(alpha_init, 0.10D0)
                  alpha_init = MIN(alpha_init, 0.95D0)
                  PNEWDT = MIN(PNEWDT, alpha_init)
               END IF
               RETURN
            END IF

C           ----- Crack activation (FI within [1.0, 1.05]) -----
            STATEV(1) = ONE

            CALL getN_PE(STATEV(5), STATEV(6), N)

C           Local crack stress at initiation (sigcr0, taucr0)
            STATEV(16:17) = MATMUL(TRANSPOSE(N), SIG)
            STATEV(13)    = STATEV(4) * XT

C           scr_old = local crack stress at initiation
            STATEV(14) = STATEV(16)
            STATEV(15) = STATEV(17)

C           Initialise crack strains to ZERO (not TINY_).
C           At initiation Ecr=0 and scr=sigcr0, placing the state
C           exactly at the top of the softening curve. This eliminates
C           the one-increment lag between STATEV(9) and STATEV(14).
            STATEV(7)  = ZERO
            STATEV(8)  = ZERO
            STATEV(9)  = ZERO
            STATEV(10) = ZERO
            STATEV(11) = ZERO
            STATEV(12) = ZERO
         END IF

      END IF    ! end pre-peak block

C     Re-read sigcr0 / taucr0 in case they were just written above
C     to store the final value when peak occurs
      sigcr0 = STATEV(16)
      taucr0 = STATEV(17)

C=======================================================================
C     POST-PEAK BLOCK  (crack is active)
C=======================================================================
      IF (STATEV(1) .GE. ONE) THEN

C        Retrieve crack history
         Ecr_MAX   = STATEV(7:8)
         Ecr_OLD   = STATEV(9:10)
         dEcr      = STATEV(11:12)
         scr_old   = STATEV(14:15)
         scr_start = scr_old   ! snapshot for Fix-GAP second Dcr re-evaluation

         CALL getN_PE(STATEV(5), STATEV(6), N)

C        STEP 7: Compute D_cr and predicted crack stress
         CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                   nu, eta, DT, Dcr, scr_old, scr, dEcr)
         STATEV(14:15) = scr

C        STEP 8: Compute crack strain increment, update Ecr history
         CALL calcEcr_PE(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N, dEcr, Ecr,
     1                   GIc, sigcr0, taucr0, Leff, nu, eta, DT,
     2                   scr_old, scr, PNEWDT)

C        FIX #6 (corrected by OS-4): Prevent crack interpenetration.
C        Clamp normal crack strain to zero for tensile-mode cracks.
C        !<OS-4> Bug fix: the previous code used -Ecr_OLD(1) but at this
C                point Ecr_OLD has already been updated inside calcEcr_PE
C                to equal Ecr (negative). The correct formula uses STATEV(9)
C                which still holds the TRUE previous crack strain.
C                Previous (wrong): dEcr(1) = -Ecr_OLD(1)
C                Correct         : dEcr(1) = ZERO - STATEV(9)
         IF (MODE .GT. ZERO) THEN
            IF (Ecr(1) .LT. ZERO) THEN
               dEcr(1)    = ZERO - STATEV(9)   !<OS-4> correct formula
               Ecr(1)     = ZERO               !<OS-4>
               Ecr_OLD(1) = ZERO               !<OS-4> keep consistent
            END IF
         END IF

         STATEV(7:8)   = Ecr_MAX
         STATEV(9:10)  = Ecr_OLD
         STATEV(11:12) = dEcr

C        FIX-GAP (Fix #9 re-enabled): Re-compute D_cr and scr using the
C        UPDATED Ecr_OLD (= Ecr_new, now in STATEV(9:10)) and the NEW dEcr.
C
C        Root cause of the gap: the first calcDcr_PE call above (Step 7)
C        used Ecr_OLD_old. With OS-1 total form: scr(1) = D_sec * Ecr_OLD_old.
C        But STATEV(9) stores Ecr_OLD_new → one-increment lag → plotted point
C        (STATEV(9), STATEV(14)) falls off the secant line → visible gap.
C
C        Fix: call calcDcr_PE again with Ecr_OLD_new.  Now:
C          Secant branch (OS-1): scr(1) = D_sec * Ecr_OLD_new
C                                → STATEV(14) consistent with STATEV(9) ✓
C          Softening branch    : scr(1) = sigcr0 + slope_soft*Ecr_OLD_new
C                                → scr at Ecr_new on the softening curve ✓
C        Also refreshes Dcr for the STEP 9 tangent computation.
         CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                   nu, eta, DT, Dcr, scr_start, scr, dEcr)
         STATEV(14:15) = scr

C        STEP 9: Compute cracked tangent D_cocr
         CALL calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)          !<Fv2-D>

         SIG = SIG + MATMUL(Dcocr, DSTRAN)
         IF (isIMPLICIT .EQ. 1) DDSDDE = Dcocr

      END IF    ! end post-peak block

      RETURN
      END SUBROUTINE sca2diso_PE

C=======================================================================
C     SUBROUTINE DEFINITIONS
C=======================================================================

C-----------------------------------------------------------------------
C     I. Plane-strain isotropic elastic stiffness matrix (4x4)
C-----------------------------------------------------------------------
      SUBROUTINE STIFFNESS_MATRIX_2D_PE(E, nu, Dco)

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: E, nu
      REAL*8, INTENT(OUT) :: Dco(4,4)
      REAL*8              :: c, G

      c   = E / ((1.0D0 - 2.0D0*nu) * (1.0D0 + nu))
      G   = E / (2.0D0 * (1.0D0 + nu))

      Dco = 0.0D0

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

C-----------------------------------------------------------------------
C     II. Crack transformation matrix N (4x2)
C-----------------------------------------------------------------------
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

C-----------------------------------------------------------------------
C     III. Wrapper for ABAQUS SPRIND
C-----------------------------------------------------------------------
      SUBROUTINE mysprind_PE(S, PS, AN, LSTR, NDI, NSHR)

      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: NDI, NSHR, LSTR
      REAL*8,  INTENT(IN)  :: S(NDI+NSHR)
      REAL*8,  INTENT(OUT) :: PS(3), AN(3,3)

      CALL SPRIND(S, PS, AN, LSTR, NDI, NSHR)

      RETURN
      END SUBROUTINE mysprind_PE

C-----------------------------------------------------------------------
C     IV. Wrapper for ABAQUS XIT
C-----------------------------------------------------------------------
      SUBROUTINE myExit_PE()
      CALL XIT
      RETURN
      END SUBROUTINE myExit_PE

C-----------------------------------------------------------------------
C     V. calcEcr_PE — compute crack strain increment and update history
C-----------------------------------------------------------------------
      SUBROUTINE calcEcr_PE(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N,
     1           dEcr, Ecr, GIc, sigcr0, taucr0, Leff, nu, eta, DT,
     2           scr_old, scr, PNEWDT)

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: Ecr_MAX(2), Ecr_OLD(2), dEcr(2)
      REAL*8, INTENT(IN)    :: DSTRAN(4), Dco(4,4), N(4,2)
      REAL*8, INTENT(OUT)   :: Ecr(2)
      REAL*8, INTENT(IN)    :: GIc, sigcr0, taucr0, Leff, nu, eta, DT
      REAL*8, INTENT(INOUT) :: scr_old(2), scr(2)
      REAL*8, INTENT(INOUT) :: PNEWDT

      REAL*8 :: Dcr(2,2)
      REAL*8 :: NT_Dco(2,4), NT_Dco_N(2,2)
      REAL*8 :: Y(2,2), invY(2,2)
      REAL*8 :: detY, alpha
      REAL*8 :: Ecr_n_ult, sigcr0_r
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20
      REAL*8, PARAMETER :: R_RES  = 1.0D-8

      LOGICAL :: in_nofail_zone
      REAL*8  :: Ecr_OLD_entry(2)

      sigcr0_r = sigcr0 * R_RES
      Ecr_n_ult = 2.0D0*GIc/(Leff*(sigcr0 - sigcr0_r))

      Ecr_OLD_entry  = Ecr_OLD
      in_nofail_zone = (Ecr_MAX(1) .LT. Ecr_n_ult)

      CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                nu, eta, DT, Dcr, scr_old, scr, dEcr)

      NT_Dco   = MATMUL(TRANSPOSE(N), Dco)
      NT_Dco_N = MATMUL(NT_Dco, N)

      Y = Dcr + NT_Dco_N

      detY = Y(1,1)*Y(2,2) - Y(1,2)*Y(2,1)
      IF (ABS(detY) .LT. SING_TOL*(ABS(Y(1,1))*ABS(Y(2,2)) + 1.0D0))
     1THEN
         WRITE(*,'(A,1PG12.4)') ' WARNING: Y near-singular, det=', detY
         PNEWDT = MIN(PNEWDT, 0.5D0)
         RETURN
      END IF

      CALL matrixInverse(Y, invY, 2, 2)

      dEcr = MATMUL(invY, MATMUL(NT_Dco, DSTRAN))

C     RELOADING OVERSHOOT (secant → softening transition):
C     With total-form traction in both softening and secant branches,
C     and the FIX-GAP second calcDcr_PE call, the traction is consistent
C     whether or not Ecr overshoots Ecr_MAX in a single step.
C     No step cut is needed here — accept the increment as-is.

C     FAILURE OVERSHOOT DETECTION (OS-2new):
C     Only cut the step if Ecr would overshoot Ecr_n_ult by more than
C     0.5%. Small overshoots are handled by the residual branch which
C     clamps scr to sigcr0_r. Only severe overshoots risk instability.
C
C     alpha floor of 0.10 prevents infinite cut loops: ABAQUS will
C     never be asked for less than 10% of the current step.           !<Fv2-C>
      IF (in_nofail_zone) THEN
         IF (dEcr(1) .GT. 0.0D0) THEN
            IF ((Ecr_OLD_entry(1)+dEcr(1)) .GT. 1.005D0*Ecr_n_ult) THEN
               alpha = (Ecr_n_ult - Ecr_OLD_entry(1)) / dEcr(1)
               alpha = MAX(alpha, 0.10D0)
               PNEWDT = MIN(PNEWDT, alpha)                       !<Fv2-C> removed PNEWDT=1.0 override
            END IF
         END IF
      END IF

      Ecr        = Ecr_OLD + dEcr
      Ecr_MAX(1) = MAX(Ecr(1), Ecr_MAX(1))
      Ecr_MAX(2) = MAX(Ecr(2), Ecr_MAX(2))
      Ecr_OLD    = Ecr

      RETURN
      END SUBROUTINE calcEcr_PE

C-----------------------------------------------------------------------
C     VI. calcDcr_PE — crack constitutive matrix (2x2)
C-----------------------------------------------------------------------
      SUBROUTINE calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0,
     1           Leff, nu, eta, DT, Dcr, scr_old, scr, dEcr)

      IMPLICIT NONE
      REAL*8, INTENT(IN)    :: Ecr_OLD(2), Ecr_MAX(2)
      REAL*8, INTENT(IN)    :: GIc, sigcr0, taucr0, Leff, nu
      REAL*8, INTENT(IN)    :: eta, DT
      REAL*8, INTENT(OUT)   :: Dcr(2,2)
      REAL*8, INTENT(IN)    :: scr_old(2), dEcr(2)
      REAL*8, INTENT(OUT)   :: scr(2)

      REAL*8 :: term, sigcr0_r, taucr0_r
      REAL*8 :: Ecr_TEMP(2)
      REAL*8 :: scr_max, Ecr_n_ult
      REAL*8 :: slope_soft, visc_term
      REAL*8, PARAMETER :: R_RES  = 1.0D-8
      REAL*8, PARAMETER :: RES_K  = 1.0D-8
      REAL*8, PARAMETER :: ECR_MIN= 1.0D-14

      IF (sigcr0 .LT. 1.0D-12) THEN
         Dcr      = 0.0D0
         Dcr(1,1) = RES_K
         Dcr(2,2) = RES_K
         scr      = scr_old
         RETURN
      END IF

      sigcr0_r  = sigcr0 * R_RES
      taucr0_r  = taucr0 * R_RES
      term      = 1.0D0 / (2.0D0*(1.0D0 + nu))

      slope_soft = -(sigcr0**2) / (2.0D0*GIc/Leff)

      Ecr_n_ult = 2.0D0*GIc / (Leff*(sigcr0 - sigcr0_r))

C     Viscous regularization: adds eta/DT to Dcr, which makes the
C     effective tangent less negative during softening. This is the
C     same approach used by ABAQUS cohesive elements. The viscous
C     term vanishes as DT->0 or eta->0.
      visc_term = 0.0D0
      IF (eta .GT. 0.0D0 .AND. DT .GT. 0.0D0) THEN
         visc_term = eta / DT
      END IF

      scr_max = sigcr0 + slope_soft * Ecr_MAX(1)
      scr_max = MAX(scr_max, 0.0D0)

      Ecr_TEMP(1) = MAX(ABS(Ecr_OLD(1)), ABS(Ecr_MAX(1)))
      Ecr_TEMP(2) = MAX(ABS(Ecr_OLD(2)), ABS(Ecr_MAX(2)))

      Dcr = 0.0D0

C     --- Normal crack stiffness D_cr(1,1) ---
      IF (Ecr_TEMP(1) .LT. Ecr_n_ult) THEN

         IF (Ecr_OLD(1) .LT. Ecr_MAX(1)) THEN
C           SECANT branch (unloading / reloading within damage surface)
            IF (ABS(Ecr_MAX(1)) .GT. ECR_MIN) THEN
               Dcr(1,1) = scr_max / Ecr_MAX(1)
            ELSE
C              Ecr_MAX ~ 0 (just initiated): use softening slope
               Dcr(1,1) = slope_soft + visc_term
            END IF
C           Total form: traction lies exactly on the secant from origin
C           to (Ecr_MAX, scr_max). No accumulation from scr_old.
            scr(1) = Dcr(1,1) * Ecr_OLD(1)
         ELSE
C           SOFTENING branch (on the damage envelope, crack propagating)
C           Add viscous regularization to the softening slope.
C           Dcr_eff = slope_soft + eta/DT.  When eta>0 this makes the
C           tangent less negative, improving Newton-Raphson convergence.
            Dcr(1,1) = slope_soft + visc_term
C           TOTAL form: scr lies exactly on the analytical softening
C           line. This eliminates the one-increment lag between
C           STATEV(14) and STATEV(9) by guaranteeing scr = f(Ecr_OLD).
C           The viscous contribution is rate-dependent: eta * dEcr/DT.
            scr(1) = sigcr0 + slope_soft * Ecr_OLD(1)
     1             + visc_term * dEcr(1)
            scr(1) = MAX(scr(1), sigcr0_r)
         END IF

      ELSE
C        RESIDUAL branch: element is fully failed
         Dcr(1,1) = RES_K
         scr(1)   = sigcr0_r
      END IF

C     --- Shear crack stiffness D_cr(2,2) ---
      IF (Ecr_TEMP(1) .LT. Ecr_n_ult) THEN
         Dcr(2,2) = Dcr(1,1) * term * 1.0D-3
      ELSE
         Dcr(2,2) = RES_K * term
      END IF

      IF (Ecr_OLD(1) .LT. Ecr_MAX(1) .AND. Ecr_TEMP(1) .LT. Ecr_n_ult)
     1THEN
         scr(2) = Dcr(2,2) * Ecr_OLD(2)
      ELSE
         scr(2) = scr_old(2) + Dcr(2,2)*dEcr(2)
      END IF

      RETURN
      END SUBROUTINE calcDcr_PE

C-----------------------------------------------------------------------
C     VII. calcDcocr_PE — cracked-continuum consistent tangent (4x4)
C
C     !<Fv2-D> PNEWDT restored as INTENT(INOUT). Near-singular tangent
C              now signals ABAQUS with MIN(PNEWDT, 0.5) and prints a
C              warning, rather than silently falling back to Dco.
C-----------------------------------------------------------------------
      SUBROUTINE calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)       !<Fv2-D>

      IMPLICIT NONE
      REAL*8, INTENT(IN)    :: Dcr(2,2), Dco(4,4), N(4,2)
      REAL*8, INTENT(OUT)   :: Dcocr(4,4)
      REAL*8, INTENT(INOUT) :: PNEWDT                            !<Fv2-D>

      REAL*8 :: term22(2,2), invT22(2,2)
      REAL*8 :: term42(4,2), term44(4,4)
      REAL*8 :: detT
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20

      term22 = Dcr + MATMUL(TRANSPOSE(N), MATMUL(Dco, N))

      detT = term22(1,1)*term22(2,2) - term22(1,2)*term22(2,1)
      IF (ABS(detT) .LT. SING_TOL*(ABS(term22(1,1))*ABS(term22(2,2))
     1    + 1.0D0)) THEN
         WRITE(*,'(A,1PG12.4)')                                  !<Fv2-D>
     1       ' WARNING: Dcocr near-singular, det=', detT
         Dcocr  = Dco
         PNEWDT = MIN(PNEWDT, 0.5D0)                            !<Fv2-D>
         RETURN
      END IF

      CALL matrixInverse(term22, invT22, 2, 2)

      term42 = MATMUL(Dco, MATMUL(N, invT22))
      term44 = MATMUL(term42, MATMUL(TRANSPOSE(N), Dco))
      Dcocr  = Dco - term44

      RETURN
      END SUBROUTINE calcDcocr_PE
