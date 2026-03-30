C***********************************************************************
C     ! Isotropic Incremental SCA 2D Plane Strain  [Enhanced v1 - Overshoot Fixes]
C     ! ----------------------------------------------------------------
C     ! Based on: claude_AFEM2D_PE.for
C     !
C     ! Overshoot fixes added in this version (see fixes_overshooting.md):
C     !
C     !  OS-1 [calcDcr_PE, secant branch]:
C     !       Replace incremental scr update with TOTAL FORM:
C     !         scr(1) = Dcr(1,1) * Ecr_OLD(1)   (was: scr_old + Dcr*dEcr_old)
C     !       This prevents the traction from inheriting an inconsistent
C     !       scr_old value when the code transitions from the softening
C     !       branch to the secant branch.
C     !
C     !  OS-2 [calcEcr_PE]:
C     !       Detect reloading overshoot — when the secant branch is active
C     !       and the computed Ecr(1) exceeds Ecr_MAX(1). Compute the
C     !       fractional step alpha at which Ecr(1) = Ecr_MAX(1) and set
C     !       PNEWDT = alpha * 0.85 so ABAQUS restarts with a smaller step.
C     !
C     !  OS-3 [REMOVED]:
C     !       The first-closure PNEWDT=0.5 was removed because it caused an
C     !       infinite ABAQUS step-cut loop: ABAQUS resets STATEV on restart,
C     !       so Ecr_OLD_save = Ecr_MAX_save on every retry → OS-3 fires again
C     !       → exponential step reduction → ABAQUS gives up. Fix-GAP below
C     !       eliminates the underlying traction discontinuity, making OS-3
C     !       redundant.
C     !
C     !  FIX-GAP [sca2diso_PE, after STATEV(9:10) update]:
C     !       Re-enable Fix #9: call calcDcr_PE a second time using the
C     !       UPDATED Ecr_OLD (= Ecr_new). With OS-1 total form in the secant
C     !       branch: scr(1) = D_sec * Ecr_new → STATEV(14) is now consistent
C     !       with STATEV(9), eliminating the one-increment lag that produced
C     !       the gap in the σ_cr vs ε_cr plot during unload/reload.
C     !
C     !  OS-4 [sca2diso_PE, interpenetration clamp]:
C     !       Fix wrong formula for dEcr when clamping Ecr(1) to zero.
C     !       Use STATEV(9) — the true previous crack strain — instead of
C     !       the already-updated Ecr_OLD.
C     !       Correct:  dEcr(1) = ZERO - STATEV(9)
C     !       Previous: dEcr(1) = -Ecr_OLD(1)   [was wrong: Ecr_OLD = Ecr here]
C     !
C     ! All other code is identical to claude_AFEM2D_PE.for.
C     ! Changed lines are marked  !<OS-n>  for traceability.
C     ! ----------------------------------------------------------------
C     ! Last Modified: 03/03/2026
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
      INTEGER     :: INDEX, I, J, K1
      REAL*8      :: E, nu, GIc, XT, G
      REAL*8      :: sigcr0, taucr0, sigcr0_r, BL
      REAL*8      :: isCRACKED, MODE
      REAL*8      :: FI_T, FI_C, FI_MAX, J_FAIL, J_MIN, PMIN, J_MAX, PMAX
      REAL*8      :: PRINS(3), T(3,3)
      REAL*8      :: Dco(NTENS,NTENS), Sco(NTENS,NTENS)

      REAL*8      :: N(4,2)
      REAL*8      :: Ecr_MAX(2), Ecr_OLD(2), Ecr_NEW(2)
      REAL*8      :: Ecr(2), dEcr(2)
      REAL*8      :: scr_old(2), scr(2)
      REAL*8      :: Dcocr(NTENS,NTENS)
      REAL*8      :: Dcr(2,2)
      REAL*8      :: scr_start(2)  ! snapshot of scr_old for Fix-GAP second Dcr call

C     --- Numerical constants ---
      REAL*8, PARAMETER :: ZERO    = 0.0D0
      REAL*8, PARAMETER :: ONE     = 1.0D0
      REAL*8, PARAMETER :: TWO     = 2.0D0
      REAL*8, PARAMETER :: TINY_   = 1.0D-8 ! seed value for crack strains
      REAL*8, PARAMETER :: RES_    = 1.0D-5 ! residual stiffness after full failure
      REAL*8, PARAMETER :: FI_CRIT = 1.0D0 ! failure index threshold

C***********************************************************************
C     READ MATERIAL PROPERTIES
C***********************************************************************
      E      = PROPS(1)      ! Young's modulus
      nu     = PROPS(2)      ! Poisson's ratio
      GIc    = PROPS(3)      ! Mode I fracture toughness
      XT     = PROPS(4)      ! Tensile (and compressive) strength
      INDEX  = INT(PROPS(5)) ! Failure criteria index

      G      = E/(TWO*(ONE + nu))
      Leff   = ONE*L !for linear elements (CPE4)

      MODE      = STATEV(2)
      sigcr0    = STATEV(19)
      taucr0    = STATEV(20)
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
         STATEV(19) = ZERO
         STATEV(20) = ZERO
         STATEV(4)  = ZERO

C        --- Linear elastic stress update ---
         SIG    = SIG + MATMUL(Dco, DSTRAN)
         DDSDDE = Dco

C        --- Principal stress computation ---
         CALL mysprind_PE(SIG, PRINS, T, 1, 3, 1)

C        Failure indices
         FI_T   = ZERO ! Tension Failure Index
         FI_C   = ZERO ! Compression Failure Index
         FI_MAX = ZERO ! Maximum Failure Index

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
            MODE   =  ONE           ! tensile
            FI_MAX = FI_T
            J_FAIL = DBLE(J_MAX)
         ELSE
            MODE   = -ONE           ! compressive
            FI_MAX = FI_C
            J_FAIL = DBLE(J_MIN)
         END IF

C        Store failure mode and crack orientation
         STATEV(2) = MODE
         STATEV(4) = FI_MAX
         STATEV(5) = T(1, INT(J_FAIL)) ! cos(theta)
         STATEV(6) = T(2, INT(J_FAIL)) ! sin(theta)

!---->   IF FAILURE CRITERIA IS SATISFIED
         IF (FI_MAX .GE. FI_CRIT) THEN
            STATEV(1) = ONE ! crack flag ON

            CALL getN_PE(STATEV(5), STATEV(6), N)

C           Local crack stress at initiation (stored as sigcr0, taucr0)
            STATEV(19:20) = MATMUL(TRANSPOSE(N), SIG)
            STATEV(13)    = STATEV(4) * XT

C           scr_old = local crack stress vector
            STATEV(14) = STATEV(16)
            STATEV(15) = STATEV(17)

C           Seed crack strains (small positive value to avoid 0)
            STATEV(7)  = TINY_
            STATEV(8)  = TINY_
            STATEV(9)  = TINY_
            STATEV(10) = TINY_
            STATEV(11) = TINY_
            STATEV(12) = TINY_
            STATEV(16) = TINY_
            STATEV(17) = TINY_
            STATEV(18) = TINY_
         END IF

      END IF    ! end pre-peak block

C     Re-read sigcr0 / taucr0 in case they were just written above
C     to store the final value when peak occurs
      sigcr0 = STATEV(19)
      taucr0 = STATEV(20)

C=======================================================================
C     POST-PEAK BLOCK  (crack is active)
C=======================================================================
      IF (STATEV(1) .GE. ONE) THEN

C        Retrieve crack history
         Ecr_MAX   = STATEV(7:8)
         Ecr_OLD   = STATEV(9:10)
         dEcr      = STATEV(11:12)
         scr_old   = STATEV(14:15)
         Ecr_NEW   = STATEV(16:17)

         scr_start = scr_old   ! snapshot for Fix-GAP second Dcr re-evaluation

         CALL getN_PE(STATEV(5), STATEV(6), N)

C        STEP 7: Compute D_cr and predicted crack stress
         CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                   nu, Dcr, scr_old, scr, dEcr)
         STATEV(14:15) = scr

C        STEP 8: Compute crack strain increment, update Ecr history
         CALL calcEcr_PE(Ecr_OLD, Ecr_MAX, Ecr_NEW, DSTRAN, Dco, N, dEcr, Ecr,
     1                   GIc, sigcr0, taucr0, Leff, nu, scr_old, scr,
     2                   PNEWDT)

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

         STATEV(7:8)    = Ecr_MAX
         STATEV(9:10)   = Ecr_NEW !updating Ecr_OLD = Ecr_NEW
         STATEV(11:12)  = dEcr
         STATEV(16:17)  = Ecr_NEW

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
C          Softening branch    : scr(1) = scr_start(1) + slope_soft * dEcr_new
C                                → scr at Ecr_new on the softening curve ✓
C        Also refreshes Dcr for the STEP 9 tangent computation.

C          CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
C      1                   nu, Dcr, scr_start, scr, dEcr)
C          STATEV(14:15) = scr

C        STEP 9: Compute cracked tangent D_cocr
         CALL calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)

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
C         Unchanged from claude_AFEM2D_PE.for
C-----------------------------------------------------------------------
      SUBROUTINE myExit_PE()
      CALL XIT
      RETURN
      END SUBROUTINE myExit_PE

C-----------------------------------------------------------------------
C     V. calcEcr_PE — compute crack strain increment and update history
C
C     Changes vs claude_AFEM2D_PE.for:
C       !<OS-2> Reloading overshoot detection + PNEWDT
C-----------------------------------------------------------------------
      SUBROUTINE calcEcr_PE(Ecr_OLD, Ecr_MAX, Ecr_NEW, DSTRAN, Dco, N,
     1           dEcr, Ecr, GIc, sigcr0, taucr0, Leff, nu,
     2           scr_old, scr, PNEWDT)

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: Ecr_MAX(2), Ecr_OLD(2), Ecr_NEW(2), dEcr(2)
      REAL*8, INTENT(IN)    :: DSTRAN(4), Dco(4,4), N(4,2)
      REAL*8, INTENT(OUT)   :: Ecr(2)
      REAL*8, INTENT(IN)    :: GIc, sigcr0, taucr0, Leff, nu
      REAL*8, INTENT(INOUT) :: scr_old(2), scr(2)
      REAL*8, INTENT(INOUT) :: PNEWDT

      REAL*8 :: Dcr(2,2)
      REAL*8 :: NT_Dco(2,4), NT_Dco_N(2,2)
      REAL*8 :: Y(2,2), invY(2,2)
      REAL*8 :: detY, alpha
      REAL*8 :: Ecr_n_ult,sigcr0_r
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20
      REAL*8, PARAMETER :: R_RES  = 1.0D-8

C     !<OS-2> Record whether we are currently in the secant branch and
C             the Ecr_MAX at entry, so we can detect a reloading overshoot
C             AFTER dEcr has been computed.
      LOGICAL :: in_secant_zone                             !<OS-2>
      LOGICAL :: in_nofail_zone

      REAL*8  :: Ecr_MAX_entry(2), Ecr_OLD_entry(2)        !<OS-2>

      sigcr0_r = sigcr0 * R_RES
      Ecr_n_ult = 2.0D0*GIc/(Leff*(sigcr0 - sigcr0_r))

      Ecr_MAX_entry  = Ecr_MAX                               !<OS-2>
      Ecr_OLD_entry  = Ecr_OLD                               !<OS-2>
      in_secant_zone = (Ecr_OLD(1) .LT. Ecr_MAX(1))        !<OS-2>
      in_nofail_zone = (Ecr_MAX(1) .LT. Ecr_n_ult)

      CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                nu, Dcr, scr_old, scr, dEcr)

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

C     !<OS-2> RELOADING OVERSHOOT DETECTION
C             If we were in the secant zone (Ecr_OLD < Ecr_MAX) and the
C             computed dEcr(1) carries Ecr(1) past the previous Ecr_MAX(1),
C             then the reloading front has hit the softening envelope
C             mid-step. The secant stiffness is no longer valid beyond
C             Ecr_MAX — the softening branch should take over.
C
C             Fix: request PNEWDT proportional to the fraction (alpha)
C             of the step that reaches Ecr_MAX. ABAQUS restarts with
C             DSTRAN' ~ alpha * DSTRAN, so Ecr stays at or below Ecr_MAX.
C
C             A safety factor of 0.85 prevents repeated boundary
C             oscillation.
      IF (in_secant_zone) THEN                                     !<OS-2>
         IF (dEcr(1) .GT. 0.0D0) THEN                             !<OS-2>
            IF ((Ecr_OLD_entry(1) + dEcr(1)) .GT. Ecr_MAX_entry(1)) THEN !<OS-2>
C              alpha = fraction of step at which Ecr(1) reaches Ecr_MAX(1)
               alpha = (Ecr_MAX_entry(1) - Ecr_OLD_entry(1))/ dEcr(1)   !<OS-2>
               alpha = MAX(alpha, 0.20D0)   ! lower bound: at least 20%  !<OS-2>
               PNEWDT = MIN(PNEWDT, alpha * 0.85D0)               !<OS-2>
C              Note: we do NOT clamp Ecr here. ABAQUS will redo this
C              increment from scratch with the smaller step, so the
C              computed dEcr (and hence Ecr_MAX update) will be correct.
                  IF (alpha .LE. 0.20D0) THEN
                        PNEWDT = 1.0d0
C                         write(*,*) "PNEWDT: ", PNEWDT
                  END IF
            END IF                                                   !<OS-2>
         END IF                                                      !<OS-2>
      END IF                                                         !<OS-2>

C-----------------------------------------------------------------------------
C     !<OS-2new> FAILURE OVERSHOOT DETECTION
C             If we were in the NO-fail zone (Ecr_MAX(1) < Ecr_n_ult) and the
C             computed dEcr(1) carries Ecr(1) past Ecr_n_ult,
C             then the loading front has hit the failure envelope
C             mid-step. The secant stiffness is no longer valid beyond
C             Ecr_n_ult — the failure branch should take over.
C
C             Fix: request PNEWDT proportional to the fraction (alpha)
C             of the step that reaches Ecr_n_ult. ABAQUS restarts with
C             DSTRAN' ~ alpha * DSTRAN, so Ecr stays at or below Ecr_n_ult.

      IF (in_nofail_zone) THEN                                     !<OS-2>
         IF (dEcr(1) .GT. 0.0D0) THEN                             !<OS-2>
            IF ((Ecr_OLD_entry(1) + dEcr(1)) .GT. Ecr_MAX_entry(1)) THEN !<OS-2>
C              alpha = fraction of step at which Ecr(1) reaches Ecr_MAX(1)
               alpha = (Ecr_MAX_entry(1) - Ecr_OLD_entry(1))/ dEcr(1)   !<OS-2>
               alpha = MAX(alpha, 0.20D0)   ! lower bound: at least 20%  !<OS-2>
               PNEWDT = MIN(PNEWDT, alpha * 0.85D0)               !<OS-2>
C              Note: we do NOT clamp Ecr here. ABAQUS will redo this
C              increment from scratch with the smaller step, so the
C              computed dEcr (and hence Ecr_MAX update) will be correct.
                  IF (alpha .LE. 0.20D0) THEN
                        PNEWDT = 1.0d0
                  END IF
            END IF                                                   !<OS-2>
         END IF                                                      !<OS-2>
      END IF 

      Ecr        = Ecr_OLD + dEcr
      Ecr_MAX(1) = MAX(Ecr(1), Ecr_MAX(1))
      Ecr_MAX(2) = MAX(Ecr(2), Ecr_MAX(2))
      Ecr_NEW    = Ecr

      RETURN
      END SUBROUTINE calcEcr_PE

C-----------------------------------------------------------------------
C     VI. calcDcr_PE — crack constitutive matrix (2x2)
C-----------------------------------------------------------------------
      SUBROUTINE calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0,
     1           Leff, nu, Dcr, scr_old, scr, dEcr)

      IMPLICIT NONE
      REAL*8, INTENT(IN)    :: Ecr_OLD(2), Ecr_MAX(2)
      REAL*8, INTENT(IN)    :: GIc, sigcr0, taucr0, Leff, nu
      REAL*8, INTENT(OUT)   :: Dcr(2,2)
      REAL*8, INTENT(IN)    :: scr_old(2), dEcr(2)
      REAL*8, INTENT(OUT)   :: scr(2)

      REAL*8 :: term, sigcr0_r, taucr0_r
      REAL*8 :: Ecr_TEMP(2)
      REAL*8 :: scr_max, Ecr_n_ult
      REAL*8 :: slope_soft
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

      scr_max = sigcr0 + slope_soft * Ecr_MAX(1)
      scr_max = MAX(scr_max, 0.0D0)

      Ecr_TEMP(1) = MAX(ABS(Ecr_OLD(1)), ABS(Ecr_MAX(1)))
      Ecr_TEMP(2) = MAX(ABS(Ecr_OLD(2)), ABS(Ecr_MAX(2)))

      Dcr = 0.0D0

C     --- Normal crack stiffness D_cr(1,1) ---
      IF (Ecr_TEMP(1) .LT. Ecr_n_ult) THEN

         IF (Ecr_OLD(1) .LT. Ecr_MAX(1)) THEN
C           SECANT branch (unloading / reloading within damage surface)
            Dcr(1,1) = scr_max / MAX(ABS(Ecr_MAX(1)), ECR_MIN)

C           For the shear component we also use the total form scaled
C           consistently with the normal component.

            scr(1) = Dcr(1,1) * Ecr_OLD(1)
         ELSE
C           SOFTENING branch (on the damage envelope, crack propagating)
            Dcr(1,1) = slope_soft

C           Incremental form is correct here because the code is
C           continuously on the softening curve; scr_old is consistent.
            scr(1) = scr_old(1) + Dcr(1,1)*dEcr(1)     ! unchanged
         END IF

      ELSE
         Dcr(1,1) = RES_K
         scr(1)   = scr_old(1) + Dcr(1,1)*dEcr(1)
      END IF

C     --- Shear crack stiffness D_cr(2,2) ---
      IF (Ecr_TEMP(1) .LT. Ecr_n_ult) THEN
         Dcr(2,2) = Dcr(1,1) * term * 1.0D-3
      ELSE
         Dcr(2,2) = RES_K * term
      END IF

      IF (Ecr_OLD(1) .LT. Ecr_MAX(1) .AND. Ecr_TEMP(1) .LT. Ecr_n_ult) THEN
         scr(2) = Dcr(2,2) * Ecr_OLD(2)
      ELSE
         scr(2) = scr_old(2) + Dcr(2,2)*dEcr(2)       ! unchanged
      END IF

      RETURN
      END SUBROUTINE calcDcr_PE

C-----------------------------------------------------------------------
C     VII. calcDcocr_PE — cracked-continuum consistent tangent (4x4)
C-----------------------------------------------------------------------
      SUBROUTINE calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)

      IMPLICIT NONE
      REAL*8, INTENT(IN)    :: Dcr(2,2), Dco(4,4), N(4,2)
      REAL*8, INTENT(OUT)   :: Dcocr(4,4)
      REAL*8, INTENT(INOUT) :: PNEWDT

      REAL*8 :: term22(2,2), invT22(2,2)
      REAL*8 :: term42(4,2), term44(4,4)
      REAL*8 :: detT
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20

      term22 = Dcr + MATMUL(TRANSPOSE(N), MATMUL(Dco, N))

      detT = term22(1,1)*term22(2,2) - term22(1,2)*term22(2,1)
      IF (ABS(detT) .LT. SING_TOL*(ABS(term22(1,1))*ABS(term22(2,2)) + 1.0D0)) THEN
         WRITE(*,'(A,1PG12.4)') ' WARNING: Dcocr Y near-singular, det=', detT
         Dcocr  = Dco
         PNEWDT = MIN(PNEWDT, 0.5D0)
         RETURN
      END IF

      CALL matrixInverse(term22, invT22, 2, 2)

      term42 = MATMUL(Dco, MATMUL(N, invT22))
      term44 = MATMUL(term42, MATMUL(TRANSPOSE(N), Dco))
      Dcocr  = Dco - term44

      RETURN
      END SUBROUTINE calcDcocr_PE
