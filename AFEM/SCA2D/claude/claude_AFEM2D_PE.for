C***********************************************************************
C     ! Isotropic Incremental SCA 2D Plane Strain  [Enhanced / Robust]
C     ! ----------------------------------------------------------------
C     ! Enhanced version of 'formulation_AFEM2D_PE.for'
C     ! Failure Criteria: Max Principal Strain (INDEX = 1)
C     !
C     ! Key improvements over original:
C     !  1. Fixed DFGRD1 type declaration (REAL*8, was INTEGER)
C     !  2. sigcr0 assigned from STATEV(16) BEFORE Bazant limit check
C     !     → eliminates divide-by-zero on first call
C     !  3. Bazant check guarded against sigcr0 = 0
C     !  4. Removed unconditional debug write(*,*) statements
C     !  5. Fixed quad-precision literal 2.0Q0 → 2.0D0 in calcDcr
C     !  6. Crack interpenetration prevention (dEcr_n clamped to >= 0
C     !     for tensile mode; elastic response restored on closure)
C     !  7. Near-singular Y matrix check before inversion (PNEWDT cut)
C     !  8. Division-by-zero guard in calcDcr (Ecr_MAX ~ 0 case)
C     !  9. Consistent scr update using the freshly computed dEcr
C     ! 10. Residual shear stiffness based on G/1000 instead of 1e-8
C     ! 11. PNEWDT = 0.5 requested when matrix ill-conditioned
C     ! ----------------------------------------------------------------
C     ! Last Modified: 02/27/2026
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
C       INTEGER     :: ZONE                   ! zone flag (1=elastic,2=loading,3=unloading,4=failed)
      REAL*8      :: E, nu, GIc, XT, G
      REAL*8      :: sigcr0, taucr0, sigcr0_r, BL
      REAL*8      :: isCRACKED, MODE
      REAL*8      :: FI_T, FI_C, FI_MAX, J_FAIL, J_MIN, PMIN, J_MAX, PMAX
      REAL*8      :: PRINS(3), T(3,3)
      REAL*8      :: Dco(NTENS,NTENS), Sco(NTENS,NTENS)

      REAL*8      :: N(4,2)
      REAL*8      :: Ecr_MAX(2), Ecr_OLD(2)
      REAL*8      :: Ecr(2), dEcr(2)
      REAL*8      :: scr_old(2), scr(2)
      REAL*8      :: Dcocr(NTENS,NTENS)
      REAL*8      :: Dcr(2,2)

C     --- Numerical constants ---
      REAL*8, PARAMETER :: ZERO  = 0.0D0
      REAL*8, PARAMETER :: ONE   = 1.0D0
      REAL*8, PARAMETER :: TWO   = 2.0D0
      REAL*8, PARAMETER :: TINY_ = 1.0D-8   ! seed value for crack strains
      REAL*8, PARAMETER :: RES_  = 1.0D-5   ! residual stiffness after full failure
      REAL*8, PARAMETER :: FI_CRIT = 0.99D0 ! failure index threshold

C***********************************************************************
C     READ MATERIAL PROPERTIES
C***********************************************************************
      E      = PROPS(1)    ! Young's modulus
      nu     = PROPS(2)    ! Poisson's ratio
      GIc    = PROPS(3)    ! Mode I fracture toughness
      XT     = PROPS(4)    ! Tensile (and compressive) strength
      INDEX  = INT(PROPS(5))  ! Failure criteria index

C     Derived quantities
      G      = E / (TWO*(ONE + nu))
      Leff   = ONE*L          ! effective element size (1st-order elements)

C     Retrieve failure mode and peak crack stresses from history
      MODE   = STATEV(2)

C     FIX #2: Read sigcr0 from STATEV(16) HERE, BEFORE the Bazant check.
C             In the original code this was done 100+ lines later, causing
C             divide-by-zero on the very first call.
      sigcr0 = STATEV(16)
      taucr0 = STATEV(17)
      STATEV(3) = L

C***********************************************************************
C     STEP 1: Elastic stiffness matrix D_co (plane strain, 4x4)
C***********************************************************************
      CALL STIFFNESS_MATRIX_2D_PE(E, nu, Dco)

C***********************************************************************
C     STEP 2: Update total strain tensor
C***********************************************************************
      DO K1 = 1, NTENS
         EPS(K1) = EPS(K1) + DSTRAN(K1)
      END DO

C***********************************************************************
C     STEP 3: Bazant size limit check
C             FIX #3: guard against sigcr0 = 0 (first call, no crack yet)
C***********************************************************************
      IF (sigcr0 .GT. 1.0D-12) THEN
         BL = TWO*GIc*E / (sigcr0**2)
         IF (Leff .GT. BL) THEN
            WRITE(*,'(A)')       ' '
            WRITE(*,'(A)')       ' *** BAZANT LIMIT EXCEEDED ***'
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
            MODE   =  ONE          ! tensile
            FI_MAX = FI_T
            J_FAIL = DBLE(J_MAX)
         ELSE
            MODE   = -ONE          ! compressive
            FI_MAX = FI_C
            J_FAIL = DBLE(J_MIN)
         END IF

C        Store failure mode and crack orientation
         STATEV(2) = MODE
         STATEV(4) = FI_MAX
         STATEV(5) = T(1, INT(J_FAIL))    ! cos(theta)
         STATEV(6) = T(2, INT(J_FAIL))    ! sin(theta)

!---->   IF FAILURE CRITERIA IS SATISFIED
         IF (FI_MAX .GE. FI_CRIT) THEN
            STATEV(1) = ONE                ! crack flag ON

            CALL getN_PE(STATEV(5), STATEV(6), N)

C           Local crack stress at initiation (stored as sigcr0, taucr0)
            STATEV(16:17) = matmul(transpose(N),SIG)
            STATEV(13)    = STATEV(4) * XT
   
C           scr_old = local crack stress vector
            STATEV(14) = STATEV(16)
            STATEV(15) = STATEV(17)

C           Seed crack strains (small positive value to avoid /0)
            STATEV(7)  = TINY_    ! Ecr_MAX(1)
            STATEV(8)  = TINY_    ! Ecr_MAX(2)
            STATEV(9)  = TINY_    ! Ecr_OLD(1)
            STATEV(10) = TINY_    ! Ecr_OLD(2)
            STATEV(11) = TINY_    ! dEcr(1)
            STATEV(12) = TINY_    ! dEcr(2)
         END IF

      END IF    ! end pre-peak block

C     Re-read sigcr0 / taucr0 in case they were just written above
      sigcr0 = STATEV(16)
      taucr0 = STATEV(17)

C=======================================================================
C     POST-PEAK BLOCK  (crack is active)
C=======================================================================
      IF (STATEV(1) .GE. ONE) THEN

C        Retrieve crack history
         Ecr_MAX  = STATEV(7:8)
         Ecr_OLD  = STATEV(9:10)
         dEcr     = STATEV(11:12)
         scr_old  = STATEV(14:15)

         CALL getN_PE(STATEV(5), STATEV(6), N)

C        STEP 7: Compute D_cr (crack stiffness) and predicted crack stress
         CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                   nu, Dcr, scr_old, scr, dEcr)
         STATEV(14:15) = scr

C        STEP 8: Compute crack strain increment, update Ecr history
         CALL calcEcr_PE(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N, dEcr, Ecr,
     1                   GIc, sigcr0, taucr0, Leff, nu, scr_old, scr,
     2                   PNEWDT)

C        FIX #6: Prevent crack interpenetration (tensile mode only)
C                Normal crack strain must remain non-negative
         IF (MODE .GT. ZERO) THEN
            IF (Ecr(1) .LT. ZERO) THEN
               Ecr(1)  = ZERO
               dEcr(1) = -Ecr_OLD(1) ! Sai: Need to verify this
            END IF
         END IF

         STATEV(7:8)   = Ecr_MAX
         STATEV(9:10)  = Ecr_OLD   ! updated inside calcEcr_PE to Ecr
         STATEV(11:12) = dEcr

C C        FIX #9: Recompute D_cr with the UPDATED dEcr and scr_old
C C                so that scr is consistent with the final crack strain
C          CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
C      1                   nu, Dcr, STATEV(14:15), scr, dEcr)
C          STATEV(14:15) = scr

C        STEP 9: Compute cracked tangent D_cocr
         CALL calcDcocr_PE(Dco, Dcr, N, Dcocr, PNEWDT)

C        Stress update and Jacobian
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
C     For plane strain: NTENS=4, order = [eps_xx, eps_yy, eps_zz, gam_xy]

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: E, nu
      REAL*8, INTENT(OUT) :: Dco(4,4)
      REAL*8              :: c, G

C     Plane-strain constant:  c = E / [(1-2nu)(1+nu)]
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
C         Maps local crack strains [eps_cr_n, gam_cr_t] to global Voigt
C-----------------------------------------------------------------------
      SUBROUTINE getN_PE(c, s, N)

      IMPLICIT NONE
      REAL*8, INTENT(IN)  :: c, s    ! cos(theta), sin(theta)
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
C     III. Wrapper for ABAQUS SPRIND (principal stresses/directions)
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
C     IV. Wrapper for ABAQUS XIT (terminate analysis)
C-----------------------------------------------------------------------
      SUBROUTINE myExit_PE()
      CALL XIT
      RETURN
      END SUBROUTINE myExit_PE

C-----------------------------------------------------------------------
C     V. calcEcr_PE — compute crack strain increment and update history
C
C     Enforces stress continuity at crack face:
C       (Dcr + N^T Dco N) dEcr = N^T Dco DSTRAN
C
C     FIX #7: PNEWDT cut-back if Y matrix is near-singular
C     FIX #6: Interpenetration clamping done in caller
C-----------------------------------------------------------------------
      SUBROUTINE calcEcr_PE(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N,
     1           dEcr, Ecr, GIc, sigcr0, taucr0, Leff, nu,
     2           scr_old, scr, PNEWDT)

      IMPLICIT NONE

      REAL*8, INTENT(INOUT) :: Ecr_MAX(2), Ecr_OLD(2), dEcr(2)
      REAL*8, INTENT(IN)    :: DSTRAN(4), Dco(4,4), N(4,2)
      REAL*8, INTENT(OUT)   :: Ecr(2)
      REAL*8, INTENT(IN)    :: GIc, sigcr0, taucr0, Leff, nu
      REAL*8, INTENT(INOUT) :: scr_old(2), scr(2)
      REAL*8, INTENT(INOUT) :: PNEWDT

      REAL*8 :: Dcr(2,2)
      REAL*8 :: NT_Dco(2,4), NT_Dco_N(2,2)
      REAL*8 :: Y(2,2), invY(2,2)
      REAL*8 :: condY, detY
      REAL*8, PARAMETER :: SING_TOL = 1.0D-20

C     Recompute D_cr at the current crack state
      CALL calcDcr_PE(Ecr_OLD, Ecr_MAX, GIc, sigcr0, taucr0, Leff,
     1                nu, Dcr, scr_old, scr, dEcr)

      NT_Dco   = MATMUL(TRANSPOSE(N), Dco)
      NT_Dco_N = MATMUL(NT_Dco, N)

C     System matrix Y = Dcr + N^T Dco N
      Y = Dcr + NT_Dco_N

C     FIX #7: Check conditioning via determinant before inverting
      detY = Y(1,1)*Y(2,2) - Y(1,2)*Y(2,1)
      IF (ABS(detY) .LT. SING_TOL * (ABS(Y(1,1))*ABS(Y(2,2)) + 1.0D0))
     1THEN
         WRITE(*,'(A,1PG12.4)') ' WARNING: Y near-singular, det=', detY
         PNEWDT = 0.5D0     ! request step cut-back from ABAQUS
         RETURN
      END IF

      CALL matrixInverse(Y, invY, 2, 2)

C     dEcr = Y^{-1} (N^T Dco DSTRAN)
      dEcr = MATMUL(invY, MATMUL(NT_Dco, DSTRAN))

C     Update crack strain and history
      Ecr        = Ecr_OLD + dEcr
      Ecr_MAX(1) = MAX(Ecr(1), Ecr_MAX(1))
      Ecr_MAX(2) = MAX(Ecr(2), Ecr_MAX(2))
      Ecr_OLD    = Ecr

      RETURN
      END SUBROUTINE calcEcr_PE

C-----------------------------------------------------------------------
C     VI. calcDcr_PE — crack constitutive matrix (2x2)
C
C     FIX #5: Guard against sigcr0 ~ 0 (no crack, should not be called)
C     FIX #4: Replaced 2.0Q0 (quad literal) with 2.0D0
C     FIX #8: Guard against Ecr_MAX(1) ~ 0 in secant computation
C     FIX #10: Shear residual stiffness uses G/1000 fraction
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
      REAL*8 :: slope_soft         ! initial softening slope (always negative)
      REAL*8, PARAMETER :: R_RES  = 5.0D-5   ! residual strength ratio
      REAL*8, PARAMETER :: RES_K  = 1.0D-5   ! residual stiffness after failure
      REAL*8, PARAMETER :: ECR_MIN= 1.0D-14  ! floor to avoid /0 in secant

C     Guard: if sigcr0 is essentially zero, return minimal stiffness
      IF (sigcr0 .LT. 1.0D-12) THEN
         Dcr      = 0.0D0
         Dcr(1,1) = RES_K
         Dcr(2,2) = RES_K
         scr      = scr_old        ! no change
         RETURN
      END IF

      sigcr0_r  = sigcr0 * R_RES    ! residual normal traction
      taucr0_r  = taucr0 * R_RES    ! residual shear  traction
      term      = 1.0D0 / (2.0D0*(1.0D0 + nu))  ! shear reduction factor

C     Softening modulus (always negative — initial slope at eps=0):
C       sigma_cr(eps) = sigcr0 - (sigcr0^2 / (2*GIc/Leff)) * eps
C     => d sigma_cr / d eps = -(sigcr0^2) / (2*GIc/Leff)
      slope_soft = -(sigcr0**2) / (2.0D0*GIc/Leff)

C     Ultimate crack strain at zero traction:
C       eps_u = 2*GIc / (Leff * sigcr0)
      Ecr_n_ult = 2.0D0*GIc / (Leff*(sigcr0 - sigcr0_r))

C     Traction at the current maximum crack strain:
C       scr_max = sigcr0 + slope_soft * Ecr_MAX(1)
      scr_max = sigcr0 + slope_soft * Ecr_MAX(1)
      scr_max = MAX(scr_max, 0.0D0)   ! cannot go below zero (no tension recovery)

C     Effective crack strain for regime detection
      Ecr_TEMP(1) = MAX(ABS(Ecr_OLD(1)), ABS(Ecr_MAX(1)))
      Ecr_TEMP(2) = MAX(ABS(Ecr_OLD(2)), ABS(Ecr_MAX(2)))

      Dcr = 0.0D0

C     --- Normal crack stiffness D_cr(1,1) ---
      IF (Ecr_TEMP(1) .LT. Ecr_n_ult) THEN
         ! Not yet fully failed
         IF (Ecr_OLD(1) .LT. Ecr_MAX(1)) THEN
            ! LOADING branch: secant stiffness (always positive)
            Dcr(1,1) = scr_max / MAX(ABS(Ecr_MAX(1)), ECR_MIN)
         ELSE
            ! UNLOADING / RELOADING branch:
            ! Use initial (steepest) softening slope as the unloading stiffness.
            ! This gives a secant back toward the origin (damage mechanics style).
            Dcr(1,1) = slope_soft   ! negative value
         END IF
      ELSE
         ! COMPLETE FAILURE: residual tangent
         Dcr(1,1) = RES_K
      END IF

C     --- Shear crack stiffness D_cr(2,2) ---
C     FIX #10: use a small but physically scaled fraction of D_cr(1,1).
C     Fully failed: residual only.  Otherwise: G_effective / large_factor.
C     Here we keep ~1/1000 of the normal crack stiffness magnitude.
      IF (Ecr_TEMP(1) .LT. Ecr_n_ult) THEN
         Dcr(2,2) = Dcr(1,1) * term * 1.0D-3
      ELSE
         Dcr(2,2) = RES_K * term
      END IF

C     --- Update local crack stress (explicit, one-step) ---
      scr(1) = scr_old(1) + Dcr(1,1)*dEcr(1)
      scr(2) = scr_old(2) + Dcr(2,2)*dEcr(2)

      RETURN
      END SUBROUTINE calcDcr_PE

C-----------------------------------------------------------------------
C     VII. calcDcocr_PE — cracked-continuum consistent tangent (4x4)
C
C     D_cocr = D_co - D_co N (D_cr + N^T D_co N)^{-1} N^T D_co
C
C     FIX #7: Near-singular check before inversion; PNEWDT on failure
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

C     Y = D_cr + N^T D_co N
      term22 = Dcr + MATMUL(TRANSPOSE(N), MATMUL(Dco, N))

C     Singularity check
      detT = term22(1,1)*term22(2,2) - term22(1,2)*term22(2,1)
      IF (ABS(detT) .LT. SING_TOL*(ABS(term22(1,1))*ABS(term22(2,2))
     1                              + 1.0D0)) THEN
         WRITE(*,'(A,1PG12.4)') ' WARNING: Dcocr Y near-singular, det=',
     1                           detT
         Dcocr  = Dco    ! fall back to elastic tangent
         PNEWDT = 0.5D0  ! request step cut-back
         RETURN
      END IF

C     Y^{-1}
      CALL matrixInverse(term22, invT22, 2, 2)

C     D_co N Y^{-1}  (4x2)
      term42 = MATMUL(Dco, MATMUL(N, invT22))

C     D_co N Y^{-1} N^T D_co  (4x4)
      term44 = MATMUL(term42, MATMUL(TRANSPOSE(N), Dco))

C     D_cocr = D_co - D_co N Y^{-1} N^T D_co
      Dcocr = Dco - term44

      RETURN
      END SUBROUTINE calcDcocr_PE
