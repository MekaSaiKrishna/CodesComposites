!===============================================================================
! Module: sca2d_formulation
! ---------------------------------------------------------------------------
! Incremental Smeared Crack Approach (SCA) for 2D isotropic materials.
! Improved version of formulation_AFEM2D.for implementing all smart changes.
!
! Smart Change #2:  Stiffness computed once per call, passed to all helpers
! Smart Change #3:  N^T*Dco*N and N^T*Dco pre-computed and shared
! Smart Change #5:  Named constants from sca2d_constants module
! Smart Change #7:  Redundant calcDcr call eliminated
! Smart Change #8:  Modules, INTENT declarations, derived types throughout
! Smart Change #9:  Improved Bazant limit check (unit 7 messaging, NOEL, PNEWDT)
!
! Failure Criteria:
!     0) Max Principal Stress (INDEX = 0)
!     1) Max Stress           (INDEX = 1)
!
! Last Updated: 02/26/2026
! meka1@purdue.edu
!===============================================================================
module sca2d_formulation
    use sca2d_constants
    use matrix_utils
    implicit none
    private
    public :: sca2diso_v2

contains

    !===========================================================================
    ! Main SCA2D Constitutive Subroutine (Improved)
    !===========================================================================
    subroutine sca2diso_v2(NTENS, NDI, NSHR, DFGRD1, NSTATV, NPROPS, NOEL, &
            stepTIME, totalTIME, DT, PNEWDT, L, PROPS, EPS, DSTRAN, SIG,   &
            STATEV, DDSDDE, DTEMP, isIMPLICIT)

        implicit none

        ! ---- Arguments with explicit INTENT (Smart Change #8) ----
        integer, intent(in)     :: NTENS, NDI, NSHR, NSTATV, NPROPS, NOEL
        integer, intent(in)     :: DFGRD1(3,3)
        real(dp), intent(in)    :: stepTIME, totalTIME, DT, L, DTEMP
        real(dp), intent(inout) :: PNEWDT
        real(dp), intent(in)    :: PROPS(NPROPS)
        real(dp), intent(inout) :: EPS(NTENS), SIG(NTENS), STATEV(NSTATV)
        real(dp), intent(in)    :: DSTRAN(NTENS)
        real(dp), intent(out)   :: DDSDDE(NTENS,NTENS)
        integer, intent(in)     :: isIMPLICIT

        ! ---- Material properties ----
        real(dp) :: E, nu, GIc, XT, sigcr0
        integer  :: INDEX

        ! ---- Local variables ----
        type(sca2d_state_t) :: state
        real(dp) :: Dco(3,3), Leff, BL
        real(dp) :: PRINS(3), T(3,3)
        real(dp) :: FI_T, FI_C, FI_MAX, PMIN, MODE
        integer  :: J_FAIL, J_MIN, I, K1

        ! ---- Crack computation variables ----
        real(dp) :: N(3,2), NT(2,3)
        real(dp) :: NT_Dco_N(2,2), NT_Dco(2,3)   ! Smart Change #3
        real(dp) :: Ecr(2), Dcr(2,2), Dcocr(3,3)
        real(dp) :: scr

        ! ---- Extract Material Properties ----
        E      = PROPS(1)
        nu     = PROPS(2)
        GIc    = PROPS(3)
        XT     = PROPS(4)
        INDEX  = int(PROPS(5))
        sigcr0 = XT

        ! ---- Unpack State Variables into derived type (Smart Change #8) ----
        call unpack_state(STATEV, NSTATV, state)
        state%elem_length = L
        Leff = 1.0d0 * L   ! Linear element (1st order)

        ! ---- Smart Change #2: Compute elastic stiffness once per call ----
        call compute_stiffness_2d_ps(E, nu, Dco)

        ! ---- Update total strain ----
        do K1 = 1, NTENS
            EPS(K1) = EPS(K1) + DSTRAN(K1)
        end do

        ! ---- Smart Change #9: Improved Bazant limit check ----
        BL = 2.0d0 * GIc * E / sigcr0**2
        if (BL < Leff) then
            ! Write to Abaqus .msg file (unit 7) with element number
            write(7,*) ' '
            write(7,*) '*** SCA2D ERROR: Bazant size-effect limit violated ***'
            write(7,*) '    Element number:             ', NOEL
            write(7,*) '    Required max element size:  ', BL
            write(7,*) '    Current element size:       ', Leff
            write(7,*) '    Action: Reduce mesh size to at most ', BL
            ! Also write to stdout for visibility
            write(*,*) ' '
            write(*,*) '*** SCA2D ERROR: Bazant limit violated at element ', NOEL
            write(*,*) '    Bazant limit = ', BL, '  Element size = ', Leff
            ! Request smaller time increment before aborting
            PNEWDT = 0.25d0
            call myExit_v2()
        end if

        ! ==================================================================
        ! STEP A: Pre-crack phase (linear elastic until failure)
        ! ==================================================================
        if (state%crack_flag < CRACK_FLAG_TOL) then

            ! Linear elastic stress update
            SIG = SIG + matmul(Dco, DSTRAN)

            ! Jacobian = elastic stiffness
            DDSDDE = Dco

            ! Compute principal stresses and directions
            call mysprind_v2(SIG, PRINS, T, 1, 2, 1)

            ! ---- Evaluate failure indices ----
            FI_T   = 0.0d0
            FI_C   = 0.0d0
            J_FAIL = 1

            ! Tensile failure index
            if (PRINS(1) > 0.0d0) FI_T = PRINS(1) / XT

            ! Compressive failure index (most negative principal stress)
            PMIN  = PRINS(1)
            J_MIN = 1
            do I = 2, 3
                if (PRINS(I) < PMIN) then
                    PMIN  = PRINS(I)
                    J_MIN = I
                end if
            end do
            if (PMIN < 0.0d0) FI_C = abs(PMIN) / XT

            ! Determine dominant failure mode
            if (FI_T >= FI_C) then
                MODE   = 1.0d0      ! Tensile mode
                FI_MAX = FI_T
                J_FAIL = 1
            else
                MODE   = -1.0d0     ! Compressive mode
                FI_MAX = FI_C
                J_FAIL = J_MIN
            end if

            ! Store failure info in state
            state%failure_mode      = MODE
            state%max_failure_index = FI_MAX
            state%crack_normal(1)   = T(1, J_FAIL)   ! cos(theta)
            state%crack_normal(2)   = T(2, J_FAIL)   ! sin(theta)

            ! Check for crack initiation
            if (FI_MAX >= 1.0d0) then
                state%crack_flag = 1.0d0
            end if

            ! Initialize crack strain history (Smart Change #5: named constant)
            state%ecr_max  = INITIAL_CRACK_STRAIN
            state%ecr_old  = INITIAL_CRACK_STRAIN
            state%decr     = INITIAL_CRACK_STRAIN
            state%scr_init = FI_MAX * XT
        end if

        ! ==================================================================
        ! STEP B: Post-crack phase (cracked material response)
        ! ==================================================================
        if (state%crack_flag > CRACK_FLAG_TOL) then

            ! Compute transformation matrix N from crack orientation
            call compute_N(state%crack_normal(1), state%crack_normal(2), N)
            NT = transpose(N)

            ! ---- Smart Change #3: Pre-compute shared matrix products ----
            ! These products are used by BOTH compute_crack_strain and
            ! compute_Dcocr, so computing them once saves redundant work.
            NT_Dco   = matmul(NT, Dco)         ! 2x3 = (2x3)(3x3)
            NT_Dco_N = matmul(NT_Dco, N)       ! 2x2 = (2x3)(3x2)

            ! ---- Smart Change #7: Compute Dcr ONCE ----
            ! The original code called calcDcr in sca2diso AND again inside
            ! calcEcr with identical inputs, producing identical results.
            ! Now it is computed once here and passed to both consumers.
            call compute_Dcr(state%ecr_old, state%ecr_max, GIc, sigcr0, &
                             Leff, nu, Dcr, state%scr_current, scr, state%decr)

            ! Compute incremental crack strain
            ! (Smart Change #1: uses analytical 2x2 inverse internally)
            ! (Smart Change #3: receives pre-computed NT_Dco_N and NT_Dco)
            ! (Smart Change #7: Dcr passed in, no internal recomputation)
            call compute_crack_strain(state%ecr_old, state%ecr_max, DSTRAN, &
                    Dcr, NT_Dco_N, NT_Dco, state%decr, Ecr)

            ! Update crack stress and stiffness in state
            state%scr_current = scr
            state%dcr_11      = Dcr(1,1)

            ! Compute degraded tangent stiffness
            ! (Smart Change #1: uses analytical 2x2 inverse internally)
            ! (Smart Change #3: receives pre-computed NT_Dco_N and NT_Dco)
            call compute_Dcocr(Dco, Dcr, N, NT_Dco_N, NT_Dco, Dcocr)

            ! Stress update with degraded stiffness
            SIG = SIG + matmul(Dcocr, DSTRAN)

            ! Jacobian update for implicit solver
            if (isIMPLICIT == 1) DDSDDE = Dcocr
        end if

        ! ---- Pack state variables back to STATEV array (Smart Change #8) ----
        call pack_state(state, STATEV, NSTATV)

    end subroutine sca2diso_v2


    !===========================================================================
    ! I) Isotropic Plane Stress Stiffness Matrix
    !    Smart Change #2: Called once per sca2diso_v2 invocation, result shared
    !===========================================================================
    subroutine compute_stiffness_2d_ps(E, nu, Dco)
        real(dp), intent(in)  :: E, nu
        real(dp), intent(out) :: Dco(3,3)
        real(dp) :: c

        c = E / (1.0d0 - nu**2)
        Dco = reshape(                                     &
            [c,      c*nu,   0.0d0,                        &
             c*nu,   c,      0.0d0,                        &
             0.0d0,  0.0d0,  c * 0.5d0 * (1.0d0 - nu)],   &
            [3, 3])
    end subroutine compute_stiffness_2d_ps


    !===========================================================================
    ! II) Transformation Matrix N (1 crack)
    !     Maps local crack tractions to global stress components
    !===========================================================================
    subroutine compute_N(c, s, N)
        real(dp), intent(in)  :: c, s
        real(dp), intent(out) :: N(3,2)

        N(1,1) = c**2
        N(1,2) = -c*s
        N(2,1) = s**2
        N(2,2) = c*s
        N(3,1) = 2.0d0*s*c
        N(3,2) = c**2 - s**2
    end subroutine compute_N


    !===========================================================================
    ! III) Principal Stress Wrapper (calls Abaqus SPRIND utility)
    !      Smart Change #8: Explicit INTENT on all arguments
    !===========================================================================
    subroutine mysprind_v2(S, PS, AN, LSTR, NDI, NSHR)
        integer, intent(in)   :: LSTR, NDI, NSHR
        real(dp), intent(in)  :: S(NDI+NSHR)
        real(dp), intent(out) :: PS(3), AN(3,3)

        call SPRIND(S, PS, AN, LSTR, NDI, NSHR)
    end subroutine mysprind_v2


    !===========================================================================
    ! IV) Graceful Exit Wrapper
    !     Smart Change #9: Diagnostics written before exit in callers
    !===========================================================================
    subroutine myExit_v2()
        call XIT
    end subroutine myExit_v2


    !===========================================================================
    ! V) Crack Stiffness Matrix Dcr
    !    Smart Change #5: All magic numbers replaced with named constants
    !
    !    Implements a linear traction-separation law:
    !      - Loading:            Secant stiffness to current softening point
    !      - Unloading/Reload:   Tangent of the softening curve
    !      - Complete failure:   Near-zero residual stiffness
    !===========================================================================
    subroutine compute_Dcr(Ecr_OLD, Ecr_MAX, GIc, sigcr0, Leff, nu, &
                           Dcr, scr_old, scr, dEcr)
        real(dp), intent(in)  :: Ecr_OLD(2), Ecr_MAX(2)
        real(dp), intent(in)  :: GIc, sigcr0, Leff, nu
        real(dp), intent(out) :: Dcr(2,2)
        real(dp), intent(in)  :: scr_old
        real(dp), intent(out) :: scr
        real(dp), intent(in)  :: dEcr(2)

        real(dp) :: sigcr0_r, shear_term, scr_max, Ecr_TEMP, ecr_critical

        ! Residual stress (Smart Change #5: named constant)
        sigcr0_r = sigcr0 * RESIDUAL_STRESS_FRACTION

        ! Shear reduction factor
        shear_term = 1.0d0 / (2.0d0 * (1.0d0 + nu))

        ! Crack stress at max crack strain (linear traction-separation law)
        scr_max = sigcr0 + (-(sigcr0)**2 / (2.0d0 * GIc / Leff)) * Ecr_MAX(1)

        ! Critical crack strain for complete failure
        ecr_critical = 2.0d0 * GIc / Leff / (sigcr0 - sigcr0_r)

        Ecr_TEMP = max(abs(Ecr_OLD(1)), abs(Ecr_MAX(1)))

        Dcr = 0.0d0

        ! Determine crack stiffness based on loading state
        if (Ecr_TEMP < ecr_critical) then
            if (Ecr_OLD(1) < Ecr_MAX(1)) then
                ! LOADING: Secant stiffness from origin to softening curve
                Dcr(1,1) = scr_max / abs(Ecr_MAX(1))
            else
                ! UNLOADING/RELOADING: Tangent of the softening law
                Dcr(1,1) = (-(sigcr0_r - sigcr0)**2 / (2.0d0 * GIc / Leff))
            end if
        else
            ! COMPLETE FAILURE: Near-zero residual (Smart Change #5)
            Dcr(1,1) = POST_FAILURE_STIFFNESS
        end if

        ! Shear crack stiffness (Smart Change #5: named constant)
        Dcr(2,2) = Dcr(1,1) * shear_term * SHEAR_RETENTION_FACTOR

        ! Update local crack stress
        scr = scr_old + Dcr(1,1) * dEcr(1)
    end subroutine compute_Dcr


    !===========================================================================
    ! VI) Incremental Crack Strain Computation
    !
    !     Smart Change #1: Uses analytical 2x2 inverse (not general LU)
    !     Smart Change #3: Receives pre-computed NT_Dco_N and NT_Dco
    !     Smart Change #7: Dcr passed in; no internal calcDcr call
    !
    !     Solves: dEcr = [Dcr + N^T*Dco*N]^(-1) * N^T*Dco*DSTRAN
    !===========================================================================
    subroutine compute_crack_strain(Ecr_OLD, Ecr_MAX, DSTRAN, &
                Dcr, NT_Dco_N, NT_Dco, dEcr, Ecr)
        real(dp), intent(inout) :: Ecr_OLD(2), Ecr_MAX(2)
        real(dp), intent(in)    :: DSTRAN(3)
        real(dp), intent(in)    :: Dcr(2,2)
        real(dp), intent(in)    :: NT_Dco_N(2,2)
        real(dp), intent(in)    :: NT_Dco(2,3)
        real(dp), intent(out)   :: dEcr(2), Ecr(2)

        real(dp) :: Y(2,2), invY(2,2)
        integer  :: flag

        ! Form the system matrix
        Y = Dcr + NT_Dco_N

        ! Smart Change #1: Analytical 2x2 inverse
        call inverse_2x2(Y, invY, flag)

        ! Crack strain increment
        dEcr = matmul(invY, matmul(NT_Dco, DSTRAN))

        ! New local crack strain
        Ecr = Ecr_OLD + dEcr

        ! Update maximum crack strain encountered (history tracking)
        Ecr_MAX = max(Ecr, Ecr_MAX)

        ! Set current as old for next increment
        Ecr_OLD = Ecr
    end subroutine compute_crack_strain


    !===========================================================================
    ! VII) Degraded (Cracked) Tangent Stiffness Dcocr
    !
    !      Smart Change #1: Uses analytical 2x2 inverse (not general LU)
    !      Smart Change #3: Receives pre-computed NT_Dco_N and NT_Dco
    !
    !      Dcocr = Dco - Dco*N*[Dcr + N^T*Dco*N]^(-1)*N^T*Dco
    !===========================================================================
    subroutine compute_Dcocr(Dco, Dcr, N, NT_Dco_N, NT_Dco, Dcocr)
        real(dp), intent(in)  :: Dco(3,3), Dcr(2,2), N(3,2)
        real(dp), intent(in)  :: NT_Dco_N(2,2)
        real(dp), intent(in)  :: NT_Dco(2,3)
        real(dp), intent(out) :: Dcocr(3,3)

        real(dp) :: Y(2,2), invY(2,2)
        real(dp) :: term32(3,2), term33(3,3)
        integer  :: flag

        ! Form the system matrix
        Y = Dcr + NT_Dco_N

        ! Smart Change #1: Analytical 2x2 inverse
        call inverse_2x2(Y, invY, flag)

        ! Dco*N*[...]^(-1)
        term32 = matmul(Dco, matmul(N, invY))

        ! Dco*N*[...]^(-1)*N^T*Dco
        term33 = matmul(term32, NT_Dco)

        ! Final degraded stiffness
        Dcocr = Dco - term33
    end subroutine compute_Dcocr

end module sca2d_formulation
