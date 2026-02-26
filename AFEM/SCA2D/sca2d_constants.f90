!===============================================================================
! Module: sca2d_constants
! ---------------------------------------------------------------------------
! Defines named constants and a derived type for SCA2D state variables.
!
! Smart Change #5: Named constants replace magic numbers throughout the code
! Smart Change #8: Derived type replaces raw STATEV array indexing
!
! Last Updated: 02/26/2026
! meka1@purdue.edu
!===============================================================================
module sca2d_constants
    implicit none

    ! Precision parameter
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! ---- Named Constants (Smart Change #5) ----

    ! Residual stress fraction: ratio of residual-to-initial crack stress
    ! Used in the traction-separation law to prevent complete zero stiffness
    real(dp), parameter :: RESIDUAL_STRESS_FRACTION = 5.0d-5

    ! Shear retention factor: scales shear crack stiffness relative to normal
    real(dp), parameter :: SHEAR_RETENTION_FACTOR = 1.0d-5

    ! Initial value assigned to crack strain components at crack initiation
    real(dp), parameter :: INITIAL_CRACK_STRAIN = 1.0d-8

    ! Near-zero stiffness assigned after complete failure
    real(dp), parameter :: POST_FAILURE_STIFFNESS = 1.0d-8

    ! Tolerance for singularity detection in 2x2 matrix inverse
    real(dp), parameter :: SINGULARITY_TOL = 1.0d-30

    ! Threshold for checking the crack flag (0.5 separates 0 from 1)
    real(dp), parameter :: CRACK_FLAG_TOL = 0.5d0

    ! ---- Validation Constants (Smart Change #6) ----

    ! Minimum number of state variables required
    integer, parameter :: N_STATEV_REQUIRED = 15

    ! Minimum number of material properties required
    integer, parameter :: N_PROPS_REQUIRED = 5

    ! Expected NTENS for 2D plane stress
    integer, parameter :: NTENS_PLANE_STRESS = 3

    ! ---- Derived Type for State Variables (Smart Change #8) ----
    ! Replaces fragile raw STATEV(i) indexing with named fields

    type :: sca2d_state_t
        real(dp) :: crack_flag         ! STATEV(1):   0=intact, 1=cracked
        real(dp) :: elem_length        ! STATEV(2):   Element characteristic length
        real(dp) :: failure_mode       ! STATEV(3):   +1=tension, -1=compression
        real(dp) :: max_failure_index  ! STATEV(4):   Max failure index at initiation
        real(dp) :: crack_normal(2)    ! STATEV(5:6): Crack normal (cos theta, sin theta)
        real(dp) :: ecr_max(2)         ! STATEV(7:8): Maximum local crack strains encountered
        real(dp) :: ecr_old(2)         ! STATEV(9:10): Previous increment crack strains
        real(dp) :: decr(2)            ! STATEV(11:12): Current crack strain increments
        real(dp) :: scr_init           ! STATEV(13):  Crack stress at initiation
        real(dp) :: scr_current        ! STATEV(14):  Current local crack stress
        real(dp) :: dcr_11             ! STATEV(15):  Current Dcr(1,1) value
    end type sca2d_state_t

contains

    !---------------------------------------------------------------------------
    ! Unpack flat STATEV array into the derived type
    !---------------------------------------------------------------------------
    subroutine unpack_state(statev, nstatv, state)
        integer, intent(in)            :: nstatv
        real(dp), intent(in)           :: statev(nstatv)
        type(sca2d_state_t), intent(out) :: state

        state%crack_flag        = statev(1)
        state%elem_length       = statev(2)
        state%failure_mode      = statev(3)
        state%max_failure_index = statev(4)
        state%crack_normal(1:2) = statev(5:6)
        state%ecr_max(1:2)      = statev(7:8)
        state%ecr_old(1:2)      = statev(9:10)
        state%decr(1:2)         = statev(11:12)
        state%scr_init          = statev(13)
        state%scr_current       = statev(14)
        state%dcr_11            = statev(15)
    end subroutine unpack_state

    !---------------------------------------------------------------------------
    ! Pack derived type back into flat STATEV array
    !---------------------------------------------------------------------------
    subroutine pack_state(state, statev, nstatv)
        type(sca2d_state_t), intent(in) :: state
        integer, intent(in)             :: nstatv
        real(dp), intent(out)           :: statev(nstatv)

        statev(1)     = state%crack_flag
        statev(2)     = state%elem_length
        statev(3)     = state%failure_mode
        statev(4)     = state%max_failure_index
        statev(5:6)   = state%crack_normal(1:2)
        statev(7:8)   = state%ecr_max(1:2)
        statev(9:10)  = state%ecr_old(1:2)
        statev(11:12) = state%decr(1:2)
        statev(13)    = state%scr_init
        statev(14)    = state%scr_current
        statev(15)    = state%dcr_11
    end subroutine pack_state

end module sca2d_constants
