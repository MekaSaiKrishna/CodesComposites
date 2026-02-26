!===============================================================================
! Improved UMAT Interface for SCA2D
! ---------------------------------------------------------------------------
! Abaqus/Standard User Material subroutine for the Incremental Smeared Crack
! Approach (SCA) on 2D isotropic materials (CPS4 - Plane Stress Elements).
!
! Smart Change #6: Input validation (NTENS, NPROPS, NSTATV, material properties)
! Smart Change #8: Uses Fortran 90 modules with explicit interfaces
!
! Compile order: sca2d_constants.f90 -> matrix_utils.f90
!                -> formulation_AFEM2D_v2.f90 -> main_SCA2D_inc_v2.f90
!
! Last Updated: 02/26/2026
! meka1@purdue.edu
!===============================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

C     Smart Change #8: USE modules for explicit interfaces
      USE sca2d_constants, ONLY: N_STATEV_REQUIRED, N_PROPS_REQUIRED,
     1                           NTENS_PLANE_STRESS
      USE sca2d_formulation, ONLY: sca2diso_v2

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C-----------------------------------------------------------------------
C     Smart Change #6: Input Validation
C-----------------------------------------------------------------------

C     Check that this is a plane stress analysis (NTENS must be 3)
      IF (NTENS .NE. NTENS_PLANE_STRESS) THEN
          WRITE(7,*) '*** SCA2D ERROR: Wrong element type ***'
          WRITE(7,*) '    Expected NTENS = ', NTENS_PLANE_STRESS,
     1               ' (plane stress)'
          WRITE(7,*) '    Got NTENS = ', NTENS
          WRITE(7,*) '    Use CPS4 or similar plane stress elements.'
          WRITE(7,*) '    Element: ', NOEL
          CALL XIT
      END IF

C     Check that enough material properties are provided
      IF (NPROPS .LT. N_PROPS_REQUIRED) THEN
          WRITE(7,*) '*** SCA2D ERROR: Insufficient properties ***'
          WRITE(7,*) '    Required: ', N_PROPS_REQUIRED
          WRITE(7,*) '    Provided: ', NPROPS
          WRITE(7,*) '    Need: E, nu, GIc, XT, INDEX'
          CALL XIT
      END IF

C     Check that enough state variables are allocated
      IF (NSTATV .LT. N_STATEV_REQUIRED) THEN
          WRITE(7,*) '*** SCA2D ERROR: Insufficient state variables ***'
          WRITE(7,*) '    Required: ', N_STATEV_REQUIRED
          WRITE(7,*) '    Allocated: ', NSTATV
          WRITE(7,*) '    Set *Depvar to at least ', N_STATEV_REQUIRED
          CALL XIT
      END IF

C     Validate material properties are physically meaningful
C     E > 0
      IF (PROPS(1) .LE. 0.0D0) THEN
          WRITE(7,*) '*** SCA2D ERROR: Invalid Youngs Modulus ***'
          WRITE(7,*) '    E = ', PROPS(1), ' (must be > 0)'
          CALL XIT
      END IF

C     0 < nu < 0.5
      IF (PROPS(2) .LE. 0.0D0 .OR. PROPS(2) .GE. 0.5D0) THEN
          WRITE(7,*) '*** SCA2D ERROR: Invalid Poissons Ratio ***'
          WRITE(7,*) '    nu = ', PROPS(2), ' (must be 0 < nu < 0.5)'
          CALL XIT
      END IF

C     GIc > 0
      IF (PROPS(3) .LE. 0.0D0) THEN
          WRITE(7,*) '*** SCA2D ERROR: Invalid Fracture Toughness ***'
          WRITE(7,*) '    GIc = ', PROPS(3), ' (must be > 0)'
          CALL XIT
      END IF

C     XT > 0
      IF (PROPS(4) .LE. 0.0D0) THEN
          WRITE(7,*) '*** SCA2D ERROR: Invalid Tensile Strength ***'
          WRITE(7,*) '    XT = ', PROPS(4), ' (must be > 0)'
          CALL XIT
      END IF

C-----------------------------------------------------------------------
C     Call the improved SCA formulation
C-----------------------------------------------------------------------
      CALL sca2diso_v2(NTENS, NDI, NSHR, DFGRD1, NSTATV, NPROPS, NOEL,
     1     TIME(1), TIME(2), DTIME, PNEWDT, CELENT, PROPS, STRAN,
     2     DSTRAN, STRESS, STATEV, DDSDDE, DTEMP, 1)

C-----------------------------------------------------------------------
      RETURN
      END
