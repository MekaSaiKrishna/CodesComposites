!***********************************************************************
! Isotropic Incremental SCA 1D 
!***********************************************************************
      SUBROUTINE sca1diso(NSTATV, NPROPS, NOEL, stepTIME, totalTIME,
     1  DT, L, PROPS, EPS, DSTRAN, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)

      IMPLICIT NONE

    ! Global Arguments
      INTEGER:: NSTATV, NPROPS, NOEL, isIMPLICIT
      REAL*8:: stepTIME, totalTIME, DT, L, DTEMP
      REAL*8:: PROPS(NPROPS), DSTRAN(6)
      REAL*8:: EPS(6), SIG(6), STATEV(NSTATV)
      REAL*8:: DDSDDE(6,6)

    ! Local Variables
      INTEGER:: INDEX, K1
      REAL*8:: Dco(6,6), Dcocr(6,6), Dcr(3,3)
      REAL*8:: PRINS(3), T(3,3), N(6,3)
      REAL*8:: E, nu, XT, GIc, sigcr0, BL
      REAL*8:: isCRACKED, MODE
      REAL*8:: Ecr_MAX(3), Ecr_OLD(3), Ecr(3), dEcr(3)

    ! 1. Extract Properties
      E      = PROPS(1) 
      nu     = PROPS(2) 
      GIc    = PROPS(3) 
      XT     = PROPS(4) 
      INDEX  = PROPS(5)
      sigcr0 = XT 

    ! 2. Initialize State Variables & Params
      isCRACKED = STATEV(1)
      BL        = (2.0d0*GIc*E/sigcr0**2) ! Bazant Limit

    ! 3. Stiffness Initialization
      CALL STIFFNESS_MATRIX(E, nu, Dco)
      EPS = EPS + DSTRAN ! Update total strain

    ! 4. Check Element Size (Regularization check)
      IF (BL < L) THEN
          WRITE(*,*) 'FATAL: BAZANT LIMIT OVERSHOOT! Element size L is too large.'
          WRITE(*,*) 'Max allowed (BL): ', BL, ' Current L: ', L
          CALL XIT()
      END IF

    ! 5. Logic for Uncracked Material
      IF (isCRACKED < 0.5d0) THEN
        
      ! Trial Stress (Elastic)
      SIG = SIG + MATMUL(Dco, DSTRAN)
      DDSDDE = Dco
  
      ! Failure Criteria Check
      IF (INDEX == 0) THEN 
      CALL mysprind(SIG, PRINS, T, 1, 3, 3)

      IF (MAXVAL(PRINS) > sigcr0) THEN
      ! Crack Initiated: Roll back stress and set flag
      SIG = SIG - MATMUL(Dco, DSTRAN)
      isCRACKED = 1.0d0
      STATEV(1) = isCRACKED

      STATEV(2) = 1.0d0 ! Mode

      ! Store orientation (Flattened T matrix)
      STATEV(9:11)  = T(1,:)
      STATEV(12:14) = T(2,:)
      STATEV(15:17) = T(3,:)
      END IF
      
      END IF
      
      END IF

    ! 6. Logic for Cracked Material
      IF (isCRACKED >= 0.99d0) THEN

      ! Retrieve state
      Ecr_MAX = STATEV(6:8)
      Ecr_OLD = STATEV(18:20)
      Ecr     = STATEV(21:23)
      T(1,:)  = STATEV(9:11)
      T(2,:)  = STATEV(12:14)
      T(3,:)  = STATEV(15:17)
  
      ! Compute N based on fixed crack orientation
      CALL getN(T, N)
  
      ! Consistently update Local Crack Strains
      CALL calcEcr(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N, dEcr, Ecr, GIc, sigcr0, L, nu)
  
      ! Update Jacobian and Stress
      CALL calcDcr(Ecr_OLD, Ecr_MAX, GIc, sigcr0, L, nu, Dcr)
      CALL calcDcocr(Dco, Dcr, N, Dcocr)
  
      SIG = SIG + MATMUL(Dcocr, DSTRAN)
      IF (isIMPLICIT == 1) DDSDDE = Dcocr
  
      ! Final State Save
      STATEV(6:8)   = Ecr_MAX
      STATEV(18:20) = Ecr_OLD
      STATEV(21:23) = Ecr
      STATEV(24)    = Dcr(1,1)
      END IF

      RETURN
      END

!-----------------------------------------------------------------------
! I) Isotropic Elastic Stiffness
!-----------------------------------------------------------------------
      SUBROUTINE STIFFNESS_MATRIX(E, nu, Dco)
      IMPLICIT NONE

      REAL*8:: E, nu
      REAL*8:: Dco(6,6)
      REAL*8:: mu, lam, twomu
  
      mu    = 0.5d0 * E / (1.0d0 + nu)
      twomu = 2.0d0 * mu
      lam   = E * nu / ((1.0d0 + nu) * (1.0d0 - 2.0d0 * nu))
  
      Dco = 0.0d0
      Dco(1,1)=lam+twomu; Dco(1,2)=lam;       Dco(1,3)=lam
      Dco(2,1)=lam;       Dco(2,2)=lam+twomu; Dco(2,3)=lam
      Dco(3,1)=lam;       Dco(3,2)=lam;       Dco(3,3)=lam+twomu
      Dco(4,4)=mu;        Dco(5,5)=mu;        Dco(6,6)=mu

      RETURN
      END
!-----------------------------------------------------------------------
! II) Transformation Matrix N (Voigt notation orientation)
!-----------------------------------------------------------------------
      SUBROUTINE getN(o,N)
      IMPLICIT NONE

      REAL*8:: o(3,3)
      REAL*8:: N(6,3)

!     ordering of principal direction: o(index of eval,1:3)   
C     o(1,:) = (/1.0d0, 0.0d0, 0.0d0/)
C     o(2,:) = (/0.0d0, 1.0d0, 0.0d0/)
C     o(3,:) = (/1.0d0, 0.0d0, 1.0d0/)
      
      N(1,1) = o(1,1)**2
      N(1,2) = 0.1D1*o(1,1)*o(2,1)
      N(1,3) = 0.1D1*o(3,1)*o(1,1)
      N(2,1) = o(1,2)**2
      N(2,2) = 0.1D1*o(1,2)*o(2,2)
      N(2,3) = 0.1D1*o(3,2)*o(1,2)
      N(3,1) = o(1,3)**2
      N(3,2) = 0.1D1*o(1,3)*o(2,3)
      N(3,3) = 0.1D1*o(3,3)*o(1,3)
      N(4,1) = 0.2D1*o(1,1)*o(1,2)
      N(4,2) = o(1,1)*o(2,2) + o(1,2)*o(2,1)
      N(4,3) = o(3,1)*o(1,2) + o(3,2)*o(1,1)
      N(5,1) = 0.2D1*o(1,3)*o(1,1)
      N(5,2) = o(1,3)*o(2,1) + o(1,1)*o(2,3)
      N(5,3) = o(3,3)*o(1,1) + o(3,1)*o(1,3)
      N(6,1) = 0.2D1*o(1,2)*o(1,3)
      N(6,2) = o(1,2)*o(2,3) + o(1,3)*o(2,2)
      N(6,3) = o(3,2)*o(1,3) + o(3,3)*o(1,2)

      RETURN
      END 

!-----------------------------------------------------------------------
! V) Calculating Incremental Crack Strain
!-----------------------------------------------------------------------
      SUBROUTINE calcEcr(Ecr_OLD, Ecr_MAX, DSTRAN, Dco, N, dEcr, Ecr, GIc, sigcr0, L, nu)
      IMPLICIT NONE

      REAL*8:: Ecr_OLD(3), Ecr_MAX(3), DSTRAN(6), Dco(6,6), N(6,3)
      REAL*8:: GIc, sigcr0, L, nu
      REAL*8:: dEcr(3), Ecr(3)
      REAL*8:: Dcr(3,3), NT(3,6), Y(3,3), invY(3,3), NT_Dco(3,6)
  
      NT = TRANSPOSE(N)
      CALL calcDcr(Ecr_OLD, Ecr_MAX, GIc, sigcr0, L, nu, Dcr)
      
      NT_Dco = MATMUL(NT, Dco)
      Y      = Dcr + MATMUL(NT_Dco, N)
      
      CALL matrixInverse(Y, invY, 3, 3)
      
      dEcr = MATMUL(invY, MATMUL(NT_Dco, DSTRAN))
      Ecr  = Ecr_OLD + dEcr

      RETURN
      END

!-----------------------------------------------------------------------
! VI) Traction-Separation Softening Matrix (Dcr)
!-----------------------------------------------------------------------
      SUBROUTINE calcDcr(Ecr_OLD, Ecr_MAX, GIc, sigcr0, L, nu, Dcr)
      IMPLICIT NONE

      REAL*8:: Ecr_OLD(3), Ecr_MAX(3), GIc, sigcr0, L, nu
      REAL*8:: Dcr(3,3)
      REAL*8:: Ecr_eff, sigcr0_r, limit
  
      sigcr0_r = sigcr0 * 5.0d-5
      limit    = (2.0d0 * GIc) / (L * (sigcr0 - sigcr0_r))
      Ecr_eff  = MAX(ABS(Ecr_OLD(1)), ABS(Ecr_MAX(1)))
      
      Dcr = 0.0d0
      IF (Ecr_eff < limit) THEN
          Dcr(1,1) = -((sigcr0_r - sigcr0)**2) / (2.0d0 * GIc / L)
      ELSE
          Dcr(1,1) = 1.0d-8 ! Fully Softened
      END IF
      
      ! Shear Coupling (Approximate)
      Dcr(2,2) = Dcr(1,1) * (1.0d0 / (2.0d0 * (1.0d0 + nu))) * 1.0d-5
      Dcr(3,3) = Dcr(2,2)

      RETURN
      END

!-----------------------------------------------------------------------
! VII) Consistent Tangent Operator (Jacobian)
!-----------------------------------------------------------------------
      SUBROUTINE calcDcocr(Dco, Dcr, N, Dcocr)
      IMPLICIT NONE

      REAL*8:: Dco(6,6), Dcr(3,3), N(6,3)
      REAL*8:: Dcocr(6,6)
      REAL*8:: invTerm33(3,3), tmp63(6,3), NT(3,6)
  
      NT = TRANSPOSE(N)
      ! Woodbury-style identity for damaged stiffness
      CALL matrixInverse(Dcr + MATMUL(NT, MATMUL(Dco, N)), invTerm33, 3, 3)
      
      tmp63 = MATMUL(Dco, MATMUL(N, invTerm33))
      Dcocr = Dco - MATMUL(tmp63, MATMUL(NT, Dco))

      RETURN
      END