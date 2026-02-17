C***********************************************************************
C     ! Isotropic Incremental SCA 1D 
C     ! ----------------------------------------------------------------
C     ! Contains the code concerned with the formulation for SCA in 
C     ! isotropic material, this is used in the 'main_UD_tot.for' file
C     ! ----------------------------------------------------------------
C     ! Last Debugged: 01/29/2026 
C     ! meka1@purdue.edu

C Failure Criteria: 
C     -0) Max Principal Stress Failure Criteria (INDEX = 0)
C     -1) Max Stress Failure Criteria           (INDEX = 1)


C***********************************************************************
C=======================================================================
C     CURRENT ISSUES 
C     1) UNLOADING/RELOADING "NOT" CONSIDERED


C***********************************************************************
      SUBROUTINE sca1diso(NSTATV, NPROPS, NOEL, stepTIME, totalTIME, 
     1 DT, L, PROPS, EPS, DSTRAN, SIG, STATEV, DDSDDE, DTEMP, isIMPLICIT)
          
      IMPLICIT NONE
      
C     STANDARD VARIABLES
      INTEGER::NSTATV,NPROPS,NOEL
      REAL*8::stepTIME,totalTIME,DT,L
      REAL*8::PROPS, SIG(6)
      REAL*8::EPS(6), DSTRAN(6)
      REAL*8::STATEV
      REAL*8::DDSDDE(6,6)      
      REAL*8::DTEMP
      INTEGER::isIMPLICIT
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)

C     DEFINED VARIABLES
      INTEGER::INDEX
      REAL*8::Dco(6,6),Sco(6,6)
      REAL*8::PRINS(3),T(3,3), N(6,3)
      REAL*8::E,nu,XT,GIc,sigcr0,sigcr0_r,term
      REAL*8::isCRACKED, MODE
      REAL*8::Ecr_MAX(3),Ecr_OLD(3)
      REAL*8::Ecr(3),dEcr(3)
      REAL*8::FI_FT
      INTEGER::ZONE

      REAL*8::EPSN, EPSF, ALF1, ALF2
      REAL*8::SIG_CORR,SIG_COMP
      REAL*8::EPS_STAR, DEPS1, DEPS2
      REAL*8::BL
      INTEGER::NTENS, K1

      REAL*8::term33(3,3),invTerm33(3,3)
      REAL*8::term63(6,3)
      REAL*8::term66(6,6)
      REAL*8::Dcocr(6,6)
      REAL*8::Dcr(3,3)

      NTENS = 6

C***********************************************************************
C       !State Variables List
C-----------------------------------------------------------------------
C       (1) Crack Flag - 0 (not failed); 1 (failed)
C       (2) MODE       - Mode of Failure
C       (3) Ecr_MAX    - stores the max strain
C       (4) ZONE       - store the loading zone in the traction separation law
C       (5) : store the stress at the maximum strain encountered (EPS_MAX)
C       (6) SLOPE      - store the slope value based on the zone 
C       (7) : stores the max stress 

C       (10)CELENT: Length of element
C-----------------------------------------------------------------------
C       !Property Variables List
C-----------------------------------------------------------------------
C       (1) Elastic Modulus
C       (2) Poisson's ratio
C       (3) Fracture Toughness-Fibre  (Tension)
C       (4) Tensile Strength      (Longitudinal)
C       (5) Failure Criteria (0 - Max Principal Stress; 1-Max Stress)
C***********************************************************************
C     !PROPERTIES LIST
      E      = props(1) ! Young's Modulus
      nu     = props(2) ! Poissons Ratio
      GIc    = props(3) ! Fracture Toughness
      XT     = props(4) ! Tensile Strength
      INDEX  = props(5) ! Failure Criteria

      sigcr0 = XT ! CRITICAL STRESS

      ! PARAMETERS
      EPSN = XT/E                    !Global Strain @ Peak
      EPSF = (2.0d0*GIc)/(XT*L)    !Global Strain @ Failure

      ALF1 = E                       ! Pre-Peak Global Slope
      ALF2 = (-1.0d0*XT)/(EPSF-EPSN) ! Post-Peak Global Softening Slope

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     ! STATE VARIABLES LIST (SDV)
      isCRACKED = STATEV(1)  ! CRACK FLAG
      MODE      = STATEV(2)  ! MODE OF FAILURE
      Ecr_MAX   = STATEV(3)  ! Stores the max local crack strain
      Ecr_OLD   = STATEV(22) ! Stores the OLD local crack strain
      ZONE      = STATEV(4)  ! Allot zone in Traction-Separation Law
      STATEV(5) = L          ! Element Length

      BL        = (2.0d0*GIc*E/sigcr0**2) !Bazant Limit for Element Size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! set a very small non-zero value for max local crack strain encountered
      statev(3) = 1.0d-8
      statev(6) = 1.0d-8
      statev(7) = 1.0d-8
      statev(8) = 1.0d-8
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C     STEP-1: Calculate Elastic Stiffness Matrix Dco    
      CALL STIFFNESS_MATRIX(E, nu, Dco)

C     STEP-2: INITIALIZE STRAIN TENSOR
      DO K1=1,NTENS
            EPS(K1) = EPS(K1) + DSTRAN(K1)
      END DO

C     STEP-3: CHECK IF ELEMENT SIZE IS SMALL ENOUGH 
      !check if element size is small enough
      if ((2.0d0*GIc*E/sigcr0**2).lt.L) then
         write(*,*) ' '
         write(*,*) 'BAZANT LIMIT OVERSHOOT!'
         write(*,*) '(2.0d0*GIc*E/sigcr0**2) = ',BL
         write(*,*) 'Current Element Size = ',L
         write(*,*) 'Try Element Size: ',BL
      call myExit()      
      end if   

C     STEP-4: ASSUME "NO" Cracks @ START
      if (isCracked.lt.0.5) then
        ! linear elasticity until crack      
        SIG = SIG + MATMUL(Dco,DSTRAN)

        ! JACOBIAN UPDATE
        DDSDDE = Dco


C     STEP-5: Check if Crack is Initiated (Failure Criteria)
      if (INDEX .eq. 0.0d0) then ! MAX PRINCIPAL STRESS FAIL CRIT
      write(*,*) "Max Principal Stress Failure Criteria"

C ****************  START FAILURE CRITERIA WRAPPER *********************
C     INPUT :  Take in STRESS STATE and STRENGTH VALUES
C     OUTPUT: 'T' matrix i.e. orientation of crack to compute N

      ! Compute PRINCIPAL STRESSES AND PRINCIPAL DIRECTIONS
      call mysprind(SIG,PRINS,T,1,3,3)

      if (dmax1(PRINS(1),PRINS(2),PRINS(3)) .gt. sigcr0) then
      SIG = SIG-matmul(Dco,DSTRAN)
      isCracked     = 1.0d0
      statev(1)     = isCracked
      statev(2)     = 1.0d0
      statev(9:11)  = T(1,:)
      statev(12:14) = T(2,:)
      statev(15:17) = T(3,:)
      endif
      endif
C ****************  END FAILURE CRITERIA WRAPPER ***********************
      end if

C-----------------------------------------------------------------------
C     STEP-6: If Crack Exists
      if (isCracked.ge.0.99d0) then
      write(*,*) "Cracked!"
      Ecr_MAX  = statev( 6: 8)  ! MAX LOCAL CRACK STRAIN ENCOUNTERED
      Ecr_OLD  = statev(18:20)  ! OLD LOCAL CRACK STRAIN 
      Ecr      = statev(21:23)  ! CURRENT LOCAL CRACK STRAIN  
      T(1,:)   = statev( 9: 11)
      T(2,:)   = statev( 12: 14)
      T(3,:)   = statev( 15: 17)

      !Compute the transformation matrix
      call getN(T,N)

C-----------------------------------------------------------------------
C     STEP-7: Compute Dcr matrix
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! HARD CODED FOR THE SAKE OF TESTING

C       sigcr0_r = sigcr0*1.0D-5
C       term     = 1.0D0/(2.0D0*(1.0D0+nu))

C       Dcr(1,1) = (-(sigcr0_r-sigcr0)**2/(2.0Q0*GIc/L))
C       Dcr(2,2) = Dcr(1,1)*term*1.0D-5
C       Dcr(3,3) = Dcr(1,1)*term*1.0D-5

      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      CALL calcDcr(Ecr_OLD,Ecr_MAX,GIc,sigcr0,L,nu,Dcr)



C-----------------------------------------------------------------------
C     STEP-8: Determine local crack strain (increment) consistently

      !INPUT : Ecr_MAX,Ecr_OLD,DSTRAN,Dco,N, GIc, sigcr0, L, nu
      !OUTPUT: Ecr, Ecr_MAX, Ecr_OLD

      CALL calcEcr(Ecr_OLD,Ecr_MAX,DSTRAN,Dco,N,dEcr,Ecr,GIc,sigcr0,L,nu)

      STATEV(6:8)  = Ecr_MAX 
      STATEV(18:20)= Ecr_OLD
      STATEV(21:23)= Ecr


      STATEV(24) = Dcr(1,1)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C     STEP-10: Compute Dcocr i.e. Jacobian
      CALL calcDcocr(Dco,Dcr,N,Dcocr)

      !_______________________
      ! DSIGMA = Dcocr*DSTRAN
      !_______________________
      !calculate the stress
      SIG = SIG + matmul(Dcocr,DSTRAN)

      ! Jacobian update
      if (isImplicit.eq.1) DDSDDE=Dcocr

      end if

      RETURN
      END
     
C=======================================================================     
C     DEFINE FUNCTIONS

C-----------------------------------------------------------------------
C     I) Computing Stiffness Matrix (ISOTROPIC)
C-----------------------------------------------------------------------
      SUBROUTINE STIFFNESS_MATRIX(E, nu, Dco)
      IMPLICIT NONE
      REAL*8::E,nu
      REAL*8::Dco(6,6)
      REAL*8::mu, twomu, lambda 

C     Isotropic Material 
      mu     = 0.5*E/( 1.0d0 + nu )
      twomu  = 2.0d0*mu
      lambda = E*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
      
      Dco=reshape(
     1 (/lambda+ twomu,  lambda, lambda, 0.0d0, 0.0d0, 0.0d0, 
     2   lambda, lambda+twomu,  lambda, 0.0d0, 0.0d0, 0.0d0, 
     3   lambda, lambda, lambda+twomu,  0.0d0, 0.0d0, 0.0d0, 
     4   0.0d0,  0.0d0,  0.0d0,  mu,    0.0d0, 0.0d0,
     5   0.0d0,  0.0d0,  0.0d0,  0.0d0, mu,    0.0d0,
     6   0.0d0,  0.0d0,  0.0d0,  0.0d0, 0.0d0, mu 
     7 /),(/6,6/)) 

      RETURN
      END

C-----------------------------------------------------------------------
C     II) Computing N matrix (1 crack)
C-----------------------------------------------------------------------
      subroutine getN(o,N)
      implicit none

!     ordering of principal direction: o(index of eval,1:3)   
C     o(1,:) = (/1.0d0, 0.0d0, 0.0d0/)
C     o(2,:) = (/0.0d0, 1.0d0, 0.0d0/)
C     o(3,:) = (/1.0d0, 0.0d0, 1.0d0/)

      real*8::o(3,3)
      real*8::N(6,3)
      
      N(1,1) = o(1,1) ** 2
      N(1,2) = 0.1D1 * o(1,1) * o(2,1)
      N(1,3) = 0.1D1 * o(3,1) * o(1,1)
      N(2,1) = o(1,2) ** 2
      N(2,2) = 0.1D1 * o(1,2) * o(2,2)
      N(2,3) = 0.1D1 * o(3,2) * o(1,2)
      N(3,1) = o(1,3) ** 2
      N(3,2) = 0.1D1 * o(1,3) * o(2,3)
      N(3,3) = 0.1D1 * o(3,3) * o(1,3)
      N(4,1) = 0.2D1 * o(1,1) * o(1,2)
      N(4,2) = o(1,1) * o(2,2) + o(1,2) * o(2,1)
      N(4,3) = o(3,1) * o(1,2) + o(3,2) * o(1,1)
      N(5,1) = 0.2D1 * o(1,3) * o(1,1)
      N(5,2) = o(1,3) * o(2,1) + o(1,1) * o(2,3)
      N(5,3) = o(3,3) * o(1,1) + o(3,1) * o(1,3)
      N(6,1) = 0.2D1 * o(1,2) * o(1,3)
      N(6,2) = o(1,2) * o(2,3) + o(1,3) * o(2,2)
      N(6,3) = o(3,2) * o(1,3) + o(3,3) * o(1,2)
      return
      end

C----------------------------------------------------------------------- 
C     (III) Finding principle directions
C----------------------------------------------------------------------- 
      subroutine mysprind(S,PS,AN,LSTR,NDI,NSHR)
      implicit none
      integer::NDI,NSHR,LSTR
      real*8::S(NDI+NSHR)
      real*8::PS(3)
      real*8::AN(3,3)
      call SPRIND(S,PS,AN,LSTR,NDI,NSHR)
      return 
      end
      
C----------------------------------------------------------------------- 
C     (IV) Exitting
C----------------------------------------------------------------------- 
      subroutine myExit()
      call XIT
      return
      end

C----------------------------------------------------------------------- 
C     (V) Calculating Incremental e_cr
C-----------------------------------------------------------------------     
      subroutine calcEcr(Ecr_OLD,Ecr_MAX,DSTRAN,Dco,N,dEcr,Ecr,
     1       GIc,sigcr0,L,nu)


      !INPUT : Ecr_OLD,Ecr_MAX,DSTRAN,Dco,N, GIc, sigcr0, L, nu
      !OUTPUT: Ecr, Ecr_MAX, Ecr_OLD
      
      implicit none

      real*8::Ecr_MAX(3),Ecr_OLD(3),DSTRAN(6)
      real*8::Dco(6,6),Dcr(3,3)
      real*8::N(6,3), NT(3,6)
      real*8::Ecr(3),dEcr(3)
      real*8::GIc,sigcr0,L,nu

      real*8::NT_Dco_N(3,3)
      real*8::NT_Dco(3,6)
      real*8::Y(3,3),invY(3,3)

      NT        = transpose(N)

C     !Calculate the Dcr matrix based on previous Ecr value and 
      !Ecr_MAX value
      call calcDcr(Ecr_OLD,Ecr_MAX,GIc,sigcr0,L,nu,Dcr)
      
      !precalculate some stuff

      !N^T*D_co*N
      NT_Dco_N=matmul(transpose(N),matmul(Dco,N))
      
      !N^T*D_co
      NT_Dco=matmul(transpose(N),Dco)

      !Dcr + NT_Dco_N
      Y = Dcr+NT_Dco_N

      ! invqY = inverse of (Dcr + NT_Dco_N)
      CALL matrixInverse(Y,invY,3,3)

      ! Increment of Local Crack Strain (Ecr)
      dEcr = matmul(invY,matmul(NT_Dco,DSTRAN))

      ! Compute the NEW Local Crack Strain
      Ecr = Ecr_OLD + dEcr

      ! UPDATE THE MAX CRACK STRAIN ENCOUNTERED
      Ecr_MAX = max(Ecr,Ecr_MAX) 

      ! SET THE COMPUTED CRACK STRAIN AS THE OLD CRACK STRAIN
      Ecr_OLD = Ecr
      
      return
      end

C----------------------------------------------------------------------- 
C     (VI) Calculating Dcr matrix
C-----------------------------------------------------------------------   
      subroutine calcDcr(Ecr_OLD,Ecr_MAX,GIc,sigcr0,L,nu,Dcr)
      implicit none

      real*8::Ecr_OLD(3),Ecr_MAX(3),GIc,sigcr0,L,Dcr(3,3),nu
      real*8::term,sigcr0_r
      real*8::Ecr_TEMP(3),Dcr_min

      sigcr0_r = sigcr0*5D-5
      term     = 1.0D0/(2.0D0*(1.0D0+nu))


      Ecr_TEMP(1) = max(abs(Ecr_OLD(1)),abs(Ecr_MAX(1)))

      Dcr      = 0.0D0

      IF (Ecr_TEMP(1).LT.(2.0Q0*GIC/L/(sigcr0-sigcr0_r))) THEN
            Dcr(1,1) = (-(sigcr0_r-sigcr0)**2/(2.0D0*GIc/L))
      ELSE
            Dcr(1,1) = 1.0D-8
      END IF
      Dcr(2,2) = Dcr(1,1)*term*1.0D-5
      Dcr(3,3) = Dcr(1,1)*term*1.0D-5

      return
      end

C----------------------------------------------------------------------- 
C     (VII) Calculating Dcocr matrix
C-----------------------------------------------------------------------

! This subroutine is not working ~ Feb17,2026
      subroutine calcDcocr(Dco,Dcr,N,Dcocr)
      implicit none

      real*8::Dcr(3,3),Dco(6,6),N(6,3),Dcocr(6,6) 
      real*8::term33(3,3),invTerm33(3,3),term63(6,3),term66(6,6)   

      ! Compute Dcocr
      term33=Dcr+matmul(transpose(N),matmul(Dco,N)) ! Dcr + N^T.Dco.N
      CALL matrixInverse(term33,invTerm33,3,3)      ! [Dcr + N^T.Dco.N]^-1
      
      term63=matmul(Dco,matmul(N,invTerm33))        ! Dco.N.[Dcr + N^T.Dco.N]^-1
      term66=matmul(term63,matmul(transpose(N),Dco))! Dco.N.([Dcr + N^T.Dco.N]^-1).N^T.Dco
      
      Dcocr=Dco-term66   ! (Dco - Dco.N.([Dcr + N^T.Dco.N]^-1).N^T.Dco)

      return
      end
























      
      