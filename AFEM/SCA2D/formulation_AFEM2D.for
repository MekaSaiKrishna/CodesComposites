C***********************************************************************
C     ! Isotropic Incremental SCA 1D 
C     ! ----------------------------------------------------------------
C     ! Contains the code concerned with the formulation for SCA in 
C     ! isotropic material, this is used in the 'main_UD_tot.for' file
C     ! 
C     ! -> UNLOADING/RELOADING ACCOUNTED
C     ! ----------------------------------------------------------------
C     ! Last Debugged: 02/17/2026 
C     ! meka1@purdue.edu

C Failure Criteria: 
C     -0) Max Principal Stress Failure Criteria (INDEX = 0)
C     -1) Max Stress Failure Criteria           (INDEX = 1)


C***********************************************************************
C=======================================================================
C     CURRENT ISSUES 


C***********************************************************************
      SUBROUTINE sca2diso(NTENS, NDI, NSHR, DFGRD1, NSTATV, NPROPS, NOEL,  
     1      stepTIME, totalTIME, DT, PNEWDT, L, PROPS, EPS, DSTRAN, SIG, STATEV, 
     2      DDSDDE, DTEMP, isIMPLICIT)
          
      IMPLICIT NONE
      
C     STANDARD VARIABLES
      INTEGER::NTENS, NDI, NSHR, DFGRD1(3,3)
      INTEGER::NSTATV,NPROPS,NOEL
      REAL*8::stepTIME,totalTIME,DT,PNEWDT,L,Leff
      REAL*8::PROPS, SIG(NTENS)
      REAL*8::EPS(NTENS), DSTRAN(NTENS)
      REAL*8::STATEV
      REAL*8::DDSDDE(NTENS,NTENS)      
      REAL*8::DTEMP
      INTEGER::isIMPLICIT
      DIMENSION STATEV(NSTATV),PROPS(NPROPS)

C     DEFINED VARIABLES
      INTEGER::INDEX
      REAL*8::Dco(NTENS,NTENS),Sco(NTENS,NTENS)
      REAL*8::PRINS(3),T(3,3)
      REAL*8::E,nu,XT,GIc,sigcr0,sigcr0_r,term
      REAL*8::isCRACKED, MODE
      REAL*8::FI_T, FI_C, FI_MAX, J_FAIL, J_MIN, PMIN
      INTEGER::I,J 
      INTEGER::ZONE

      REAL*8::BL
      INTEGER::K1

      REAL*8::N(3,2)
      REAL*8::Ecr_MAX(2),Ecr_OLD(2)
      REAL*8::Ecr(2),dEcr(2)
      REAL*8::scr_old(2),scr(2)

      REAL*8::Dcocr(2,2)
      REAL*8::Dcr(2,2)

C***********************************************************************
C       !State Variables List
C-----------------------------------------------------------------------
C       (1) Crack Flag - 0 (not failed); 1 (failed)

C-----------------------------------------------------------------------
C       !Property Variables List
C-----------------------------------------------------------------------
C       (1) Elastic Modulus
C       (2) Poisson's ratio
C       (3) Fracture Toughness-Fibre  (Tension)
C       (4) Tensile Strength          (Longitudinal)
C       (5) Failure Criteria          (0 - Max Principal Stress; 1-Max Stress)
C***********************************************************************
C     !PROPERTIES LIST
      E      = props(1) ! Young's Modulus
      nu     = props(2) ! Poissons Ratio
      GIc    = props(3) ! Fracture Toughness
      XT     = props(4) ! Tensile Strength
      INDEX  = props(5) ! Failure Criteria

      sigcr0 = XT ! CRITICAL STRESS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     ! STATE VARIABLES LIST (SDV)
      isCRACKED = STATEV(1)  ! CRACK FLAG
C       MODE      = STATEV(2)  ! MODE OF FAILURE

      STATEV(2) = L          ! Element Length

      Leff = 1.0d0*L ! IF LINEAR ELEMENT    (1st Order)
C       Leff = 2.0d0*L ! IF QUADRATIC ELEMENT (2nd Order)

      BL        = (2.0d0*GIc*E/sigcr0**2) !Bazant Limit for Element Size

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

C     STEP-1: Calculate Elastic Stiffness Matrix Dco    
      CALL STIFFNESS_MATRIX_2D_PS(E, nu, Dco)

C     STEP-2: INITIALIZE STRAIN TENSOR
      DO K1=1,NTENS
            EPS(K1) = EPS(K1) + DSTRAN(K1)
      END DO

C     STEP-3: CHECK IF ELEMENT SIZE IS SMALL ENOUGH 
      !check if element size is small enough
      if ((2.0d0*GIc*E/sigcr0**2).lt.Leff) then
         write(*,*) ' '
         write(*,*) 'BAZANT LIMIT OVERSHOOT!'
         write(*,*) '(2.0d0*GIc*E/sigcr0**2) = ',BL
         write(*,*) 'Current Element Size = ',Leff
         write(*,*) 'Try Element Size: ',BL
      call myExit()      
      end if   


C     STEP-4: ASSUME "NO" Cracks @ START
      if (isCracked.lt.0.5) then

      ! linear elasticity until crack      
      SIG = SIG + MATMUL(Dco,DSTRAN)

      ! JACOBIAN UPDATE
      DDSDDE = 1.0d0*Dco

C     STEP-5: Check if Crack is Initiated (Failure Criteria)

C ****************  START FAILURE CRITERIA WRAPPER *********************
C     INPUT :  Take in STRESS STATE and STRENGTH VALUES
C     OUTPUT: 'T' matrix i.e. orientation of crack to compute N

      ! Compute PRINCIPAL STRESSES AND PRINCIPAL DIRECTIONS
      call mysprind(SIG,PRINS,T,1,2,1)

      ! 1. Determine failure index and tracking index J_FAIL
      ! Initialize
      FI_T   = 0.0D0 ! Tension Failure Index
      FI_C   = 0.0D0 ! Compression Failure Index
      J_FAIL = 1     ! Default to first principal direction

      ! Tensile Check
      IF (PRINS(1) .GT. 0.0D0) FI_T = PRINS(1)/XT

      ! Compressive Check (finding the most negative)
      PMIN  = PRINS(1)
      J_MIN = 1

      DO I = 2, 3
         IF (PRINS(I) .LT. PMIN) THEN
            PMIN = PRINS(I)
            J_MIN = I
         END IF
      END DO

      IF (PMIN .LT. 0.0D0) FI_C = ABS(PMIN)/XT

      ! 2. Identify the Dominant Mode and store the Vector
      IF (FI_T .GE. FI_C) THEN
         ! TENSILE MODE
         MODE   = 1.0d0
         FI_MAX = FI_T
         J_FAIL = 1  ! Direction of max tensile stress
      ELSE
         ! COMPRESSIVE MODE
         MODE   = -1.0d0
         FI_MAX = FI_C
         J_FAIL = J_MIN ! Direction of max compressive stress
      END IF

      ! 3. Store in State Variables (STATEV)
      ! STATEV(10): Max Failure Index
      ! STATEV(11): X-component of failure vector
      ! STATEV(12): Y-component of failure vector
      ! Z - COMPONENT OF FAILURE VECTOR is 0 for IN-PLANE FAILURE

      STATEV(3)  = MODE   ! +1 (Tensile), -1 (Compressive)
      STATEV(4) = FI_MAX
      STATEV(5) = T(1, J_FAIL) !--> cos(theta)
      STATEV(6) = T(2, J_FAIL) !--> sin(theta)


      IF (FI_MAX .GE. 1.0D0) THEN
      STATEV(1) = 1.0D0  ! Set the failure flag
      ! Optional: Print a message to the .log file for debugging
      ! WRITE(7,*) 'Failure initiated at Element:', NOEL
C       statev(14) = FI_MAX*XT
      
      !Compute the transformation matri
      call getN(STATEV(5),STATEV(6),N)
      
C       scr_old = matmul(transpose(N),SIG) 
      statev(14:15) = matmul(transpose(N),SIG) 
      
      END IF


C C ****************  END FAILURE CRITERIA WRAPPER ***********************
      ! Max Local Crack Strain (initial)
      statev(7)  = 1.0d-8
      statev(8)  = 1.0d-8

      ! OLD Local Crack Strain (initial)      
      statev(9)  = 1.0d-8
      statev(10) = 1.0d-8
      
      ! CURRENT INCREMENT IN Local Crack Strain
      statev(11) = 1.0d-8
      statev(12) = 1.0d-8

      ! local crack stress
      statev(13) = FI_MAX*XT !initialization to STRENGTH

      end if

C     STEP-6: If Crack Exists
      IF (STATEV(1) .GT. 0.5D0) THEN       

      Ecr_MAX  = statev(7:8)   ! MAX LOCAL CRACK STRAIN ENCOUNTERED
      Ecr_OLD  = statev(9:10)  ! OLD LOCAL CRACK STRAIN 
      dEcr     = statev(11:12) ! CURRENT INCREMENT IN LOCAL CRACK STRAIN  


C     STEP-7: Compute Dcr matrix
      scr_old = statev(14:15)

      CALL calcDcr(Ecr_OLD,Ecr_MAX,GIc,sigcr0,Leff,nu,Dcr,scr_old,scr,dEcr)
      statev(14:15) = scr


C     STEP-8: Determine local crack strain (increment) consistently

      !INPUT : Ecr_MAX,Ecr_OLD,DSTRAN,Dco,N, GIc, sigcr0, L, nu
      !OUTPUT: Ecr, Ecr_MAX, Ecr_OLD

      CALL calcEcr(Ecr_OLD,Ecr_MAX,DSTRAN,Dco,N,dEcr,Ecr,GIc,sigcr0,Leff,nu,scr_old,scr)

      STATEV(7:8)   = Ecr_MAX 
      STATEV(9:10)  = Ecr_OLD
      STATEV(11:12) = dEcr

C       STATEV(15)   = Dcr(1,1)

C     STEP-9: Compute Dcocr i.e. Jacobian
      CALL calcDcocr(Dco,Dcr,N,Dcocr)

      !calculate the stress
      SIG = SIG + matmul(Dcocr,DSTRAN)

      ! Jacobian update
      if (isImplicit.eq.1) DDSDDE=Dcocr
      
      END IF

      RETURN
      END
     
C=======================================================================     
C     DEFINE SUBROUTINE FUNCTIONS

C-----------------------------------------------------------------------
C     I) Computing Stiffness Matrix (ISOTROPIC)
C-----------------------------------------------------------------------
      SUBROUTINE STIFFNESS_MATRIX_2D_PS(E, nu, Dco)
      IMPLICIT NONE
      REAL*8::E,nu
      REAL*8::Dco(3,3)
      REAL*8::const

C     Isotropic Material 
      const = E/(1.0d0 - nu**2.0d0)

      Dco=reshape(
     1 (/const*1.0d0,  const*nu, 0.0d0, 
     2   const*nu, const*1.0d0, 0.0d0, 
     3   0.0d0, 0.0d0, const*0.50d0*(1-nu)
     4 /),(/3,3/))

      RETURN
      END

C-----------------------------------------------------------------------
C     II) Computing N matrix (1 crack)
C-----------------------------------------------------------------------
      subroutine getN(c,s,N)
      implicit none

      real*8::c, s
      real*8::N(3,2)
      
      N(1,1) = c**2
      N(1,2) = -1.0D0*c*s
      N(2,1) = s**2
      N(2,2) = 1.0D0*c*s
      N(3,1) = 2.0d0*s*c
      N(3,2) = (c**2) - (s**2)
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
     1       GIc,sigcr0,Leff,nu,scr_old,scr)


      !INPUT : Ecr_OLD,Ecr_MAX,DSTRAN,Dco,N, GIc, sigcr0, L, nu
      !OUTPUT: Ecr, Ecr_MAX, Ecr_OLD
      
      implicit none

      real*8::Ecr_MAX(2),Ecr_OLD(2),DSTRAN(3)
      real*8::Dco(3,3),Dcr(2,2)
      real*8::N(3,2), NT(2,3)
      real*8::Ecr(2),dEcr(2)
      real*8::GIc,sigcr0,Leff,nu

      real*8::NT_Dco_N(2,2)
      real*8::NT_Dco(2,3)
      real*8::Y(2,2),invY(2,2)

      real*8::scr_old(2),scr(2) 

      NT = transpose(N)

C     !Calculate the Dcr matrix based on previous Ecr value and 
      !Ecr_MAX value
      call calcDcr(Ecr_OLD,Ecr_MAX,GIc,sigcr0,Leff,nu,Dcr,scr_old,scr,dEcr)
      
      !precalculate some stuff
      NT_Dco_N=matmul(transpose(N),matmul(Dco,N))
      
      NT_Dco=matmul(transpose(N),Dco)

      Y = Dcr+NT_Dco_N

      ! invqY = inverse of (Dcr + NT_Dco_N)
      CALL matrixInverse(Y,invY,2,2)

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
      subroutine calcDcr(Ecr_OLD,Ecr_MAX,GIc,sigcr0,Leff,nu,Dcr,scr_old,scr,dEcr)
      implicit none

      real*8::Ecr_OLD(2),Ecr_MAX(2),GIc,sigcr0,Leff,Dcr(2,2),nu
      real*8::term,sigcr0_r
      real*8::Ecr_TEMP(2),Dcr_min

      real*8::scr_old(2),scr(2),dEcr(2)
      real*8::scr_max !local crack stress at max local crack strain

      sigcr0_r = sigcr0*5D-5
      term     = 1.0D0/(2.0D0*(1.0D0+nu))


      !local crack stress @ max local crack strain (LINEAR_TRACTION_SEPARATION_LAW)
      scr_max = sigcr0 + (-(sigcr0)**2/(2.0D0*GIc/Leff))*Ecr_MAX(1)


      Ecr_TEMP(1) = max(abs(Ecr_OLD(1)),abs(Ecr_MAX(1)))

      Dcr      = 0.0D0

      ! BEFORE COMPLETE FAILURE
      IF (Ecr_TEMP(1).LT.(2.0Q0*GIc/Leff/(sigcr0-sigcr0_r))) THEN
            ! LOADING 
            IF (Ecr_OLD(1) .LT. Ecr_MAX(1)) THEN
                  Dcr(1,1) =  scr_max/abs(Ecr_MAX(1)) 
            ! UNLOADING/RELOADING
            ELSE
                  Dcr(1,1) = (-(sigcr0_r-sigcr0)**2/(2.0D0*GIc/Leff))
            END IF
      ! AFTER COMPLETE FAILURE
      ELSE
            Dcr(1,1) = 1.0D-8
      END IF

      Dcr(2,2) = Dcr(1,1)*term*1.0D-5

      scr = scr_old + matmul(Dcr,dEcr)
      return
      end

C----------------------------------------------------------------------- 
C     (VII) Calculating Dcocr matrix
C-----------------------------------------------------------------------
      subroutine calcDcocr(Dco,Dcr,N,Dcocr)
      implicit none

      real*8::Dcr(2,2),Dco(3,3),N(3,2),Dcocr(3,3) 
      real*8::term22(2,2),invTerm22(2,2)
      real*8::term32(3,2),term33(3,3)   

      term22=Dcr+matmul(transpose(N),matmul(Dco,N)) ! Dcr + N^T.Dco.N
      CALL matrixInverse(term22,invTerm22,2,2)      ! [Dcr + N^T.Dco.N]^-1
      
      term32=matmul(Dco,matmul(N,invTerm22))        ! Dco.N.[Dcr + N^T.Dco.N]^-1
      term33=matmul(term32,matmul(transpose(N),Dco))! Dco.N.([Dcr + N^T.Dco.N]^-1).N^T.Dco
      
      Dcocr=Dco-term33   ! (Dco - Dco.N.([Dcr + N^T.Dco.N]^-1).N^T.Dco)

      return
      end
























      
      