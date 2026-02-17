C     ! Isotropic Incremental SCA 1D  (Damping==0)
C     ! -----------------------------------------------------
C     ! Interface to use smeared crack approach within Abaqus/Standard
C     ! -----------------------------------------------------
C     ! Last Debugged: 01/29/2026  
C     ! meka1@purdue.edu

C     Description: AFME is to improve the convergence and substitute with
C                  Newton-Raphson Method

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)



C     TOTAL FORMULATION
C     if using meka_formulation_AFEM4.for or meka_formulation_AFEM6.for
C       DIMENSION EPS(6)
      
C       EPS = STRAN+DSTRAN
     
C       call  sca3dortho(NSTATV, NPROPS, NOEL, TIME(1), TIME(2),
C      1     DTIME, CELENT, PROPS, EPS, STRESS, STATEV, DDSDDE, DTEMP,1)  


C-----------------------------------------------------------------------
C     INCREMENTAL FORMULATION
C C     if using meka_formulation_AFEM5.for  or meka_formulation_AFEM6.for
C       call  sca3dortho(NSTATV, NPROPS, NOEL, TIME(1), TIME(2),
C      1     DTIME, CELENT, PROPS, STRAN, DSTRAN, STRESS, STATEV, 
C      2     DDSDDE, DTEMP,1)  


C     if using meka_formulation_AFEM5.for 
      call  sca1diso(NSTATV, NPROPS, NOEL, TIME(1), TIME(2),
     1     DTIME, CELENT, PROPS, STRAN, DSTRAN, STRESS, STATEV, 
     2     DDSDDE, DTEMP,1)  

C-----------------------------------------------------------------------
      RETURN
      END


C     SHEAR 23 not accounted
C     Crack Planes are Determined "SIMPLY"
C       include 'meka_formulationLatest.for'

C     CRACK PLANES for 22,12,13 are based on DZ code
C     MAX STRESS FAILURE CRTIERIA
C       include 'meka_formulation_MAXS.for'
      include 'meka_formulation_AFEM5.for'
C       include 'standard_support.for'
      include 'matrix_inverse.for'