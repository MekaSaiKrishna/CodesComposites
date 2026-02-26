C     ! Isotropic Incremental SCA [1 crack]  (Damping==0)

C     ! FOR 2D MODELS (CPS4 - PLANE STRESS ELEMENTS)
C     ! -----------------------------------------------------
C     ! Interface to use smeared crack approach within Abaqus/Standard
C     ! -----------------------------------------------------
C     ! Last Debugged: 02/25/2026  
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


C-----------------------------------------------------------------------
C     INCREMENTAL FORMULATION
C     if using meka_formulation_AFEM5.for 
      call  sca2diso(NTENS, NDI, NSHR, DFGRD1, NSTATV, NPROPS, NOEL, 
     1     TIME(1), TIME(2), DTIME, PNEWDT, CELENT, PROPS, STRAN, DSTRAN, STRESS,  
     2     STATEV, DDSDDE, DTEMP,1)  

C-----------------------------------------------------------------------
      RETURN
      END

C     MAX PRINCIPAL STRESS FAILURE CRTIERIA
      include 'formulation_AFEM2D.for'
      include 'matrix_inverse.for'