! things to address:
! - compression
! - crack closing
! - automatically find element size
! - include check for element size

      subroutine sca2d(nstatev, nprops,noel,npt,layer,
     *     stepTime, totalTime, dt, L,
     *     props, eps, sig, statev, ddsdde, isImplicit)


      
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !Engineering Ordering
      !eps11,eps22,gam12
      
      implicit none
      
      integer::nstatev,nprops,noel,layer,npt,isImplicit
      real*8::stepTime,totalTime,dt,L
      real*8::props(nprops)
      real*8::eps(3)
      real*8::sig(3)
      real*8::statev(nstatev)
      real*8::ddsdde(3,3)
      
      
      real*8::E11
      real*8::E22
      real*8::G12
      real*8::nu12
      real*8::sigcrf0, sigcrf
      real*8::sigcrm0, sigcrm
      real*8::taucrm0, taucrm
      real*8::GICf
      real*8::GICm
      real*8::GIICm
      real*8::Sco(3,3)
      real*8::Dco(3,3)
      real*8::Dcr(4,4)
      real*8::Dda(4,4) !damping matrix
      real*8::Dcocr(3,3)
      real*8::invDcrDco(4,4)
      real*8::epscr(4)
      real*8::epscrold(4)
      real*8::epscrf
      real*8::epscrm
      real*8::gamcrm
      real*8::epscrmax(4)
      real*8::epscrfmax
      real*8::epscrmmax
      real*8::gamcrmmax
      real*8::gamcrfmax !*** The maximum crack strain in history
      real*8::Ecrf, Ecrf0
      real*8::Ecrm, Ecrm0
      real*8::Gcrm, Gcrm0
      real*8::sigOld(3)
      real*8::l1 !characteristic length in the 1-direction
      real*8::l2 !characteristic length in the 2-direction
      real*8::term33(4,3)
      real*8::modeIIadvance
      real*8::dGIf,dGIm,dGIIm
      real*8::damp
      real*8::GGICf !*** Shear modulus after fiber breakage
      real*8::N1(3,4), N2(3,4) !*** The transformation law N
      real*8::alpha  !*** The factor that used to deduct GGICf
      real*8::term11(4,4),term22(4,3)
      
      !Property list
      !( 1) elestic modulus in fiber direction
      !( 2) elastic modulus perpendicular to fiber direction
      !( 3) shear modulus
      !( 4) Poissons ratio
      !( 5) critical stress in fiber direction
      !( 6) critical stress perpendicular to fiber direction
      !( 7) critical stress in shear
      !( 8) fracture energy mode I in fiber direction
      !( 9) fracture energy mode I in matrix direction
      !(10) fracture energy mode II in matrix direction
      !(11) Shear modulus after fiber breakage !***
      
      !State variable list
      !( 1) characteristic length in the 1 direction
      !( 2) characteristic length in the 2 direction
      !( 3) maximum crack strain 11
      !( 4) maximum crack strain 21 (engineering) !*** Fiber shear
      !( 5) maximum crack strain 22
      !( 6) maximum crack strain 12 (engineering)
      !( 7) GIf/GICf
      !( 8) GIm/GICm
      !( 9) GIIm/GIICm
      !( 10) Status for element deletion
      !(11) last crack strain 11
      !(12) last crack strain 21 (enginnering) !***
      !(13) last crack strain 22
      !(14) last crack strain 12 (engineering)

      !Define the transformation law
      N1 = reshape((/1.0d0, 0.0d0, 0.0d0 ,
     *              0.0d0, 0.0d0, 0.0d0 ,
     *              0.0d0, 1.0d0, 0.0d0 ,
     *              0.0d0, 0.0d0, 1.0d0/),(/3,4/)) ! Matrix crack

      N2 = reshape((/1.0d0, 0.0d0, 0.0d0 ,
     *              0.0d0, 0.0d0, 1.0d0 ,
     *              0.0d0, 1.0d0, 0.0d0 ,
     *              0.0d0, 0.0d0, 0.0d0/),(/3,4/)) ! Fiber crack



      !get properties
      E11       =props(1)
      E22       =props(2)
      G12       =props(3)
      nu12      =props(4)
      sigcrf0   =props(5)
      sigcrm0   =props(6)
      taucrm0   =props(7)
      GICf      =props(8)
      GICm      =props(9)
      GIICm     =props(10)
      GGICf     =props(11)
      alpha     =props(12)
      damp      =0.01d0 !!MAKE VAR            

      ! construct (continuum) compliance matrix
      Sco=reshape((/
     *     1.0d0/E11, -nu12/E11, 0.0d0 ,
     *     -nu12/E11, 1.0d0/E22, 0.0d0 ,
     *         0.0d0,     0.0d0, 1.0d0/G12 /),(/3,3/))
      
      
      
      ! invert for stiffness matrix
      call inverse(Sco,3,Dco)      
  
      !Deterime element size
c      l1=statev(1)
c      l2=statev(2)

       l1 = L
       l2 = L
      
      !set material point to not deleted
      statev(10)=1.0d0
      
      !check if elements small enough
      if ((2.0d0*GICf*E11/sigcrf0**2.lt.l1).or.
     *    (2.0d0*GICm*E22/sigcrm0**2.lt.l2).or.
     *    (2.0d0*GIICm*G12/taucrm0**2.lt.l2)) then
      write(*,*) ' '
      write(*,*) '2.0d0*GICf*E11/sigcrf0**2,l1'
      write(*,*) 2.0d0*GICf*E11/sigcrf0**2,l1
      write(*,*) '2.0d0*GICm*E22/sigcrm0**2,l2'
      write(*,*) 2.0d0*GICm*E22/sigcrm0**2,l2
      write(*,*) '2.0d0*GIICm*G12/taucrm0**2,l2'
      write(*,*) 2.0d0*GIICm*G12/taucrm0**2,l2
      write(*,*)
      write(*,*) 'Elements too large for given G vals'
      write(*,*) ' '
      call myExit()
      endif
      
!      if ((l1.le.1e-9).or.(l2.le.1e-9)) then
!      write(*,*) ' '
!      write(*,*) 'No characteristic length prescribed'
!      write(*,*)
!      call myExit()
!      endif
             
      
      !obtain maximum cracking strain that the material has seen 
      !at any point intime
      epscrfmax=statev(3)
      gamcrfmax=statev(4)
      epscrmmax=statev(5)
      gamcrmmax=statev(6)
      
      epscr=0.0d0
      
      !if not crack present:
      if ((epscrfmax.lt.1.d-16).and.(epscrmmax.lt.1d-16).and.
     *    (gamcrmmax.lt.1.d-16).and.(gamcrfmax.lt.1.d-16)) then
      ! linear elasticity: snew=D*Eps
      sig=matmul(Dco,eps)
      if (isImplicit.eq.1) ddsdde=Dco
      !     if new stress larger max stress:
      if((sig(1).gt.sigcrf0).or.
     *   (sig(2).gt.sigcrm0).or.
     *   (dabs(sig(3)).gt.taucrm0)) then
      !     initialize crack strain
      epscr=1.0d-14
      epscrfmax=1.0d-14
      gamcrfmax=1.0d-14
      epscrmmax=1.0d-14
      gamcrmmax=1.0d-14
      endif
!     else
      endif      

      if (((epscrmmax.gt.1d-16).or.(gamcrmmax.gt.1.d-16))
     *  .and.(epscrfmax.le.1.d-8)) then
C       write(*,*) 'Matrix crack!!!!: '
C       write(*,*) epscrfmax,gamcrfmax,epscrmmax,gamcrmmax
      !crack band
      epscrold=(/statev(11),statev(12),statev(13),statev(14)/) !(/epscrf,gamcrf,epscrm,gamcrm/)
      epscrmax=(/epscrfmax,gamcrfmax,epscrmmax,gamcrmmax/)
      !form the damping matrix
      Dda=reshape((/damp,0.0d0,0.0d0,0.0d0,
     # 0.0d0,damp,0.0d0,0.0d0,
     # 0.0d0,0.0d0,damp,0.0d0,
     # 0.0d0,0.0d0,0.0d0,damp/),(/4,4/))
      !find the current crack strain (i.e at the end of the increment)      
      call calcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,
     # GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2,dt,GGICf,N1,alpha)
      
      !calculate the new crack modulus Dcr
      call calcDcr(epscr,epscrmax,GICm,GICf,GIICm,
     #       sigcrf0,sigcrm0,taucrm0,l1,l2,Dcr,GGICf,alpha)
      
      !calculate combined stiffness of crack and continuum
      term11 = MATMUL(MATMUL(TRANSPOSE(N1),Dco),N1)
      term22 = MATMUL(TRANSPOSE(N1),Dco)
      call inverse(Dcr+term11+1.0d0/dt*Dda,4,invDcrDco)  
      term33=matmul(invDcrDco,term22)
      
      !calculate new stress
      Dcocr=Dco-matmul(TRANSPOSE(term22),term33)
      sig=matmul(Dcocr,eps)
     #   -matmul(matmul(matmul(TRANSPOSE(term22),invDcrDco),
     #   1.0d0/dt*Dda),epscrold)

      if (isImplicit.eq.1) ddsdde=Dcocr
      

           
!      !calculate the energy dissipated so far
!      dGIf =dmax1(0.0d0,depscr(1)*l1)*(2*sig(1)-dsig(1))*0.5d0
!      dGIm =dmax1(0.0d0,depscr(2)*l2)*(2*sig(2)-dsig(2))*0.5d0
!      dGIIm=dabs(       depscr(3)*l2 *(2*sig(3)-dsig(3))*0.5d0)
!     #                                          *modeIIadvance
!      statev(7)=statev(7)+dGIf/GICf
!      statev(8)=statev(8)+dGIm/GICm
!      statev(9)=statev(9)+dGIIm/GIICm
!**********Under this line is the Fiber tension failure mode************ 
      else if ((epscrfmax.gt.1.d-8)) then
C       write(*,*) 'Fiber crack!!!!: '
C       write(*,*) epscrfmax,gamcrfmax,epscrmmax,gamcrmmax

      !crack band
      epscrold=(/statev(11),statev(12),statev(13),statev(14)/) !(/epscrf,gamcrf,epscrm,gamcrm/)
      epscrmax=(/epscrfmax,gamcrfmax,epscrmmax,gamcrmmax/)
      !form the damping matrix
      Dda=reshape((/damp,0.0d0,0.0d0,0.0d0,
     # 0.0d0,damp,0.0d0,0.0d0,
     # 0.0d0,0.0d0,damp,0.0d0,
     # 0.0d0,0.0d0,0.0d0,damp/),(/4,4/))
      !find the current crack strain (i.e at the end of the increment)      
      call calcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,
     # GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2,dt,GGICf,N2,alpha)
      
      !calculate the new crack modulus Dcr
      call calcDcr(epscr,epscrmax,GICm,GICf,GIICm,
     #       sigcrf0,sigcrm0,taucrm0,l1,l2,Dcr,GGICf,alpha)
      
      !calculate combined stiffness of crack and continuum
      term11 = MATMUL(MATMUL(TRANSPOSE(N2),Dco),N2)
      term22 = MATMUL(TRANSPOSE(N2),Dco)
      call inverse(Dcr+term11+1.0d0/dt*Dda,4,invDcrDco)  
      term33=matmul(invDcrDco,term22)
      
      !calculate new stress
      Dcocr=Dco-matmul(TRANSPOSE(term22),term33)
      sig=matmul(Dcocr,eps)
     #   -matmul(matmul(matmul(TRANSPOSE(term22),invDcrDco),
     #   1.0d0/dt*Dda),epscrold)

      if (isImplicit.eq.1) ddsdde=Dcocr
      

           
!      !calculate the energy dissipated so far
!      dGIf =dmax1(0.0d0,depscr(1)*l1)*(2*sig(1)-dsig(1))*0.5d0
!      dGIm =dmax1(0.0d0,depscr(2)*l2)*(2*sig(2)-dsig(2))*0.5d0
!      dGIIm=dabs(       depscr(3)*l2 *(2*sig(3)-dsig(3))*0.5d0)
!     #                                          *modeIIadvance
!      statev(7)=statev(7)+dGIf/GICf
!      statev(8)=statev(8)+dGIm/GICm
!      statev(9)=statev(9)+dGIIm/GIICm
!**********Above this line is the Fiber tension failure mode************       

      endif


      !update the maximum crack strain       
      statev(3)=dmax1(epscr(1),epscrfmax)
      statev(4)=dmax1(dabs(epscr(2)),gamcrfmax)
      statev(5)=dmax1(epscr(3),epscrmmax)
      statev(6)=dmax1(dabs(epscr(4)),gamcrmmax) !!?? absolute values??
      
      statev(11)=epscr(1)
      statev(12)=epscr(2)
      statev(13)=epscr(3)
      statev(14)=epscr(4)
      
      if ((statev(7).gt.0.95).or.
     #    (statev(8).gt.0.95).or.
     #    (statev(9).gt.0.95)) then
           statev(10)=0.0d0
      endif
      
      return
      end

!-----------------------------------------------------------------------
! subroutine to find the crack strain at the end of the increment
! using Newton's method (Broyden's method does not work due to too large
! variation in the Jacobian)
 
       subroutine calcEpscr(epscr,epscrmax,epscrold,eps,Dco,Dda,
     # GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2,dt,GGICf,N,alpha)
      
      implicit none
      real*8::epscr(4),epscrmax(4),Dco(3,3),Dda(4,4),eps(3),epscrold(4) 
      real*8::GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2,dt
      real*8::GGICf,N(3,4),alpha
      integer::imax,i,k


      
      ! everything internal is handled as quad precision due to 
      ! very small initial cracks
      real*16::J(4,4),invJ(4,4),Ftol,delta
      real*16::x(4),dx(4,1),Fo(4),Fn(4),xPlusDx(4),x0(4) 
      real*16::qepscrmax(4),qDco(3,3),qDcr(4,4), qeps(3),qepscrold(4)
      real*16::qDda(4,4)
      real*16::qGICm,qGICf,qGIICm,qsigcrf0,qsigcrm0,qtaucrm0,ql1,ql2,q
      real*16::qdti
      real*16::qGGICf
      real*16::qN(3,4)
      real*16::qalpha
      real*16::term11(4,4) !*** Nt*Dco*N
      real*16::term22(4,3) !*** Nt*Dco
            
      qepscrmax  = epscrmax
      qepscrold  = epscrold
      qeps       = eps
      qDco       = Dco
      qDda       = Dda
      qGICm      = GICm
      qGICf      = GICf
      qGIICm     = GIICm
      qsigcrf0   = sigcrf0
      qsigcrm0   = sigcrm0
      qtaucrm0   = taucrm0
      ql1        = l1
      ql2        = l2
      qdti       = 1.0d0/dt
      qGGICf     = GGICf
      qN         = N
      qalpha     = alpha
           
      imax=50  !maximum number of iterations
      Ftol=1.0d-4 !tolerance on the force
      
      !initial guess of crackstrain
      x=0.0d0
      
      !iterate:
      do i=1,imax
      !calc current force
       call qcalcDcr(x,qepscrmax,qGICm,qGICf,qGIICm,
     #       qsigcrf0,qsigcrm0,qtaucrm0,ql1,ql2,qDcr,qGGICf,qalpha)
      term11 = MATMUL(MATMUL(TRANSPOSE(qN),qDco),qN)
      term22 = MATMUL(TRANSPOSE(qN),qDco)
      Fo=matmul(qDcr+term11+qdti*qDda,x)
     #  -matmul(term22,qeps)-matmul(qdti*qDda,qepscrold)   
      if (qabs(Fo(1))+qabs(Fo(2))+qabs(Fo(3))+qabs(Fo(4)).lt.Ftol) exit 
      !calc jacobian
          do k=1,4
          delta=qsign(1.0q-16,x(k))
          xPlusDx=x
          xPlusDx(k)=x(k)+delta
       call qcalcDcr(xPlusDx,qepscrmax,qGICm,qGICf,qGIICm,
     #       qsigcrf0,qsigcrm0,qtaucrm0,ql1,ql2,qDcr,qGGICf,qalpha)
            Fn=matmul(qDcr+term11+qdti*qDda,xPlusDx)
     #        -matmul(term22,qeps)-matmul(qdti*qDda,qepscrold)
          J(:,k)=(Fn-Fo)/delta
          end do
      !dx=-matmul(J^-1,F)
      call qinverse(J,4,invJ) !q~

      dx(:,1)=-matmul(invJ,Fo)
      !x=x+dx
      x=x+dx(:,1)
      !check for convergence
      end do

      if (i.ge.imax) then
      !A little bit of diagnostics
      write(*,*)
      write(*,*) 'Newtons`s method failed to converge in max number'
      write(*,*) 'of increments'
      write(*,*) 'Dco='
      write(*,*) Dco
      write(*,*) 'Dda='
      write(*,*) Dda
      write(*,*) 'eps='
      write(*,*) eps
      write(*,*) 'epscrmax='
      write(*,*) epscrmax
      write(*,*) 'epscrold'
      write(*,*) epscrold
      write(*,*) 'GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2,dt'
      write(*,*) GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2,dt
      write(*,*) 'The final values of the iteration'
      write(*,*) 'epscr='
      write(*,*) dble(x)
      write(*,*) 'F'
      write(*,*) dble(Fo)
 
      
      call myExit()
      endif 
      
      epscr=dble(x)

      
      return
      end

!----------------------------------------------------------------------------


      
! ----------------------------------------------------------------------------      
! function to find the crack slope
! using exponentials 
      subroutine calcDcr(epscr,epscrmax,GICm,GICf,GIICm,
     #       sigcrf0,sigcrm0,taucrm0,l1,l2,Dcr,GGICf,alpha)
      implicit none
      
      real*8::GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2
      real*8::epscr(4),epscrmax(4),Dcr(4,4)
      real*8::GGICf,alpha
      real*8::epscr_temp(4)
      real*8::sigcrf0_r,sigcrm0_r,taucrm0_r

      Dcr=0.0d0
      epscr_temp(1)=max(dabs(epscr(1)),dabs(epscrmax(1)))
      epscr_temp(2)=max(dabs(epscr(2)),dabs(epscrmax(2)))
      epscr_temp(3)=max(dabs(epscr(3)),dabs(epscrmax(3)))
      epscr_temp(4)=max(dabs(epscr(4)),dabs(epscrmax(4)))

      sigcrf0_r = sigcrf0*1d-10
      sigcrm0_r = sigcrm0*1d-10
      taucrm0_r = taucrm0*1d-10

      !fiber mode I
      IF (epscr_temp(1).LT.(2.0d0*GICf/l1/(sigcrf0-sigcrf0_r))) THEN
       Dcr(1,1) = (-(sigcrf0_r-sigcrf0)**2/(2.0d0*GICf/l1)*
     1 epscr_temp(1)+sigcrf0)/epscr_temp(1)
      ELSE
       Dcr(1,1) = sigcrf0_r/epscr_temp(1)
      END IF 
      
      !Fiber mode II. However, it is not real mode II.
      if (dabs(epscr(1)).lt.1.d-16) then
        Dcr(2,2)=1.d10*GGICf
      else
        Dcr(2,2)=alpha*GGICf
      end if

      !matrix mode I
      IF (epscr_temp(3).LT.(2.0d0*GICm/l2/(sigcrm0-sigcrm0_r))) THEN
       Dcr(3,3) = (-(sigcrm0_r-sigcrm0)**2/(2.0d0*GICm/l2)*
     1 epscr_temp(3)+sigcrm0)/epscr_temp(3)
      ELSE
       Dcr(3,3) = sigcrm0_r/epscr_temp(3)
      END IF 
      
      !matrix mode II
      !the slope vs. crack strain is symmetric with respect to the y-axis
      !therefore only absolute values are of interest
      IF (epscr_temp(4).LT.(2.0d0*GIICm/l2/(taucrm0-taucrm0_r))) THEN
       Dcr(4,4) = (-(taucrm0_r-taucrm0)**2/(2.0d0*GIICm/l2)*
     1 epscr_temp(4)+taucrm0)/epscr_temp(4)
      ELSE
       Dcr(4,4) = taucrm0_r/epscr_temp(4)
      END IF 
      
!*************This is for DCB********************
      !Matrix mode II. However, it is not real mode II.
C       if (dabs(epscrmax(3)).GE.(2.0d0*GICm/l2/(sigcrm0-sigcrm0_r))) then
      if (dabs(epscrmax(3)).GE.(1.d-8)) then
        Dcr(4,4)=alpha*Dcr(3,3)/(2.0d0*(1.0d0+0.25))
        Dcr(1,1)=alpha*alpha*Dcr(3,3)
      end if
      return
      end
      
      !---------------------------
      !same subroutine in quads
      subroutine qcalcDcr(epscr,epscrmax,GICm,GICf,GIICm,
     #       sigcrf0,sigcrm0,taucrm0,l1,l2,Dcr,GGICf,alpha)
      implicit none
      
      real*16::GICm,GICf,GIICm,sigcrf0,sigcrm0,taucrm0,l1,l2
      real*16::epscr(4),epscrmax(4),Dcr(4,4)
      real*16::GGICf,alpha
      real*16::epscr_temp(4)
      real*16::sigcrf0_r,sigcrm0_r,taucrm0_r

      Dcr=0.0d0
      epscr_temp(1)=max(qabs(epscr(1)),qabs(epscrmax(1)))
      epscr_temp(2)=max(qabs(epscr(2)),qabs(epscrmax(2)))
      epscr_temp(3)=max(qabs(epscr(3)),qabs(epscrmax(3)))
      epscr_temp(4)=max(qabs(epscr(4)),qabs(epscrmax(4)))

      sigcrf0_r = sigcrf0*1d-10
      sigcrm0_r = sigcrm0*1d-10
      taucrm0_r = taucrm0*1d-10

      !fiber mode I
      IF (epscr_temp(1).LT.(2.0Q0*GICf/l1/(sigcrf0-sigcrf0_r))) THEN
       Dcr(1,1) = (-(sigcrf0_r-sigcrf0)**2/(2.0Q0*GICf/l1)*
     1 epscr_temp(1)+sigcrf0)/epscr_temp(1)
      ELSE
       Dcr(1,1) = sigcrf0_r/epscr_temp(1)
      END IF 

      !Fiber mode II. However, it is not real mode II.
      if (qabs(epscr(1)).le.1.d-8) then
        Dcr(2,2)=1.d10*GGICf
      else
        Dcr(2,2)=alpha*GGICf
      end if
      
      !matrix mode I
      IF (epscr_temp(3).LT.(2.0d0*GICm/l2/(sigcrm0-sigcrm0_r))) THEN
       Dcr(3,3) = (-(sigcrm0_r-sigcrm0)**2/(2.0d0*GICm/l2)*
     1 epscr_temp(3)+sigcrm0)/epscr_temp(3)
      ELSE
       Dcr(3,3) = sigcrm0_r/epscr_temp(3)
      END IF 
      
      !matrix mode II
      !the slope vs. crack strain is symmetric with respect to the y-axis
      !therefore only absolute values are of interest
      IF (epscr_temp(4).LT.(2.0d0*GIICm/l2/(taucrm0-taucrm0_r))) THEN
       Dcr(4,4) = (-(taucrm0_r-taucrm0)**2/(2.0d0*GIICm/l2)*
     1 epscr_temp(4)+taucrm0)/epscr_temp(4)
      ELSE
       Dcr(4,4) = taucrm0_r/epscr_temp(4)
      END IF 

!*************This is for DCB********************
      !Matrix mode II. However, it is not real mode II.
C       if (qabs(epscrmax(3)).GE.(2.0d0*GICm/l2/(sigcrm0-sigcrm0_r))) then
      if (qabs(epscrmax(3)).GE.(1.d-8)) then
        Dcr(4,4)=alpha*Dcr(3,3)/(2.0d0*(1.0d0+0.25))
        Dcr(1,1)=alpha*alpha*Dcr(3,3)
      end if   
      return
      end

      
! invert 3x3 matrix directly
! call inverse(A,invA)
      include 'myaux.for'


