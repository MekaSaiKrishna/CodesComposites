C ! invert a 3x3 matrix directly

C       subroutine inverseMatrix3x3(J,invJ)
      
C       implicit none
      
C       real*8,dimension(3,3)::J,invJ
C       real*8::detJ
      
C       detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
C      1    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
C      2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      
C       invJ(1,:)=(/ (J(2,2)*J(3,3) - J(2,3)*J(3,2)), 
C      1            -(J(1,2)*J(3,3) - J(1,3)*J(3,2)),  
C      2             (J(1,2)*J(2,3) - J(1,3)*J(2,2)) /)
C       invJ(2,:)=(/-(J(2,1)*J(3,3) - J(2,3)*J(3,1)),  
C      1             (J(1,1)*J(3,3) - J(1,3)*J(3,1)), 
C      2            -(J(1,1)*J(2,3) - J(1,3)*J(2,1)) /)
C       invJ(3,:)=(/ (J(2,1)*J(3,2) - J(2,2)*J(3,1)), 
C      1            -(J(1,1)*J(3,2) - J(1,2)*J(3,1)),  
C      2             (J(1,1)*J(2,2) - J(1,2)*J(2,1)) /)
      
C       invJ=invJ/detJ
      
C       return
C       end
      
C ! invert a 3x3 matrix directly, quadruple precision
C       subroutine qinverseMatrix3x3(J,invJ)
      
C       implicit none
      
C       real*16,dimension(3,3)::J,invJ
C       real*16::detJ
      
C       detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
C      1    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
C      2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      
C       invJ(1,:)=(/ (J(2,2)*J(3,3) - J(2,3)*J(3,2)), 
C      1            -(J(1,2)*J(3,3) - J(1,3)*J(3,2)),  
C      2             (J(1,2)*J(2,3) - J(1,3)*J(2,2)) /)
C       invJ(2,:)=(/-(J(2,1)*J(3,3) - J(2,3)*J(3,1)),  
C      1             (J(1,1)*J(3,3) - J(1,3)*J(3,1)), 
C      2            -(J(1,1)*J(2,3) - J(1,3)*J(2,1)) /)
C       invJ(3,:)=(/ (J(2,1)*J(3,2) - J(2,2)*J(3,1)), 
C      1            -(J(1,1)*J(3,2) - J(1,2)*J(3,1)),  
C      2             (J(1,1)*J(2,2) - J(1,2)*J(2,1)) /)
      
C       invJ=invJ/detJ
      
C       return
C       end
 
C ! invert a 4x4 matrix directly
C       subroutine inverseMatrix4x4(J,invJ)
      
C       implicit none
      
C       real*8,dimension(4,4)::J,invJ
C       real*8::detJ
      
C       invJ=0.0d0
C       detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
C      1    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
C      2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      
C       invJ(1,1:3)=(/ (J(2,2)*J(3,3) - J(2,3)*J(3,2)), 
C      1            -(J(1,2)*J(3,3) - J(1,3)*J(3,2)),  
C      2             (J(1,2)*J(2,3) - J(1,3)*J(2,2)) /)
C       invJ(2,1:3)=(/-(J(2,1)*J(3,3) - J(2,3)*J(3,1)),  
C      1             (J(1,1)*J(3,3) - J(1,3)*J(3,1)), 
C      2            -(J(1,1)*J(2,3) - J(1,3)*J(2,1)) /)
C       invJ(3,1:3)=(/ (J(2,1)*J(3,2) - J(2,2)*J(3,1)), 
C      1            -(J(1,1)*J(3,2) - J(1,2)*J(3,1)),  
C      2             (J(1,1)*J(2,2) - J(1,2)*J(2,1)) /)

C       invJ=invJ/detJ
C       invJ(4,4)=1.0d0/J(4,4)
      
C       return
C       end


C ! invert a 3x3 matrix directly, quadruple precision
C       subroutine qinverseMatrix4x4(J,invJ)
      
C       implicit none
      
C       real*16,dimension(4,4)::J,invJ
C       real*16::detJ
      
C       invJ=0.0d0
C       detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))
C      1    -J(1,2)*(J(2,1)*J(3,3)-J(2,3)*J(3,1))
C      2    +J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))

C       invJ(1,1:3)=(/ (J(2,2)*J(3,3) - J(2,3)*J(3,2)), 
C      1            -(J(1,2)*J(3,3) - J(1,3)*J(3,2)),  
C      2             (J(1,2)*J(2,3) - J(1,3)*J(2,2)) /)
C       invJ(2,1:3)=(/-(J(2,1)*J(3,3) - J(2,3)*J(3,1)),  
C      1             (J(1,1)*J(3,3) - J(1,3)*J(3,1)), 
C      2            -(J(1,1)*J(2,3) - J(1,3)*J(2,1)) /)
C       invJ(3,1:3)=(/ (J(2,1)*J(3,2) - J(2,2)*J(3,1)), 
C      1            -(J(1,1)*J(3,2) - J(1,2)*J(3,1)),  
C      2             (J(1,1)*J(2,2) - J(1,2)*J(2,1)) /)

C       invJ=invJ/detJ
C       invJ(4,4)=1.0d0/J(4,4)
      
C       return
C       end


      subroutine inverse(ori_m,N,inv_m)
      implicit none
      integer indx(N),I,J,K,N
      real*8 ori_m(N,N),inv_m(N,N),B(N,N),
     1 ori_mt(N,N),inv_mt(N,N),ztm,critical
C
c      tell whether it is zero matrix
      critical = 1e-20
      ztm = 0.
      do 101 i = 1,3
        do 102 j = 1,3
          ztm = ztm + abs(ori_m(i,j))
102     continue
101   continue
      
      if (ztm < critical) then
        print *,"Zero Matrix"
        stop
      endif 
c
      do  103 I = 1, N
        do 104 J = 1, N
          B(I,J) = 0.0
          if(I==J)then
            B(I,J) = 1.
          endif
104     continue
103   continue
C
      do 105 I = 1,N
        do 106 j = 1,N
        ori_mt(I,J) = ori_m(i,j)
106     continue 
105   continue

      do 107 I = 1,N
        do 108 J = 1,N
          ori_m(I,J) = ori_mt(J,I)
108     continue
107   continue

      CALL PIGS(ori_m,N,indx)
C
      do 109 I = 1, N-1
        do 110 J = I+1, N
          do 111 K = 1, N
            B(indx(J),K) = B(indx(J),K)
     *                    -ori_m(indx(J),I)*B(indx(I),K)
111       continue
110     continue
109   continue
C
      do 112 I = 1, N
        inv_m(N,I) = B(indx(N),I)/ori_m(indx(N),N)
        do 113 J = N-1, 1, -1
          inv_m(J,I) = B(indx(J),I)
          do 114 K = J+1, N
            inv_m(J,I) = inv_m(J,I)-ori_m(indx(J),K)*inv_m(K,I)
114       continue
          inv_m(J,I) =  inv_m(J,I)/ori_m(indx(J),J)
113     continue
112   continue

      do 115 I = 1,N
        do 116 j = 1,N
        inv_mt(I,J) = inv_m(i,j)
116     continue 
115   continue

      do 117 I = 1,N
        do 118 J = 1,N
          inv_m(I,J) = inv_mt(J,I)
118     continue
117   continue
C
      RETURN
      END
c----------------------------------------------------
      subroutine PIGS(ori_m,N,indx)
      implicit none      
      integer indx(N),I,J,K,N
      real*8 ori_m(N,N),C(N),critical,C1,PI1,PI,ITMP,PJ
C
      critical = 1e-20
      do  119 I = 1, N
        indx(I) = I
 119   continue
C
        do 120 I = 1, N
          C1= 0.0
          do 121 J = 1, N
            C1 = MAX(C1,ABS(ori_m(I,J)))
 121      continue
          C(I) = C1
 120     continue
C
      do 122 J = 1, N-1
        PI1 = 0.0
        do 123 I = J, N
          if (abs(C(indx(I)))<critical) then
            print *,"Singular Matrix1"
            stop
          endif
          PI = ABS(ori_m(indx(I),J))/C(indx(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
 123    continue
C
        ITMP    = indx(J)
        indx(J) = indx(K)
        indx(K) = ITMP
        do 124 I = J+1, N
          if (abs(ori_m(indx(J),J))<critical) then
            print *,"Singular Matrix2"
            stop
          endif
          PJ  = ori_m(indx(I),J)/ori_m(indx(J),J)
C
          ori_m(indx(I),J) = PJ
C
          do 125 K = J+1, N
            ori_m(indx(I),K) = ori_m(indx(I),K)-PJ*ori_m(indx(J),K)
 125      continue
 124   continue
 122  continue
C
      RETURN
      END
C-----------------------------------------------------------------------

      subroutine qinverse(ori_m,N,inv_m)
      implicit none
      integer indx(N),I,J,K,N
      real*16 ori_m(N,N),inv_m(N,N),B(N,N),
     1 ori_mt(N,N),inv_mt(N,N),ztm,critical
C
c      tell whether it is zero matrix
      critical = 1Q-2000
      ztm = 0.Q0
      do 101 i = 1,3
        do 102 j = 1,3
          ztm = ztm + abs(ori_m(i,j))
102     continue
101   continue
      
      if (ztm < critical) then
        print *,"Zero Matrix"
        stop
      endif 
c
      do  103 I = 1, N
        do 104 J = 1, N
          B(I,J) = 0.0
          if(I==J)then
            B(I,J) = 1.
          endif
104     continue
103   continue
C
      do 105 I = 1,N
        do 106 j = 1,N
        ori_mt(I,J) = ori_m(i,j)
106     continue 
105   continue

      do 107 I = 1,N
        do 108 J = 1,N
          ori_m(I,J) = ori_mt(J,I)
108     continue
107   continue

      CALL qPIGS(ori_m,N,indx)
C
      do 109 I = 1, N-1
        do 110 J = I+1, N
          do 111 K = 1, N
            B(indx(J),K) = B(indx(J),K)
     *                    -ori_m(indx(J),I)*B(indx(I),K)
111       continue
110     continue
109   continue
C
      do 112 I = 1, N
        inv_m(N,I) = B(indx(N),I)/ori_m(indx(N),N)
        do 113 J = N-1, 1, -1
          inv_m(J,I) = B(indx(J),I)
          do 114 K = J+1, N
            inv_m(J,I) = inv_m(J,I)-ori_m(indx(J),K)*inv_m(K,I)
114       continue
          inv_m(J,I) =  inv_m(J,I)/ori_m(indx(J),J)
113     continue
112   continue

      do 115 I = 1,N
        do 116 j = 1,N
        inv_mt(I,J) = inv_m(i,j)
116     continue 
115   continue

      do 117 I = 1,N
        do 118 J = 1,N
          inv_m(I,J) = inv_mt(J,I)
118     continue
117   continue
C
      RETURN
      END
c----------------------------------------------------
      subroutine qPIGS(ori_m,N,indx)
      implicit none      
      integer indx(N),I,J,K,N
      real*16 ori_m(N,N),C(N),critical,C1,PI1,PI,ITMP,PJ
C
      critical = 1Q-2000
      do  119 I = 1, N
        indx(I) = I
 119   continue
C
        do 120 I = 1, N
          C1= 0.0
          do 121 J = 1, N
            C1 = MAX(C1,ABS(ori_m(I,J)))
 121      continue
          C(I) = C1
 120     continue
C
      do 122 J = 1, N-1
        PI1 = 0.0
        do 123 I = J, N
          if (abs(C(indx(I)))<critical) then
            print *,"Singular Matrix1"
            stop
          endif
          PI = ABS(ori_m(indx(I),J))/C(indx(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ELSE
          ENDIF
 123    continue
C
        ITMP    = indx(J)
        indx(J) = indx(K)
        indx(K) = ITMP
        do 124 I = J+1, N
          if (abs(ori_m(indx(J),J))<critical) then
            print *,"Singular Matrix2"
            stop
          endif
          PJ  = ori_m(indx(I),J)/ori_m(indx(J),J)
C
          ori_m(indx(I),J) = PJ
C
          do 125 K = J+1, N
            ori_m(indx(I),K) = ori_m(indx(I),K)-PJ*ori_m(indx(J),K)
 125      continue
 124   continue
 122  continue
C
      RETURN
      END