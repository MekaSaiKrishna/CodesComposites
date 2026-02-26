!===============================================================================
! Module: matrix_utils
! ---------------------------------------------------------------------------
! Linear algebra utilities for the SCA2D formulation.
!
! Smart Change #1:  Analytical 2x2 matrix inverse (no LU overhead)
! Smart Change #4:  Quad-precision duplicates removed (only double precision)
! Smart Change #10: Allocatable arrays in LU decomposition (no fixed NMAX)
!
! Based on: Numerical Recipes in Fortran (Crout's method with partial pivoting)
! Reference: http://www.mpi-hd.mpg.de/astrophysik/HEA/internal/Numerical_Recipes/f2-3.pdf
!
! Last Updated: 02/26/2026
! meka1@purdue.edu
!===============================================================================
module matrix_utils
    implicit none
    private
    public :: inverse_2x2, matrix_inverse_general, solve_lin_sys_lu

contains

    !---------------------------------------------------------------------------
    ! Smart Change #1: Analytical 2x2 matrix inverse
    !
    ! For a 2x2 matrix A = [a b; c d], the inverse is:
    !   A^(-1) = (1/det) * [d, -b; -c, a]   where det = ad - bc
    !
    ! This replaces the general LU-based matrixInverse for the 2x2 case,
    ! eliminating decomposition, pivoting, and back-substitution overhead.
    !---------------------------------------------------------------------------
    subroutine inverse_2x2(A, invA, flag)
        real(8), intent(in)  :: A(2,2)
        real(8), intent(out) :: invA(2,2)
        integer, intent(out) :: flag
        real(8) :: det

        det = A(1,1)*A(2,2) - A(1,2)*A(2,1)

        if (abs(det) < 1.0d-30) then
            flag = 0
            invA = 0.0d0
            return
        end if

        flag = 1
        invA(1,1) =  A(2,2) / det
        invA(1,2) = -A(1,2) / det
        invA(2,1) = -A(2,1) / det
        invA(2,2) =  A(1,1) / det
    end subroutine inverse_2x2

    !---------------------------------------------------------------------------
    ! General matrix inverse using LU decomposition
    ! Smart Change #10: No fixed NMAX limit; uses allocatable arrays
    !
    ! Input:  a(np,np) - matrix to invert (DESTROYED on output)
    !         n        - logical dimension
    !         np       - physical dimension
    ! Output: y(np,np) - inverse of original a
    !---------------------------------------------------------------------------
    subroutine matrix_inverse_general(a, y, n, np)
        integer, intent(in)    :: n, np
        real(8), intent(inout) :: a(np,np)
        real(8), intent(out)   :: y(np,np)

        integer, allocatable :: indx(:)
        integer :: flag, i, j
        real(8) :: d

        allocate(indx(n))

        ! Set up identity matrix
        y = 0.0d0
        do i = 1, n
            y(i,i) = 1.0d0
        end do

        ! Decompose once
        call ludcmp_v2(a, n, np, indx, d, flag)

        ! Solve for each column of the inverse
        if (flag == 1) then
            do j = 1, n
                call lubksb_v2(a, n, np, indx, y(1,j))
            end do
        end if

        deallocate(indx)
    end subroutine matrix_inverse_general

    !---------------------------------------------------------------------------
    ! Solve Ax=b via LU decomposition (double precision only)
    ! Smart Change #4: Single version only; quad-precision duplicate removed
    !
    ! Input:  A(n,n) - coefficient matrix (DESTROYED on output)
    !         b(n)   - right-hand side vector
    !         n      - system size
    ! Output: b(n)   - solution vector (overwrites RHS)
    !         flag   - 1=success, 0=singular matrix
    !---------------------------------------------------------------------------
    subroutine solve_lin_sys_lu(A, b, n, flag)
        integer, intent(in)    :: n
        real(8), intent(inout) :: A(n,n), b(n)
        integer, intent(out)   :: flag

        integer, allocatable :: indx(:)
        real(8) :: d

        allocate(indx(n))

        call ludcmp_v2(A, n, n, indx, d, flag)
        if (flag == 1) then
            call lubksb_v2(A, n, n, indx, b)
        end if

        deallocate(indx)
    end subroutine solve_lin_sys_lu

    !---------------------------------------------------------------------------
    ! LU Decomposition using Crout's method with partial pivoting
    ! Smart Change #10: Uses allocatable vv(:) instead of fixed vv(NMAX)
    !
    ! Replaces a(1:n,1:n) in-place with its LU decomposition.
    ! indx(1:n) records the row permutation from partial pivoting.
    ! d = +/-1 records parity of row swaps.
    ! flag = 1 on success, 0 if matrix is singular.
    !---------------------------------------------------------------------------
    subroutine ludcmp_v2(a, n, np, indx, d, flag)
        integer, intent(in)    :: n, np
        real(8), intent(inout) :: a(np,np)
        integer, intent(out)   :: indx(n)
        real(8), intent(out)   :: d
        integer, intent(out)   :: flag

        real(8), parameter :: TINY = 1.0d-20
        real(8), allocatable :: vv(:)
        integer :: i, imax, j, k
        real(8) :: aamax, dum, sumval

        allocate(vv(n))

        flag = 1
        d = 1.0d0
        imax = 1

        ! Get implicit scaling of each row
        do i = 1, n
            aamax = 0.0d0
            do j = 1, n
                if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))
            end do
            if (aamax == 0.0d0) then
                flag = 0
                deallocate(vv)
                return
            end if
            vv(i) = 1.0d0 / aamax
        end do

        ! Crout's method with partial pivoting
        do j = 1, n
            ! Upper triangular part
            do i = 1, j-1
                sumval = a(i,j)
                do k = 1, i-1
                    sumval = sumval - a(i,k)*a(k,j)
                end do
                a(i,j) = sumval
            end do

            ! Lower triangular part + pivot search
            aamax = 0.0d0
            do i = j, n
                sumval = a(i,j)
                do k = 1, j-1
                    sumval = sumval - a(i,k)*a(k,j)
                end do
                a(i,j) = sumval
                dum = vv(i)*abs(sumval)
                if (dum >= aamax) then
                    imax = i
                    aamax = dum
                end if
            end do

            ! Row interchange if needed
            if (j /= imax) then
                do k = 1, n
                    dum = a(imax,k)
                    a(imax,k) = a(j,k)
                    a(j,k) = dum
                end do
                d = -d
                vv(imax) = vv(j)
            end if

            indx(j) = imax
            if (a(j,j) == 0.0d0) a(j,j) = TINY

            ! Divide by pivot
            if (j /= n) then
                dum = 1.0d0 / a(j,j)
                do i = j+1, n
                    a(i,j) = a(i,j)*dum
                end do
            end if
        end do

        deallocate(vv)
    end subroutine ludcmp_v2

    !---------------------------------------------------------------------------
    ! Forward/Back Substitution for an LU-decomposed system
    !
    ! Given the LU decomposition from ludcmp_v2, solves L*U*x = P*b.
    ! Solution overwrites b. Optimized to skip leading zeros in b.
    !---------------------------------------------------------------------------
    subroutine lubksb_v2(a, n, np, indx, b)
        integer, intent(in)    :: n, np
        real(8), intent(in)    :: a(np,np)
        integer, intent(in)    :: indx(n)
        real(8), intent(inout) :: b(n)

        integer :: i, ii, j, ll
        real(8) :: sumval

        ii = 0

        ! Forward substitution (unscramble permutation)
        do i = 1, n
            ll = indx(i)
            sumval = b(ll)
            b(ll) = b(i)
            if (ii /= 0) then
                do j = ii, i-1
                    sumval = sumval - a(i,j)*b(j)
                end do
            else if (sumval /= 0.0d0) then
                ii = i
            end if
            b(i) = sumval
        end do

        ! Back substitution
        do i = n, 1, -1
            sumval = b(i)
            do j = i+1, n
                sumval = sumval - a(i,j)*b(j)
            end do
            b(i) = sumval / a(i,i)
        end do
    end subroutine lubksb_v2

end module matrix_utils
