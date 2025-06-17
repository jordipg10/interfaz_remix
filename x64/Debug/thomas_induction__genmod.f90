        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:42 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE THOMAS_INDUCTION__genmod
          INTERFACE 
            SUBROUTINE THOMAS_INDUCTION(A,B,I,X_I,TOL,X_K)
              USE MATRICES_M
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B(:)
              INTEGER(KIND=4), INTENT(IN) :: I
              REAL(KIND=8), INTENT(IN) :: X_I
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: X_K(:)
            END SUBROUTINE THOMAS_INDUCTION
          END INTERFACE 
        END MODULE THOMAS_INDUCTION__genmod
