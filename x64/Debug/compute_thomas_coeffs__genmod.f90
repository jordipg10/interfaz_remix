        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:38 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_THOMAS_COEFFS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_THOMAS_COEFFS(A,B,TOL,C_TILDE,D_TILDE)
              USE METODOS_SIST_LIN_M
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: C_TILDE(:)
              REAL(KIND=8), INTENT(OUT) :: D_TILDE(:)
            END SUBROUTINE COMPUTE_THOMAS_COEFFS
          END INTERFACE 
        END MODULE COMPUTE_THOMAS_COEFFS__genmod
