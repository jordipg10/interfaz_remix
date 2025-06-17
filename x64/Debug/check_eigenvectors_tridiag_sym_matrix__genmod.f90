        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:59 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX__genmod
          INTERFACE 
            SUBROUTINE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX(THIS,      &
     &TOLERANCE)
              USE MATRICES_M, ONLY :                                    &
     &          TRIDIAG_SYM_MATRIX_C
              CLASS (TRIDIAG_SYM_MATRIX_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
            END SUBROUTINE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX
          END INTERFACE 
        END MODULE CHECK_EIGENVECTORS_TRIDIAG_SYM_MATRIX__genmod
