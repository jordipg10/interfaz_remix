        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_INVERSE_TRIDIAG_MATRIX__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_INVERSE_TRIDIAG_MATRIX(THIS,TOL,INV_MAT)
              USE METODOS_SIST_LIN_M, ONLY :                            &
     &          THOMAS,                                                 &
     &          TRIDIAG_MATRIX_C
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: INV_MAT(:,:)
            END SUBROUTINE COMPUTE_INVERSE_TRIDIAG_MATRIX
          END INTERFACE 
        END MODULE COMPUTE_INVERSE_TRIDIAG_MATRIX__genmod
