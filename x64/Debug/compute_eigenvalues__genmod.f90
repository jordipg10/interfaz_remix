        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:28 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_EIGENVALUES__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_EIGENVALUES(THIS)
              USE MATRICES_M, ONLY :                                    &
     &          TRIDIAG_SYM_MATRIX_C,                                   &
     &          ARRAY_C
              CLASS (ARRAY_C) :: THIS
            END SUBROUTINE COMPUTE_EIGENVALUES
          END INTERFACE 
        END MODULE COMPUTE_EIGENVALUES__genmod
