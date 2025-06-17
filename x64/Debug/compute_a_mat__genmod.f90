        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:44 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_A_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_A_MAT(THIS,THETA,E_MAT)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          TRIDIAG_MATRIX_C
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: E_MAT
            END SUBROUTINE COMPUTE_A_MAT
          END INTERFACE 
        END MODULE COMPUTE_A_MAT__genmod
