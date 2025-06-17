        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_A_MAT_ODE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_A_MAT_ODE(THIS,A_MAT)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          TRIDIAG_MATRIX_C,                                       &
     &          DIAG_MATRIX_C,                                          &
     &          PROD_TRIDIAG_DIAG_MAT,                                  &
     &          PROD_DIAG_TRIDIAG_MAT
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              TYPE (TRIDIAG_MATRIX_C), INTENT(OUT) :: A_MAT
            END SUBROUTINE COMPUTE_A_MAT_ODE
          END INTERFACE 
        END MODULE COMPUTE_A_MAT_ODE__genmod
