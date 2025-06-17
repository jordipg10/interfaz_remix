        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:50 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LUMPED_A_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LUMPED_A_MAT(THIS,A_MAT_LUMPED)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          DIAG_MATRIX_C
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              TYPE (DIAG_MATRIX_C), INTENT(OUT) :: A_MAT_LUMPED
            END SUBROUTINE COMPUTE_LUMPED_A_MAT
          END INTERFACE 
        END MODULE COMPUTE_LUMPED_A_MAT__genmod
