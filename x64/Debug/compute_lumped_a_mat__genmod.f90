        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 10:33:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LUMPED_A_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LUMPED_A_MAT(THIS,A_MAT_LUMPED)
              USE TIME_DISCR_M
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          DIAG_MATRIX_C
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              TYPE (DIAG_MATRIX_C), INTENT(OUT) :: A_MAT_LUMPED
            END SUBROUTINE COMPUTE_LUMPED_A_MAT
          END INTERFACE 
        END MODULE COMPUTE_LUMPED_A_MAT__genmod
