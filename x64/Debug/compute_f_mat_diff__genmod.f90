        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:39 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_F_MAT_DIFF__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_F_MAT_DIFF(THIS)
              USE CHAR_PARAMS_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          SPATIAL_DISCR_RAD_C,                                    &
     &          DIAG_MATRIX_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
            END SUBROUTINE COMPUTE_F_MAT_DIFF
          END INTERFACE 
        END MODULE COMPUTE_F_MAT_DIFF__genmod
