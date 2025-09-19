        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 10:32:05 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_TRANS_MAT_DIFF__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_TRANS_MAT_DIFF(THIS)
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE VECTORS_M
              USE MATRICES_M
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C
              CLASS (DIFFUSION_1D_C) :: THIS
            END SUBROUTINE COMPUTE_TRANS_MAT_DIFF
          END INTERFACE 
        END MODULE COMPUTE_TRANS_MAT_DIFF__genmod
