        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:11 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_TRANS_MAT_DIFF__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_TRANS_MAT_DIFF(THIS)
              USE BCS_M
              USE SPATIAL_DISCR_RAD_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C
              CLASS (DIFFUSION_1D_C) :: THIS
            END SUBROUTINE COMPUTE_TRANS_MAT_DIFF
          END INTERFACE 
        END MODULE COMPUTE_TRANS_MAT_DIFF__genmod
