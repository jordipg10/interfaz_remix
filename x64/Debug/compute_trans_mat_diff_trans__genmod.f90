        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 10:33:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_TRANS_MAT_DIFF_TRANS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_TRANS_MAT_DIFF_TRANS(THIS)
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE CONC_M
              USE TIME_DISCR_M
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE DIFFUSION_TRANSIENT_M
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
            END SUBROUTINE COMPUTE_TRANS_MAT_DIFF_TRANS
          END INTERFACE 
        END MODULE COMPUTE_TRANS_MAT_DIFF_TRANS__genmod
