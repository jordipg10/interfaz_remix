        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 10:32:29 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_DIFFUSION_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_DIFFUSION_TRANSIENT(THIS)
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M, ONLY :                            &
     &          STAB_PARAMS_DIFF_C,                                     &
     &          STAB_PARAMS_C
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M, ONLY :                          &
     &          DIFF_PROPS_HETEROG_C,                                   &
     &          PROPS_C
              USE CONC_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_HETEROG_C,                                   &
     &          TIME_DISCR_C
              USE TIME_FCT_M
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
            END SUBROUTINE INITIALISE_DIFFUSION_TRANSIENT
          END INTERFACE 
        END MODULE INITIALISE_DIFFUSION_TRANSIENT__genmod
