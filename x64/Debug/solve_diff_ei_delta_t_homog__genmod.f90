        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 15:24:59 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_DIFF_EI_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_DIFF_EI_DELTA_T_HOMOG(THIS,THETA,TIME_OUT, &
     &OUTPUT)
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
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HOMOG_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_DIFF_EI_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE SOLVE_DIFF_EI_DELTA_T_HOMOG__genmod
