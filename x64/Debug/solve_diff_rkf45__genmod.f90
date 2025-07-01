        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 15:00:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_DIFF_RKF45__genmod
          INTERFACE 
            SUBROUTINE SOLVE_DIFF_RKF45(THIS,DELTA_T_INIT,TOLERANCE)
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE DIFFUSION_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE VECTORS_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HOMOG_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: DELTA_T_INIT
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
            END SUBROUTINE SOLVE_DIFF_RKF45
          END INTERFACE 
        END MODULE SOLVE_DIFF_RKF45__genmod
