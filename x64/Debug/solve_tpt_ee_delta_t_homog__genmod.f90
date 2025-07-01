        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 15:00:34 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_TPT_EE_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_TPT_EE_DELTA_T_HOMOG(THIS,TIME_OUT,OUTPUT)
              USE TRANSPORT_STAB_PARAMS_M
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE DIFFUSION_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          TRIDIAG_MATRIX_C,                                       &
     &          MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_RECHARGE,        &
     &          MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_DISCHARGE,       &
     &          MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_EVAP,            &
     &          MASS_BALANCE_ERROR_ADE_TRANS_PMF_RECHARGE,              &
     &          MASS_BALANCE_ERROR_ADE_TRANS_PMF_DISCHARGE,             &
     &          MASS_BALANCE_ERROR_ADE_TRANS_PMF_EVAP
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_TPT_EE_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE SOLVE_TPT_EE_DELTA_T_HOMOG__genmod
