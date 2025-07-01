        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:39 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_EVAP__genmod
          INTERFACE 
            FUNCTION MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_EVAP(THIS,  &
     &CONC_OLD,CONC_NEW,DELTA_T,DELTA_X) RESULT(MASS_BAL_ERR)
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
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:)
              REAL(KIND=8), INTENT(IN) :: CONC_NEW(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: DELTA_X
              REAL(KIND=8) :: MASS_BAL_ERR
            END FUNCTION MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_EVAP
          END INTERFACE 
        END MODULE MASS_BALANCE_ERROR_ADE_TRANS_DIRICHLET_EVAP__genmod
