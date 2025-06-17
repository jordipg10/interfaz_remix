        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:47 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE__genmod
          INTERFACE 
            FUNCTION MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE(   &
     &THIS) RESULT(MASS_BAL_ERR)
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE TRANSPORT_M, ONLY :                                   &
     &          TRANSPORT_1D_C,                                         &
     &          MESH_1D_EULER_HOMOG_C
              CLASS (TRANSPORT_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8) :: MASS_BAL_ERR
            END FUNCTION MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE
          END INTERFACE 
        END MODULE                                                      &
     &MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE__genmod
