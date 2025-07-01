        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 15:00:34 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_TPT_1D_STAT__genmod
          INTERFACE 
            SUBROUTINE SOLVE_TPT_1D_STAT(THIS)
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE VECTORS_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE DIFFUSION_M
              USE TRANSPORT_M, ONLY :                                   &
     &          TRANSPORT_1D_C,                                         &
     &          PDE_1D_C,                                               &
     &          DIFFUSION_1D_C,                                         &
     &          MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE
              CLASS (TRANSPORT_1D_C) :: THIS
            END SUBROUTINE SOLVE_TPT_1D_STAT
          END INTERFACE 
        END MODULE SOLVE_TPT_1D_STAT__genmod
