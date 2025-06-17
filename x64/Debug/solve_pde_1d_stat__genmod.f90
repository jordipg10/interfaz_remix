        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:37 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_1D_STAT__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_1D_STAT(THIS)
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE TRANSPORT_M, ONLY :                                   &
     &          TRANSPORT_1D_C,                                         &
     &          PDE_1D_C,                                               &
     &          DIFFUSION_1D_C,                                         &
     &          MASS_BALANCE_ERROR_ADE_STAT_DIRICHLET_DISCHARGE,        &
     &          TRIDIAG_MATRIX_C,                                       &
     &          THOMAS
              CLASS (PDE_1D_C) :: THIS
            END SUBROUTINE SOLVE_PDE_1D_STAT
          END INTERFACE 
        END MODULE SOLVE_PDE_1D_STAT__genmod
