        !COMPILER-GENERATED INTERFACE MODULE: Mon Sep  8 16:52:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_1D_STAT__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_1D_STAT(THIS)
              USE VECTORS_M
              USE MATRICES_M
              USE TIME_FCT_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              CLASS (PDE_1D_C) :: THIS
            END SUBROUTINE SOLVE_PDE_1D_STAT
          END INTERFACE 
        END MODULE SOLVE_PDE_1D_STAT__genmod
