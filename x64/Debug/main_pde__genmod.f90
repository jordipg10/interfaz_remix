        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:53 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAIN_PDE__genmod
          INTERFACE 
            SUBROUTINE MAIN_PDE(THIS,ROOT)
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M, ONLY :                                         &
     &          PDE_1D_C
              CLASS (PDE_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE MAIN_PDE
          END INTERFACE 
        END MODULE MAIN_PDE__genmod
