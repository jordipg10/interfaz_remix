        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 10:31:47 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MAIN_PDE__genmod
          INTERFACE 
            SUBROUTINE MAIN_PDE(THIS,ROOT)
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M, ONLY :                                         &
     &          PDE_1D_C
              CLASS (PDE_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE MAIN_PDE
          END INTERFACE 
        END MODULE MAIN_PDE__genmod
