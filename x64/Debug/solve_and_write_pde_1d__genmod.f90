        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 15:15:00 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_AND_WRITE_PDE_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_AND_WRITE_PDE_1D(THIS,ROOT,TIME_OUT)
              USE TRANSPORT_TRANSIENT_M
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              CLASS (PDE_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
            END SUBROUTINE SOLVE_AND_WRITE_PDE_1D
          END INTERFACE 
        END MODULE SOLVE_AND_WRITE_PDE_1D__genmod
