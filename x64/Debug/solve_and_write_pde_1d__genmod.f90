        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:40 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_AND_WRITE_PDE_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_AND_WRITE_PDE_1D(THIS,TIME_OUT)
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE TRANSPORT_TRANSIENT_M
              CLASS (PDE_1D_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
            END SUBROUTINE SOLVE_AND_WRITE_PDE_1D
          END INTERFACE 
        END MODULE SOLVE_AND_WRITE_PDE_1D__genmod
