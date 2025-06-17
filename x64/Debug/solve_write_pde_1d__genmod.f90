        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:39 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_WRITE_PDE_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_WRITE_PDE_1D(THIS,TIME_OUT)
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          PDE_1D_C
              CLASS (PDE_1D_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
            END SUBROUTINE SOLVE_WRITE_PDE_1D
          END INTERFACE 
        END MODULE SOLVE_WRITE_PDE_1D__genmod
