        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:27 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_RKF45__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_RKF45(THIS,DELTA_T_INIT,TOLERANCE)
              USE VECTORS_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HOMOG_C
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: DELTA_T_INIT
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
            END SUBROUTINE SOLVE_PDE_RKF45
          END INTERFACE 
        END MODULE SOLVE_PDE_RKF45__genmod
