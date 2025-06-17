        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_B_ODE__genmod
          INTERFACE 
            FUNCTION COMPUTE_B_ODE(THIS) RESULT(B)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          DIAG_MATRIX_C
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8) ,ALLOCATABLE :: B(:)
            END FUNCTION COMPUTE_B_ODE
          END INTERFACE 
        END MODULE COMPUTE_B_ODE__genmod
