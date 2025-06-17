        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:48 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_K_RKF45__genmod
          INTERFACE 
            FUNCTION COMPUTE_K_RKF45(THIS,DELTA_T,CONC_RK4) RESULT(K)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          PROD_MAT_VEC
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: CONC_RK4(:)
              REAL(KIND=8) ,ALLOCATABLE :: K(:,:)
            END FUNCTION COMPUTE_K_RKF45
          END INTERFACE 
        END MODULE COMPUTE_K_RKF45__genmod
