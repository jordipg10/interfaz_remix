        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:44 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE UPDATE_TIME_STEP_RKF45_NEW__genmod
          INTERFACE 
            SUBROUTINE UPDATE_TIME_STEP_RKF45_NEW(DELTA_T_OLD,TOLERANCE,&
     &CONC_RK4,CONC_RK5,DELTA_T_NEW)
              REAL(KIND=8), INTENT(IN) :: DELTA_T_OLD
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
              REAL(KIND=8), INTENT(IN) :: CONC_RK4(:)
              REAL(KIND=8), INTENT(IN) :: CONC_RK5(:)
              REAL(KIND=8), INTENT(OUT) :: DELTA_T_NEW
            END SUBROUTINE UPDATE_TIME_STEP_RKF45_NEW
          END INTERFACE 
        END MODULE UPDATE_TIME_STEP_RKF45_NEW__genmod
