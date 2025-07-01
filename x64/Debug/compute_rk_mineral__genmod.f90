        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RK_MINERAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RK_MINERAL(THIS,ACT_CAT,SATURATION,      &
     &REACT_SURF,TEMP,RK)
              USE KIN_PARAMS_M
              USE KIN_MINERAL_PARAMS_M
              USE REACTION_M
              USE KIN_REACTION_M
              USE KIN_MINERAL_M, ONLY :                                 &
     &          KIN_MINERAL_C
              CLASS (KIN_MINERAL_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: ACT_CAT(:)
              REAL(KIND=8), INTENT(IN) :: SATURATION
              REAL(KIND=8), INTENT(IN) :: REACT_SURF
              REAL(KIND=8), INTENT(IN) :: TEMP
              REAL(KIND=8), INTENT(OUT) :: RK
            END SUBROUTINE COMPUTE_RK_MINERAL
          END INTERFACE 
        END MODULE COMPUTE_RK_MINERAL__genmod
