        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:11 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DRK_DC_MINERAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DRK_DC_MINERAL(THIS,CONC_SP,ACT_SP,      &
     &LOG_ACT_COEFFS_SP,ACT_CAT,SATURATION,REACT_SURF,TEMP,DRK_DC)
              USE KIN_PARAMS_M
              USE KIN_MINERAL_PARAMS_M
              USE REACTION_M
              USE KIN_REACTION_M
              USE KIN_MINERAL_M, ONLY :                                 &
     &          KIN_MINERAL_C
              CLASS (KIN_MINERAL_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC_SP(:)
              REAL(KIND=8), INTENT(IN) :: ACT_SP(:)
              REAL(KIND=8), INTENT(IN) :: LOG_ACT_COEFFS_SP(:)
              REAL(KIND=8), INTENT(IN) :: ACT_CAT(:)
              REAL(KIND=8), INTENT(IN) :: SATURATION
              REAL(KIND=8), INTENT(IN) :: REACT_SURF
              REAL(KIND=8), INTENT(IN) :: TEMP
              REAL(KIND=8), INTENT(OUT) :: DRK_DC(:)
            END SUBROUTINE COMPUTE_DRK_DC_MINERAL
          END INTERFACE 
        END MODULE COMPUTE_DRK_DC_MINERAL__genmod
