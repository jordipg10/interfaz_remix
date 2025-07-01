        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RK_MONOD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RK_MONOD(THIS,CONC,RK)
              USE KIN_PARAMS_M
              USE MONOD_PARAMS_M
              USE REACTION_M
              USE KIN_REACTION_M
              USE REDOX_KIN_REACTION_M, ONLY :                          &
     &          REDOX_KIN_C
              CLASS (REDOX_KIN_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: RK
            END SUBROUTINE COMPUTE_RK_MONOD
          END INTERFACE 
        END MODULE COMPUTE_RK_MONOD__genmod
