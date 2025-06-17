        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:57 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_ACT_COEFF_GAINES_THOMAS__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_ACT_COEFF_GAINES_THOMAS(THIS,CATION, &
     &CEC,LOG_ACT_COEFF)
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE GAINES_THOMAS_M, ONLY :                               &
     &          GAINES_THOMAS_C,                                        &
     &          SPECIES_C
              CLASS (GAINES_THOMAS_C) :: THIS
              CLASS (SPECIES_C), INTENT(IN) :: CATION
              REAL(KIND=8), INTENT(IN) :: CEC
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFF
            END SUBROUTINE COMPUTE_LOG_ACT_COEFF_GAINES_THOMAS
          END INTERFACE 
        END MODULE COMPUTE_LOG_ACT_COEFF_GAINES_THOMAS__genmod
