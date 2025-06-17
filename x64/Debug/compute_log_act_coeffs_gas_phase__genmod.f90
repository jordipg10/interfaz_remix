        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:58 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_ACT_COEFFS_GAS_PHASE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_ACT_COEFFS_GAS_PHASE(THIS,           &
     &LOG_ACT_COEFFS)
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE GAS_M
              USE PHASE_M
              USE GAS_PHASE_M
              CLASS (GAS_PHASE_C) :: THIS
              REAL(KIND=8), INTENT(OUT) :: LOG_ACT_COEFFS(:)
            END SUBROUTINE COMPUTE_LOG_ACT_COEFFS_GAS_PHASE
          END INTERFACE 
        END MODULE COMPUTE_LOG_ACT_COEFFS_GAS_PHASE__genmod
