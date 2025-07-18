        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:56 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE(THIS,   &
     &OUT_PROD,CONC,LOG_JACOBIAN_ACT_COEFFS)
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE AQ_SPECIES_M
              USE PHASE_M
              USE AQ_PHASE_M, ONLY :                                    &
     &          AQ_PHASE_C
              CLASS (AQ_PHASE_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: OUT_PROD(:,:)
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: LOG_JACOBIAN_ACT_COEFFS(:,:)
            END SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE
          END INTERFACE 
        END MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_AQ_PHASE__genmod
