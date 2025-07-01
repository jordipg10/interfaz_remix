        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:51 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_ACT_COEFFS_SURF_COMPL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_ACT_COEFFS_SURF_COMPL(THIS,CONC,     &
     &LOG_ACTS)
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE SOLID_M
              USE PHASE_M
              USE SURF_COMPL_M
              CLASS (SURFACE_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: LOG_ACTS(:)
            END SUBROUTINE COMPUTE_LOG_ACT_COEFFS_SURF_COMPL
          END INTERFACE 
        END MODULE COMPUTE_LOG_ACT_COEFFS_SURF_COMPL__genmod
