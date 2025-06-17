        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:54 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_SURF_COMPL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_SURF_COMPL(THIS, &
     &IONIC_ACT,LOG_ACT_COEFFS,CONC,LOG_JACOBIAN_ACT_COEFFS)
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE SOLID_M
              USE PHASE_M
              USE SURF_COMPL_M
              CLASS (SURFACE_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: IONIC_ACT
              REAL(KIND=8), INTENT(IN) :: LOG_ACT_COEFFS(:)
              REAL(KIND=8), INTENT(IN) :: CONC(:)
              REAL(KIND=8), INTENT(OUT) :: LOG_JACOBIAN_ACT_COEFFS(:,:)
            END SUBROUTINE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_SURF_COMPL
          END INTERFACE 
        END MODULE COMPUTE_LOG_JACOBIAN_ACT_COEFFS_SURF_COMPL__genmod
