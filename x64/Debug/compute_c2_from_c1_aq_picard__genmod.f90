        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_C2_FROM_C1_AQ_PICARD__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_C2_FROM_C1_AQ_PICARD(THIS,C2_INIT,C2,    &
     &NITER,CV_FLAG)
              USE VECTORS_M
              USE MATRICES_M
              USE GAS_CHEMISTRY_M
              USE MINERAL_ZONE_M
              USE CV_PARAMS_M
              USE MONOD_PARAMS_M
              USE REDOX_KIN_REACTION_M
              USE KIN_PARAMS_M
              USE KIN_MINERAL_PARAMS_M
              USE KIN_MINERAL_M
              USE LIN_KIN_REACTION_M
              USE KIN_REACTION_M
              USE SPECIATION_ALGEBRA_M
              USE REACTION_M
              USE EQ_REACTION_M
              USE MINERAL_M
              USE AQ_SPECIES_M
              USE AQ_PHASE_M
              USE CHEM_SYSTEM_M
              USE GAS_M
              USE GAS_PHASE_M
              USE EXCH_SITES_CONV_M
              USE SOLID_M
              USE PHASE_M
              USE SURF_COMPL_M
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE REACTIVE_ZONE_LAGR_M
              USE SOLID_CHEMISTRY_M
              USE PARAMS_AQ_SOL_M
              USE LOCAL_CHEMISTRY_M
              USE AQUEOUS_CHEMISTRY_M, ONLY :                           &
     &          AQUEOUS_CHEMISTRY_C,                                    &
     &          INF_NORM_VEC_REAL,                                      &
     &          LU_LIN_SYST
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C2_INIT(:)
              REAL(KIND=8), INTENT(OUT) :: C2(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE COMPUTE_C2_FROM_C1_AQ_PICARD
          END INTERFACE 
        END MODULE COMPUTE_C2_FROM_C1_AQ_PICARD__genmod
