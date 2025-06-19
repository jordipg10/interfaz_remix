        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 19 17:19:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WATER_MIXING_ITER_EE_KIN_IDEAL_LUMP__genmod
          INTERFACE 
            SUBROUTINE WATER_MIXING_ITER_EE_KIN_IDEAL_LUMP(THIS,C1_OLD, &
     &C_TILDE,DELTA_T,THETA,CONC_NC)
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
     &          REACTION_ITERATION_EE_KIN_LUMP,                         &
     &          COMPUTE_MOLALITIES,                                     &
     &          COMPUTE_IONIC_ACT,                                      &
     &          COMPUTE_LOG_ACT_COEFF_WAT,                              &
     &          COMPUTE_ACTIVITIES_AQ,                                  &
     &          COMPUTE_PH,                                             &
     &          COMPUTE_SALINITY,                                       &
     &          COMPUTE_MOLARITIES
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(OUT) :: CONC_NC(:)
            END SUBROUTINE WATER_MIXING_ITER_EE_KIN_IDEAL_LUMP
          END INTERFACE 
        END MODULE WATER_MIXING_ITER_EE_KIN_IDEAL_LUMP__genmod
