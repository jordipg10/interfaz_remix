        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 18 20:04:01 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WATER_MIXING_ITER_EI_KIN_ANAL_IDEAL_OPT2_LUMP__genmod
          INTERFACE 
            SUBROUTINE WATER_MIXING_ITER_EI_KIN_ANAL_IDEAL_OPT2_LUMP(   &
     &THIS,C1_OLD,C_TILDE,DELTA_T,THETA,CONC_NC)
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
     &          INITIALISE_ITERATIVE_METHOD,                            &
     &          NEWTON_EI_RK_KIN_AQ_ANAL_IDEAL_OPT2_LUMP
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(INOUT) :: CONC_NC(:)
            END SUBROUTINE WATER_MIXING_ITER_EI_KIN_ANAL_IDEAL_OPT2_LUMP
          END INTERFACE 
        END MODULE WATER_MIXING_ITER_EI_KIN_ANAL_IDEAL_OPT2_LUMP__genmod
