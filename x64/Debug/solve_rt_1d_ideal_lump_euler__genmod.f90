        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 17:53:49 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_RT_1D_IDEAL_LUMP_EULER__genmod
          INTERFACE 
            SUBROUTINE SOLVE_RT_1D_IDEAL_LUMP_EULER(THIS,ROOT)
              USE CHEM_OUT_OPTIONS_M
              USE VECTORS_M
              USE MATRICES_M
              USE GAS_CHEMISTRY_M
              USE TARGET_M
              USE MINERAL_ZONE_M
              USE CV_PARAMS_M
              USE REACTIVE_ZONE_LAGR_M
              USE SOLID_CHEMISTRY_M
              USE PARAMS_AQ_SOL_M
              USE LOCAL_CHEMISTRY_M
              USE AQUEOUS_CHEMISTRY_M, ONLY :                           &
     &          AQUEOUS_CHEMISTRY_C,                                    &
     &          WATER_MIXING_ITER_EE_EQ_KIN_IDEAL_LUMP,                 &
     &          WATER_MIXING_ITER_EE_KIN_IDEAL_LUMP,                    &
     &          WATER_MIXING_ITER_EI_EQ_KIN_ANAL_IDEAL_OPT2_LUMP,       &
     &          WATER_MIXING_ITER_EI_KIN_ANAL_IDEAL_OPT2_LUMP,          &
     &          MIXING_ITER_COMP_EXCH_IDEAL,                            &
     &          MIXING_ITER_COMP_IDEAL_LUMP,                            &
     &          COMPUTE_RK_TILDE_EXPL,                                  &
     &          COMPUTE_C_TILDE
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
              USE EXCH_SITES_CONV_M
              USE SURF_COMPL_M
              USE SOLID_M
              USE MINERAL_M
              USE GAS_M
              USE GAS_PHASE_M
              USE AQ_SPECIES_M
              USE PHASE_M
              USE AQ_PHASE_M
              USE PARAMS_SPEC_VOL_M
              USE PARAMS_ACT_COEFF_M
              USE SPECIES_M
              USE CHEM_SYSTEM_M
              USE RT_1D_M, ONLY :                                       &
     &          RT_1D_TRANSIENT_C,                                      &
     &          RT_1D_STAT_C,                                           &
     &          RT_1D_C
              CLASS (RT_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE SOLVE_RT_1D_IDEAL_LUMP_EULER
          END INTERFACE 
        END MODULE SOLVE_RT_1D_IDEAL_LUMP_EULER__genmod
