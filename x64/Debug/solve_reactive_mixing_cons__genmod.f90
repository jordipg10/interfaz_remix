        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:09 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_REACTIVE_MIXING_CONS__genmod
          INTERFACE 
            SUBROUTINE SOLVE_REACTIVE_MIXING_CONS(THIS,ROOT,            &
     &MIXING_RATIOS_CONC,MIXING_RATIOS_RK_INIT,MIXING_WATERS_INDICES,   &
     &MIXING_WATERS_INDICES_DOM,TIME_DISCR_TPT,INT_METHOD_CHEM_REACTS,  &
     &MIXING_RATIOS_RK)
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_C
              USE CHEM_OUT_OPTIONS_M
              USE VECTORS_M
              USE MATRICES_M
              USE GAS_CHEMISTRY_M
              USE MINERAL_ZONE_M
              USE CV_PARAMS_M
              USE REACTIVE_ZONE_LAGR_M
              USE SOLID_CHEMISTRY_M
              USE PARAMS_AQ_SOL_M
              USE LOCAL_CHEMISTRY_M
              USE AQUEOUS_CHEMISTRY_M, ONLY :                           &
     &          AQUEOUS_CHEMISTRY_C,                                    &
     &          REAL_ARRAY_C,                                           &
     &          INT_ARRAY_C,                                            &
     &          MIXING_ITER_COMP_IDEAL,                                 &
     &          WATER_MIXING_ITER_EE_EQ_KIN_IDEAL,                      &
     &          WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL,                &
     &          WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL_OPT2,           &
     &          MIXING_ITER_COMP_EXCH_IDEAL,                            &
     &          INF_NORM_VEC_REAL,                                      &
     &          WATER_MIXING_ITER_EE_KIN_IDEAL,                         &
     &          WATER_MIXING_ITER_EFI_KIN_ANAL_IDEAL,                   &
     &          COMPUTE_C_TILDE,                                        &
     &          COMPUTE_RK_TILDE_EXPL,                                  &
     &          COMPUTE_RK_TILDE_IMPL_OPT1,                             &
     &          WATER_MIXING_ITER_EI_EQ_KIN_ANAL_IDEAL_OPT2,            &
     &          WATER_MIXING_ITER_EI_KIN_ANAL_IDEAL_OPT2,               &
     &          COMPUTE_RK_TILDE_IMPL_OPT2,                             &
     &          COMPUTE_RK_TILDE_IMPL_OPT3,                             &
     &          COMPUTE_RK_TILDE_IMPL_OPT4
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
              USE CHEMISTRY_LAGR_M, ONLY :                              &
     &          CHEMISTRY_C
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              CLASS (REAL_ARRAY_C), INTENT(IN) :: MIXING_RATIOS_CONC
              CLASS (REAL_ARRAY_C), INTENT(IN) :: MIXING_RATIOS_RK_INIT
              CLASS (INT_ARRAY_C), INTENT(IN) :: MIXING_WATERS_INDICES
              CLASS (INT_ARRAY_C), INTENT(IN) ::                        &
     &MIXING_WATERS_INDICES_DOM
              CLASS (TIME_DISCR_C), INTENT(IN) :: TIME_DISCR_TPT
              INTEGER(KIND=4), INTENT(IN) :: INT_METHOD_CHEM_REACTS
              CLASS (REAL_ARRAY_C), INTENT(INOUT) :: MIXING_RATIOS_RK
            END SUBROUTINE SOLVE_REACTIVE_MIXING_CONS
          END INTERFACE 
        END MODULE SOLVE_REACTIVE_MIXING_CONS__genmod
