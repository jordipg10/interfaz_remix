        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:00 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_REACTIVE_MIXING_BCS_DEP_T__genmod
          INTERFACE 
            SUBROUTINE SOLVE_REACTIVE_MIXING_BCS_DEP_T(THIS,ROOT,UNIT,  &
     &MIXING_RATIOS,MIXING_WATERS_INDICES,TIME_DISCR_TPT,               &
     &INT_METHOD_CHEM_REACTS,SPATIAL_DISCR_TPT,D,Q,PHI,ANAL_SOL_COMP)
              USE SPATIAL_DISCR_M, ONLY :                               &
     &          SPATIAL_DISCR_C
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
     &          REACTIVE_ZONE_C,                                        &
     &          AQ_PHASE_C,                                             &
     &          GET_C1_AQ,                                              &
     &          GET_C1,                                                 &
     &          WATER_MIXING_ITER_EE_EQ_KIN_IDEAL,                      &
     &          WATER_MIXING_ITER_EFI_EQ_KIN_ANAL_IDEAL_OPT2,           &
     &          COMPUTE_C_TILDE,                                        &
     &          MIXING_ITER_COMP_IDEAL,                                 &
     &          MIXING_ITER_COMP_EXCH_IDEAL,                            &
     &          REAL_ARRAY_C,                                           &
     &          INT_ARRAY_C,                                            &
     &          SOLID_CHEMISTRY_C
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
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              CLASS (REAL_ARRAY_C), INTENT(IN) :: MIXING_RATIOS
              CLASS (INT_ARRAY_C), INTENT(IN) :: MIXING_WATERS_INDICES
              CLASS (TIME_DISCR_C), INTENT(IN) :: TIME_DISCR_TPT
              INTEGER(KIND=4), INTENT(IN) :: INT_METHOD_CHEM_REACTS
              CLASS (SPATIAL_DISCR_C), INTENT(IN) :: SPATIAL_DISCR_TPT
              REAL(KIND=8), INTENT(IN) :: D
              REAL(KIND=8), INTENT(IN) :: Q
              REAL(KIND=8), INTENT(IN) :: PHI
              REAL(KIND=8) :: ANAL_SOL_COMP
              EXTERNAL ANAL_SOL_COMP
            END SUBROUTINE SOLVE_REACTIVE_MIXING_BCS_DEP_T
          END INTERFACE 
        END MODULE SOLVE_REACTIVE_MIXING_BCS_DEP_T__genmod
