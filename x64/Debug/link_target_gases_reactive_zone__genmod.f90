        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:58 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LINK_TARGET_GASES_REACTIVE_ZONE__genmod
          INTERFACE 
            SUBROUTINE LINK_TARGET_GASES_REACTIVE_ZONE(THIS,I,          &
     &DOM_INDICES,EXT_INDICES)
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
              USE AQUEOUS_CHEMISTRY_M
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
              CLASS (CHEMISTRY_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: I
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: DOM_INDICES(:&
     &)
              INTEGER(KIND=4) ,ALLOCATABLE, INTENT(OUT) :: EXT_INDICES(:&
     &)
            END SUBROUTINE LINK_TARGET_GASES_REACTIVE_ZONE
          END INTERFACE 
        END MODULE LINK_TARGET_GASES_REACTIVE_ZONE__genmod
