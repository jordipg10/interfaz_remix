        !COMPILER-GENERATED INTERFACE MODULE: Thu Jun 19 17:20:10 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_CHEMISTRY__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_CHEMISTRY(THIS,ROOT,PATH_DB,          &
     &UNIT_CHEM_SYST_FILE,UNIT_LOC_CHEM_FILE,                           &
     &UNIT_TARGET_WATERS_INIT_FILE,UNIT_OUTPUT_FILE)
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
              CLASS (CHEMISTRY_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              CHARACTER(*), INTENT(IN) :: PATH_DB
              INTEGER(KIND=4), INTENT(IN) :: UNIT_CHEM_SYST_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_LOC_CHEM_FILE
              INTEGER(KIND=4), INTENT(IN) ::                            &
     &UNIT_TARGET_WATERS_INIT_FILE
              INTEGER(KIND=4), INTENT(IN) :: UNIT_OUTPUT_FILE
            END SUBROUTINE INITIALISE_CHEMISTRY
          END INTERFACE 
        END MODULE INITIALISE_CHEMISTRY__genmod
