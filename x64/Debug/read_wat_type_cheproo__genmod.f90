        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:35 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_WAT_TYPE_CHEPROO__genmod
          INTERFACE 
            SUBROUTINE READ_WAT_TYPE_CHEPROO(THIS,N_P_AQ,NUM_CSTR,MODEL,&
     &JAC_OPT,UNIT,NITER,CV_FLAG,SURF_CHEM)
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
     &          SOLID_CHEMISTRY_C,                                      &
     &          REACTIVE_ZONE_C,                                        &
     &          GAS_CHEMISTRY_C,                                        &
     &          GAS_C,                                                  &
     &          MINERAL_C,                                              &
     &          MINERAL_ZONE_C,                                         &
     &          AQ_SPECIES_C,                                           &
     &          SPECIES_C,                                              &
     &          GAS_PHASE_C,                                            &
     &          CAT_EXCH_C
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: N_P_AQ
              INTEGER(KIND=4), INTENT(IN) :: NUM_CSTR
              INTEGER(KIND=4), INTENT(IN) :: MODEL
              INTEGER(KIND=4), INTENT(IN) :: JAC_OPT
              INTEGER(KIND=4), INTENT(IN) :: UNIT
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
              TYPE (SOLID_CHEMISTRY_C) ,OPTIONAL, INTENT(IN) ::         &
     &SURF_CHEM
            END SUBROUTINE READ_WAT_TYPE_CHEPROO
          END INTERFACE 
        END MODULE READ_WAT_TYPE_CHEPROO__genmod
