        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELIMINATE_NON_FLOWING_SPECIES__genmod
          INTERFACE 
            SUBROUTINE ELIMINATE_NON_FLOWING_SPECIES(THIS,NF_IND)
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
              CLASS (REACTIVE_ZONE_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: NF_IND(:)
            END SUBROUTINE ELIMINATE_NON_FLOWING_SPECIES
          END INTERFACE 
        END MODULE ELIMINATE_NON_FLOWING_SPECIES__genmod
