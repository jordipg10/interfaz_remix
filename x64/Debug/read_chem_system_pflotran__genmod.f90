        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:22 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_CHEM_SYSTEM_PFLOTRAN__genmod
          INTERFACE 
            SUBROUTINE READ_CHEM_SYSTEM_PFLOTRAN(THIS,PATH,UNIT)
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
              USE CHEM_SYSTEM_M, ONLY :                                 &
     &          CHEM_SYSTEM_C,                                          &
     &          SPECIES_C,                                              &
     &          AQ_SPECIES_C,                                           &
     &          MINERAL_C,                                              &
     &          GAS_C,                                                  &
     &          SURFACE_C,                                              &
     &          KIN_PARAMS_C,                                           &
     &          KIN_REACTION_C,                                         &
     &          KIN_REACTION_POLY_C,                                    &
     &          EQ_REACTION_C
              CLASS (CHEM_SYSTEM_C) :: THIS
              CHARACTER(*), INTENT(IN) :: PATH
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE READ_CHEM_SYSTEM_PFLOTRAN
          END INTERFACE 
        END MODULE READ_CHEM_SYSTEM_PFLOTRAN__genmod
