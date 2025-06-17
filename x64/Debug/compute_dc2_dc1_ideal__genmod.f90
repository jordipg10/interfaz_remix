        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:45 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_DC2_DC1_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_DC2_DC1_IDEAL(THIS,C1,C2,DC2_DC1)
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
     &          AQUEOUS_CHEMISTRY_C
              CLASS (AQUEOUS_CHEMISTRY_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1(:)
              REAL(KIND=8), INTENT(IN) :: C2(:)
              REAL(KIND=8), INTENT(OUT) :: DC2_DC1(:,:)
            END SUBROUTINE COMPUTE_DC2_DC1_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_DC2_DC1_IDEAL__genmod
