        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:59 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RES_JAC_RES_ANAL_IDEAL__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RES_JAC_RES_ANAL_IDEAL(THIS,INDICES_ICON,&
     &N_ICON,INDICES_CONSTRAINS,CTOT,DC2AQ_DC1,RES,JAC_RES)
              USE VECTORS_M
              USE MATRICES_M, ONLY :                                    &
     &          INT_ARRAY_C
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
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              CLASS (INT_ARRAY_C), INTENT(IN) :: INDICES_ICON
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:,:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              REAL(KIND=8), INTENT(IN) :: DC2AQ_DC1(:,:)
              REAL(KIND=8), INTENT(OUT) :: RES(:)
              REAL(KIND=8), INTENT(OUT) :: JAC_RES(:,:)
            END SUBROUTINE COMPUTE_RES_JAC_RES_ANAL_IDEAL
          END INTERFACE 
        END MODULE COMPUTE_RES_JAC_RES_ANAL_IDEAL__genmod
