        !COMPILER-GENERATED INTERFACE MODULE: Wed Jun 18 20:03:51 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_RT_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_RT_1D(THIS,ROOT,PATH_PY)
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
              USE CHEMISTRY_LAGR_M
              USE RT_1D_M, ONLY :                                       &
     &          RT_1D_C,                                                &
     &          RT_1D_TRANSIENT_C
              CLASS (RT_1D_C), INTENT(IN) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              CHARACTER(*) ,OPTIONAL, INTENT(IN) :: PATH_PY
            END SUBROUTINE WRITE_RT_1D
          END INTERFACE 
        END MODULE WRITE_RT_1D__genmod
