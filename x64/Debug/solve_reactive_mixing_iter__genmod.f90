        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:36 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_REACTIVE_MIXING_ITER__genmod
          INTERFACE 
            SUBROUTINE SOLVE_REACTIVE_MIXING_ITER(THIS,C1_OLD,          &
     &MIXING_RATIOS,CONC_OLD,POROSITY,DELTA_T,SOLVER)
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
     &          INF_NORM_VEC_REAL,                                      &
     &          COMPUTE_C_TILDE,                                        &
     &          REACTIVE_ZONE_C,                                        &
     &          AQ_PHASE_C,                                             &
     &          SOLID_CHEMISTRY_C
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
              REAL(KIND=8), INTENT(IN) :: MIXING_RATIOS(:)
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:,:)
              REAL(KIND=8), INTENT(IN) :: POROSITY
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              INTERFACE 
                SUBROUTINE SOLVER(THIS,C1_OLD,C2NC,C_TILDE,CONC_NC,     &
     &DELTA_T)
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
                  USE AQUEOUS_CHEMISTRY_M, ONLY :                       &
     &              AQUEOUS_CHEMISTRY_C,                                &
     &              INF_NORM_VEC_REAL,                                  &
     &              COMPUTE_C_TILDE,                                    &
     &              REACTIVE_ZONE_C,                                    &
     &              AQ_PHASE_C,                                         &
     &              SOLID_CHEMISTRY_C
                  CLASS (AQUEOUS_CHEMISTRY_C), INTENT(INOUT) :: THIS
                  REAL(KIND=8), INTENT(IN) :: C1_OLD(:)
                  REAL(KIND=8), INTENT(IN) :: C2NC(:)
                  REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
                  REAL(KIND=8), INTENT(INOUT) :: CONC_NC(:)
                  REAL(KIND=8), INTENT(IN) :: DELTA_T
                END SUBROUTINE SOLVER
              END INTERFACE 
            END SUBROUTINE SOLVE_REACTIVE_MIXING_ITER
          END INTERFACE 
        END MODULE SOLVE_REACTIVE_MIXING_ITER__genmod
