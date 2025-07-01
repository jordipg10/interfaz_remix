        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:17 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_CONC_ANAL__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_CONC_ANAL(THIS,ICON,N_ICON,           &
     &INDICES_CONSTRAINS,CTOT,NITER,CV_FLAG)
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
     &          INT_ARRAY_C,                                            &
     &          LU_LIN_SYST,                                            &
     &          OUTER_PROD_VEC,                                         &
     &          INF_NORM_VEC_REAL
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: N_ICON(:)
              INTEGER(KIND=4), INTENT(IN) :: INDICES_CONSTRAINS(:,:)
              REAL(KIND=8), INTENT(IN) :: CTOT(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE INITIALISE_CONC_ANAL
          END INTERFACE 
        END MODULE INITIALISE_CONC_ANAL__genmod
