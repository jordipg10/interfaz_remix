        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:26 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NEWTON_EI_RK_KIN_AQ_ANAL_IDEAL_OPT2__genmod
          INTERFACE 
            SUBROUTINE NEWTON_EI_RK_KIN_AQ_ANAL_IDEAL_OPT2(THIS,C_TILDE,&
     &RK_OLD,RK_TILDE,MIX_RATIO_RK,DELTA_T,THETA,CONC_NC,NITER,CV_FLAG)
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
     &          LU_LIN_SYST
              CLASS (AQUEOUS_CHEMISTRY_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: C_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: RK_OLD(:)
              REAL(KIND=8), INTENT(IN) :: RK_TILDE(:)
              REAL(KIND=8), INTENT(IN) :: MIX_RATIO_RK
              REAL(KIND=8), INTENT(IN) :: DELTA_T
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(INOUT) :: CONC_NC(:)
              INTEGER(KIND=4), INTENT(OUT) :: NITER
              LOGICAL(KIND=4), INTENT(OUT) :: CV_FLAG
            END SUBROUTINE NEWTON_EI_RK_KIN_AQ_ANAL_IDEAL_OPT2
          END INTERFACE 
        END MODULE NEWTON_EI_RK_KIN_AQ_ANAL_IDEAL_OPT2__genmod
