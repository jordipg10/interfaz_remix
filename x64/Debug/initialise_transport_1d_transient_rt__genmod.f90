        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 15:14:56 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_TRANSPORT_1D_TRANSIENT_RT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT_RT(THIS,ROOT,  &
     &MESH_TYPE)
              USE TRANSPORT_STAB_PARAMS_M, ONLY :                       &
     &          STAB_PARAMS_TPT_C
              USE CONC_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_HETEROG_C,                                   &
     &          TIME_DISCR_C
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              INTEGER(KIND=4), INTENT(IN) :: MESH_TYPE
            END SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT_RT
          END INTERFACE 
        END MODULE INITIALISE_TRANSPORT_1D_TRANSIENT_RT__genmod
