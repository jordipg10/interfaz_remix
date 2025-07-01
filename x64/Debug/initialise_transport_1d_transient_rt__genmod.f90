        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:58:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_TRANSPORT_1D_TRANSIENT_RT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT_RT(THIS,ROOT)
              USE TRANSPORT_STAB_PARAMS_M, ONLY :                       &
     &          STAB_PARAMS_TPT_C
              USE DIFFUSION_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_HETEROG_C,                                   &
     &          TIME_DISCR_C
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE CHAR_PARAMS_M, ONLY :                                 &
     &          CHAR_PARAMS_C
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT_RT
          END INTERFACE 
        END MODULE INITIALISE_TRANSPORT_1D_TRANSIENT_RT__genmod
