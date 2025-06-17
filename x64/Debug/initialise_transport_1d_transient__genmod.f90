        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:23 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_TRANSPORT_1D_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT(THIS,ROOT)
              USE TRANSPORT_STAB_PARAMS_M
              USE CHAR_PARAMS_M, ONLY :                                 &
     &          CHAR_PARAMS_TPT_C
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE INITIALISE_TRANSPORT_1D_TRANSIENT
          END INTERFACE 
        END MODULE INITIALISE_TRANSPORT_1D_TRANSIENT__genmod
