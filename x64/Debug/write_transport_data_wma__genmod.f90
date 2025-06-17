        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:39 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_DATA_WMA__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_DATA_WMA(THIS,UNIT)
              USE TRANSPORT_STAB_PARAMS_M
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE CHAR_PARAMS_M
              USE BCS_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              INTEGER(KIND=4), INTENT(IN) :: UNIT
            END SUBROUTINE WRITE_TRANSPORT_DATA_WMA
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_DATA_WMA__genmod
