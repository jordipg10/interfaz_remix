        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:59 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_1D_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_1D_TRANSIENT(THIS,TIME_OUT,OUTPUT&
     &)
              USE TRANSPORT_STAB_PARAMS_M
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE CHAR_PARAMS_M
              USE BCS_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_TRANSPORT_1D_TRANSIENT
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_1D_TRANSIENT__genmod
