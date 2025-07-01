        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:59 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_E_MAT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_E_MAT(THIS,E_MAT,K)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_HETEROG_C,                                   &
     &          TRIDIAG_MATRIX_C
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              TYPE (TRIDIAG_MATRIX_C), INTENT(OUT) :: E_MAT
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_E_MAT
          END INTERFACE 
        END MODULE COMPUTE_E_MAT__genmod
