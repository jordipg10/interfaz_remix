        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:26 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_B_VEC_LIN_SYST__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_B_VEC_LIN_SYST(THIS,THETA,CONC_OLD,B_VEC,&
     &K)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          DIAG_MATRIX_C
              CLASS (PDE_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: CONC_OLD(:)
              REAL(KIND=8), INTENT(OUT) :: B_VEC(:)
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_B_VEC_LIN_SYST
          END INTERFACE 
        END MODULE COMPUTE_B_VEC_LIN_SYST__genmod
