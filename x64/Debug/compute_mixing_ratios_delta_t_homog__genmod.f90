        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG(THIS,        &
     &A_MAT_LUMPED)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C,                               &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          DIAG_MATRIX_C,                                          &
     &          TRIDIAG_MATRIX_C,                                       &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          PROD_TRIDIAG_MAT_MAT,                                   &
     &          COMPUTE_INVERSE_TRIDIAG_MATRIX
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              TYPE (DIAG_MATRIX_C) ,OPTIONAL, INTENT(OUT) ::            &
     &A_MAT_LUMPED
            END SUBROUTINE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE COMPUTE_MIXING_RATIOS_DELTA_T_HOMOG__genmod
