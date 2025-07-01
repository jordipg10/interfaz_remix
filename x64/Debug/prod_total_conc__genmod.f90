        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:08 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_TOTAL_CONC__genmod
          INTERFACE 
            SUBROUTINE PROD_TOTAL_CONC(THIS,A_MAT,TIME)
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE DIFFUSION_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TRIDIAG_MATRIX_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
              CLASS (TRIDIAG_MATRIX_C), INTENT(IN) :: A_MAT
              REAL(KIND=8) ,OPTIONAL, INTENT(IN) :: TIME
            END SUBROUTINE PROD_TOTAL_CONC
          END INTERFACE 
        END MODULE PROD_TOTAL_CONC__genmod
