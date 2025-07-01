        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 15:00:13 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_DIFF_EE_DELTA_T_HETEROG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_DIFF_EE_DELTA_T_HETEROG(THIS,TIME_OUT,     &
     &OUTPUT)
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE DIFFUSION_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE VECTORS_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C,                               &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HETEROG_C,                                   &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          TRIDIAG_MATRIX_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_DIFF_EE_DELTA_T_HETEROG
          END INTERFACE 
        END MODULE SOLVE_DIFF_EE_DELTA_T_HETEROG__genmod
