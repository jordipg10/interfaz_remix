        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 15:00:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_DIFF_EI_DELTA_T_HOMOG_MRMT__genmod
          INTERFACE 
            SUBROUTINE SOLVE_DIFF_EI_DELTA_T_HOMOG_MRMT(THIS,THETA,     &
     &TIME_OUT,OUTPUT)
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_MODEL_M
              USE MRMT_M, ONLY :                                        &
     &          MRMT_C,                                                 &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TRIDIAG_MATRIX_C,                                       &
     &          THOMAS
              CLASS (MRMT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_DIFF_EI_DELTA_T_HOMOG_MRMT
          END INTERFACE 
        END MODULE SOLVE_DIFF_EI_DELTA_T_HOMOG_MRMT__genmod
