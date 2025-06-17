        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_EI_DELTA_T_HOMOG__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_EI_DELTA_T_HOMOG(THIS,THETA,TIME_OUT,  &
     &OUTPUT)
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C,                               &
     &          DIFFUSION_1D_TRANSIENT_C,                               &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          PDE_1D_TRANSIENT_C,                                     &
     &          TRIDIAG_MATRIX_C,                                       &
     &          THOMAS
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: THETA
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_PDE_EI_DELTA_T_HOMOG
          END INTERFACE 
        END MODULE SOLVE_PDE_EI_DELTA_T_HOMOG__genmod
