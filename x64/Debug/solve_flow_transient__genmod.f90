        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 15:00:12 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_FLOW_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE SOLVE_FLOW_TRANSIENT(THIS,TIME_OUT,OUTPUT)
              USE STABILITY_PARAMETERS_M
              USE STAB_PARAMS_FLOW_M
              USE PROPERTIES_M
              USE FLOW_PROPS_HETEROG_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_C,                                           &
     &          TIME_DISCR_HETEROG_C
              USE VECTORS_M
              USE MATRICES_M, ONLY :                                    &
     &          TRIDIAG_MATRIX_C
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE FLOW_TRANSIENT_M, ONLY :                              &
     &          FLOW_TRANSIENT_C
              CLASS (FLOW_TRANSIENT_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_FLOW_TRANSIENT
          END INTERFACE 
        END MODULE SOLVE_FLOW_TRANSIENT__genmod
