        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 15:15:03 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_1D_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_1D_TRANSIENT(THIS,ROOT,TIME_OUT, &
     &OUTPUT)
              USE TRANSPORT_STAB_PARAMS_M
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE CONC_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C), INTENT(IN) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_TRANSPORT_1D_TRANSIENT
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_1D_TRANSIENT__genmod
