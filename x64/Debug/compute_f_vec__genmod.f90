        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 19 10:31:52 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_F_VEC__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_F_VEC(THIS,K)
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C,                                     &
     &          TIME_DISCR_HETEROG_C
              USE TIME_FCT_M
              USE BCS_M
              USE TARGET_M
              USE SPATIAL_DISCR_M
              USE PDE_TRANSIENT_M, ONLY :                               &
     &          PDE_1D_TRANSIENT_C
              CLASS (PDE_1D_TRANSIENT_C) :: THIS
              INTEGER(KIND=4) ,OPTIONAL, INTENT(IN) :: K
            END SUBROUTINE COMPUTE_F_VEC
          END INTERFACE 
        END MODULE COMPUTE_F_VEC__genmod
