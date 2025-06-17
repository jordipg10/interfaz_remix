        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:25 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_TRANS_MAT_TPT__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_TRANS_MAT_TPT(THIS)
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE BCS_M
              USE SPATIAL_DISCR_1D_M
              USE TRANSPORT_M, ONLY :                                   &
     &          TRANSPORT_1D_C
              CLASS (TRANSPORT_1D_C) :: THIS
            END SUBROUTINE COMPUTE_TRANS_MAT_TPT
          END INTERFACE 
        END MODULE COMPUTE_TRANS_MAT_TPT__genmod
