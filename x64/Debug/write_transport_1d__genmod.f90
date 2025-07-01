        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:57 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_TRANSPORT_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_TRANSPORT_1D(THIS)
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE VECTORS_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE DIFFUSION_M
              USE TRANSPORT_M, ONLY :                                   &
     &          TRANSPORT_1D_C,                                         &
     &          MESH_1D_EULER_HOMOG_C
              CLASS (TRANSPORT_1D_C), INTENT(IN) :: THIS
            END SUBROUTINE WRITE_TRANSPORT_1D
          END INTERFACE 
        END MODULE WRITE_TRANSPORT_1D__genmod
