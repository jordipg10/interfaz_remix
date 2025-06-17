        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:20 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_TRANSPORT_1D__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_TRANSPORT_1D(THIS)
              USE TRANSPORT_PROPERTIES_HETEROG_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE TRANSPORT_M, ONLY :                                   &
     &          TRANSPORT_1D_C,                                         &
     &          TPT_PROPS_HETEROG_C,                                    &
     &          MESH_1D_EULER_HOMOG_C,                                  &
     &          MESH_1D_EULER_HETEROG_C,                                &
     &          BCS_T,                                                  &
     &          SPATIAL_DISCR_C
              CLASS (TRANSPORT_1D_C) :: THIS
            END SUBROUTINE INITIALISE_TRANSPORT_1D
          END INTERFACE 
        END MODULE INITIALISE_TRANSPORT_1D__genmod
