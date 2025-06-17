        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:46 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_DIFFUSION_1D__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_DIFFUSION_1D(THIS)
              USE PROPERTIES_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C,                                         &
     &          DIFF_PROPS_HETEROG_C,                                   &
     &          SPATIAL_DISCR_C,                                        &
     &          MESH_1D_EULER_HOMOG_C,                                  &
     &          MESH_1D_EULER_HETEROG_C,                                &
     &          SPATIAL_DISCR_RAD_C,                                    &
     &          BCS_T
              CLASS (DIFFUSION_1D_C) :: THIS
            END SUBROUTINE INITIALISE_DIFFUSION_1D
          END INTERFACE 
        END MODULE INITIALISE_DIFFUSION_1D__genmod
