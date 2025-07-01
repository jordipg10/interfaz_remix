        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:58:13 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_DIFFUSION_1D__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_DIFFUSION_1D(THIS,ROOT)
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M, ONLY :                          &
     &          DIFF_PROPS_HETEROG_C
              USE VECTORS_M
              USE MATRICES_M
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE SPATIAL_DISCR_M
              USE PDE_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C
              CLASS (DIFFUSION_1D_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE INITIALISE_DIFFUSION_1D
          END INTERFACE 
        END MODULE INITIALISE_DIFFUSION_1D__genmod
