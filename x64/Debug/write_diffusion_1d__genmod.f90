        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:37 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_DIFFUSION_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_DIFFUSION_1D(THIS,TIME_OUT,OUTPUT)
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE VECTORS_M
              USE MATRICES_M
              USE BCS_M
              USE SPATIAL_DISCR_1D_M
              USE PDE_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C
              CLASS (DIFFUSION_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_DIFFUSION_1D
          END INTERFACE 
        END MODULE WRITE_DIFFUSION_1D__genmod
