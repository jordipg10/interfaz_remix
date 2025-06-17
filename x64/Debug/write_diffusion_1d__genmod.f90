        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:51 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_DIFFUSION_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_DIFFUSION_1D(THIS,TIME_OUT,OUTPUT)
              USE BCS_M
              USE SPATIAL_DISCR_1D_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C
              CLASS (DIFFUSION_1D_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_DIFFUSION_1D
          END INTERFACE 
        END MODULE WRITE_DIFFUSION_1D__genmod
