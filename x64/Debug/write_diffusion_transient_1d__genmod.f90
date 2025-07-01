        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:53:15 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_DIFFUSION_TRANSIENT_1D__genmod
          INTERFACE 
            SUBROUTINE WRITE_DIFFUSION_TRANSIENT_1D(THIS,TIME_OUT,OUTPUT&
     &)
              USE STABILITY_PARAMETERS_M
              USE DIFF_STAB_PARAMS_M
              USE PROPERTIES_M
              USE DIFF_PROPS_HETEROG_M
              USE DIFFUSION_M
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M
              USE VECTORS_M
              USE BCS_M
              USE SPATIAL_DISCR_RAD_M
              USE PDE_M
              USE PDE_TRANSIENT_M
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C
              CLASS (DIFFUSION_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_DIFFUSION_TRANSIENT_1D
          END INTERFACE 
        END MODULE WRITE_DIFFUSION_TRANSIENT_1D__genmod
