        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:40 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_DIFFUSION_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE WRITE_DIFFUSION_TRANSIENT(THIS,TIME_OUT,OUTPUT)
              USE VECTORS_M
              USE CHAR_PARAMS_M
              USE BCS_M
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C
              CLASS (DIFFUSION_1D_TRANSIENT_C), INTENT(IN) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(IN) :: OUTPUT(:,:)
            END SUBROUTINE WRITE_DIFFUSION_TRANSIENT
          END INTERFACE 
        END MODULE WRITE_DIFFUSION_TRANSIENT__genmod
