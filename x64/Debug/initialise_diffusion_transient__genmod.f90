        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:33 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INITIALISE_DIFFUSION_TRANSIENT__genmod
          INTERFACE 
            SUBROUTINE INITIALISE_DIFFUSION_TRANSIENT(THIS)
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE CHAR_PARAMS_M, ONLY :                                 &
     &          CHAR_PARAMS_DIFF_C
              USE DIFFUSION_TRANSIENT_M, ONLY :                         &
     &          DIFFUSION_1D_TRANSIENT_C
              CLASS (DIFFUSION_1D_TRANSIENT_C) :: THIS
            END SUBROUTINE INITIALISE_DIFFUSION_TRANSIENT
          END INTERFACE 
        END MODULE INITIALISE_DIFFUSION_TRANSIENT__genmod
