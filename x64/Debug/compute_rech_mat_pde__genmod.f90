        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:43 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_RECH_MAT_PDE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_RECH_MAT_PDE(THIS)
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C,                               &
     &          DIFFUSION_1D_TRANSIENT_C
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE PDE_M, ONLY :                                         &
     &          PDE_1D_C
              CLASS (PDE_1D_C) :: THIS
            END SUBROUTINE COMPUTE_RECH_MAT_PDE
          END INTERFACE 
        END MODULE COMPUTE_RECH_MAT_PDE__genmod
