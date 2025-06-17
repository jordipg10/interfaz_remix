        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:42 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_PDE_1D__genmod
          INTERFACE 
            SUBROUTINE SOLVE_PDE_1D(THIS,TIME_OUT,OUTPUT)
              USE MATRICES_M, ONLY :                                    &
     &          TRIDIAG_MATRIX_C
              USE BCS_M
              USE SPATIAL_DISCR_M
              USE DIFFUSION_M, ONLY :                                   &
     &          DIFFUSION_1D_C
              CLASS (PDE_1D_C) :: THIS
              REAL(KIND=8), INTENT(IN) :: TIME_OUT(:)
              REAL(KIND=8), INTENT(OUT) :: OUTPUT(:,:)
            END SUBROUTINE SOLVE_PDE_1D
          END INTERFACE 
        END MODULE SOLVE_PDE_1D__genmod
