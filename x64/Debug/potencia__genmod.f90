        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:10 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE POTENCIA__genmod
          INTERFACE 
            SUBROUTINE POTENCIA(A,Z0,TOLERANCE,RHO,NITER)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(INOUT) :: Z0(:)
              REAL(KIND=8), INTENT(IN) :: TOLERANCE
              REAL(KIND=8), INTENT(OUT) :: RHO
              INTEGER(KIND=4), INTENT(OUT) :: NITER
            END SUBROUTINE POTENCIA
          END INTERFACE 
        END MODULE POTENCIA__genmod
