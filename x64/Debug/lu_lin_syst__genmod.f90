        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:10 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LU_LIN_SYST__genmod
          INTERFACE 
            SUBROUTINE LU_LIN_SYST(A,B,TOL,X)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: X(:)
            END SUBROUTINE LU_LIN_SYST
          END INTERFACE 
        END MODULE LU_LIN_SYST__genmod
