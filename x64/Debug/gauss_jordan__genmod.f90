        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:36 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAUSS_JORDAN__genmod
          INTERFACE 
            SUBROUTINE GAUSS_JORDAN(A,B,TOL,X,ERROR)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: X(:)
              INTEGER(KIND=4) ,OPTIONAL, INTENT(OUT) :: ERROR
            END SUBROUTINE GAUSS_JORDAN
          END INTERFACE 
        END MODULE GAUSS_JORDAN__genmod
