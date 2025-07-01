        !COMPILER-GENERATED INTERFACE MODULE: Tue Jul  1 14:52:15 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE LU__genmod
          INTERFACE 
            SUBROUTINE LU(A,TOL,L,U,ERROR)
              REAL(KIND=8), INTENT(IN) :: A(:,:)
              REAL(KIND=8), INTENT(IN) :: TOL
              REAL(KIND=8), INTENT(OUT) :: L(:,:)
              REAL(KIND=8), INTENT(OUT) :: U(:,:)
              LOGICAL(KIND=4), INTENT(OUT) :: ERROR
            END SUBROUTINE LU
          END INTERFACE 
        END MODULE LU__genmod
