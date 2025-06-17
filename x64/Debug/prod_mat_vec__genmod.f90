        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:56 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE PROD_MAT_VEC__genmod
          INTERFACE 
            FUNCTION PROD_MAT_VEC(A,B) RESULT(X)
              USE MATRICES_M, ONLY :                                    &
     &          ARRAY_C,                                                &
     &          TRIDIAG_MATRIX_C
              CLASS (ARRAY_C), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B(:)
              REAL(KIND=8) ,ALLOCATABLE :: X(:)
            END FUNCTION PROD_MAT_VEC
          END INTERFACE 
        END MODULE PROD_MAT_VEC__genmod
