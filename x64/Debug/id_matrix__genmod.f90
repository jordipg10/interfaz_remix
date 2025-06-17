        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:06 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ID_MATRIX__genmod
          INTERFACE 
            FUNCTION ID_MATRIX(N)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8) :: ID_MATRIX(N,N)
            END FUNCTION ID_MATRIX
          END INTERFACE 
        END MODULE ID_MATRIX__genmod
