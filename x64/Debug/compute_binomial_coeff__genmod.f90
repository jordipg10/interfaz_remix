        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:24:53 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_BINOMIAL_COEFF__genmod
          INTERFACE 
            FUNCTION COMPUTE_BINOMIAL_COEFF(M,N) RESULT(BIN_COEFF)
              INTEGER(KIND=4), INTENT(IN) :: M
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4) :: BIN_COEFF
            END FUNCTION COMPUTE_BINOMIAL_COEFF
          END INTERFACE 
        END MODULE COMPUTE_BINOMIAL_COEFF__genmod
