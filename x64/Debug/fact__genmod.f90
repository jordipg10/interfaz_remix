        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:23:15 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE FACT__genmod
          INTERFACE 
            RECURSIVE FUNCTION FACT(N) RESULT(RES)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4) :: RES
            END FUNCTION FACT
          END INTERFACE 
        END MODULE FACT__genmod
