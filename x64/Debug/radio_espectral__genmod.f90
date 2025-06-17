        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:49:07 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE RADIO_ESPECTRAL__genmod
          INTERFACE 
            FUNCTION RADIO_ESPECTRAL(LAMBDA) RESULT(RHO)
              REAL(KIND=8), INTENT(IN) :: LAMBDA(:)
              REAL(KIND=8) :: RHO
            END FUNCTION RADIO_ESPECTRAL
          END INTERFACE 
        END MODULE RADIO_ESPECTRAL__genmod
