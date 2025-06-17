        !COMPILER-GENERATED INTERFACE MODULE: Tue Jun 17 18:50:18 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE READ_TRANSPORT_DATA_WMA__genmod
          INTERFACE 
            SUBROUTINE READ_TRANSPORT_DATA_WMA(THIS,ROOT)
              USE TRANSPORT_STAB_PARAMS_M
              USE TRANSPORT_PROPERTIES_HETEROG_M, ONLY :                &
     &          TPT_PROPS_HETEROG_C
              USE CHAR_PARAMS_M
              USE TIME_DISCR_M, ONLY :                                  &
     &          TIME_DISCR_HOMOG_C
              USE BCS_M, ONLY :                                         &
     &          BCS_T
              USE SPATIAL_DISCR_1D_M, ONLY :                            &
     &          MESH_1D_EULER_HOMOG_C
              USE TRANSPORT_TRANSIENT_M, ONLY :                         &
     &          TRANSPORT_1D_TRANSIENT_C
              CLASS (TRANSPORT_1D_TRANSIENT_C) :: THIS
              CHARACTER(*), INTENT(IN) :: ROOT
            END SUBROUTINE READ_TRANSPORT_DATA_WMA
          END INTERFACE 
        END MODULE READ_TRANSPORT_DATA_WMA__genmod
