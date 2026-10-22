!> @file chemistry_interface_m.f90
!> @brief Public interface module for the reactive-mixing interface procedures.
!> @details Groups the explicit interfaces of the five reactive-mixing interface
!>          procedures used by the interactive driver (main_interfaz.f90) and
!>          re-exports the chemistry container type chemistry_c from chemistry_m,
!>          so that a client only needs to `use chemistry_interface_m`.
!>
!>          The interface procedures themselves are implemented as standalone
!>          external subroutines (interfaz_esp_arch.f90, interfaz_comps_arch.f90).
!>          This module provides their explicit interfaces so they can be called
!>          safely and referenced through procedure pointers, e.g.
!>          `procedure(interfaz_esp_arch), pointer :: p_interfaz`.
!>
!>          Interface summary:
!>            - interfaz_esp_arch            : kinetic-only (no equilibrium reactions).
!>            - interfaz_comps_arch_eq       : equilibrium-only, rows=components, cols=targets.
!>            - interfaz_comps_arch_eq_T     : equilibrium-only, transposed (rows=targets, cols=components).
!>            - interfaz_comps_arch_eq_kin   : equilibrium + kinetic, rows=components, cols=targets.
!>            - interfaz_comps_arch_eq_kin_T : equilibrium + kinetic, transposed (rows=targets, cols=components).
module chemistry_interface_m
    !> Re-export the chemistry container type so clients only need this module.
    use chemistry_m, only: chemistry_c
    implicit none
    private
    !> Public API: chemistry class plus the five reactive-mixing interfaces.
    public :: chemistry_c
    public :: interfaz_esp_arch
    public :: interfaz_comps_arch_eq
    public :: interfaz_comps_arch_eq_T
    public :: interfaz_comps_arch_eq_kin
    public :: interfaz_comps_arch_eq_kin_T
    !> Driver helper procedures used by main_interfaz.
    public :: read_chemistry_interface
    public :: get_num_aq_comps_tar_wat

    !> Explicit interfaces for the external reactive-mixing interface subroutines.
    !> All five share an identical argument list.
    interface
        !> @brief Reactive-mixing interface for a kinetic-only system (no equilibrium reactions).
        !> @param[in,out] this         Chemistry object.
        !> @param[in]     path         Path for input and output files.
        !> @param[in]     num_aq_comps Number of aqueous components (= aqueous species here).
        !> @param[in]     file_in      Input file with concentrations after conservative transport.
        !> @param[in]     Delta_t      Time step.
        !> @param[in]     file_out     Output file with concentrations after reactive mixing.
        subroutine interfaz_esp_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: num_aq_comps
            character(len=*), intent(in) :: file_in
            real(kind=8), intent(in) :: Delta_t
            character(len=*), intent(in) :: file_out
        end subroutine interfaz_esp_arch

        !> @brief Reactive-mixing interface for an equilibrium-only system.
        !>        Input file: rows=components, columns=targets.
        !> @param[in,out] this         Chemistry object.
        !> @param[in]     path         Path for input and output files.
        !> @param[in]     num_aq_comps Number of aqueous components.
        !> @param[in]     file_in      Input file with component concentrations after conservative transport.
        !> @param[in]     Delta_t      Time step.
        !> @param[in]     file_out     Output file with concentrations after reactive mixing.
        subroutine interfaz_comps_arch_eq(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: num_aq_comps
            character(len=*), intent(in) :: file_in
            real(kind=8), intent(in) :: Delta_t
            character(len=*), intent(in) :: file_out
        end subroutine interfaz_comps_arch_eq

        !> @brief Reactive-mixing interface for an equilibrium-only system (transposed input).
        !>        Input file: rows=targets, columns=components.
        !> @param[in,out] this         Chemistry object.
        !> @param[in]     path         Path for input and output files.
        !> @param[in]     num_aq_comps Number of aqueous components.
        !> @param[in]     file_in      Input file with component concentrations after conservative transport.
        !> @param[in]     Delta_t      Time step.
        !> @param[in]     file_out     Output file with concentrations after reactive mixing.
        subroutine interfaz_comps_arch_eq_T(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: num_aq_comps
            character(len=*), intent(in) :: file_in
            real(kind=8), intent(in) :: Delta_t
            character(len=*), intent(in) :: file_out
        end subroutine interfaz_comps_arch_eq_T

        !> @brief Reactive-mixing interface for an equilibrium + kinetic system.
        !>        Input file: rows=components, columns=targets.
        !> @param[in,out] this         Chemistry object.
        !> @param[in]     path         Path for input and output files.
        !> @param[in]     num_aq_comps Number of aqueous components.
        !> @param[in]     file_in      Input file with component concentrations after conservative transport.
        !> @param[in]     Delta_t      Time step.
        !> @param[in]     file_out     Output file with concentrations and kinetic rates after reactive mixing.
        subroutine interfaz_comps_arch_eq_kin(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: num_aq_comps
            character(len=*), intent(in) :: file_in
            real(kind=8), intent(in) :: Delta_t
            character(len=*), intent(in) :: file_out
        end subroutine interfaz_comps_arch_eq_kin

        !> @brief Reactive-mixing interface for an equilibrium + kinetic system (transposed input).
        !>        Input file: rows=targets, columns=components.
        !> @param[in,out] this         Chemistry object.
        !> @param[in]     path         Path for input and output files.
        !> @param[in]     num_aq_comps Number of aqueous components.
        !> @param[in]     file_in      Input file with component concentrations after conservative transport.
        !> @param[in]     Delta_t      Time step.
        !> @param[in]     file_out     Output file with concentrations and kinetic rates after reactive mixing.
        subroutine interfaz_comps_arch_eq_kin_T(this,path,num_aq_comps,file_in,Delta_t,file_out)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: path
            integer(kind=4), intent(in) :: num_aq_comps
            character(len=*), intent(in) :: file_in
            real(kind=8), intent(in) :: Delta_t
            character(len=*), intent(in) :: file_out
        end subroutine interfaz_comps_arch_eq_kin_T

        !> @brief Read all chemistry input data for the reactive-mixing interface executable.
        !> @details Standalone external subroutine (read_chemistry_interface.f90) that
        !>          populates the chemistry object with the chemical system, initial and
        !>          boundary conditions, and reactive zones, then links the target waters
        !>          to their reactive zones.
        !> @param[in,out] this    Chemistry object to populate.
        !> @param[in]     root    Root of the input/output file names.
        !> @param[in]     path_pb Path for the problem input/output files.
        !> @param[in]     path_DB Path for the chemical databases.
        !> @param[in]     num_tar Expected number of target (domain) waters from the mesh.
        subroutine read_chemistry_interface(this,root,path_pb,path_DB,num_tar)
            import chemistry_c
            implicit none
            class(chemistry_c) :: this
            character(len=*), intent(in) :: root
            character(len=*), intent(in) :: path_pb
            character(len=*), intent(in) :: path_DB
            integer(kind=4), intent(in) :: num_tar
        end subroutine read_chemistry_interface
    end interface
contains
    !> @brief Number of aqueous components of the first target water.
    !> @details Returns the number of aqueous primary species (components) of the
    !>          first target water, assuming all target waters share the same count.
    !>          Used by main_interfaz to size the u_tilde input file (one entry per
    !>          aqueous component per target).
    !> @param[in] this Chemistry object.
    !> @return num_aq_comps Number of aqueous components of the first target water.
    function get_num_aq_comps_tar_wat(this) result(num_aq_comps)
        implicit none
        class(chemistry_c), intent(in) :: this
        integer(kind=4) :: num_aq_comps
        num_aq_comps = this%waters(this%tar_wat_indices(1))%solid_chemistry%reactive_zone% &
            speciation_alg%num_aq_prim_species
    end function get_num_aq_comps_tar_wat
end module chemistry_interface_m
