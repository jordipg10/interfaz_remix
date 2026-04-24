!> \file initialise_chemistry.f90
!> \brief Top-level initialisation of the chemistry object.
!> \details
!> Reads chemical system, initial/boundary conditions, and reactive
!> zone definitions, then allocates and populates all chemistry
!> data structures.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup io
!> \see chemistry_m, read_chemistry
subroutine initialise_chemistry(this,root,path_pb,path_DB)
    use chemistry_m, only: chemistry_c
    implicit none
    class(chemistry_c) :: this
    character(len=*), intent(in) :: root
    character(len=*), intent(in) :: path_pb
    character(len=*), intent(in) :: path_DB
    !integer(kind=4), intent(in) :: unit_chem_syst_file
    !character(len=*), intent(in) :: chem_syst_file
    !integer(kind=4), intent(in) :: unit_loc_chem_file
    !character(len=*), intent(in) :: loc_chem_file
    !integer(kind=4), intent(in) :: unit_target_waters_init_file
    !character(len=*), intent(in) :: target_waters_init_file
    !integer(kind=4), intent(in) :: unit_output_file
    !character(len=*), intent(in) :: output_file
!**************************************************************************************************
!> Read chemistry
    call this%read_chemistry(root,path_pb,path_DB)
    !print *, this%waters(19)%aq_phase%num_species
end subroutine

