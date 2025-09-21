!> Lectura quimica
subroutine read_chemistry(this,root,path_DB)
    use chemistry_m,   only: chemistry_c
    implicit none
    class(chemistry_c) :: this
    character(len=*), intent(in) :: root
    character(len=*), intent(in) :: path_DB
    !integer(kind=4), intent(in) :: unit_chem_opts_file
    !integer(kind=4), intent(in) :: unit_chem_syst_file
    !!character(len=*), intent(in) :: chem_syst_file
    !integer(kind=4), intent(in) :: unit_loc_chem_file
    !!character(len=*), intent(in) :: loc_chem_file
    !integer(kind=4), intent(in) :: unit_target_waters_file
    !!character(len=*), intent(in) :: target_waters_init_file
    !integer(kind=4), intent(in) :: unit_output_file
    !character(len=*), intent(in) :: output_file
    
    integer(kind=4) :: i
    integer(kind=4) :: unit_chem_syst !> chemical system file unit
    integer(kind=4) :: unit_chem_opts !> chemical options file unit
    integer(kind=4) :: unit_loc_chem !> local chemistry file unit
    integer(kind=4) :: unit_tw !> target waters file unit
    integer(kind=4) :: unit_output !> output file unit
    
    !> File units (arbitrary)
    unit_chem_syst=1
    unit_chem_opts=5
    unit_loc_chem=2
    unit_tw=3
    !unit_discr=10
    !unit_tpt=11
    unit_output=4
    !unit_res=6

    call this%read_chem_opts(root,unit_chem_opts) !> read chemical options
    
    if (this%read_opt==1) then !> CHEPROO
        call this%read_chemistry_CHEPROO(root,path_DB,unit_chem_syst,unit_loc_chem,unit_tw, &
            unit_output)
    else if (this%read_opt==2) then !> PHREEQC
        !call this%read_chemistry_PHREEQC()
        error stop "PHREEQC data input not fully implemented yet"
    else if (this%read_opt==3) then !> PFLOTRAN
        !call this%read_chemistry_PFLOTRAN()
        error stop "PFLOTRAN data input not fully implemented yet"
    else
        error stop "This data input option is not implemented yet"
    end if
!> Autentica chapuza
    if (this%num_reactive_zones==1) then
        !do i=1,this%num_target_solids
        !    !call this%target_solids_init(i)%set_reactive_zone(this%reactive_zones(1))
        !    !call this%target_solids(i)%set_reactive_zone(this%reactive_zones(1))
        !end do
        !do i=1,this%num_target_gases
        !    call this%target_gases(i)%set_reactive_zone(this%reactive_zones(1))
        !end do
    else if (this%num_reactive_zones==3) then
        !print *, this%num_target_solids/2
        !do i=1,this%num_target_solids/2
        !    !call this%target_solids_init(i)%set_reactive_zone(this%reactive_zones(1))
        !    call this%target_solids(i)%set_reactive_zone(this%reactive_zones(1))
        !end do
        !!call this%target_solids_init(this%num_target_solids/2+1)%set_reactive_zone(this%reactive_zones(2))
        !call this%target_solids(this%num_target_solids/2+1)%set_reactive_zone(this%reactive_zones(2))
        !do i=this%num_target_solids/2+2,this%num_target_solids
        !    !call this%target_solids_init(i)%set_reactive_zone(this%reactive_zones(3))
        !    call this%target_solids(i)%set_reactive_zone(this%reactive_zones(3))
        !end do
        !do i=1,this%num_target_gases
        !    call this%target_gases(i)%set_reactive_zone(this%reactive_zones(1))
        !end do
    end if
    !print *, this%target_waters(17)%aq_phase%num_species
end subroutine