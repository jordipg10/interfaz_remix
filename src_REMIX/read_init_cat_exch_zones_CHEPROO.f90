subroutine read_init_cat_exch_zones_CHEPROO(this,unit,init_cat_exch_zones)
use chemistry_m, only: chemistry_c, solid_chemistry_c, reactive_zone_c, &
    mineral_zone_c
    use Gaines_Thomas_m, only: Gaines_Thomas_c
    implicit none
    class(chemistry_c) :: this
    integer(kind=4), intent(in) :: unit !> file
    type(solid_chemistry_c), intent(out), allocatable :: init_cat_exch_zones(:)
    !type(reactive_zone_c), intent(inout), allocatable, optional :: reactive_zones(:)
    
    integer(kind=4) :: i,j,k,ndtype,idtype,num_gas_rz,icon,num_ads_zones,num_mins_glob,num_rz
    integer(kind=4), allocatable :: valences(:)
    character(len=256) :: str,constrain,label
    real(kind=8) :: c_int,c_ext,spec_sorb_surf,CEC
    type(reactive_zone_c) :: react_zone
    type(mineral_zone_c) :: min_zone
    type(reactive_zone_c), allocatable :: aux_react_zones(:)
    type(Gaines_Thomas_c), target :: Gaines_Thomas

    read(unit,*) ndtype !> number of surface adsorption zones
    
    allocate(init_cat_exch_zones(ndtype))

!> Vamos a asumir que los complejos de superficie est�n en el mismo orden que en el sistema quimico
    num_ads_zones=0 !> counter adsorption zones
    do
        read(unit,*) idtype
        if (idtype<0 .or. idtype>ndtype) error stop
        read(unit,*) str, c_int
        read(unit,*) str, c_ext
        read(unit,*) str, spec_sorb_surf
        read(unit,*) str, CEC
        call react_zone%set_chem_syst_react_zone(this%chem_syst)
        call react_zone%set_cat_exch_zone()
        react_zone%cat_exch_zone%convention=>Gaines_Thomas !> by default
        call react_zone%set_num_solids()
        call react_zone%set_non_flowing_species()
        call react_zone%set_CV_params(this%CV_params)
        init_cat_exch_zones(idtype)%num_solids=react_zone%cat_exch_zone%num_surf_compl
        call init_cat_exch_zones(idtype)%set_reactive_zone(react_zone)
        call min_zone%set_chem_syst_min_zone(this%chem_syst) !> set chemical system for mineral zone by default
        call init_cat_exch_zones(idtype)%set_mineral_zone(min_zone) !> mineral zone by default
        call init_cat_exch_zones(idtype)%allocate_conc_solids()
        call init_cat_exch_zones(idtype)%set_CEC(CEC) !> cation exchange capacity
        call init_cat_exch_zones(idtype)%allocate_equivalents()
        call init_cat_exch_zones(idtype)%allocate_log_act_coeffs_solid_chem()
        valences=this%chem_syst%aq_phase%get_valences()
        call init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%compute_log_act_coeffs_ads_cats(valences(&
            init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%exch_cat_indices),CEC,&
            init_cat_exch_zones(idtype)%log_act_coeffs(2:init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%num_surf_compl))
        call init_cat_exch_zones(idtype)%allocate_activities()
        call init_cat_exch_zones(idtype)%set_conc_free_site() !> 'x-'
        call init_cat_exch_zones(idtype)%compute_num_solids_solid_chem() !> 
        !init_cat_exch_zones(idtype)%activities(1)=1d0-SUM(init_cat_exch_zones(idtype)%activities(2:init_cat_exch_zones(idtype)%reactive_zone%cat_exch_zone%num_surf_compl))
        num_ads_zones=num_ads_zones+1
        if (num_ads_zones==ndtype) exit
    end do
    !if (present(reactive_zones)) then
        if (allocated(this%reactive_zones)) then
            num_gas_rz=size(this%reactive_zones)
            allocate(aux_react_zones(num_gas_rz))
            do i=1,num_gas_rz
                call aux_react_zones(i)%assign_react_zone(this%reactive_zones(i))
            end do
            num_rz=num_gas_rz+ndtype*(1+num_gas_rz)
            !deallocate(this%reactive_zones)
            !allocate(this%reactive_zones(num_rz))
            call this%allocate_reactive_zones(num_rz)
            do i=1,num_gas_rz
                call this%reactive_zones(i)%assign_react_zone(aux_react_zones(i))
            end do
            do i=1,ndtype
                call this%reactive_zones(num_gas_rz+i)%set_chem_syst_react_zone(this%chem_syst)
                call this%reactive_zones(num_gas_rz+i)%set_cat_exch_zone(init_cat_exch_zones(i)%reactive_zone%cat_exch_zone)
                !call init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(num_gas_rz+i))
            end do
            do i=1,num_gas_rz
                do j=1,ndtype
                    call this%reactive_zones(num_gas_rz+i*ndtype+j)%set_chem_syst_react_zone(this%chem_syst)
                    call this%reactive_zones(num_gas_rz+i*ndtype+j)%set_gas_phase(this%reactive_zones(i)%gas_phase)
                    call this%reactive_zones(num_gas_rz+i*ndtype+j)%set_cat_exch_zone(&
                    init_cat_exch_zones(i)%reactive_zone%cat_exch_zone)
                end do
            end do
        else
            !allocate(this%reactive_zones(ndtype))
            call this%allocate_reactive_zones(ndtype)
            do i=1,ndtype
                call this%reactive_zones(i)%assign_react_zone(init_cat_exch_zones(i)%reactive_zone)
                !call init_cat_exch_zones(i)%set_reactive_zone(this%reactive_zones(i))
            end do
        end if
    !end if
end subroutine