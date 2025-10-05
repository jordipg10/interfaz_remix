!< Mineral kinetic reaction subclass
module kin_mineral_m
    use kin_mineral_params_m, only: kin_mineral_params_c
    use kin_reaction_m, only: kin_reaction_c, write_reaction_sup
    use mineral_m, only: mineral_c
    use aq_phase_m, only: aq_phase_c, aq_species_c
    implicit none
    save
    type, public, extends(kin_reaction_c) :: kin_mineral_c
        type(kin_mineral_params_c) :: params !> parameters to compute reaction rate
    contains
    !> Set
        procedure, public :: set_mineral_params
        procedure, public :: set_indices_aq_phase_min
    !> Compute
        procedure, public :: compute_rk_mineral
        procedure, public :: compute_drk_dc_mineral
    !> Write
        procedure, public :: write_params=>write_min_params
        procedure, public :: write_reaction=>write_min_reaction
    end type
        
    interface        
        subroutine compute_rk_mineral(this,act_cat,saturation,react_surf,temp,rk)
            import kin_mineral_c
            implicit none
            class(kin_mineral_c), intent(in) :: this
            real(kind=8), intent(in) :: act_cat(:)
            real(kind=8), intent(in) :: saturation !> (chapuza)
            real(kind=8), intent(in) :: react_surf
            real(kind=8), intent(in) :: temp
            real(kind=8), intent(out) :: rk
        end subroutine
                
        subroutine compute_drk_dc_mineral(this,conc_sp,act_sp,log_act_coeffs_sp,act_cat,saturation,react_surf,temp,drk_dc)
            import kin_mineral_c
            implicit none
            class(kin_mineral_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_sp(:) !> concentrations aqueous species reaction
            real(kind=8), intent(in) :: act_sp(:) !> activities aqueous species reaction
            real(kind=8), intent(in) :: log_act_coeffs_sp(:) !> activity coefficients aqueous species reaction
            real(kind=8), intent(in) :: act_cat(:) !> activities catalysers
            real(kind=8), intent(in) :: saturation
            real(kind=8), intent(in) :: react_surf
            real(kind=8), intent(in) :: temp !> Kelvin
            real(kind=8), intent(out) :: drk_dc(:) !> must be already allocated
        end subroutine
                
    end interface
    
    contains
        subroutine set_mineral_params(this,mineral_params)
            implicit none
            class(kin_mineral_c) :: this
            type(kin_mineral_params_c), intent(in) :: mineral_params
            this%params=mineral_params
        end subroutine
                       
       subroutine set_indices_aq_phase_min(this,aq_phase)
            implicit none
            class(kin_mineral_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase
            
            integer(kind=4) :: i,aq_species_ind
            type(aq_species_c) :: DOC
            logical :: flag
            
            call this%allocate_indices_aq_phase_kin_react(this%num_species-1)
            do i=1,this%num_aq_rk
                call aq_phase%is_species_in_aq_phase(this%species(i),flag,aq_species_ind)
                if (flag.eqv..true.) then
                    this%indices_aq_phase(i)=aq_species_ind
                else
                    error stop "Mineral reactant is not in aqueous phase"
                end if
            end do
       end subroutine
       
       subroutine write_min_params(this,unit)
            implicit none
            class(kin_mineral_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            write(unit,*) this%params%act_energy
            write(unit,*) this%params%num_par_reacts
            write(unit,*) this%params%k
       end subroutine
       
       subroutine write_min_reaction(this,unit)
            implicit none
            class(kin_mineral_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            call write_reaction_sup(this,unit)
            call this%write_params(unit)
        end subroutine
end module