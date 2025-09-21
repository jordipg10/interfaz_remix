!> Redox kinetic reaction module
module redox_kin_reaction_m
    use kin_reaction_m, only: kin_reaction_c, write_reaction_sup
    use Monod_params_m, only: Monod_params_c
    use aq_phase_m, only : aq_phase_c, aq_species_c
    implicit none
    save
    type, public, extends(kin_reaction_c) :: redox_kin_c
        type(Monod_params_c) :: params !> kinetic parameters Monod
        real(kind=8) :: yield=0d0 !> yield factor in cell synthesis (QUIZAS MEJOR HACER UNA SUBCLASE PARA ESTO)
    contains
    !> Set
        procedure, public :: set_Monod_params
        procedure, public :: set_yield
    !> Compute
        procedure, public :: compute_rk_drk_dc_Monod
        procedure, public :: compute_drk_dc_Monod
        procedure, public :: compute_rk_Monod
        procedure, public :: compute_biomass_growth
        procedure, public :: compute_rk_biomass
    !> Rearrange
        procedure, public :: rearrange_indices_aq_phase_Monod
    !> Check
        !procedure, public :: check_indices_aq_phase_Monod
    !> Write
        procedure, public :: write_params=>write_Monod_params
        procedure, public :: write_reaction=>write_redox_reaction
    end type
    
    !type, public :: monod_type
    !>    PetscInt :: id
    !>    character(len=MAXWORDLENGTH) :: species_name
    !>    PetscReal :: half_saturation_constant
    !>    PetscReal :: threshold_concentration
    !>    type(monod_type), pointer :: next
    !>  end type monod_type
    
    interface
     
        subroutine compute_rk_Monod(this,conc,rk)
            import redox_kin_c
            implicit none
            class(redox_kin_c), intent(in) :: this
            real(kind=8), intent(in) :: conc(:)
            real(kind=8), intent(out) :: rk
        end subroutine

        
        subroutine compute_rk_drk_dc_Monod(this,conc,rk,drk_dc)
            import redox_kin_c
            implicit none
            class(redox_kin_c), intent(in) :: this
            real(kind=8), intent(in) :: conc(:)
            real(kind=8), intent(out) :: rk
            real(kind=8), intent(out) :: drk_dc(:)
        end subroutine
        
        subroutine compute_drk_dc_Monod(this,conc,rk,drk_dc)
            import redox_kin_c
            implicit none
            class(redox_kin_c), intent(in) :: this
            real(kind=8), intent(in) :: conc(:)
            real(kind=8), intent(in) :: rk
            real(kind=8), intent(out) :: drk_dc(:)
        end subroutine
        
      
    end interface
    
    contains

        
        subroutine set_Monod_params(this,Monod_params)
            implicit none
            class(redox_kin_c) :: this
            class(Monod_params_c), intent(in) :: Monod_params
            this%params=Monod_params
        end subroutine
        
        
        
        subroutine rearrange_indices_aq_phase_Monod(this,aq_phase_old,aq_phase_new)
            implicit none
            class(redox_kin_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase_old
            class(aq_phase_c), intent(in) :: aq_phase_new
            
            integer(kind=4) :: i,ind_new
            type(aq_species_c) :: DOC
            logical :: flag
            
            do i=1,this%params%num_terms
                call aq_phase_new%is_species_in_aq_phase(aq_phase_old%aq_species(this%indices_aq_phase(i)),flag,ind_new)
                if (flag.eqv..true.) then
                    this%indices_aq_phase(i)=ind_new
                else
                    error stop "Inhibitor/electron acceptor/electron donor is not in aqueous phase"
                end if
            end do
        end subroutine
        
        subroutine write_Monod_params(this,unit)
            implicit none
            class(redox_kin_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            write(unit,*) this%params%k_inh
        end subroutine
        
        subroutine write_redox_reaction(this,unit)
            implicit none
            class(redox_kin_c) :: this
            integer(kind=4), intent(in) :: unit !> file unit
            call write_reaction_sup(this,unit)
            call this%write_params(unit)
        end subroutine
        
        subroutine set_yield(this,yield)
        class(redox_kin_c) :: this !> redox kinetic reaction object
        real(kind=8), intent(in) :: yield !> yield factor in cell synthesis
        if (yield<0d0) error stop "Yield factor must be positive"
        this%yield=yield
        end subroutine
        
        subroutine compute_biomass_growth(this,conc,vol_frac_bf,porosity,rk_g)
        !< Computes the growth rate of the biomass
        !> We follow equation (6) from the paper of Wang et al (WRR, 2024)
        class(redox_kin_c) :: this
        real(kind=8), intent(in) :: conc(:) !> conc=[conc_microbe,conc_inh,conc_acceptor,conc_donor]
        real(kind=8), intent(in) :: vol_frac_bf !> volume fraction of biofilm
        real(kind=8), intent(in) :: porosity !> porosity
        real(kind=8), intent(out) :: rk_g !> growth rate
        rk_g=this%params%rate_cst*(1d0-vol_frac_bf/porosity)*conc(1)*this%params%k_inh(1)*conc(3)*conc(4)/((this%params%k_M(2)+&
            conc(4))*(this%params%k_M(1)+conc(3))*(this%params%k_inh(1)+conc(2)))
        end subroutine
        
        subroutine compute_rk_biomass(this,yield,rk_g)
        class(redox_kin_c) :: this
        real(kind=8), intent(in) :: yield !> yield factor in cell synthesis (NO NECESARIO QUIZAS)
        real(kind=8), intent(in) :: rk_g !> growth rate of the biomass
        real(kind=8) :: rk !> reaction rate of the biomass
        rk=rk_g/yield
        end subroutine
        
        !subroutine compute_biomass_decay(this,yield,rk_g)
        !class(redox_kin_c) :: this
        !real(kind=8), intent(in) :: yield !> yield factor in cell synthesis (NO NECESARIO QUIZAS)
        !real(kind=8), intent(in) :: rk_g !> growth rate of the biomass
        !real(kind=8) :: rk !> reaction rate of the biomass
        !rk=rk_g/yield
        !end subroutine
        
end module