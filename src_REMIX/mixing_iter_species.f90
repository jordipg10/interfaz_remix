!> Sets aqueous species concentrations after conservative mixing iteration
!! Also computes remaining aqueous state variables
!> We assume all species are aqueous
subroutine mixing_iter_species(this,c1_old,c2nc_ig,c_tilde,conc_nc,conc_comp,porosity,Delta_t)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:)
    real(kind=8), intent(in) :: c2nc_ig(:) 
    real(kind=8), intent(in) :: c_tilde(:)
    real(kind=8), intent(out) :: conc_nc(:)
    real(kind=8), intent(out) :: conc_comp(:) !> concentration components
    real(kind=8), intent(in), optional :: porosity !> (NOT NECESSARY)
    real(kind=8), intent(in), optional :: Delta_t !> time step (NOT NECESSARY)
!> Process    
    !> We set concentrations after mixing
        call this%set_concentrations(c_tilde)
        conc_nc=THIS%get_conc_nc()
        conc_comp=conc_nc
    !> We change units
        call this%compute_molalities()
    !> We compute ionic activity
        call this%compute_ionic_act() 
    !> We compute log_10 activity coefficients
        call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs)
        call this%compute_log_act_coeff_wat()
    !> We compute activities
        call this%compute_activities_aq()
    !> We compute other aqueous chemistry properties
        call this%compute_pH()
        call this%compute_salinity()
        call this%compute_molarities()
end subroutine