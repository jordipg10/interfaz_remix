!> Sets aqueous species concentrations after conservative mixing iteration
!! Also computes remaining aqueous state variables
!> We assume all species are aqueous
subroutine mixing_iter_species(this,c1_old,c_hat, &
    mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:)
    !real(kind=8), intent(in) :: c2v_ig(:) 
    real(kind=8), intent(in) :: c_hat(:)
    real(kind=8), intent(in) :: mix_ratio_r !> porosity (NOT NECESSARY)
    real(kind=8), intent(in) :: Delta_t !> time step (NOT NECESSARY)
    real(kind=8), intent(in) :: theta !> (NOT NECESSARY)
    real(kind=8), intent(inout) :: conc_nc(:)
    real(kind=8), intent(inout) :: conc_comp(:) !> concentration components
    
    
!> Process    
    !> We set concentrations after mixing
        call this%set_concentrations(c_hat)
        conc_nc=THIS%get_conc_nc()
        conc_comp=conc_nc
    !> We change units
        call this%compute_molalities()
    !> We compute ionic activity
        call this%compute_ionic_strength() 
    !> We compute log_10 activity coefficients
        call this%compute_log_act_coeffs_aq_chem()
        call this%compute_log_act_coeff_wat()
    !> We compute activities
        call this%compute_activities_aq()
    !> We compute other aqueous chemistry properties
        call this%compute_pH()
        call this%compute_salinity()
        call this%compute_molarities()
end subroutine