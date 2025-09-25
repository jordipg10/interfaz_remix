!> Computes aqueous species concentrations after iteration of WMA-EE method in kinetic chemical system
!! Aqueous chemistry object is not associated to any solid chemistry
subroutine water_mixing_iter_EE_kin(this,c1_old,c2nc_ig,c_tilde,rk_tilde,Delta_t,conc_nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c,  inf_norm_vec_real, reaction_iteration_EE_kin, &
    compute_molalities, compute_ionic_act, compute_log_act_coeff_wat, compute_activities_aq, compute_pH, &
    compute_salinity, compute_molarities
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations at previous time step
    real(kind=8), intent(in) :: c2nc_ig(:) !> secondary variable activity species concentrations after mixing
    real(kind=8), intent(in) :: c_tilde(:) !> aqueous species concentrations after mixing
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations after iteration
    !real(kind=8), intent(out) :: conc_comp(:) 
    ! real(kind=8), intent(in), optional :: porosity !> NOT NECESSARY
    ! real(kind=8), intent(in), optional :: Delta_t !> time step (NOT NECESSARY)
!> Variables
    real(kind=8), allocatable :: conc_react(:)
!> Pre-process
    allocate(conc_react(size(c_tilde)))
!> Aqueous concentrations
    call this%reaction_iteration_EE_kin(Delta_t,rk_tilde,conc_react) !> chemical part
    conc_nc=c_tilde+conc_react !> we sum both parts
    !conc_comp=conc_nc
    this%concentrations=conc_nc !> all species are aqueous
!> We compute aqueous chemistry attributes
    call this%compute_molalities()
    call this%compute_ionic_act()
    call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs)
    call this%compute_log_act_coeff_wat()
    call this%compute_activities_aq()
    call this%compute_pH()
    call this%compute_salinity()
    call this%compute_molarities()
!> Post-process
    deallocate(conc_react)
end subroutine