!> Computes variable activity species concentrations after time iteration of WMA-EE method in kinetic chemical system
!! We assume concentrations are expressed in molalities
subroutine reactive_mixing_iter_EE_kin(this,c1_old,c_hat,mix_ratio_r,Delta_t,theta,conc_nc,conc_comp)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: c1_old(:) !> primary species concentrations before previous time step
    !real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    !real(kind=8), intent(in) :: R_tilde(:) !> reaction contributions after mixing
    real(kind=8), intent(in) :: c_hat(:) !> mixing contribution to variable activity species concentrations
    real(kind=8), intent(in) :: mix_ratio_r !> mixing ratio for reaction contributions
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity species concentrations after time step
    real(kind=8), intent(inout) :: conc_comp(:) !> component concentrations after time step
    !real(kind=8), intent(in) :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: conc_react(:) !> chemical part of concentrations
    real(kind=8), allocatable :: rk_new(:) !> chemical part of concentrations
!> Pre-process
    allocate(conc_react(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    !allocate(rk_new(this%indices_rk%num_cols))
    ! !> Update old kinetic reaction rates
    ! call this%update_rk_old() !> we update old values of kinetic reaction rates
    ! call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> Aqueous concentrations
    !call this%compute_rk_new(rk_new) !> we compute chemical part of concentrations
    call this%compute_react_term_EE_kin(Delta_t,mix_ratio_r,conc_react) !> we compute chemical part of concentrations
    !conc_react=Delta_t*mix_ratio_r*rk_new !> we compute chemical part of concentrations
    conc_nc=c_hat+conc_react !> we sum chemical and mixing parts
    conc_comp=conc_nc !> we compute component concentrations (trivial)
!> We update aqueous chemistry attributes
    call this%set_conc_var_act_species(conc_nc) !> chapuza
    call this%compute_log_act_coeffs_aq_chem()
    call this%compute_activities_aq()
    call this%compute_pH()
    call this%compute_salinity()
!> Post-process
    deallocate(conc_react)
end subroutine