!> Computes variable activity species concentrations after time iteration of WMA-EE method in kinetic chemical system
!! We assume concentrations are expressed in molalities
subroutine water_mixing_iter_EE_kin_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c,  inf_norm_vec_real, reaction_iteration_EE_kin, &
        compute_molalities, compute_ionic_act, compute_log_act_coeff_wat, compute_activities_aq, compute_pH, &
        compute_salinity, compute_molarities
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c1_old(:) !> variable activity concentrations after previous time step
    !real(kind=8), intent(in) :: c2nc_ig(:)
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates contribution after mixing
    real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target (NOT USED)
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations after current time step (must be already allocated)
    !real(kind=8), intent(out) :: conc_comp(:)
    ! real(kind=8), intent(in), optional :: porosity !>  porosity
    ! real(kind=8), intent(in), optional :: Delta_t !> time step
!> Variables
    real(kind=8), allocatable :: conc_react(:) !> chemical part of concentrations
!> Pre-process
    allocate(conc_react(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
    ! !> Update old kinetic reaction rates
    ! call this%update_rk_old() !> we update old values of kinetic reaction rates
    ! call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
!> Aqueous concentrations
    call this%reaction_iteration_EE_kin(Delta_t,rk_tilde,conc_react) !> we compute chemical part of concentrations
    conc_nc=c_tilde+conc_react !> we sum chemical and mixing parts
!> We update aqueous chemistry attributes
    call this%set_conc_var_act_species(conc_nc) !> chapuza
    call this%set_act_aq_species()
    call this%compute_pH()
    call this%compute_salinity()
!> Post-process
    deallocate(conc_react)
end subroutine