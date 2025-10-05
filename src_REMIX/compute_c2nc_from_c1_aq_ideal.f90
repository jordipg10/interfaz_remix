!> Computes concentration of secondary species with variable activity from concentration of primary species explicitly using mass action law
!!  We assume the activity coefficients are ideal
!! This subroutine is meant to be used only when all primary species are aqueous
subroutine compute_c2nc_from_c1_aq_ideal(this,c2nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c
    implicit none
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(out) :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
    
    integer(kind=4) :: n_p_aq
    real(kind=8), allocatable :: log_c2nc(:)
    
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species

    !call this%compute_molalities() !> we change units to use the mass action law (the salinity is already computed)

    log_c2nc=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log10(this%concentrations(&
        this%ind_var_act_species(1:n_p_aq))))+this%solid_chemistry%reactive_zone%speciation_alg%logK_star
    c2nc=10**log_c2nc

    call this%set_conc_sec_var_act_species(c2nc) !> we update secondary variable activity species concentrations

    call this%compute_activities_aq() !> we update activities of aqueous species
    call this%compute_salinity() !> we update salinity to compute molarities
    !call this%compute_molarities() !> we compute molarities

 end subroutine 