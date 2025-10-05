!> Computes concentration of secondary species from concentration of primary species using Picard method and mass action law
!! We assume the initial guess of secondary aqueous concentrations is already set in the aquoeus chemnistry object
!! This subroutine is meant to be used only when all primary species are aqueous
subroutine compute_c2_from_c1_aq_Picard(this,c2_init,c2,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst

    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(in) :: c2_init(:) !> chapuza (dim=n_eq)
    real(kind=8), intent(out) :: c2(:) !> chapuza (dim=n_eq)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag
!> Variables
    integer(kind=4) :: n_p,n_sec_aq,n_e,n_p_aq
    real(kind=8), allocatable :: log_gamma1_old(:),log_gamma2_old(:),log_c2_new(:),c2_old(:),c2_new(:)
!> Pre-process
    CV_flag=.false.
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_e=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
    n_sec_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
    
    allocate(log_gamma1_old(n_p),log_gamma2_old(n_e),log_c2_new(n_e))
        
    niter=0
    
    !c1=1d0 !> chapuza
    !c1(1:n_p_aq)=this%concentrations(1:n_p_aq)
    c2_old=c2_init
    log_gamma1_old=0d0 !> chapuza
    log_gamma2_old=0d0 !> chapuza
!> Process
    do
        niter=niter+1 !> we update number of iterations
        !call this%compute_molalities()
        call this%compute_ionic_act() !> we compute ionic activity
        !print *, this%ionic_act
        !call this%compute_salinity()
        !call this%compute_molarities()
        call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs) !> we compute log activity coefficients aqueous species
        call this%compute_activities_aq()
        call this%compute_log_act_coeff_wat()
        
        log_gamma1_old(1:n_p)=this%log_act_coeffs(1:n_p) !> we set log activity coefficients aqueous primary species
        log_gamma2_old(1:n_sec_aq)=this%log_act_coeffs(n_p+1:this%aq_phase%num_species) !> we set log activity coefficients aqueous secondary species
        !log_gamma2_old(n_sec_aq+1:n_e)=LOG10(act_ph)-LOG10(c2_old(n_sec_aq+1:n_e)) !> chapuza
        log_c2_new=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_gamma1_old+&
        log10(THIS%concentrations(1:n_p)))+this%solid_chemistry%reactive_zone%speciation_alg%logK_tilde-log_gamma2_old !> mass action law
        c2_new=10**log_c2_new
        
        call this%set_conc_sec_aq_species(c2_new(1:n_sec_aq)) !> we update aqueous secondary species
        if (inf_norm_vec_real((c2_new-c2_old)/c2_old)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then
            CV_flag=.true.
            exit !> CV reached
        else if (niter==this%solid_chemistry%reactive_zone%CV_params%niter_max) then
            print *, inf_norm_vec_real((c2_new-c2_old)/c2_old)
            error stop "Too many Picard iterations in speciation"
        else
            c2_old=c2_new
        end if
    end do
    c2=c2_new
    call this%compute_ionic_act() !> we compute ionic activity
    call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs) !> we compute log activity coefficients aqueous species
    call this%compute_activities_aq()
    call this%compute_log_act_coeff_wat()
end subroutine compute_c2_from_c1_aq_Picard