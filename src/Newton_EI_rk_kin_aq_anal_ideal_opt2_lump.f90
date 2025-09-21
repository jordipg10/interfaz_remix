!> Performs Newton method in reactive mixing iteration using Euler fully implicit in chemical reactions
!! The chemical system has only kinetic reactions

!> Computes concentration of variable activity species in target j at time step k+1
        
!> The Jacobians are computed analytically
    
!> We assume all primary species are aqueous

subroutine Newton_EI_rk_kin_aq_anal_ideal_opt2_lump(this,c_tilde,rk_old,Delta_t,theta,conc_nc,niter,CV_flag)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    real(kind=8), intent(in) :: rk_old(:) !> kinetic reaction rates at previous time step
    !real(kind=8), intent(in) :: Rk_tilde(:) !> kinetic reaction amount contributions after mixing
    !real(kind=8), intent(in) :: Rk_tilde_up(:) !> kinetic reaction amount contributions after mixing
    !real(kind=8), intent(inout) :: Rk_tilde_down(:) !> kinetic reaction amount contributions after mixing
    !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
    !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (will be updated)
    integer(kind=4), intent(out) :: niter !> number of iterations
    logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
!> Variables
    real(kind=8), allocatable :: rk(:) !> kinetic reaction rates
    real(kind=8), allocatable :: drk_dc(:,:) !> Jacobian of kinetic reactions
    real(kind=8), allocatable :: dfk_dc(:,:) !> Jacobian Newton residual
    real(kind=8), allocatable :: conc_nc_old(:) !> old variable activity species concentrations
    real(kind=8), allocatable :: rk_avg(:) !> averaged kinetic reaction rates
    real(kind=8), allocatable :: Delta_c(:) !> primary concentrations difference
    real(kind=8), allocatable :: c_rk_tilde(:) !> reaction part of variable activity concentrations after mixing
    real(kind=8), allocatable :: fk(:) !> Newton residual
    integer(kind=4) :: n_p !> number of primary species
    integer(kind=4) :: n_nc !> number of variable activity species
    integer(kind=4) :: sp_ind !> species index
    logical :: zero_flag !> flag for zero concentration
    real(kind=8) :: mix_ratio_Rk_new !> new mixing ratio
    real(kind=8), parameter :: alpha=1.1 !> mixing ratios factor
!> Pre-process
    !mix_ratio_Rk_new=mix_ratio_Rk !> we initialise new mixing ratio (chapuza)
    niter=0 !> we initialise number of iterations
    CV_flag=.false. !> we initialise convergence flag
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    n_nc=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    allocate(dfk_dc(n_nc,n_nc),Delta_c(n_nc),&
        drk_dc(this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts,n_nc),&
        rk(this%solid_chemistry%mineral_zone%num_minerals_kin+&
        this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts),&
        c_rk_tilde(n_nc),fk(n_nc))
    drk_dc=0d0 !> we initialise Jacobian of kinetic reactions
    dfk_dc=0d0 !> we initialise Jacobian of Newton residual
    conc_nc_old=conc_nc !> we store old variable activity species concentrations
!> Process
    !> We compute reaction part of concentrations after mixing (llamalo de otra forma)
        !c_rk_tilde=Delta_t*rk_tilde
    !> We get old kinetic reaction rates
        !rk_old=this%get_rk()
    !> Newton loop
        do 
            niter=niter+1 !> we update number of iterations
            if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
                print *, inf_norm_vec_real(fk)
                print *, "Too many iterations in subroutine Newton_EI_rk_kin_aq_anal_ideal_opt2_lump"
                exit
            end if
            !conc_nc=this%get_conc_nc()
        !> We compute kinetic reaction rates and its Jacobian analitically
            call this%compute_rk_Jac_rk_anal(rk,drk_dc)
        !> We compute the ponderated average of kinetic reaction rates
            rk_avg=theta*rk+(1d0-theta)*rk_old
        !> Newton residual
            fk=conc_nc-c_tilde-Delta_t*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_avg)
        !> Check convergence
            if (inf_norm_vec_real(fk)<this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> CV reached
                CV_flag=.true.
                call this%compute_Rk_mean(theta,Delta_t) !> we update mean reaction amounts
                exit
            else
                call this%compute_dfk_dc_aq_EfI_lump(drk_dc,Delta_t,theta,dfk_dc) !> computes Jacobian of Newton resiudal
                call LU_lin_syst(dfk_dc,-fk,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c) !> solves linear system dfk_dc1*Delta_c1=-fk, where c1_new=c1_old+Delta_c1
                if (inf_norm_vec_real(Delta_c/conc_nc)<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2)&
                    then !> we check relative tolerance
                    print *, "Newton solution not accurate enough"
                    exit
                else
                    call this%update_conc_nc(conc_nc,Delta_c,zero_flag,sp_ind) !> updates conc_nc and Delta_c
                    if (zero_flag) then
                        return
                    end if
                    !call this%compute_molalities() !> we change units to compute activities
                    call this%set_act_aq_species() !> we update aqueous activities (chapuza)
                    !call this%compute_salinity() !> we compute salinity to change back to molarities
                    !call this%compute_molarities() !> we change units to compute reaction rates
                end if
            end if
        end do
!> Post-process
        deallocate(fk,dfk_dc,Delta_c,drk_dc,rk)
end subroutine