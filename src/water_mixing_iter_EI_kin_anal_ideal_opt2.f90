!> Computes variable activity species concentrations after iteration of WMA using Euler implicit in chemical reactions for kinetic system
!! We assume all primary species are aqueous
!! The Jacobians are computed analytically
subroutine water_mixing_iter_EI_kin_anal_ideal_opt2(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, initialise_iterative_method, &
        Newton_EfI_rk_kin_aq_anal_ideal
    implicit none
!> Arguments
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object at current time step
    real(kind=8), intent(in) :: c1_old(:) !> aqueous variable activity concentrations at previous time step
    real(kind=8), intent(in) :: c_tilde(:) !> variable activity species concentrations after mixing
    !real(kind=8), intent(in) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(in) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
    real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    real(kind=8), intent(inout) :: mix_ratio_Rk !> mixing ratio of kinetic reaction amount in this target
    !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(out) :: conc_nc(:) !> variable activity species concentrations (already allocated)
    !real(kind=8), intent(in), optional :: porosity !> at this target
!> Variables
    real(kind=8), allocatable :: conc_nc_old(:) !> var act concentrations at beginning of current time step
    !real(kind=8), allocatable :: R_tilde_up(:) !> reaction amount contributions after mixing
    !real(kind=8), allocatable :: R_tilde_down(:) !> reaction amount contributions after mixing
    real(kind=8), allocatable :: R_tilde(:) !> reaction amount contributions after mixing
    real(kind=8), allocatable :: rk_old(:) !> old kinetic reaction rates
    !integer(kind=4) :: n_nc !> number of variable activity species
    !integer(kind=4) :: i !> counter of species
    integer(kind=4) :: k_div !> counter of time step divisions
    integer(kind=4) :: k !> counter of time steps
    integer(kind=4) :: niter !> number of Newton iterations
    real(kind=8) :: mu !> Newton initialisation parameter
    real(kind=8) :: Delta_t_bis !> reduced time step
    !real(kind=8) :: mix_ratio_Rk_new !> new mixing ratio
    logical :: CV_flag !> convergence flag
    !real(kind=8), parameter :: alpha=1.05 !> mixing ratio factor
!> Pre-process
    !n_nc=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
    !mix_ratio_Rk_new=mix_ratio_Rk !> we initialise new mixing ratio (chapuza)
    !> We get var act species concentrations
        !conc_nc_old=this%get_conc_nc()
    !> Newton initialisation parameter
        mu=0d0 
!> Process
        k=0 !> we initialise time step counter
        k_div=0 !> we initialise time step division counter
        Delta_t_bis=Delta_t !> we initialise auxiliary time step
        !R_tilde_up=Delta_t*rk_tilde_up !> we compute reaction part of concentrations after mixing
        !R_tilde_down=Delta_t*rk_tilde_down !> we compute reaction part of concentrations after mixing
        R_tilde=Delta_t*rk_tilde !> we compute reaction part of concentrations after mixing
        ! !> We check negative concentrations after mixing
        ! i=1
        ! do
        !     !if (c_tilde(i)+R_tilde_up(i)<=0d0) then
        !         !error stop "Negative concentration after mixing upstream waters in subroutine water_mixing_iter_EI_kin_aq_anal_ideal_opt2"
        !     !if (c_tilde(i)+R_tilde_up(i)+R_tilde_down(i)<=0d0 .and. R_tilde_down(i)<0d0) then
        !     if (c_tilde(i)+R_tilde(i)<=0d0) then
        !         print *, "Warning: negative concentration after mixing in subroutine water_mixing_iter_EI_kin_aq_anal_ideal_opt2"
        !         print *, "Species index: ", i
        !         !R_tilde_down(i)=R_tilde_down(i)/2d0 !> autentica chapuza
        !         !R_tilde_down=R_tilde_down*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new) !> we reduce contributions of reaction amounts in this target (chapuza)
        !         !R_tilde_up=R_tilde_up*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new) !> we reduce contributions of reaction amounts in this target (chapuza)
        !         R_tilde=R_tilde*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new) !> we reduce contributions of reaction amounts in this target (chapuza)
        !         mix_ratio_Rk_new=mix_ratio_Rk_new*alpha !> we increase mixing ratio in this target (chapuza)
        !         print *, mix_ratio_Rk_new*Delta_t
        !         if (Delta_t*mix_ratio_Rk_new>1d0) then
        !             !> chapuza
        !             mix_ratio_Rk_new=1d0/Delta_t
        !             !R_tilde_down=0d0
        !             !R_tilde_up=0d0
        !             R_tilde=0d0
        !         end if
        !     else if (i<n_nc) then
        !         i=i+1
        !     else
        !         exit
        !     end if
        ! end do
    !> We get old kinetic reaction rates
        rk_old=this%get_rk()
    !!> Update old kinetic reaction rates
        !call this%update_rk_old() !> we update old values of kinetic reaction rates
        !call this%solid_chemistry%update_rk_old() !> we update old values of kinetic reaction rates
    !> Initial guess variabel activity concentrations
        !call initialise_iterative_method(c1_old,conc_nc_old,mu,conc_nc)
        call this%initialise_conc_nc_iterative_method_ideal(mu,conc_nc) !> we initialise Newton method
        !call this%set_conc_var_act_species(conc_nc) !> chapuza
        !call this%set_conc_var_act_species(c_tilde+R_tilde+Delta_t*mix_ratio_Rk*matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,rk_old)) !> chapuza
        call this%set_act_aq_species() !> chapuza
        !conc_nc=this%get_conc_nc() !> we get variable activity species concentrations
        do
            do !> loop until convergence is reached
            !> We apply Newton method to compute aqueous concentrations
                call this%Newton_EI_rk_kin_aq_anal_ideal_opt2(c_tilde,rk_old,R_tilde,mix_ratio_Rk,Delta_t_bis,theta,conc_nc,&
                    niter,CV_flag)
            !> We check convergence
                if (CV_flag.eqv..false.) then !> no CV
                    if (mu<1d0) then
                        mu=mu+0.25 !> we increase Newton initialisation parameter
                        !call initialise_iterative_method(c1_old,c1,mu,conc_nc)
                        !call this%set_conc_var_act_species(conc_nc) !> chapuza
                        !call this%set_act_aq_species() !> chapuza
                    else
                        mu=0d0
                        k_div=k_div+1
                        if (k_div>this%solid_chemistry%reactive_zone%CV_params%k_div_max) then
                            print *, "Too many time step divisions"
                            R_tilde=0d0 !> we set reaction contributions to zero
                            mix_ratio_Rk=1d0 !> we set mixing ratio to 1
                        end if
                        Delta_t_bis=Delta_t_bis/2d0 !> we reduce time step
                    end if
                    !call initialise_iterative_method(c1_old,conc_nc_old,mu,conc_nc) !> chapuza
                    call this%initialise_conc_nc_iterative_method_ideal(mu,conc_nc) !> we reinitialise Newton method
                    !call this%set_conc_var_act_species(conc_nc) !> chapuza
                    call this%set_act_aq_species() !> chapuza
                    !conc_nc=this%get_conc_nc() !> we get variable activity species concentrations
                else 
                    k=k+1
                    exit
                end if
            end do
            if (abs(2d0**(k_div)-k)<this%solid_chemistry%reactive_zone%CV_params%zero) then
                exit
            end if
        end do
!< Post-process
    !call this%set_conc_var_act_species(conc_nc)
        deallocate(R_tilde,rk_old)
end subroutine