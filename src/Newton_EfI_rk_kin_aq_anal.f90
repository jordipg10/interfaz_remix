!!> Performs Newton method in reactive mixing iteration using Euler fully implicit in chemical reactions
!!! We assume the chemical system is kinetic
!
!!> Computes concentration of aqueous species in target j at time step k+1
!    
!!> We assume all species are aqueous
!    
!!> The Jacobians are computed analytically
!
!subroutine Newton_EfI_rk_kin_aq_anal(this,c_tilde,porosity,Delta_t,niter,CV_flag)
!    use aqueous_chemistry_m, only: aqueous_chemistry_c, inf_norm_vec_real, LU_lin_syst
!    implicit none
!!> Arguments
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: c_tilde(:)
!    real(kind=8), intent(in) :: porosity !> in target j
!    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
!    integer(kind=4), intent(out) :: niter !> number of iterations
!    logical, intent(out) :: CV_flag !> TRUE: CV, FALSE: no CV
!!> Variables
!    real(kind=8), allocatable :: rk(:) !> Jacobian of kinetic reactions
!    real(kind=8), allocatable :: drk_dc(:,:) !> Jacobian of kinetic reactions
!    real(kind=8), allocatable :: dfk_dc(:,:) !> Jacobian Newton residual
!    real(kind=8), allocatable :: Delta_c_aq(:) !> aqueous concentrations increment
!    real(kind=8), allocatable :: fk(:) !> Newton residual
!    integer(kind=4) :: n_aq !> number of aqueous species
!    integer(kind=4) :: niter_Picard !> number of iterations in Picard method
!!> Pre-process
!    niter=0
!    CV_flag=.false.
!    n_aq=this%aq_phase%num_species
!    allocate(dfk_dc(n_aq,n_aq),Delta_c_aq(n_aq),drk_dc(this%solid_chemistry%reactive_zone%chem_syst%num_kin_reacts,n_aq))
!!> Process
!    !> We start Newton loop
!        do
!        !> we update number of iterations
!            niter=niter+1 
!            if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
!                print *, "Too many Newton iterations"
!                exit
!            end if
!        !> we compute kinetic reaction rates and its Jacobian analytically
!            call this%compute_rk_Jac_rk_anal(rk,drk_dc)
!         !> Newton residual
!            fk=this%concentrations-c_tilde-(Delta_t/porosity)*matmul(transpose(this%solid_chemistry%reactive_zone%chem_syst%Sk),rk)
!        !> We check convergence
!            if (inf_norm_vec_real(fk)<this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> CV reached
!                CV_flag=.true.
!                exit
!            else
!            !> compute Jacobian of Newton residual
!                call this%compute_dfk_dc_aq_EfI(drk_dc,porosity,Delta_t,dfk_dc) 
!            !> solve linear system dfk_dc*Delta_c_aq=-fk using LU decomposition
!                call LU_lin_syst(dfk_dc,-fk,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c_aq) 
!                if (inf_norm_vec_real(Delta_c_aq)<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then 
!                    print *, "Newton solution not accurate enough"
!                    exit
!                else
!                    call this%update_conc_aq_species(Delta_c_aq) !> updates concentration aqueous species and Delta_c_aq
!                end if
!            end if
!        end do
!    !> We compute aqueous chemistry attributes
!        call this%compute_ionic_act()
!        call this%aq_phase%compute_log_act_coeffs_aq_phase(this%ionic_act,this%params_aq_sol,this%log_act_coeffs)
!        call this%compute_log_act_coeff_wat()
!        call this%compute_activities_aq()
!        call this%compute_pH()
!        call this%compute_salinity()
!!> Post-process
!    deallocate(dfk_dc,Delta_c_aq,drk_dc,fk)
!end subroutine