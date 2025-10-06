!!> Performs Newton method in reactive mixing iteration using Euler fully implicit in chemical reactions
!!! The chemical system has equilibrium and kinetic reactions
!
!!> Computes concentration of variable activity species in target j at time step k+1
!    
!!> We assume all primary species are aqueous
!    
!!> The Jacobians are computed analytically
!
!subroutine Newton_EfI_rk_eq_kin_aq_anal(this,c2nc_ig,c_tilde,porosity,Delta_t,conc_nc,niter,CV_flag)
!    use aqueous_chemistry_m, only: aqueous_chemistry_c
!    use metodos_sist_lin_m
!    use vectors_m
!    implicit none
!!> Arguments
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
!    real(kind=8), intent(in) :: c_tilde(:)
!    real(kind=8), intent(in) :: porosity !> in target j
!    real(kind=8), intent(in) :: Delta_t !> (k+1)-th time step
!    real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
!    integer(kind=4), intent(out) :: niter !> number of iterations
!    logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
!!> Variables
!    real(kind=8), allocatable :: rk(:) !> Jacobian of kinetic reactions
!    real(kind=8), allocatable :: drk_dc(:,:) !> Jacobian of kinetic reactions
!    real(kind=8), allocatable :: dfk_dc1(:,:) !> Jacobian Newton residual
!    real(kind=8), allocatable :: c2nc_old(:) !> secondary variable activity concentrations previous iteration
!    real(kind=8), allocatable :: c2nc_new(:) !> secondary variable activity concentrations next iteration
!    real(kind=8), allocatable :: Delta_c1(:) !> primary concentrations difference
!    real(kind=8), allocatable :: u_tilde(:) !> component concentrations after mixing
!    real(kind=8), allocatable :: fk(:) !> Newton residual
!    integer(kind=4) :: n_p !> number of primary species
!    integer(kind=4) :: n_nc !> number of variable activity species
!    integer(kind=4) :: niter_Picard !> number of iterations in Picard method
!    logical :: CV_flag_Picard !> convergence flag Picard
!!> Pre-process
!    niter=0
!    CV_flag=.false.
!    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
!    n_nc=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species
!    allocate(c2nc_new(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions),dfk_dc1(n_p,n_p),Delta_c1(n_p),&
!    drk_dc(this%solid_chemistry%reactive_zone%chem_syst%num_kin_reacts,n_nc))
!    c2nc_old=c2nc_ig
!    drk_dc=0d0
!!> Process
!    !> We compute component concentrations after mixing
!        u_tilde=this%compute_u_tilde(c_tilde)
!    !> Newton loop
!        do 
!            niter=niter+1 !> we update number of iterations
!            if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
!                print *, "Too many iterations in subroutine Newton_EfI_rk_eq_kin_aq_anal"
!                exit
!            end if
!            conc_nc(1:n_p)=this%concentrations(1:n_p) !> chapuza
!        !> Compute concentration of secondary variable activity species from concentration of primary species using mass action law
!            call this%compute_c2nc_from_c1_aq_Picard(c2nc_old,c2nc_new,niter_Picard,CV_flag_Picard)
!            conc_nc(n_p+1:n_nc)=c2nc_new !> chapuza
!        !> We compute kinetic reaction rates and its Jacobian analitically
!            call this%compute_rk_Jac_rk_anal(rk,drk_dc)
!        !> Newton residual
!            fk=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)-u_tilde-(Delta_t/porosity)*&
!            matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,this%rk) 
!        !> Check convergence
!            if (inf_norm_vec_real(fk)<this%solid_chemistry%reactive_zone%CV_params%abs_tol) then !> CV reached
!                CV_flag=.true.
!                exit
!            else
!                call this%compute_dfk_dc1_aq_EfI(conc_nc(n_p+1:n_nc),drk_dc,porosity,Delta_t,dfk_dc1) !> computes Jacobian of Newton resiudal
!                call LU_lin_syst(dfk_dc1,-fk,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1) !> solves linear system dfk_dc1*Delta_c1=-fk, where c1_new=c1_old+Delta_c1
!                !call Gauss_Jordan(dfk_dc1,-fk,this%solid_chemistry%reactive_zone%CV_params%zero,Delta_c1) !> solves linear system dfk_dc1*Delta_c1=-fk, where c1_new=c1_old+Delta_c1
!                if (inf_norm_vec_real(Delta_c1/this%concentrations(1:n_p))<this%solid_chemistry%reactive_zone%CV_params%rel_tol)&
!                 then !> we check relative tolerance
!                    print *, "Newton solution not accurate enough"
!                    exit
!                else
!                    call this%update_conc_aq_prim_species(Delta_c1) !> updates c1 and Delta_c1
!                    c2nc_old=c2nc_new
!                end if
!            end if
!        end do
!    !> We update aqueous chemistry attributes
!        call this%compute_pH()
!!> Post-process
!    deallocate(fk,dfk_dc1,Delta_c1,drk_dc,u_tilde)
!end subroutine