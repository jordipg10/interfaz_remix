!!> Performs Picard method in reactive mixing iteration for Euler fully implicit in chemical reactions
!!! The chemical system has equilibrium and kinetic reactions
!
!!> Computes concentration of variable activity species in target j at time step k+1
!    
!!> We assume all variable activity species are aqueous
!    
!subroutine Picard_EfI_rk_eq_kin_aq_chem(this,mixing_ratios,porosity,Delta_t,mixing_waters,niter,CV_flag)
!    use aqueous_chemistry_m
!    use metodos_sist_lin_m
!    use vectors_m
!    implicit none
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: mixing_ratios(:) !> the first element corresponds to argument "this", and the rest correspond to argument "mixing waters" (in the same order)
!    real(kind=8), intent(in) :: porosity !> at target j
!    real(kind=8), intent(in) :: Delta_t !> time step
!    class(aqueous_chemistry_c), intent(in) :: mixing_waters(:) !> including external and boundary waters
!    integer(kind=4), intent(out) :: niter !> number of iterations
!    logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
!            
!    real(kind=8), allocatable :: c1(:),c_nc_new(:),dc2nc_dc1(:,:),rk(:),Um_Skt_prod(:,:),fk(:),dfk_dc1(:,:),dfk_dc(:,:),um_tilde(:),Delta_c1(:),Delta_c(:),conc_kin(:),gamma_1(:),gamma_2(:),conc_Sk(:),c_tilde(:),conc(:)
!    integer(kind=4) :: i,n_p,n_p_aq,err,rk_ind,l,niter_spec
!    integer(kind=4), allocatable :: kin_ind(:)
!    
!    niter=0
!    CV_flag=.false.
!    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
!    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
!
!    allocate(Delta_c1(n_p_aq),c_tilde(n_p_aq))
!    
!    call this%compute_c_tilde(mixing_ratios,mixing_waters,c_tilde)
!    do 
!        niter=niter+1 !> we update number of iterations
!        if (niter>this%solid_chemistry%reactive_zone%CV_params%niter_max) then
!            print *, "Too many Picard iterations"
!            error stop
!        end if
!        call this%compute_c2nc_from_c1_aq_Picard(niter_spec,CV_flag) !> computes concentration of secondary non constant activity species from concentration of primary species using mass action law
!        call this%compute_rk_aq_chem()
!        c_nc_new=c_tilde+Delta_t*matmul(transpose(this%solid_chemistry%reactive_zone%chem_syst%Sk(:,1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)),this%rk)
!        Delta_c1=c_nc_new(1:n_p_aq)-this%concentrations(1:n_p_aq)
!        if (inf_norm_vec_real(Delta_c1/this%concentrations(1:n_p_aq))<this%solid_chemistry%reactive_zone%CV_params%rel_tol) then !> CV reached
!            CV_flag=.true.
!            exit
!        else
!            call this%update_conc_aq_prim_species(Delta_c1) !> updates c1 and Delta_c1
!        end if
!    end do
!end subroutine