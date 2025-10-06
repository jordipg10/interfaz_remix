!!> Computes r_eq(u_j)
!subroutine compute_r_eq_aq_chem_De_Simoni(this,grad_conc_comp,dispersion,porosity)
!    use aqueous_chemistry_m
!    use BCs_m
!    implicit none
!    class(aqueous_chemistry_c) :: this
!    real(kind=8), intent(in) :: grad_conc_comp(:) ! du_j/dx
!    real(kind=8), intent(in) :: dispersion
!    real(kind=8), intent(in) :: porosity
!    !real(kind=8), intent(out) :: r_eq(:) ! (tiene que estar alocatado)
!            
!    integer(kind=4) :: i,n,niter,rk_ind,l
!    
!    !allocate(grad_conc_comp(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
!        
!    !call this%compute_grad_conc_comp(Delta_x,BCs,grad_conc_comp)
!    
!    this%r_eq=0d0
!    
!    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
!        if (this%solid_chemistry%reactive_zone%eq_reactions(i)%react_type==2 .and. this%solid_chemistry%reactive_zone%eq_reactions(i)%num_species==3) then
!            call this%solid_chemistry%reactive_zone%eq_reactions(i)%react_rate_bin_syst_eq_1D(this%conc_comp(i),grad_conc_comp(i),dispersion,porosity,this%r_eq(i))
!        end if
!    end do
!end subroutine