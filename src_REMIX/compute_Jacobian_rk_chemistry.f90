!! Computes Jacobian of rk(c_j)
!subroutine compute_Jacobian_rk_chemistry(this,drk_dc)
!    use chemistry_m
!    implicit none
!    class(chemistry_c), intent(in) :: this
!    !real(kind=8), intent(in) :: conc(:) ! c_j or [c1,c2]_j
!    !real(kind=8), intent(in) :: rk(:) ! rk_j
!    real(kind=8), intent(out) :: drk_dc(:,:) ! Jacobian (must be allocated)
!    
!    integer(kind=4) :: i,k,n,m,l,rk_ind,cat_ind,p,inh_ind,DOC_ind
!    real(kind=8), allocatable :: conc(:),conc_rk_min
!    
!    !if (this%solid_chemistry%reactive_zone%num_eq_reactions
!    
!    drk_dc=0d0
!    !do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
!    !    call this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%compute_drk_dc(conc,rk(i),drk_dc(i,:))
!    !end do
!    do i=1,this%num_loc_min_kin_reacts
!        !call get species_ind
!        ! call get_act_cat
!        call this%loc_min_kin_reacts(i)%compute_drk_dc_min(drk_dc(i,:))
!    end do
!    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
!        !call this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%compute_drk_dc(conc,rk(i),drk_dc(i,:))
!    end do
!end subroutine