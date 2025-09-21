!> Computes mineral kinetic reaction rate gradient
subroutine compute_drk_dc_mineral(this,conc_sp,act_sp,log_act_coeffs_sp,act_cat,saturation,react_surf,temp,drk_dc)
    use kin_mineral_m, only: kin_mineral_c
    implicit none
    class(kin_mineral_c), intent(in) :: this !> mineral kinetic reaction object
    real(kind=8), intent(in) :: conc_sp(:) !> concentrations aqueous species reaction
    real(kind=8), intent(in) :: act_sp(:) !> activities aqueous species reaction
    real(kind=8), intent(in) :: log_act_coeffs_sp(:) !> log10 activity coefficients aqueous species reaction
    real(kind=8), intent(in) :: act_cat(:) !> activities catalysers
    real(kind=8), intent(in) :: saturation !> saturation
    real(kind=8), intent(in) :: react_surf !> reactive surface
    real(kind=8), intent(in) :: temp !> temperature [K]
    real(kind=8), intent(out) :: drk_dc(:) !> must be already allocated
    
    integer(kind=4) :: n_sp,i,l,j,k,zeta
    real(kind=8) :: rk_j,sum,prod,prod_sat
    real(kind=8), parameter :: R=8.31446261815324 !> J/(mol*K)
    
    drk_dc=0d0 !> initialize gradient
    n_sp=size(act_sp) !> number of species
    do l=1,n_sp !> loop over species
        sum=0d0 !> initialize sum
        prod_sat=(10**log_act_coeffs_sp(l))*this%eq_cst !> initialize saturation product
        do j=1,l-1 !> loop over species except l-th species
            prod_sat=prod_sat*act_sp(j) !> product of activities species
        end do
        do j=l+1,n_sp !> loop over species except l-th species
            prod_sat=prod_sat*act_sp(j) !> product of activities species
        end do
        do j=1,this%params%num_par_reacts !> loop over parallel reactions
            prod=1d0 !> initialize product
            do i=1,this%params%num_cat !> loop over catalysers
                prod=prod*act_cat(i)**this%params%p(j,i) !> product of activities catalysers
            end do
            sum=sum+prod*this%params%theta(j)*(conc_sp(l)**(this%params%theta(j)-1d0))*(prod_sat**this%params%theta(j))*&
                this%params%k(j)*this%params%eta(j)*((&
                saturation**this%params%theta(j))-1d0)**(this%params%eta(j)-1d0) !> sum of parallel reactions
        end do
        drk_dc(l)=sum*(-this%stoichiometry(l)) !> gradient of reaction rate
    end do
!> Flag for dissolution or precipitation
    if (saturation<1d0) then !> dissolution
        zeta=-1
    else if (saturation>1d0) then !> precipitation
        zeta=1
    end if
    drk_dc=drk_dc*react_surf*exp(-this%params%act_energy/(R*temp)) !> multiply by reactive surface and Arrhenius
end subroutine