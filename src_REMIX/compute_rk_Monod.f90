!> Computes Monod reaction rate
subroutine compute_rk_Monod(this,conc,rk)
    use redox_kin_reaction_m, only: redox_kin_c
    implicit none
    class(redox_kin_c), intent(in) :: this
    real(kind=8), intent(in) :: conc(:) !> conc=[conc_inh,conc_acceptor,conc_donor]
    real(kind=8), intent(out) :: rk !> reaction rate
    
    integer(kind=4) :: n,m,j,k,l
    real(kind=8) :: prod_cat,prod_inh
    real(kind=8), allocatable :: conc_inh(:),conc_M(:)
    
    prod_cat=1d0 !> product catalyser terms
    prod_inh=1d0 !> product inhibitor terms
    
    rk=this%params%rate_cst !> rate constant
    do j=1,this%params%n_inh
        prod_inh=prod_inh*this%params%k_inh(j)/(this%params%k_inh(j)+conc(j))
    end do
    rk=rk*prod_inh !> inhibition factors
    do j=1,2
        prod_cat=prod_cat*conc(this%params%n_inh+j)/(this%params%k_M(j)+conc(this%params%n_inh+j))
    end do
    rk=rk*prod_cat !> electron acceptor & donor
    !rk=rk*(1d0-conc(this%params%n_t+1)/this%params%cb_max) !> logistic factor
end subroutine