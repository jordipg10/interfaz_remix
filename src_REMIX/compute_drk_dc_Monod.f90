!> Computes Monod kinetic reaction rate and its gradient
!!> drk_dc=rk*dlogT/dc
!!>        =rk*k/(c*(k+c)) if electron acceptor/donor
!!>        =-rk/(k+c) if inhibitor
subroutine compute_drk_dc_Monod(this,conc,rk,drk_dc)
    use redox_kin_reaction_m, only: redox_kin_c
    implicit none
    class(redox_kin_c), intent(in) :: this !> redox kinetic reaction
    real(kind=8), intent(in) :: conc(:) !> conc=[conc_inh,conc_acceptor,conc_donor]
    real(kind=8), intent(in) :: rk !> kinetic reaction rate
    real(kind=8), intent(out) :: drk_dc(:) !> gradient of rk (must be already allocated)
    
    integer(kind=4) :: i,n,m,j,k,l
    real(kind=8) :: prod_inh
    real(kind=8), allocatable :: conc_inh(:),conc_M(:),prod_grad(:)
    
    !prod_cat=1d0 !> product catalyser terms
    !prod_inh=1d0 !> product inhibitor terms
    allocate(prod_grad(2)) !> rk/c_M (chapuza para el gradiente)
    
    !rk=this%params%rate_cst !> rate constant
    !prod_grad=this%params%rate_cst !> rate constant
!> inhibition factors
    !do j=1,this%params%n_inh
    !    rk=rk*this%params%k_inh(j)/(this%params%k_inh(j)+conc(j))
    !end do
    !prod_grad=rk
    !rk=rk*prod_inh 
!> electron acceptor & donor
    !do j=1,2
    !    rk=rk*conc(this%params%n_inh+j)/(this%params%k_M(j)+conc(this%params%n_inh+j))
    !end do
!> Chapuza
    prod_grad(1)=prod_grad(1)*conc(this%params%n_inh+2)/((this%params%k_M(1)+conc(this%params%n_inh+1))*(this%params%k_M(2)+&
        conc(this%params%n_inh+2)))
    prod_grad(2)=prod_grad(2)*conc(this%params%n_inh+1)/((this%params%k_M(1)+conc(this%params%n_inh+1))*(this%params%k_M(2)+&
        conc(this%params%n_inh+2)))
    !rk=rk*prod_cat
    !rk=rk*(1d0-conc(this%params%n_t+1)/this%params%cb_max) !> logistic factor
!> Gradient
    do j=1,this%params%n_inh
        drk_dc(j)=-rk/(this%params%k_inh(j)+conc(j))
    end do
    do j=1,2
        drk_dc(this%params%n_inh+j)=prod_grad(j)*this%params%k_M(j)/(this%params%k_M(j)+conc(this%params%n_inh+j))
        !drk_dc(this%params%n_inh+j)=rk*this%params%k_M(j)/(this%params%k_M(j)+conc(this%params%n_inh+j))**2
    end do
end subroutine