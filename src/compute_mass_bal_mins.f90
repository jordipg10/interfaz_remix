!> Computes volumetric fractions of minerals after an iteration of Euler explicit method
subroutine compute_mass_bal_mins(this,Delta_t) 
    use solid_chemistry_m, only: solid_chemistry_c
    implicit none
!> Arguments
    class(solid_chemistry_c) :: this
    !real(kind=8), intent(in) :: r_e(:) !> equilibrium reaction rates of mineral reactions
    real(kind=8), intent(in) :: Delta_t !> time step
!> Variables
    integer(kind=4) :: i !> counter minerals
!> Process
    do i=1,this%mineral_zone%num_minerals_kin
        this%vol_fracts(i)=this%vol_fracts(i)+Delta_t*this%mineral_zone%chem_syst%minerals(&
            this%mineral_zone%ind_min_chem_syst(i))%mineral%mol_vol*this%rk(i)
    end do
    do i=1,this%reactive_zone%num_minerals
        this%vol_fracts(this%mineral_zone%num_minerals_kin+i)=this%vol_fracts(this%mineral_zone%num_minerals_kin+i)+&
            Delta_t*this%reactive_zone%chem_syst%minerals(this%reactive_zone%ind_mins_chem_syst(i))%mineral%mol_vol*this%r_eq(i)
    end do
end subroutine