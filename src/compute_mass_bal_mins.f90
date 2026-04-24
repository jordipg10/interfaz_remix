!> \file compute_mass_bal_mins.f90
!> \brief Updates mineral volumetric fractions using mass balance after one time step of Euler explicit method
!> \details This subroutine computes the updated volumetric fractions of minerals based on mass balance
!> considerations after one iteration of the explicit Euler method for reactive transport. The volumetric
!> fraction update accounts for:
!> - Kinetic minerals: minerals undergoing dissolution/precipitation via kinetic reactions
!> - Equilibrium minerals: minerals controlled by equilibrium reactions
!>
!> The volumetric fraction update is based on the following mass balance:
!> \f[
!> \phi_i^{n+1} = \phi_i^n + \Delta t \cdot V_{\text{mol},i} \cdot r_{k,i}
!> \f]
!> where:
!> - \f$\phi_i\f$ is the volumetric fraction of mineral i [-]
!> - \f$V_{\text{mol},i}\f$ is the molar volume of mineral i [L³/M]
!> - \f$r_{k,i}\f$ is the reaction rate (kinetic or equilibrium) [M/(L³·T)]
!> - \f$\Delta t\f$ is the time step [T]
!>
!> The positive reaction rate indicates precipitation (increasing volumetric fraction), while negative
!> rate indicates dissolution (decreasing volumetric fraction).
!>
!> \param[in,out] this Solid chemistry object containing mineral volumetric fractions, reaction rates, and mineral properties
!> \param[in] Delta_t Time step for the explicit Euler update [T]

subroutine compute_mass_bal_mins(this,Delta_t) 
    use solid_chemistry_m, only: solid_chemistry_c !> Import solid chemistry class containing mineral data and reaction rates
    implicit none !> Enforce explicit variable declarations
!> Arguments
    class(solid_chemistry_c) :: this !> Solid chemistry object containing mineral zone, reactive zone, volumetric fractions, and reaction rates [-]
    real(kind=8), intent(in) :: Delta_t !> Time step for explicit Euler update [T]
!> Variables
    integer(kind=4) :: i !> Loop counter for minerals (kinetic and equilibrium) [-]
!> Process
!> Update volumetric fractions for kinetic minerals (dissolution/precipitation governed by kinetic rate laws)
    do i=1,this%mineral_zone%num_minerals_kin !> Loop over all kinetic minerals in mineral zone
        this%vol_fracts_mins(i)=this%vol_fracts_mins(i)+Delta_t*this%mineral_zone%chem_syst%minerals(& !> Update volumetric fraction: φᵢ = φᵢ + Δt·Vₘₒₗ·rₖ
            this%mineral_zone%ind_min_chem_syst(i))%mineral%mol_vol*this%rk_mean(i) !> Multiply time step by molar volume [L³/M] and mean kinetic reaction rate [M/(L³·T)]
        if (this%vol_fracts_mins(i)<0d0) then !> Ensure volumetric fraction does not become negative
            this%vol_fracts_mins(i)=0d0 !> Set volumetric fraction to zero if negative
        end if
    end do !> End loop over kinetic minerals
!> Update volumetric fractions for equilibrium minerals (minerals controlled by equilibrium reactions)
    do i=1,this%reactive_zone%num_minerals !> Loop over all equilibrium minerals in reactive zone
        this%vol_fracts_mins(this%mineral_zone%num_minerals_kin+i)=this%vol_fracts_mins(this%mineral_zone%num_minerals_kin+i)+& !> Update volumetric fraction at offset index (after kinetic minerals): φᵢ = φᵢ + Δt·Vₘₒₗ·rₑ
            Delta_t*this%reactive_zone%chem_syst%minerals(this%reactive_zone%ind_mins_chem_syst(i))%mineral%mol_vol*this%re_mean(i) !> Multiply time step by molar volume [L³/M] and mean equilibrium reaction rate [M/(L³·T)]
        if (this%vol_fracts_mins(this%mineral_zone%num_minerals_kin+i)<0d0) then !> Ensure volumetric fraction does not become negative
            this%vol_fracts_mins(this%mineral_zone%num_minerals_kin+i)=0d0 !> Set volumetric fraction to zero if negative
        end if
    end do !> End loop over equilibrium minerals
end subroutine