!> \file compute_mass_bal_mins.f90
!> \brief Updates mineral-to-pore volume ratios using mass balance after one time step of Euler explicit method
!> \details This subroutine computes the updated mineral-to-pore volume ratios based on mass balance
!> considerations after one iteration of the explicit Euler method for reactive transport. The update accounts for:
!> - Kinetic minerals: minerals undergoing dissolution/precipitation via kinetic reactions
!> - Equilibrium minerals: minerals controlled by equilibrium reactions
!>
!> The update is based on the following mass balance:
!> \f[
!> \psi_i^{n+1} = \psi_i^n + \Delta t \cdot V_{\text{mol},i} \cdot r_{i}
!> \f]
!> where:
!> - \f$\psi_i = V_{\text{min},i}/V_{\text{pore}}\f$ is the mineral-to-pore volume ratio of mineral i [-]; it is
!>   non-negative but NOT bounded by 1 (mineral volume may exceed pore volume)
!> - \f$V_{\text{mol},i}\f$ is the molar volume of mineral i [L³/M]
!> - \f$r_{i}\f$ is the reaction rate per unit pore volume (kinetic or equilibrium) [M/(L³·T)]
!> - \f$\Delta t\f$ is the time step [T]
!>
!> The reaction mixing ratio \f$\lambda_r\f$ is NOT applied here: `rxn_rate_eq` already holds the local
!> rate \f$R_e/\Delta t\f$ with \f$R_e=\rho_e/\lambda_r\f$ (see compute_Re), so the mineral in the cell
!> gains/loses the full reaction amount; \f$\lambda_r\f$ only scales its footprint on the mixed aqueous phase.
!>
!> The positive reaction rate indicates precipitation (increasing ratio), while negative
!> rate indicates dissolution (decreasing ratio).
!>
!> \param[in,out] this Solid chemistry object containing mineral-to-pore volume ratios, reaction rates, and mineral properties
!> \param[in] Delta_t Time step for the explicit Euler update [T]

subroutine compute_mass_bal_mins(this,Delta_t)
    use solid_chemistry_m, only: solid_chemistry_c !> Import solid chemistry class containing mineral data and reaction rates
    implicit none !> Enforce explicit variable declarations
!> Arguments
    class(solid_chemistry_c) :: this !> Solid chemistry object containing mineral zone, reactive zone, volumetric fractions, and reaction rates [-]
    real(kind=8), intent(in) :: Delta_t !> Time step for explicit Euler update [T]
    !real(kind=8), intent(in) :: lambda_r !> reaction mixing ratio; scales the eq amount to the aqueous contribution lambda_r*Re
!> Variables
    integer(kind=4) :: i !> Loop counter for minerals (kinetic and equilibrium) [-]
!> Process
!> Update volumetric fractions for kinetic minerals (dissolution/precipitation governed by kinetic rate laws)
    do i=1,this%mineral_zone%num_minerals_kin !> Loop over all kinetic minerals in mineral zone
        this%min_pore_vol_ratios(i)=this%min_pore_vol_ratios(i)+Delta_t*this%mineral_zone%chem_syst%minerals(& !> Update volumetric fraction: φᵢ = φᵢ + Δt·Vₘₒₗ·rₖ
            this%mineral_zone%ind_min_chem_syst(i))%mineral%mol_vol*this%rxn_rate_kin_avg(i) !> Multiply time step by molar volume [L³/M] and mean kinetic reaction rate [M/(L³·T)]
        if (this%min_pore_vol_ratios(i)<0d0) then !> Ensure volumetric fraction does not become negative
            this%min_pore_vol_ratios(i)=0d0 !> Set volumetric fraction to zero if negative
            write(0,*) "Mineral-to-pore volume ratio of kinetic mineral became negative: ", &
                trim(this%mineral_zone%chem_syst%minerals(this%mineral_zone%ind_min_chem_syst(i))%mineral%name)
            error stop "Mineral-to-pore volume ratio of kinetic mineral became negative"
        end if
    end do !> End loop over kinetic minerals
!> Update volumetric fractions for equilibrium minerals (minerals controlled by equilibrium reactions)
    do i=1,this%reactive_zone%num_minerals !> Loop over all equilibrium minerals in reactive zone
        this%min_pore_vol_ratios(this%mineral_zone%num_minerals_kin+i)=this%min_pore_vol_ratios(this%mineral_zone%num_minerals_kin+i)+& !> Update volumetric fraction at offset index (after kinetic minerals): φᵢ = φᵢ + Δt·Vₘₒₗ·rₑ
            Delta_t*this%reactive_zone%chem_syst%minerals(this%reactive_zone%ind_mins_chem_syst(i))%mineral%mol_vol*&
            this%rxn_rate_eq(i) !> molar volume [L³/M] × eq reaction rate [M/(L³·T)] × time step × reaction mixing ratio lr
        if (this%min_pore_vol_ratios(this%mineral_zone%num_minerals_kin+i)<0d0) then !> Ensure volumetric fraction does not become negative
            this%min_pore_vol_ratios(this%mineral_zone%num_minerals_kin+i)=0d0 !> Set volumetric fraction to zero if negative
            write(0,*) "Mineral-to-pore volume ratio of equilibrium mineral became negative: ", &
                trim(this%reactive_zone%chem_syst%minerals(this%reactive_zone%ind_mins_chem_syst(i))%mineral%name)
            error stop "Mineral-to-pore volume ratio of equilibrium mineral became negative"
        end if
    end do !> End loop over equilibrium minerals
end subroutine