!> \file gas_chemistry_m.f90
!> \brief Module for gas-phase chemistry, including partial pressures, activities, and equilibrium/kinetic routines.
!>
!> This module defines the `gas_chemistry_c` type and its associated procedures for managing gas-phase chemical processes in reactive transport simulations. It supports partial pressure calculation, activity coefficients, equilibrium and kinetic routines, and coupling to reactive zones. The module provides memory management, property calculation, iterative solvers, and assignment routines for gas-phase chemical state.
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2024

!> This class contains the partial pressures (activities) of gases
module gas_chemistry_m
    use local_chemistry_m, only: local_chemistry_c
    use reactive_zone_m, only: reactive_zone_c
    implicit none
    save
    private
    !public :: gas_chemistry_c
    !************************************************************************************************************************************************************************************************!
    type, public, extends(local_chemistry_c) :: gas_chemistry_c !> gas chemistry subclass
        class(reactive_zone_c), pointer :: reactive_zone !> reactive zone containing gas phase
        integer(kind=4), allocatable :: ind_gases_eq_cst_act(:) !> indices of gases in equilibrium with constant activity in chemical system gas phase
        integer(kind=4) :: num_gases_eq_act_fct_time=0 !> number of gases in equilibrium with activity as a function of time in chemical system gas phase
        integer(kind=4), allocatable :: ind_gases_eq_act_fct_time(:) !> indices of gases in equilibrium with activity as a function of time in chemical system gas phase
    contains
    !> Set
        procedure :: set_reactive_zone
        procedure :: set_indices_gases
        procedure :: set_concentrations=>set_conc_gases
    !> Allocate
        procedure :: allocate_partial_pressures
        procedure :: allocate_conc_gases
        procedure :: allocate_log_act_coeffs_gases
        procedure :: allocate_ind_gases_eq_cst_act
        procedure :: allocate_ind_gases_eq_act_fct_time
        procedure :: allocate_reaction_rates_gas_species_chem
    !> Compute
        procedure :: compute_log_act_coeffs_gases
        procedure :: compute_log_act_coeffs_gas
        procedure :: compute_partial_pressures
        procedure :: compute_conc_gases_ideal
        procedure :: compute_conc_gases_iter
        procedure :: compute_vol_gas_species_conc
        procedure :: compute_vol_gas_act_coeffs
        procedure :: compute_pressure
    !> Update
        procedure :: update_conc_gases
    !> Assign
        procedure :: copy_gas_chemistry
    end type
    
    interface
        !subroutine read_gas_species_chem_init(this,filename,gas_zones,line,num_tar)
        !>    import gas_zone_c
        !>    import gas_chemistry_c
        !>    implicit none
        !>    class(gas_chemistry_c) :: this
        !>    character(len=*), intent(in) :: filename
        !>    class(gas_zone_c), intent(in) :: gas_zones(:)
        !>    integer(kind=4), intent(inout) :: line
        !>    integer(kind=4), intent(out) :: num_tar
        !end subroutine
        
        !subroutine compute_conc_gases_iter(this,r_vec,Delta_t)
        !    import gas_chemistry_c
        !    implicit none
        !    class(gas_chemistry_c) :: this
        !    real(kind=8), intent(in) :: r_vec(:) !> reaction rates
        !    !real(kind=8), intent(in) :: porosity
        !    real(kind=8), intent(in) :: Delta_t !> time step
        !end subroutine
        
        
    end interface
    
    
    
    contains
        
        !subroutine set_partial_pressures(this,partial_pressures)
        !>    implicit none
        !>    class(gas_chemistry_c) :: this
        !>    real(kind=8), intent(in) :: partial_pressures(:)
        !>    !if (this%gas_zone%num_species/=size(partial_pressures)) error stop "Dimension error in set_partial_pressures"
        !>    this%activities=partial_pressures
        !end subroutine
        
        subroutine allocate_partial_pressures(this) !< units are atm
            implicit none
            class(gas_chemistry_c) :: this
!> \brief Allocate partial pressures array for gas species.
!> \details Allocates memory for storing partial pressures (activities) of all gas species in the reactive zone. Units are atm.
!> \param this Gas chemistry object.
            allocate(this%activities(this%reactive_zone%gas_phase%num_species))
        end subroutine
        
        subroutine allocate_conc_gases(this) !> units are moles
            implicit none
            class(gas_chemistry_c) :: this
            if (allocated(this%concentrations)) then
                deallocate(this%concentrations)
            end if
!> \brief Allocate concentrations array for gas species.
!> \details Allocates and initializes the concentrations array for all gas species in the reactive zone. Units are moles.
!> \param this Gas chemistry object.
            allocate(this%concentrations(this%reactive_zone%gas_phase%num_species))
            this%concentrations=0d0
            this%conc_old=this%concentrations
            this%conc_old_old=this%conc_old
        end subroutine

        subroutine allocate_log_act_coeffs_gases(this) !> <units are atm>
            implicit none
            class(gas_chemistry_c) :: this
            allocate(this%log_act_coeffs(this%reactive_zone%gas_phase%num_species),&
                this%log_Jacobian_act_coeffs(this%reactive_zone%gas_phase%num_species,this%reactive_zone%gas_phase%num_species))
            this%log_act_coeffs=0d0
            this%log_Jacobian_act_coeffs=0d0
!> \brief Allocate log activity coefficients and Jacobian matrix for gases.
!> \details Allocates and initializes arrays for log activity coefficients and their derivatives for all gas species. Units are atm.
!> \param this Gas chemistry object.
        end subroutine
        !
        !subroutine set_gas_zone(this,gas_zone)
        !>    implicit none
        !>    class(gas_chemistry_c) :: this
        !>    class(gas_zone_c), intent(in), target :: gas_zone
        !>    this%gas_zone=>gas_zone
        !end subroutine
        
        subroutine set_reactive_zone(this,reactive_zone)
            implicit none
            class(gas_chemistry_c) :: this
            class(reactive_zone_c), intent(in), target :: reactive_zone
            if (.not. associated(reactive_zone%chem_syst)) then
                error stop "Reactive zone has no chemical system assigned"
            else 
                this%reactive_zone=>reactive_zone
            end if
!> \brief Associate gas chemistry object with a reactive zone.
!> \details Sets pointer to the reactive zone, validating the chemical system association.
!> \param this Gas chemistry object.
!> \param reactive_zone Reactive zone object to associate.
        end subroutine
        
        subroutine set_indices_gases(this)
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i,j,k 
            j=0
            k=0
            do i=1,this%reactive_zone%gas_phase%num_species
                if (.not. this%reactive_zone%gas_phase%gases(i)%cst_act_flag) then
                    j=j+1
                    this%var_act_species_indices(j)=i
                else
                    k=k+1
                    this%cst_act_species_indices(k)=i
                end if
            end do
!> \brief Set indices for gas species classification (variable vs constant activity).
!> \details Classifies gas species as variable or constant activity and stores their indices for later use.
!> \param this Gas chemistry object.
        end subroutine
        
       subroutine compute_log_act_coeffs_gases(this)
!> \brief Compute log activity coefficients for gas species.
!> \details Calculates log activity coefficients for all gas species using the ideal gas law and current temperature/volume.
!> \param this Gas chemistry object.
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            do i=1,this%reactive_zone%gas_phase%num_species
                this%log_act_coeffs(i)=log10(R*this%temp/this%volume)
            end do
       end subroutine

       subroutine compute_log_act_coeffs_gas(this,ind_gas)
!> \brief Compute log activity coefficient for a single gas.
!> \details Calculates log activity coefficient for a single gas using the ideal gas law and current temperature/volume.
!> \param this Gas chemistry object.
            implicit none
            class(gas_chemistry_c) :: this !> gas chemistry object
            integer(kind=4), intent(in) :: ind_gas !> index of the gas species
            
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            
            this%log_act_coeffs(ind_gas)=log10(R*this%temp/this%volume)
       end subroutine
       
       subroutine compute_partial_pressures(this)
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            do i=1,this%reactive_zone%gas_phase%num_species
                this%activities(i)=this%concentrations(i)*R*this%temp/this%volume
            end do
!> \brief Compute partial pressures for gas species from concentrations.
!> \details Calculates partial pressures (activities) for all gas species using the ideal gas law and current temperature/volume.
!> \param this Gas chemistry object.
       end subroutine
       
       subroutine compute_conc_gases_ideal(this) !> ideal gas equation
       !> Concentrations are expressed in moles
            implicit none
            class(gas_chemistry_c) :: this
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            do i=1,this%reactive_zone%gas_phase%num_species
                this%concentrations(i)=this%activities(i)*this%volume/(this%temp*R)
            end do
!> \brief Compute concentrations for gas species using ideal gas law.
!> \details Calculates concentrations for all gas species from their activities using the ideal gas law.
!> \param this Gas chemistry object.
       end subroutine
       
       subroutine compute_conc_gases_iter(this,Delta_t,wat_vol,r_aq) !> gas conservation equation
       !> Concentrations are expressed in moles
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(in) :: porosity !> porosity
            real(kind=8), intent(in) :: wat_vol !> water volume
            real(kind=8), intent(in), optional :: r_aq(:)
            integer(kind=4) :: i
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            
            if (present(r_aq)) then
                if (this%reactive_zone%gas_phase%num_gases_kin>0) then
                    do i=1,this%reactive_zone%gas_phase%num_cst_act_species
                        this%concentrations(this%cst_act_species_indices(i))=this%concentrations(this%cst_act_species_indices(i))+&
                            Delta_t*dot_product(this%reactive_zone%chem_syst%stoich_mat_gas(:,this%cst_act_species_indices(i)),&
                            r_aq)*wat_vol
                    end do
                    if (this%concentrations(this%cst_act_species_indices(i))<0d0) then
                        print *, "Negative concentration in gas chemistry"
                        print *, "Gas: ", this%reactive_zone%gas_phase%gases(this%cst_act_species_indices(i))%name
                    end if
                else if (this%reactive_zone%gas_phase%num_gases_eq>0) then
                    do i=1,this%reactive_zone%gas_phase%num_cst_act_species
                        !this%concentrations(this%cst_act_species_indices(i))=this%concentrations(this%cst_act_species_indices(i))+&
                        !    Delta_t*dot_product(this%reactive_zone%chem_syst%stoich_mat_gas(:,this%cst_act_species_indices(i)),&
                        !    [this%re_mean(i),r_aq])*wat_vol
                        this%concentrations(this%cst_act_species_indices(i))=this%concentrations(this%cst_act_species_indices(i))+&
                            Delta_t*this%re_mean(i)*wat_vol
                        if (this%concentrations(this%cst_act_species_indices(i))<0d0) then
                            print *, "Negative concentration in gas chemistry"
                            print *, "Gas: ", this%reactive_zone%gas_phase%gases(this%cst_act_species_indices(i))%name
                        end if
                    end do
                end if
            end if
!> \brief Compute concentrations for gas species using conservation equation and optional reaction rates.
!> \details Updates concentrations for gas species based on time step, water volume, and optional reaction rates. Handles both kinetic and equilibrium gases.
!> \param this Gas chemistry object.
!> \param Delta_t Time step size.
!> \param wat_vol Water volume [L].
!> \param r_aq Optional array of reaction rates
       end subroutine
       
       subroutine compute_vol_gas_species_conc(this)
       !> Computes total volume of gas from concentrations
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            this%volume=sum(this%concentrations)*R*this%temp/this%pressure !> we choose first gas species in gas phase arbitrarily
!> \brief Compute total volume of gas from concentrations.
!> \details Calculates the total gas volume from the sum of concentrations using the ideal gas law.
!> \param this Gas chemistry object.
       end subroutine
       
       subroutine compute_vol_gas_act_coeffs(this)
       !> Computes total volume of gas from activity coefficients
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), parameter :: R=0.08205746 !> [atm*L/mol*K]
            this%volume=R*this%temp/(10**this%log_act_coeffs(1)) !> we choose first gas species in gas phase arbitrarily
!> \brief Compute total volume of gas from activity coefficients.
!> \details Calculates the total gas volume from activity coefficients using the ideal gas law.
!> \param this Gas chemistry object.
       end subroutine
       
       subroutine compute_pressure(this)
       !> Computes total pressure of gas
            implicit none
            class(gas_chemistry_c) :: this
            this%pressure=sum(this%activities) !> we sum the partial pressures
!> \brief Compute total pressure of gas phase.
!> \details Sums the partial pressures of all gas species to obtain the total pressure.
!> \param this Gas chemistry object.
       end subroutine
       
       subroutine update_conc_gases(this,conc_gases)
       !> Updates concentration of gases
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_gases(:)
            this%concentrations=conc_gases
!> \brief Update concentrations of gases.
!> \details Updates the concentrations array for all gas species.
!> \param this Gas chemistry object.
!> \param conc_gases Array of new concentrations [mol].
       end subroutine
       
       subroutine set_conc_gases(this,conc)
       !> Sets concentration of gases
            implicit none
            class(gas_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
            if (size(conc)/=this%reactive_zone%gas_phase%num_species) error stop "Dimension error in set_conc_gases"
            this%concentrations=conc
!> \brief Set concentrations of gases with validation.
!> \details Sets the concentrations array for all gas species, validating the input size.
!> \param this Gas chemistry object.
!> \param conc Array of concentrations to set [mol].
       end subroutine
       
        !subroutine initialise_log_act_coeffs_gas_species_chem(this)
        !    implicit none
        !    class(gas_chemistry_c) :: this
        !    allocate(this%log_act_coeffs(this%reactive_zone%gas_phase%num_species))
        !    this%log_act_coeffs=0d0
        !end subroutine
        !
        !subroutine initialise_log_Jacobian_act_coeffs_gas_species_chem(this)
        !    implicit none
        !    class(gas_chemistry_c) :: this
        !    allocate(this%log_Jacobian_act_coeffs(this%reactive_zone%gas_phase%num_species,this%reactive_zone%gas_phase%num_species))
        !    this%log_Jacobian_act_coeffs=0d0
        !end subroutine
       
       subroutine allocate_ind_gases_eq_cst_act(this)
           class(gas_chemistry_c) :: this
!> \brief Allocate indices for gases in equilibrium with constant activity.
!> \details Allocates the index array for gases in equilibrium with constant activity.
!> \param this Gas chemistry object.
           allocate(this%ind_gases_eq_cst_act(this%reactive_zone%gas_phase%num_gases_eq_cst_act))
       end subroutine
       
       subroutine allocate_ind_gases_eq_act_fct_time(this,num_gases_eq_act_fct_time)
       class(gas_chemistry_c) :: this
       integer(kind=4), intent(in), optional :: num_gases_eq_act_fct_time
!> \brief Allocate indices for gases in equilibrium with activity as a function of time.
!> \details Allocates the index array for gases in equilibrium with activity as a function of time, with optional input for number of such gases.
!> \param this Gas chemistry object.
!> \param num_gases_eq_act_fct_time Optional number of gases with time-dependent activity.
       if (present(num_gases_eq_act_fct_time)) then
           if (this%num_gases_eq_act_fct_time<0) then
               error stop "Dimension error in allocate_ind_gases_eq_act_fct_time"
           else 
               this%num_gases_eq_act_fct_time=num_gases_eq_act_fct_time
            end if
       end if
       allocate(this%ind_gases_eq_act_fct_time(this%num_gases_eq_act_fct_time))
       end subroutine
       
       subroutine allocate_reaction_rates_gas_species_chem(this)
           class(gas_chemistry_c) :: this
           if (allocated(this%re_mean) .and. allocated(this%Re)) then
               deallocate(this%re_mean,this%Re)
           end if
!> \brief Allocate reaction rates arrays for gas-phase reactions.
!> \details Allocates and initializes arrays for mean and instantaneous reaction rates for equilibrium gases.
!> \param this Gas chemistry object.
           allocate(this%re_mean(this%reactive_zone%gas_phase%num_gases_eq))
           allocate(this%Re(this%reactive_zone%gas_phase%num_gases_eq))
           this%re_mean=0d0 !> by default
           this%Re=0d0 !> by default
       end subroutine

       subroutine copy_gas_chemistry(this,other)
           implicit none
           class(gas_chemistry_c) :: this
           class(gas_chemistry_c), intent(in) :: other
!> \brief Deep copy assignment of gas chemistry object.
!> \details Copies all properties and associations from another gas chemistry object, including reactive zone and index arrays.
!> \param this Gas chemistry object.
!> \param other Source gas chemistry object to copy from.
           !call this%set_id(other%id)
        !    this%temp=other%temp
        !    this%pressure=other%pressure
        !    this%volume=other%volume
           call this%copy_local_chemistry(other)
           if (associated(other%reactive_zone)) then
               call this%set_reactive_zone(other%reactive_zone)
           end if
           this%num_gases_eq_act_fct_time=other%num_gases_eq_act_fct_time
           if (allocated(other%ind_gases_eq_cst_act)) then
               this%ind_gases_eq_cst_act=other%ind_gases_eq_cst_act
           end if
           if (allocated(other%ind_gases_eq_act_fct_time)) then
               this%ind_gases_eq_act_fct_time=other%ind_gases_eq_act_fct_time
           end if
       end subroutine
end module