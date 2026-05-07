!> This subroutine is the interface to solve a reactive mixing iteration for aqueous species concentrations
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It assumes that there are no equilibrium reactions
!> It reads the concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
subroutine interfaz_esp_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous species (= number of aqueous components)
    character(len=*), intent(in) :: file_in !> name of file containing aqueous species concentrations after solving conservative transport
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing aqueous species concentrations after solving reactive mixing
!> Variables
    integer(kind=4) :: i,j,k,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: n_aq_kin,n_min_kin,n_gas_kin !> per-water kinetic reaction counts
    integer(kind=4) :: n_kin_tot !> total kinetic rate rows written per water
    real(kind=8), allocatable :: c_tilde(:,:) !> concentration aqueous species after solving conservative transport
    real(kind=8), allocatable :: c_react(:) !> reaction contribution to aqueous species concentrations after solving reactive mixing
    real(kind=8), allocatable :: c_new(:,:) !> concentration aqueous species after solving reactive mixing
    real(kind=8), allocatable :: rk_out(:,:) !> kinetic reaction rates per target water [n_kin_tot, num_target_waters]
!> Pre-process
    allocate(c_tilde(num_aq_comps,this%num_target_waters),c_new(num_aq_comps,this%num_target_waters))
    allocate(c_react(num_aq_comps))
    !> Determine reaction-rate sizes from the first target water (assumed uniform across targets).
    associate(tw0 => this%waters(this%tar_wat_indices(1)))
        n_aq_kin = tw0%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts + &
                   tw0%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        n_min_kin= tw0%solid_chemistry%mineral_zone%num_minerals_kin
        n_gas_kin= tw0%solid_chemistry%reactive_zone%chem_syst%num_gas_kin_reacts
    end associate
    n_kin_tot = n_aq_kin + n_min_kin + n_gas_kin
    allocate(rk_out(max(n_kin_tot,1),this%num_target_waters))
    rk_out=0d0
!> Process
    !> We read the aqueous species concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (c_tilde(i,j), j=1,this%num_target_waters) !> we read one row (species) at a time
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters !> loop over target waters
        tw_idx = this%tar_wat_indices(j)
        associate(tw => this%waters(tw_idx))
        call tw%compute_react_term_EE_kin(Delta_t,1.0d0,c_react) !> chemical part of aqueous species concentrations
        c_new(:,j)=c_tilde(:,j)+c_react !> we sum transport and reaction parts
        !> Kinetic reaction rates effectively used in this iteration (Euler-explicit uses rk_old).
        !> Order: aqueous (linear+redox) | mineral | gas
        if (n_kin_tot > 0) then
            k=0
            do i=1,n_aq_kin
                rk_out(k+i,j) = tw%rk_old(i)
            end do
            k=k+n_aq_kin
            do i=1,n_min_kin
                rk_out(k+i,j) = tw%solid_chemistry%rk_old(i)
            end do
            !> Gas kinetic rates not yet implemented in gas_chemistry_c
        end if
        end associate
    end do
!> Post-process
    !> We write the aqueous species concentrations after solving reactive mixing
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Aqueous species concentrations after solving reactive mixing iteration (rows: species, columns: & 
        targets):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (c_new(i,j), j=1,this%num_target_waters)
    end do
    write(2,*)
    write(2,*) "Kinetic reaction rates (rows: reactions [aqueous (linear+redox) | minerals], & 
        columns: targets):"
    write(2,*)
    if (n_kin_tot > 0) then
        do i=1,n_kin_tot
            write(2,"(2x,*(ES15.5))") (rk_out(i,j), j=1,this%num_target_waters)
        end do
    else
        write(2,*) " (no kinetic reactions)"
    end if
    close(2)
    deallocate(c_tilde,c_react,c_new,rk_out)
end subroutine
