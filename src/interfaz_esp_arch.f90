!> This subroutine is the interface to solve a reactive mixing iteration for aqueous species concentrations
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> For each target water it solves the kinetic system implicitly with Newton-Raphson
!> (theta=1 backward Euler, lumped kinetic mixing ratio=1), falling back to explicit Euler
!> for any cell whose implicit Newton solve does not converge. The interface performs a single
!> step from an initial state whose kinetic reaction rates are zero (rk_old=0), hence the fully
!> implicit weighting; for the same reason the explicit fallback carries no kinetic reaction
!> It assumes that there are no equilibrium reactions
!> It reads the concentrations after conservative transport in a file, and writes the concentrations
!> and the end-of-step kinetic reaction rates after reactive mixing in a different file
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
    integer(kind=4) :: i,j,tw_idx !> loop variables / shorthand for target water index
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_p !> number of primary species (= variable activity species, kinetic-only)
    integer(kind=4) :: n_aq_kin,n_min_kin,n_gas_kin !> per-water kinetic reaction counts
    integer(kind=4) :: n_kin_tot !> total kinetic rate rows written per water
    logical :: CV_flag !> TRUE if Newton converges, FALSE otherwise
    real(kind=8), allocatable :: c_tilde(:,:) !> concentration aqueous species after solving conservative transport
    real(kind=8), allocatable :: c_react(:) !> reaction contribution to aqueous species concentrations (explicit fallback)
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
        n_p      = tw0%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
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
        !> --------------------------------------------------------------------
        !> Step 1: fully implicit reactive mixing (Newton-Raphson, theta=1 backward
        !> Euler, lumped ratio=1). For a kinetic-only system the variable-activity
        !> species are all primary, so c_hat = c_tilde directly (no speciation). If the
        !> implicit Newton does not converge (stiff cell), fall back to the explicit (EE)
        !> update for that cell so the run always completes.
        !> --------------------------------------------------------------------
        block
            real(kind=8), allocatable :: c_hat(:)  !< mixed primary species (c_tilde + non-aqueous primaries from current state)
            real(kind=8), allocatable :: conc_p(:) !< primary species concentrations: initial guess on entry, solution on exit
            real(kind=8), parameter :: theta=1d0   !< reaction time weighting: fully implicit (no previous-step rates exist in a single-step call)
            allocate(c_hat(n_p),conc_p(n_p))
            !> Start from the current primary state (covers any non-aqueous primary species), then
            !> overwrite the aqueous part with the conservatively transported concentrations.
            c_hat = tw%get_c1()
            c_hat(1:num_aq_comps) = c_tilde(:,j)
            conc_p = c_hat !< initial guess for the implicit Newton solver
            call tw%Newton_EI_kin_anal_ideal_opt2(c_hat,1d0,Delta_t,theta,conc_p,niter,CV_flag)
            if (CV_flag) then
                c_new(:,j)=conc_p(1:num_aq_comps)
            else
                !> Explicit-Euler fallback for this stiff cell. With the zero previous-step rates of the
                !> initial state it reduces to c_new = c_tilde (no kinetic reaction in this cell), so warn.
                print *, "WARNING: target water ", tw_idx, ": implicit kinetic Newton solve did not converge; ", &
                    "explicit fallback with zero previous-step rates (no kinetic reaction in this step)"
                call tw%compute_react_term_EE_kin(Delta_t,1.0d0,c_react)
                c_new(:,j)=c_tilde(:,j)+c_react
                do i=1,num_aq_comps
                    if (c_new(i,j)<0d0) c_new(i,j)=0d0
                end do
                conc_p(1:num_aq_comps)=c_new(:,j)
                call tw%set_conc_prim_species(conc_p)
                call tw%set_act_aq_species()
            end if
            deallocate(c_hat,conc_p)
        end block
        !> --------------------------------------------------------------------
        !> Step 2: kinetic reaction rates evaluated at the end-of-step (converged)
        !> concentrations. Order: aqueous (linear+redox) | mineral | gas.
        !> --------------------------------------------------------------------
        if (n_aq_kin+n_min_kin > 0) then
            block
                real(kind=8), allocatable :: rk_new_cell(:) !< end-of-step kinetic rates for this cell
                allocate(rk_new_cell(n_aq_kin+n_min_kin))
                call tw%compute_rxn_rate_kin(rk_new_cell)
                rk_out(1:n_aq_kin+n_min_kin,j) = rk_new_cell
            end block
        end if
        !> Gas kinetic rates not yet implemented in gas_chemistry_c
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
