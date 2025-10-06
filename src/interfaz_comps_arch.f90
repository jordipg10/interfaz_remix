!> This subroutine is the interface to solve a reactive mixing iteration for components
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It reads the concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
subroutine interfaz_comps_arch(this,path,num_aq_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_aq_comps !> number of aqueous components
    character(len=*), intent(in) :: file_in !> name of file containing component concentrations after solving conservative transport
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing component concentrations after solving conservative transport
!> Variables
    integer(kind=4) :: i,j !> loop variables
    integer(kind=4) :: niter !> number of iterations in Newton algorithm
    integer(kind=4) :: n_nc !> number of variable activity species
    logical :: CV_flag !> TRUE if converges, FALSE otherwise
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport
    real(kind=8), allocatable :: u_react(:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: u_new(:,:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: conc_nc(:,:) !> concentrations of variable activity species after solving reactive mixing
!> Pre-process
    allocate(u_tilde(num_aq_comps,this%num_target_waters_dom),u_new(num_aq_comps,this%num_target_waters_dom))
    allocate(u_react(num_aq_comps))
    n_nc=this%get_num_aq_var_act_species_dom() !> we get number of variable activity species in the domain
    allocate(conc_nc(n_nc,this%num_target_waters_dom)) !> we allocate conc_nc for all target waters in the domain
!> Process
    !> We read the component concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (u_tilde(i,j), j=1,this%num_target_waters_dom) !> we read one row (component) at a time
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters_dom !> loop over target waters in the domain
        call this%target_waters(this%dom_tar_wat_indices(j))%reaction_iteration_EE_eq_kin_lump(Delta_t,u_react) !> chemical part of components
        u_new(:,j)=u_tilde(:,j)+u_react !> we sum transport and reaction parts
        call this%target_waters(this%dom_tar_wat_indices(j))%compute_c_nc_from_u_aq_Newton_ideal(u_new(:,j),conc_nc(:,j),niter,&
            CV_flag) !> we compute concentrations of variable activity species from component concentrations
        if (.not.CV_flag) then
            print *, "Target water index: ", this%dom_tar_wat_indices(j)
            print *, "No convergence in speciation after reactive mixing iteration"
            print *, "Try reducing the time step"
            error stop
        end if
    end do
!> Post-process
    !> We write the component concentrations after solving reactive mixing
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Component concentrations after solving reactive mixing iteration (rows: components, columns: & 
        target waters in the domain):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (u_new(i,j), j=1,this%num_target_waters_dom)
    end do
    write(2,*)
    write(2,*) "Concentrations of variable activity species after solving reactive mixing iteration (rows: species, & 
        columns: target waters in the domain):"
    write(2,*)
    do i=1,n_nc
        write(2,"(2x,*(ES15.5))") (conc_nc(i,j), j=1,this%num_target_waters_dom)
    end do
    close(2)
    deallocate(u_tilde,u_react,u_new,conc_nc)
end subroutine