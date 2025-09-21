!> This subroutine is the interface to solve a reactive mixing iteration for components
!> It has to be called only once in each time step
!> It is supposed to be called after the conservative transport has been solved
!> It assumes there are both equilibrium and kinetic reactions
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It reads the concentrations after conservative transport in a file, and writes the concentrations after reactive mixing in a different file
subroutine interfaz_comps_arch(this,path,num_comps,file_in,Delta_t,file_out)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    character(len=*), intent(in) :: path !> path for input and output files
    integer(kind=4), intent(in) :: num_comps !> number of components
    character(len=*), intent(in) :: file_in !> name of file containing component concentrations after solving conservative transport
    !integer(kind=4), intent(in) :: unit_in !> file unit
    real(kind=8), intent(in) :: Delta_t !> time step
    character(len=*), intent(in) :: file_out !> name of file containing component concentrations after solving conservative transport
    !integer(kind=4), intent(in) :: unit_out !> file unit
!> Variables
    integer(kind=4) :: i,j !> loop variables
    real(kind=8), allocatable :: u_tilde(:,:) !> concentration components after solving conservative transport
    real(kind=8), allocatable :: u_react(:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: u_new(:,:) !> concentration components after solving reactive mixing
!> Pre-process
    allocate(u_tilde(num_comps,this%num_target_waters_dom),u_new(num_comps,this%num_target_waters_dom))
    allocate(u_react(num_comps))
!> Process
    !> We read the component concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_comps
        read(1,*) (u_tilde(i,j), j=1,this%num_target_waters_dom)
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters_dom
        call this%target_waters(this%dom_tar_wat_indices(j))%reaction_iteration_EE_eq_kin_lump(Delta_t,u_react) !> chemical part of components
        u_new(:,j)=u_tilde(:,j)+u_react !> we sum transport and reaction parts
        !print *, u_react
    end do
    !> We write the component concentrations after solving reactive mixing
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    do i=1,num_comps
        write(2,"(1x,*(ES15.5))") (u_new(i,j), j=1,this%num_target_waters_dom)
    end do
    close(2)
!> Post-process
    deallocate(u_tilde,u_react,u_new)
    !print *, this%num_target_waters_dom
end subroutine