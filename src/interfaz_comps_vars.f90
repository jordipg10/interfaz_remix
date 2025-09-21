!> This subroutine is the interface to solve a reactive mixing iteration for components in a target water
!> It has to be called for each target water in the domain at every time step
!> It is supposed to be called after the conservative transport has been solved
!> It assumes there are both equilibrium and kinetic reactions
!> It uses Euler explicit and applies lumping to the kinetic mixing ratios
!> It receives the concentrations after conservative transport as an argument, and returns the concentrations after reactive mixing as an argument as well
subroutine interfaz_comps_vars(this,u_tilde,Delta_t,u_new)
    use chemistry_m, only: chemistry_c
    implicit none
!> Arguments
    class(chemistry_c) :: this !> chemistry object
    !integer(kind=4), intent(in) :: num_comps !> number of components
    real(kind=8), intent(in) :: u_tilde(:,:) !> concentrations after solving conservative transport
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: u_new(:,:) !> concentrations after solving reactive mixing (must be allocated)
!> Variables
    integer(kind=4) :: j !> loop index
    integer(kind=4) :: num_comps !> number of components
    real(kind=8), allocatable :: u_react(:) !> concentration components due to reactions
!> Pre-Process
    allocate(u_react(num_comps))
!> Process
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters_dom
        call this%target_waters(this%dom_tar_wat_indices(j))%reaction_iteration_EE_eq_kin_lump(Delta_t,u_react) !> chemical part of components
        u_new(:,j)=u_tilde(:,j)+u_react !> we sum transport and reaction parts
    end do
!> Post-process
    deallocate(u_react)
end subroutine