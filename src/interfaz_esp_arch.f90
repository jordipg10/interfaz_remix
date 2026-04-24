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
    integer(kind=4) :: i,j !> loop variables
    real(kind=8), allocatable :: c_tilde(:,:) !> concentration aqueous species after solving conservative transport
    real(kind=8), allocatable :: c_react(:) !> concentration components after solving reactive mixing
    real(kind=8), allocatable :: c_new(:,:) !> concentration components after solving reactive mixing
!> Pre-process
    allocate(c_tilde(num_aq_comps,this%num_target_waters),c_new(num_aq_comps,this%num_target_waters))
    allocate(c_react(num_aq_comps))
!> Process
    !> We read the aqueous species concentrations after solving conservative transport
    open(unit=1,file=path//file_in,status='old',action='read')
    do i=1,num_aq_comps
        read(1,*) (c_tilde(i,j), j=1,this%num_target_waters) !> we read one row (species) at a time
    end do
    close(1)
    !> We solve reactive mixing for each target water
    do j=1,this%num_target_waters !> loop over target waters in the domain
        call this%waters(this%tar_wat_indices(j))%compute_react_term_EE_kin(Delta_t,1.0d0,c_react) !> chemical part of aqueous species concentrations
        c_new(:,j)=c_tilde(:,j)+c_react !> we sum transport and reaction parts
    end do
!> Post-process
    !> We write the aqueous species concentrations after solving reactive mixing
    open(unit=2,file=path//file_out,status='unknown',form='formatted')
    write(2,*) "Aqueous species concentrations after solving reactive mixing iteration (rows: species, columns: & 
        target waters in the domain):"
    write(2,*)
    do i=1,num_aq_comps
        write(2,"(2x,*(ES15.5))") (c_new(i,j), j=1,this%num_target_waters)
    end do
    close(2)
    deallocate(c_tilde,c_react,c_new)
end subroutine
