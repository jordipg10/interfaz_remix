!> Writes data and results of 1D diffusion equation
subroutine write_diffusion_1D(this,Time_out,output)
    use diffusion_m, only: diffusion_1D_c
    use spatial_discr_1D_m
    implicit none
    !> Variables
    class(diffusion_1D_c), intent(in) :: this !> 1D diffusion object
    real(kind=8), intent(in) :: Time_out(:)
    real(kind=8), intent(in) :: output(:,:)

    real(kind=8), allocatable :: props_mat(:,:),Delta_t(:),stab_params(:)
    integer(kind=4) :: i,j,k
    
    open(unit=1,file='diffusion_1D.out',status='unknown')
    write(1,"(2x,'Equation:',5x,'0 = T*c + g',/)")
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c)
        write(1,"(2x,'Length of domain:',F15.5/)") mesh%measure
        write(1,"(2x,'Number of cells:',I5/)") this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,'Mesh size:',F15.5/)") mesh%Delta_x
    end select
    if (this%BCs%BCs_label(1) == 1 .and. this%BCs%BCs_label(2) == 1) then
        write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',10x,*(F15.5),/)") this%BCs%conc_inf, this%BCs%conc_out
    else if (this%BCs%BCs_label(1) == 2 .and. this%BCs%BCs_label(2) == 2) then
        write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
    end if
    if (this%spatial_discr%scheme == 1) then
        write(1,"(2x,'Scheme:',10x,'CFD',/)")
    else if (this%spatial_discr%scheme == 2) then
        write(1,"(2x,'Scheme:',10x,'IFD',/)")
    else if (this%spatial_discr%scheme == 3) then
        write(1,"(2x,'Scheme:',10x,'Upwind',/)")
    end if
    write(1,"(2x,'Properties:'/)")
    write(1,"(10x,'Dispersion:'/)")
    write(1,"(/,2x,'Transition matrix T (with BCs):'/)") 
    write(1,"(17x,2F15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
    do i=2,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag-1
        write(1,"(2x,3F15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
    end do
    write(1,"(2x,2F15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1),this%trans_mat%diag(this%spatial_discr%Num_targets)
    write(1,"(/,2x,'Source term g:'/)")
    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,F15.5)") this%source_term_PDE(i)
    end do
    write(1,"(/,2x,'Cell'/)")
    do i=1,this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,I5,ES15.5)") i, this%conc(i)
    end do
    rewind(1)
    close(1)
end subroutine