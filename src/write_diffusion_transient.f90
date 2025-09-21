!!> Writes data and results of 1D transient diffusion equation
!subroutine write_diffusion_transient(this,Time_out,output)
!    use diffusion_transient_m, only: diffusion_1D_transient_c
!    use spatial_discr_rad_m
!    use vectors_m
!    use matrices_m, only: tridiag_matrix_c
!    implicit none
!    !> Variables
!    class(diffusion_1D_transient_c), intent(in) :: this !> diffusion object
!    real(kind=8), intent(in) :: Time_out(:) !> dimensionless time
!    real(kind=8), intent(in) :: output(:,:) !> dimensionless concentration
!
!    integer(kind=4) :: Num_output,n
!    real(kind=8) :: sum_flux,sum_MRMT,c_m,imm_por
!    real(kind=8), allocatable :: Delta_t(:),Delta_r(:),bd_flux(:),F_im(:)
!    integer(kind=4) :: i,j,k
!    type(tridiag_matrix_c) :: A_mat
!    real(kind=8), parameter :: tol=1d-12
!    character(len=256) :: file_out
!    
!    n=this%spatial_discr%Num_targets
!    Num_output=size(Time_out)
!    write(file_out,"('diffusion_transient_dim',I1,'_n',I2,'.out')") this%spatial_discr%get_dim(), n
!    open(unit=1,file=file_out,status='unknown')
!    if (this%dimensionless.eqv..true.) then
!        write(1,"(2x,'Dimensionless equation:',5x,'F*dc_D/dt_D = T*c_D + g_D',/)")
!        select type (mesh=>this%spatial_discr)
!        type is (spatial_discr_rad_c)
!            Delta_r=this%spatial_discr%get_Cell_size()
!            write(1,"(2x,'Dimension:',I5/)") mesh%dim
!            if (mesh%dim == 1) then
!                write(1,"(2x,'Length of domain:',F15.5/)") mesh%measure
!            else
!                write(1,"(2x,'Radius:',F15.5/)") mesh%radius
!            end if
!        end select
!        write(1,"(2x,'Number of cells:',I5/)") n
!        write(1,"(2x,'Dimensionless mesh:'/)")
!        do i=1,n
!            write(1,"(2x,ES15.5)") Delta_r(i)
!        end do
!        if (this%BCs%BCs_label(1) == 1 .and. this%BCs%BCs_label(2) == 1) then
!            write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',/)")
!        else if (this%BCs%BCs_label(1) == 2 .and. this%BCs%BCs_label(2) == 2) then
!            write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
!        end if
!        if (this%spatial_discr%scheme == 1) then
!            write(1,"(2x,'Scheme:',10x,'CFD'/)")
!        else if (this%spatial_discr%scheme == 2) then
!            write(1,"(2x,'Scheme:',10x,'IFD'/)")
!        end if
!        write(1,"(2x,'Properties:'/)")
!        write(1,"(10x,'Porosity:',ES15.5,10x,'Dispersion:',ES15.5/)") this%diff_props_heterog%porosity(1), &
!        this%diff_props_heterog%dispersion(1)
!        if (this%sol_method == 1) then
!            write(1,"(2x,'Method:',10x,'Numerical in space and time',/)")
!            write(1,"(2x,'Dimensionless time step:'/)")
!            write(1,"(2x,ES15.5/)") this%time_discr%get_Delta_t()
!        else if (this%sol_method == 2) then
!            write(1,"(2x,'Method:',10x,'Eigendecomposition',/)") 
!        end if
!        write(1,"(2x,'Dimensionless final time:'/)")
!        write(1,"(2x,ES15.5/)") this%time_discr%Final_time
!        if (this%sol_method.eq.1) then
!            if (this%time_discr%int_method == 1) then
!                write(1,"(2x,'Integration method:',10x,'Euler explicit'/)")
!            else if (this%time_discr%int_method == 2 .or. this%time_discr%int_method == 3) then
!                write(1,"(2x,'Integration method:',10x,'Euler fully implicit'/)")
!            else if (this%time_discr%int_method == 4) then
!                write(1,"(2x,'Integration method:',10x,'Crank-Nicolson'/)")
!            end if
!        end if
!        select type (mesh=>this%spatial_discr)
!        type is (spatial_discr_rad_c)
!            write(1,"(2x,'Characteristic parameters:'/)")
!            write(1,"(10x,'t_c:',ES15.5,10x,'r_c:',ES15.5/)") this%char_params%char_time, this%char_params%char_length
!        end select
!    
!        write(1,"(2x,'F:'/)") 
!        do i=1,n
!            write(1,"(2x,ES15.5)") this%F_mat%diag(i)
!        end do
!    
!        write(1,"(/,2x,'Transition matrix (with BCs):'/)")
!        write(1,"(17x,2ES15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
!        do i=2,n-1
!            write(1,"(2x,3ES15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
!        end do
!        write(1,"(2x,2ES15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1), &
!        this%trans_mat%diag(this%spatial_discr%Num_targets)
!    
!        if (this%sol_method.eq.2) then
!            call this%compute_A_mat_ODE(A_mat)
!            call A_mat%compute_eigenvalues()
!            call A_mat%compute_eigenvectors()
!            allocate(A_mat%eigenvectors(n,n))
!            call A_mat%check_eigenvectors_tridiag_sym_matrix(tol)
!             write(1,"(/,2x,'A:'/)") 
!             write(1,"(17x,2ES15.5)") A_mat%diag(1), A_mat%sub(1)
!             do i=2,n-1
!                write(1,"(2x,3ES15.5)") A_mat%sub(i-1), A_mat%diag(i), A_mat%sub(i)
!             end do
!             write(1,"(2x,2ES15.5/)") A_mat%sub(n-1), A_mat%diag(n)
!             write(1,"(/,2x,'Eigenvalues of A:'/)")
!             do i=1,n
!                write(1,"(2x,ES15.5)") A_mat%eigenvalues(i)
!             end do
!             write(1,"(/,2x,'Eigenvectors of A (por columnas):'/)")
!             do i=1,n
!                write(1,"(2x,*(ES15.5))") (A_mat%eigenvectors(i,j), j=1,n)
!             end do
!        end if
!        write(1,"(/,2x,'Dimensionless source term:'/)")
!        do i=1,n
!            write(1,"(2x,ES15.5)") this%source_term_PDE(i)
!        end do
!        write(1,"(/,2x,'Initial dimensionless concentration:'/)")
!        do i=1,n
!            write(1,"(2x,F15.5)") this%conc_init(i)
!        end do
!   
!        if (this%sol_method.eq.1) then
!            write(1,"(/,2x,'Cell',*(ES20.5)/)") (Time_out(k), k=1,Num_output)
!            do i=1,n
!                write(1,"(2x,I4,*(F20.5))") i,(output(i,k), k=1,Num_output)
!            end do
!        else
!            call this%prod_total_conc(A_mat)
!            write(1,"(/,2x,'Cell',*(ES15.5),15x,'Limit'/)") (Time_out(k), k=1,Num_output)
!            do i=1,n
!                write(1,"(2x,I4,*(F15.5),F15.5)") i,(output(i,k), k=1,Num_output), this%conc(i)
!            end do
!            allocate(bd_flux(Num_output))
!            do k=1,Num_output
!                sum_flux=0d0
!                do j=1,n
!                    sum_flux=sum_flux+(A_mat%eigenvectors(n,j)**2)*exp(-A_mat%eigenvalues(j)*Time_out(k))/A_mat%eigenvalues(j)
!                end do
!                bd_flux(k)=(this%BCs%conc_inf*this%diff_props_heterog%dispersion(1)/this%char_params%char_length)*sum_flux*4d0/(&
!                Delta_r(n)**2)
!            end do
!            write(1,"(/,2x,'Dimensionless boundary flux:',/)")
!            write(1,"(6x,*(ES15.5)/)") (bd_flux(k)*this%char_params%char_length/(this%BCs%conc_inf*&
!            this%diff_props_heterog%dispersion(1)), k=1,Num_output)
!            c_m=this%BCs%conc_inf
!            allocate(F_im(Num_output))
!            do k=1,Num_output
!                sum_MRMT=0d0
!                do j=1,n
!                    sum_MRMT=sum_MRMT+((A_mat%eigenvectors(n,j)**2)/A_mat%eigenvalues(j))*exp(-A_mat%eigenvalues(j)*Time_out(k))
!                end do
!                F_im(k)=sum_MRMT*this%BCs%conc_inf*4d0*this%diff_props_heterog%dispersion(1)/(this%char_params%char_length*&
!                Delta_r(n)**2)
!            end do
!            if (inf_norm_vec_real(F_im-bd_flux)>=tol) error stop "Fluxes are not equal"
!            write(1,"(/,2x,'Alphas:'/)")
!            do i=1,n
!                write(1,"(ES15.5)") A_mat%eigenvalues(i)/this%char_params%char_time
!            end do
!            imm_por=4d0*this%diff_props_heterog%porosity(1)*this%char_params%char_length*sum((A_mat%eigenvectors(n,:)**2)/(&
!            A_mat%eigenvalues**2))/(Delta_r(n)**2)
!            !write(1,"(/,2x,'Porosidad inmovil para que probabilidades sumen 1:'/)")
!            !write(1,"(ES15.5)") imm_por
!            write(1,"(/,2x,'Probabilities:'/)")
!            do i=1,n
!                write(1,"(ES15.5)") 4d0*this%diff_props_heterog%porosity(1)*this%char_params%char_length*(&
!                A_mat%eigenvectors(n,i)**2)/(imm_por*(Delta_r(n)**2)*(A_mat%eigenvalues(i)**2))
!            end do
!        end if
!    end if
!    rewind(1)
!    close(1)
!end subroutine