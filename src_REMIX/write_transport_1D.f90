!> Writes data and results of 1D transport equation
subroutine write_transport_1D(this)
    use transport_m, only: transport_1D_c, mesh_1D_Euler_homog_c
    implicit none
    !> Variables
    class(transport_1D_c), intent(in) :: this !> 1D transport object

    real(kind=8) :: Delta_x,x,D
    real(kind=8), parameter :: pi=4d0*atan(1d0)
    real(kind=8), allocatable :: anal_sol_yo(:),der_anal_sol_yo(:),adv_term(:),adv_term_anal(:)
    integer(kind=4) :: i,n,n_flux
    character(len=256) :: file_out
    
    n=this%spatial_discr%Num_targets-this%spatial_discr%targets_flag !> number of cells
    
    if (this%spatial_discr%scheme.eq.1) then
        if (this%dimless.eqv..true.) then
            write(file_out,"('transport_1D_CFDS_adim.out')")
        else
            write(file_out,"('transport_1D_CFDS.out')")
        end if
    else if (this%spatial_discr%scheme.eq.2) then
        if (this%dimless.eqv..true.) then
            write(file_out,"('transport_1D_IFDS_adim.out')")
        else
            write(file_out,"('transport_1D_IFDS.out')")
        end if
    end if
    open(unit=1,file=file_out,status='unknown')
    write(1,"(2x,'Equation:',5x,'0 = T*c + g',/)")
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c)
        Delta_x=mesh%Delta_x
        write(1,"(2x,'Length of domain:',F15.5/)") mesh%measure
        write(1,"(2x,'Number of cells:',I5/)") this%spatial_discr%Num_targets-this%spatial_discr%targets_flag
        write(1,"(2x,'Mesh size:',F15.5/)") mesh%Delta_x
    end select
    if (this%BCs%BCs_label(1).eq.1 .and. this%BCs%BCs_label(2).eq.1) then
        write(1,"(2x,'Boundary conditions:',10x,'Dirichlet',10x,2F15.5)") this%BCs%conc_inf, this%BCs%conc_out
        
    else if (this%BCs%BCs_label(1).eq.2 .and. this%BCs%BCs_label(2).eq.2) then
        write(1,"(2x,'Boundary conditions:',10x,'Neumann homogeneous',/)")
    end if
    if (this%spatial_discr%scheme.eq.1) then
        write(1,"(/,2x,'Scheme:',10x,'CFD',/)")
    else if (this%spatial_discr%scheme.eq.2) then
        write(1,"(/,2x,'Scheme:',10x,'IFD',/)")
    else if (this%spatial_discr%scheme.eq.3) then
        write(1,"(/,2x,'Scheme:',10x,'Upwind',/)")
    end if

    if (this%dimless.eqv..false.) then
        write(1,"(2x,'Properties:'/)")
        write(1,"(10x,'Dispersion:',/)")
        do i=1,n
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%dispersion(i)
        end do
        write(1,"(/,10x,'Flux interfaces:'/)")
        n_flux=size(this%tpt_props_heterog%flux_int)
        do i=1,n_flux
            write(1,"(20x,ES15.5)") this%tpt_props_heterog%flux_int(i)
        end do
    end if
      
    write(1,"(/,2x,'Transition matrix T (with BCs):'/)") 
    write(1,"(17x,2F15.5)") this%trans_mat%diag(1), this%trans_mat%super(1)    
    do i=2,n-1
        write(1,"(2x,3F15.5)") this%trans_mat%sub(i-1), this%trans_mat%diag(i), this%trans_mat%super(i)
    end do
    write(1,"(2x,2F15.5/)") this%trans_mat%sub(this%spatial_discr%Num_targets-1),this%trans_mat%diag(this%spatial_discr%Num_targets)
    write(1,"(/,2x,'Source term g:'/)")
    do i=1,n
        write(1,"(2x,F15.5)") this%source_term_PDE(i)
    end do
   
    write(1,"(/,2x,'Cell'/)")
    do i=1,n
        write(1,"(2x,I5,ES15.5)") i, this%conc(i)
    end do

    ! allocate(anal_sol_yo(n),der_anal_sol_yo(n))
    ! if (this%tpt_props_heterog%source_term_order.eq.0) then !> Solucion analitica para q(x)=-x, c(0)=1, c(L)=0
    !     write(1,"(/,2x,'Analytical solution & derivative:'/)")
    !     do i=1,n
    !         x=(2*i-1)*Delta_x/2d0
    !         anal_sol_yo(i)=anal_sol_tpt_1D_stat_flujo_lin(this,x)
    !         der_anal_sol_yo(i)=der_anal_sol_tpt_1D_stat_flujo_lin(this,x)
    !         write(1,"(2x,ES15.5,10x,ES15.5,10x,ES15.5)") x, anal_sol_yo(i), der_anal_sol_yo(i)
    !     end do
    !     write(1,"(/,2x,'L^2 error advection term:'/)")
    !     allocate(adv_term_anal(n),adv_term(n))
    !     if (this%spatial_discr%scheme == 1) then !> CFDS
    !         adv_term_anal(1)=-this%tpt_props_heterog%flux(1)*der_anal_sol_yo(1)
    !         adv_term(1)=-this%tpt_props_heterog%flux(1)*(this%conc(2)-2d0*this%BCs%conc_inf+this%conc(1))/(2d0*Delta_x)
    !         do i=2,n-1
    !             adv_term_anal(i)=-this%tpt_props_heterog%flux(i)*der_anal_sol_yo(i)
    !             adv_term(i)=-this%tpt_props_heterog%flux(i)*(this%conc(i+1)-this%conc(i-1))/(2d0*Delta_x)
    !         end do
    !         adv_term_anal(n)=-this%tpt_props_heterog%flux(n)*der_anal_sol_yo(n)
    !         adv_term(n)=-this%tpt_props_heterog%flux(n)*(2d0*this%BCs%conc_out-this%conc(n)-this%conc(n-1))/(2d0*Delta_x)
    !     else if (this%spatial_discr%scheme == 2) then !> IFDS
    !         adv_term(1)=-(this%tpt_props_heterog%flux(2)*(this%conc(2)-this%conc(1))+this%BCs%flux_inf*(2d0*this%conc(1)-&
    !         2d0*this%BCs%conc_inf))/(2d0*Delta_x)
    !         adv_term_anal(1)=-0.5*(this%tpt_props_heterog%flux(1)+this%tpt_props_heterog%flux(2))*der_anal_sol_yo(1)
    !         do i=2,n-1
    !             adv_term_anal(i)=-0.5*(this%tpt_props_heterog%flux(i)+this%tpt_props_heterog%flux(i+1))*der_anal_sol_yo(i)
    !             adv_term(i)=-0.5*(this%tpt_props_heterog%flux(i+1)*(this%conc(i+1)-this%conc(i))+(this%tpt_props_heterog%flux(i)*&
    !             (this%conc(i)-this%conc(i-1))))/Delta_x
    !         end do
    !         adv_term_anal(n)=-0.5*(this%tpt_props_heterog%flux(n)+this%tpt_props_heterog%flux(n+1))*der_anal_sol_yo(n)
    !         adv_term(n)=-(this%BCs%flux_out*(2d0*this%BCs%conc_out-2d0*this%conc(n))+this%tpt_props_heterog%flux(n)*(this%conc(n)-&
    !         this%conc(n-1)))/(2d0*Delta_x)
    !     end if
    !     write(1,*) p_norm_vec(adv_term_anal-adv_term,2)
    !     write(1,"(/,2x,'L^2 error concentrations:'/)")
    !     write(1,*) p_norm_vec(this%conc-anal_sol_yo,2)
    !     write(1,"(/,2x,'Mean L^2 error concentrations:'/)")
    !     write(1,*) p_norm_vec(this%conc-anal_sol_yo,2)/sqrt(n*1d0)
    ! else if (this%tpt_props_heterog%source_term_order == 1) then
    !     write(1,"(/,2x,'Derivative analytical solution:'/)")
    !     do i=1,n
    !         x=(2*i-1)*Delta_x/2d0
    !         der_anal_sol_yo(i)=der_anal_sol_tpt_1D_stat_flujo_cuad(this,x)
    !         write(1,"(2x,ES15.5)") der_anal_sol_yo(i)
    !     end do
    ! end if
    rewind(1)
    close(1)
end subroutine