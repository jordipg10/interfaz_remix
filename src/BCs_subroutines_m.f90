!> This module contains the subroutines that impose BCs in transition matrix and source term of a PDE
module BCs_subroutines_m
    use PDE_m, only: PDE_1D_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c
    use transport_m, only: transport_1D_c
    use transport_transient_m, only: transport_1D_transient_c, diffusion_1D_transient_c
    use flow_transient_m, only: flow_transient_c
    implicit none
    save
    contains
         
         subroutine Dirichlet_BCs_PDE(this)
         !> Imposes Dirichlet boundary conditions in transition matrix & source term
            implicit none
            class(PDE_1D_c) :: this
            
            integer(kind=4) :: n,opcion,n_flux
            real(kind=8) :: r_n,r_n_12
            
            n=this%spatial_discr%Num_targets
            
            select type (this)
            type is (transport_1D_c)
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%targets_flag.eq.1) then !> edge centred finite differences
                        this%trans_mat%diag(1)=0d0
                        this%trans_mat%super(1)=0d0
                        this%trans_mat%diag(n)=0d0
                        this%trans_mat%sub(n-1)=0d0
                        !this%source_term_PDE(1)=this%BCs%conc_inf
                        !this%source_term_PDE(1)=this%BCs%conc_out
                    else if (mesh%targets_flag.eq.0 .and. mesh%scheme.eq.1) then !> cell centred finite differences, traditional scheme
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) - &
                        3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) + &
                        this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x) - &
                        3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x) + &
                        this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%bd_mat(1)=(this%tpt_props_heterog%flux_cent(1)+ &
                            2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        this%bd_mat(2)=-(this%tpt_props_heterog%flux_cent(n)-&
                            2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n)+this%BCs%conc_out*this%bd_mat(2)
                    else if (mesh%targets_flag.eq.0 .and. mesh%scheme.eq.2) then !> cell-centred FD, Petchame & Carrera (2024)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux_int(1)/mesh%Delta_x + &
                            this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(&
                            mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux_int(n+1)/mesh%Delta_x - &
                            this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(&
                            mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        !this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+&
                        !    this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        !this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n+1)-&
                        !    2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                        this%bd_mat(1)=(this%tpt_props_heterog%flux_int(1)+&
                            this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%bd_mat(2)=-(this%tpt_props_heterog%flux_int(n+1)-&
                            2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        this%source_term_PDE(n)=this%source_term_PDE(n)+this%BCs%conc_out*this%bd_mat(2)
                    else
                        error stop "BCs not implemented yet"
                    end if
                end select
            type is (transport_1D_transient_c)
                select type (mesh=>this%spatial_discr) 
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%targets_flag.eq.1) then !> ECFD
                        this%trans_mat%diag(1)=0d0
                        this%trans_mat%super(1)=0d0
                        this%trans_mat%diag(n)=0d0
                        this%trans_mat%sub(n-1)=0d0
                        !this%source_term_PDE(1)=this%conc_init(1)
                        !this%source_term_PDE(n)=this%conc_init(n)
                    else if (mesh%targets_flag.eq.0 .and. mesh%scheme.eq.1) then !> CCFD traditional
                        this%trans_mat%diag(1)=this%trans_mat%diag(1) - this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) - &
                            3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n) + this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x) - &
                            3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        !this%source_term_PDE(1)=this%source_term_PDE(1) + this%BCs%conc_inf*(this%tpt_props_heterog%flux(1) + &
                        !    2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        !this%source_term_PDE(n)=this%source_term_PDE(n) - this%BCs%conc_out*(this%tpt_props_heterog%flux(n) - &
                        !    2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                        this%bd_mat(1)=(this%tpt_props_heterog%flux_cent(1) + &
                            2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1) + this%BCs%conc_inf*this%bd_mat(1)
                        this%bd_mat(2)=-(this%tpt_props_heterog%flux_cent(n) - &
                            2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n) + this%BCs%conc_out*this%bd_mat(2)
                    else if (mesh%targets_flag.eq.0 .and. mesh%scheme.eq.2) then !> CCFD proposed
                        this%trans_mat%diag(1)=this%trans_mat%diag(1) - this%tpt_props_heterog%flux_int(1)/mesh%Delta_x + &
                            this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(&
                            mesh%Delta_x**2)
                        !this%trans_mat%super(1)=-this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) + &
                        !    this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n) + this%tpt_props_heterog%flux_int(n+1)/mesh%Delta_x - &
                            this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n+1)/(&
                            mesh%Delta_x**2)
                        !this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x) + &
                        !    this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%bd_mat(1)=(this%tpt_props_heterog%flux_int(1) + &
                            2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1) + this%BCs%conc_inf*this%bd_mat(1)
                        this%bd_mat(2)=-(this%tpt_props_heterog%flux_int(n+1) - &
                            2d0*this%tpt_props_heterog%dispersion(n+1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(n)=this%source_term_PDE(n) + this%BCs%conc_out*this%bd_mat(2)
                        !this%source_term_PDE(1)=this%source_term_PDE(1) + this%BCs%conc_inf*(this%tpt_props_heterog%flux(1) + &
                        !    2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        !this%source_term_PDE(n)=this%source_term_PDE(n) - this%BCs%conc_out*(this%tpt_props_heterog%flux(n+1) - &
                        !    2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else
                        error stop "BCs not implemented yet"
                    end if
                end select
            !type is (diffusion_1D_transient_c)
            !    select type (mesh=>this%spatial_discr)
            !    type is (spatial_discr_rad_c)
            !        if (this%dimless.eqv..true.) then
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)-2d0/mesh%Delta_r(n) !> dimensionless
            !            this%source_term_PDE(n)=2d0/mesh%Delta_r(n) !> dimensionless
            !            this%bd_mat(2)=2d0/mesh%Delta_r(n) !> dimensionless
            !        end if
            !    end select
            type is (flow_transient_c)
                select type (mesh=>this%spatial_discr)
                type is (spatial_discr_rad_c)
                    if (this%dimless) then
                        !> Inflow
                        this%trans_mat%diag(1)=-4d0*(mesh%r_min_D**(mesh%dim-1)/mesh%Delta_r_D(1)+((mesh%r_min_D+mesh%Delta_r_D(1))&
                            **(mesh%dim-1))/(sum(mesh%Delta_r_D(1:2)))) !> dimensionless
                        this%bd_mat(1)=4d0*mesh%r_min_D**(mesh%dim-1)/mesh%Delta_r_D(1) !> dimensionless
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%bd_mat(1)*this%BCs%head_inf_D !> dimensionless
                        !> Outflow
                        this%trans_mat%diag(n)=-4d0*(mesh%r_max_D**(mesh%dim-1)/mesh%Delta_r_D(n)+((mesh%r_max_D-mesh%Delta_r_D(n))&
                            **(mesh%dim-1))/(sum(mesh%Delta_r_D(n-1:n)))) !> dimensionless
                        this%bd_mat(2)=4d0*mesh%r_max_D**(mesh%dim-1)/mesh%Delta_r_D(n) !> dimensionless
                        this%source_term_PDE(n)=this%source_term_PDE(n)+this%bd_mat(2)*this%BCs%head_out_D !> dimensionless
                    end if
                end select
            end select
         end subroutine
         
         subroutine Dirichlet_Neumann_BCs_PDE(this)
         !> Imposes Dirichlet & Neumann boundary conditions in transition matrix & source term
            implicit none
            class(PDE_1D_c) :: this
            
            integer(kind=4) :: n,opcion,n_flux
            real(kind=8) :: r_n,r_n_12
            
            n=this%spatial_discr%Num_targets
            
            select type (this)
            !type is (transport_1D_c)
            !    select type (mesh=>this%spatial_discr)
            !    type is (mesh_1D_Euler_homog_c)
            !        if (mesh%targets_flag.eqv.1) then
            !            this%trans_mat%diag(1)=0d0
            !            this%trans_mat%super(1)=0d0
            !            this%trans_mat%diag(n)=0d0
            !            this%trans_mat%sub(n-1)=0d0
            !            !this%source_term_PDE(1)=this%BCs%conc_inf
            !            !this%source_term_PDE(1)=this%BCs%conc_out
            !        else if (mesh%targets_flag.eqv.0 .and. mesh%scheme.eqv.1) then !> CFDS
            !            this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
            !            this%trans_mat%super(1)=-this%tpt_props_heterog%flux(1)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
            !            this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
            !            this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
            !            this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n)-2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
            !        else if (mesh%targets_flag.eqv.0 .and. mesh%scheme.eqv.2) then !> IFDS
            !            this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux(1)/mesh%Delta_x + this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
            !            this%trans_mat%super(1)=-this%tpt_props_heterog%flux(2)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)+this%tpt_props_heterog%flux(n+1)/mesh%Delta_x - this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
            !            this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux(n)/(2d0*mesh%Delta_x) + this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
            !            this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*(this%tpt_props_heterog%flux(1)+2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
            !            this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n+1)-2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
            !        else
            !            error stop "BCs not implemented yet"
            !        end if
            !    end select
            type is (transport_1D_transient_c)
                select type (mesh=>this%spatial_discr) 
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%targets_flag.eq.1) then !> ECFDS
                        this%trans_mat%diag(1)=0d0
                        this%trans_mat%super(1)=0d0
                        this%trans_mat%sub(n-1)=2d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    else if (mesh%targets_flag.eq.0 .and. mesh%scheme.eq.1) then !> CCFD traditional
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) - &
                        3d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        this%bd_mat(1)=(this%tpt_props_heterog%flux_cent(1)+ &
                            2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        !this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n)- &
                        !2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else if (mesh%targets_flag.eq.0 .and. mesh%scheme.eq.2) then !> CCFD proposed
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%tpt_props_heterog%flux_int(1)/mesh%Delta_x + &
                            this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(1)/(&
                            mesh%Delta_x**2)
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) + &
                            this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        this%bd_mat(1)=(this%tpt_props_heterog%flux_int(1)+ &
                            2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_x)/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        !this%source_term_PDE(n)=this%source_term_PDE(n)-this%BCs%conc_out*(this%tpt_props_heterog%flux(n+1)- &
                        !2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                    else
                        error stop "BCs not implemented yet"
                    end if
                end select
            !type is (diffusion_1D_transient_c)
            !    select type (mesh=>this%spatial_discr)
            !    type is (spatial_discr_rad_c)
            !        if (this%dimensionless.eqv..true.) then
            !            this%trans_mat%diag(n)=this%trans_mat%diag(n)-2d0/mesh%Delta_r(n) !> dimensionless
            !            this%source_term_PDE(n)=2d0/mesh%Delta_r(n) !> dimensionless
            !        end if
            !    end select
            end select
         end subroutine
        
         subroutine Neumann_homog_BCs(this)
         !> Imposes Neumann homogeneous boundary conditions in transition matrix
            implicit none
            class(PDE_1D_c) :: this
            
            integer(kind=4) :: n
            
            n=this%spatial_discr%Num_targets

            select type (this)
            type is (transport_1D_transient_c)
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (mesh%scheme.eq.1 .and. mesh%targets_flag==0) then !> CCFD traditional
                        this%trans_mat%super(1)= - this%tpt_props_heterog%flux_cent(1)/(2d0*mesh%Delta_x) +  & 
                        this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1) = this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2) + &
                        this%tpt_props_heterog%flux_cent(n)/(2d0*mesh%Delta_x)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%trans_mat%super(1)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    else if (mesh%scheme.eq.1 .and. mesh%targets_flag==1) then !> ECFD traditional
                        this%trans_mat%super(1)=2d0*this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=2d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%trans_mat%super(1)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    else if (mesh%scheme.eq.2) then !> proposed FD scheme
                        this%trans_mat%super(1)= - this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x) +  & 
                        this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1) = this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2) + &
                        this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)-this%trans_mat%super(1)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    end if
                end select
            end select
         end subroutine 
         
         subroutine Robin_Neumann_homog_BCs(this)
            !> Imposes Robin BC inflow & Neumann homogeneous BC outflow in transition matrix & sink/source term
            implicit none
            class(PDE_1D_c) :: this
            
            integer(kind=4) :: n
            real(kind=8) :: q_cent_1,q_cent_n,q_inf,q_out,D_1,D_n,c_inf,a,b,alpha,beta,gamma
            
            n=this%spatial_discr%Num_targets
            
            select type (this)
            type is (transport_1D_transient_c)
                q_cent_1=this%tpt_props_heterog%flux_cent(1)
                !q_cent_2=this%tpt_props_heterog%flux_cent(2)
                q_cent_n=this%tpt_props_heterog%flux_cent(n)
                q_inf=this%BCs%flux_inf
                q_out=this%BCs%flux_out
                D_1=this%tpt_props_heterog%dispersion(1)
                D_n=this%tpt_props_heterog%dispersion(n)
                c_inf=this%BCs%conc_inf
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (this%tpt_props_heterog%homog_flag.eqv..true. .and. mesh%scheme<3 .and. mesh%targets_flag==0) then !> CCFD homogeneous
                        this%trans_mat%super(1)=-q_cent_1/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)= D_n/(mesh%Delta_x**2) + q_cent_n/(2d0*mesh%Delta_x)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)+(2d0*D_1-q_inf*mesh%Delta_x)*(q_cent_1*mesh%Delta_x+2*D_1)/(&
                            2*q_inf*mesh%Delta_x**3+4*D_1*mesh%Delta_x**2) - 2*D_1/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf*(q_1*mesh%Delta_x+2*D_1)/(&
                        !    q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
                        this%bd_mat(1)=q_inf*(q_cent_1*mesh%Delta_x+2*D_1)/(&
                            q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                    else if ((mesh%scheme.eq.2 .and. mesh%targets_flag.eq.0) .AND. (this%tpt_props_heterog%homog_flag.eqv..false.))&
                        then !> CCFD proposed scheme, heterogeneous
                        this%trans_mat%super(1)=-this%tpt_props_heterog%flux_int(2)/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)= this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) + D_n/(mesh%Delta_x**2)
                        !this%trans_mat%diag(1)=this%trans_mat%diag(1)-(q_inf**2+D_1*(3*q_inf*mesh%Delta_x+2*D_1)/(mesh%Delta_x**2))/(q_inf*mesh%Delta_x+2*D_1) + q_2/(2*mesh%Delta_x)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf/mesh%Delta_x + this%tpt_props_heterog%flux_int(2)/(&
                            2*mesh%Delta_x) - &
                            D_1/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf/mesh%Delta_x
                        this%bd_mat(1)=q_inf/mesh%Delta_x
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                    else if (mesh%targets_flag.eq.1 .AND. this%tpt_props_heterog%homog_flag.eqv..true.) then !> edge centred homogeneous
                        this%trans_mat%super(1)=2d0*D_1/(mesh%Delta_x**2)
                        this%trans_mat%sub(n-1)=2d0*D_n/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf*(q_cent_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x) - &
                            2*D_1/(mesh%Delta_x**2)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        !this%source_term_PDE(1)=this%source_term_PDE(1) + q_inf*(q_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
                        this%bd_mat(1)=q_inf*(q_cent_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                    end if
                type is (spatial_discr_rad_c)
                    if (this%dimless.eqv..true.) then
                        ! this%trans_mat%diag(n)=this%trans_mat%diag(n)-2d0/mesh%Delta_r(n) !> dimensionless
                        ! this%source_term_PDE(n)=2d0/mesh%Delta_r(n) !> dimensionless
                        ! this%bd_mat(2)=2d0/mesh%Delta_r(n) !> dimensionless
                    else if ((mesh%scheme.eq.2 .and. mesh%targets_flag.eq.0) .AND. (this%tpt_props_heterog%homog_flag.eqv..false.))&
                        then !> proposed scheme, cell centred, heterogeneous properties
                        a=mesh%Delta_r(1)
                        b=0.5d0*(mesh%Delta_r(1)+mesh%Delta_r(2))
                        alpha=2d0*D_1+q_inf*a
                        beta=2d0*D_1+q_inf*b
                        gamma=2d0*D_1-q_inf*a
                        this%trans_mat%diag(1)=this%trans_mat%diag(1) + gamma*(1d0+beta/alpha)/(2d0*a*(a+b)) - 2d0*D_1/(a*b) - &
                            0.5*(q_inf/a - this%tpt_props_heterog%flux_int(2)/b)
                        this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                        this%bd_mat(1)=q_inf*(1d0+beta/alpha)/(a+b)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        !this%trans_mat%diag(n)=this%trans_mat%diag(n)-2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_r(1) !> 
                        !this%bd_mat(1)=2d0*this%tpt_props_heterog%dispersion(1)/mesh%Delta_r(1) !> 
                        !this%source_term_PDE(1)=this%BCs%conc_inf*this%bd_mat(1)
                    end if
                end select
            end select
         end subroutine
         
         subroutine Robin_Dirichlet_BCs(this)
         !> Imposes Robin BC inflow & Dirichlet BC outflow in transition matrix & sink/source term
            implicit none
            class(PDE_1D_c) :: this
            
            integer(kind=4) :: n
            real(kind=8) :: q_cent_1,q_cent_2,q_cent_n,q_inf,q_out,D_1,D_n,c_inf
            
            n=this%spatial_discr%Num_targets
            
            select type (this)
            type is (transport_1D_transient_c)
                q_cent_1=this%tpt_props_heterog%flux_cent(1)
                q_cent_2=this%tpt_props_heterog%flux_cent(2)
                q_cent_n=this%tpt_props_heterog%flux_cent(n)
                q_inf=this%BCs%flux_inf
                q_out=this%BCs%flux_out
                D_1=this%tpt_props_heterog%dispersion(1)
                D_n=this%tpt_props_heterog%dispersion(n)
                c_inf=this%BCs%conc_inf
                select type (mesh=>this%spatial_discr)
                type is (mesh_1D_Euler_homog_c)
                    if (this%tpt_props_heterog%homog_flag.eqv. .true. .and. mesh%targets_flag==0) then !> CCFD homogeneous
                        this%trans_mat%super(1)=-q_cent_1/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1)+(2d0*D_1-q_inf*mesh%Delta_x)*(q_cent_1*mesh%Delta_x+2*D_1)/(&
                            2*q_inf*mesh%Delta_x**3+4*D_1*mesh%Delta_x**2) - 2*D_1/(mesh%Delta_x**2)
                        this%bd_mat(1)=q_inf*(q_cent_1*mesh%Delta_x+2*D_1)/(&
                            q_inf*mesh%Delta_x**2+2*D_1*mesh%Delta_x)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        if (mesh%scheme.eq.1) then !> traditional scheme
                            this%trans_mat%diag(n)=this%trans_mat%diag(n) + q_cent_n/(2d0*mesh%Delta_x) - &
                                3d0*this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                            this%trans_mat%sub(n-1)=q_cent_n/(2d0*mesh%Delta_x) + &
                                this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                            this%bd_mat(2)=-(q_cent_n - &
                                2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                            this%source_term_PDE(n)=this%source_term_PDE(n) + this%BCs%conc_out*this%bd_mat(2)
                        else if (mesh%scheme .eq. 2) then !> proposed scheme
                            this%trans_mat%diag(n)=this%trans_mat%diag(n) + q_out/mesh%Delta_x - &
                                this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) - 3d0*this%tpt_props_heterog%dispersion(n)/(&
                                mesh%Delta_x**2)
                            this%trans_mat%sub(n-1)=this%tpt_props_heterog%flux_int(n)/(2d0*mesh%Delta_x) + &
                                this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2)
                            this%bd_mat(2)=-(q_out - &
                                2d0*this%tpt_props_heterog%dispersion(n)/mesh%Delta_x)/mesh%Delta_x
                            this%source_term_PDE(n)=this%source_term_PDE(n) + this%BCs%conc_out*this%bd_mat(2)
                        end if
                    !else if ((mesh%scheme.eq.2 .and. mesh%targets_flag.eq.0) .AND. (this%tpt_props_heterog%homog_flag.eqv..false.))&
                    !    then !> proposed scheme
                    !    this%trans_mat%super(1)=-q_2/(2d0*mesh%Delta_x)+D_1/(mesh%Delta_x**2)
                    !    !this%trans_mat%sub(n-1)= q_n/(2d0*mesh%Delta_x) + D_n/(mesh%Delta_x**2)
                    !    !this%trans_mat%diag(1)=this%trans_mat%diag(1)-(q_inf**2+D_1*(3*q_inf*mesh%Delta_x+2*D_1)/(mesh%Delta_x**2))/(q_inf*mesh%Delta_x+2*D_1) + q_2/(2*mesh%Delta_x)
                    !    this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf/mesh%Delta_x + q_2/(2*mesh%Delta_x) - &
                    !        D_1/(mesh%Delta_x**2)
                    !    !this%trans_mat%diag(n)=this%trans_mat%diag(n)-this%trans_mat%sub(n-1)
                    !    !this%source_term_PDE(1)=this%source_term_PDE(1)+q_inf*c_inf/mesh%Delta_x
                    !    this%bd_mat(1)=q_inf/mesh%Delta_x
                    !    this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                    else if (mesh%targets_flag.eq.1 .AND. this%tpt_props_heterog%homog_flag.eqv..true.) then !> edge centred homogeneous
                        this%trans_mat%super(1)=2d0*D_1/(mesh%Delta_x**2)
                        this%trans_mat%diag(1)=this%trans_mat%diag(1) - q_inf*(q_cent_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x) - &
                            2*D_1/(mesh%Delta_x**2)
                        this%bd_mat(1)=q_inf*(q_cent_1*mesh%Delta_x+2d0*D_1)/(D_1*mesh%Delta_x)
                        this%source_term_PDE(1)=this%source_term_PDE(1)+this%BCs%conc_inf*this%bd_mat(1)
                        this%trans_mat%diag(n)=0d0
                        this%trans_mat%sub(n-1)=0d0
                    end if
                end select
            type is (flow_transient_c)
                select type (mesh=>this%spatial_discr)
                type is (spatial_discr_rad_c)
                    if (this%dimless.eqv..true.) then
                        !> Inflow
                        !this%trans_mat%diag(1)=this%trans_mat%diag(1)-4d0*mesh%r_min/mesh%Delta_r(1) !> dimensionless
                        !this%bd_mat(1)=4d0*mesh%r_min/mesh%Delta_r(1) !> dimensionless
                        !this%source_term_PDE(1)=this%bd_mat(1)*this%BCs%head_inf !> dimensionless
                        !!> Outflow
                        !this%trans_mat%diag(n)=this%trans_mat%diag(n)-4d0*mesh%r_max/mesh%Delta_r(n) !> dimensionless
                        !this%bd_mat(2)=4d0*mesh%r_max/mesh%Delta_r(n) !> dimensionless
                        !this%source_term_PDE(n)=this%bd_mat(2)*this%BCs%head_out !> dimensionless
                    end if
                end select
            end select
         end subroutine
end module