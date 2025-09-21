subroutine compute_trans_mat_tpt_transient(this)
!> T: transition matrix (tridiagonal, negative stoich_mat_react_zonemi-definite)
!> rows sum = 0 if r=0
!> F*dc/dt=T*c+g
    use transport_transient_m, only: transport_1D_transient_c
    use spatial_discr_1D_m, only: mesh_1D_Euler_homog_c, mesh_1D_Euler_heterog_c
    use spatial_discr_rad_m, only: spatial_discr_rad_c
    use transport_properties_heterog_m
    implicit none
    
    class(transport_1D_transient_c) :: this
    
    real(kind=8) :: sign_flux !> sign of flux
    real(kind=8) :: r_i_12 !> radius for radial symmetry
    real(kind=8) :: a,b !> auxiliary variables for distances between points
    integer(kind=4) :: i,n
    
    n=this%spatial_discr%Num_targets

    if (this%spatial_discr%adapt_ref.eq.1) then
        deallocate(this%trans_mat%sub,this%trans_mat%diag,this%trans_mat%super)
        call this%allocate_trans_mat()
    end if
    
    select type (mesh=>this%spatial_discr)
    type is (mesh_1D_Euler_homog_c) !> homogeneous mesh
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then !> dimensionless, cell centred
            !if (mesh%scheme.eq.1) then !> traditional centered finite differences
            if (mesh%scheme<3 .and. mesh%scheme>0) then !> centered finite differences
                this%trans_mat%sub=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x)
                this%trans_mat%super=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x)
            !else if (mesh%scheme.eq.2) then !> proposed by Petchamé & Carrera (2024)
            !    this%trans_mat%super(1)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x)
            !    do i=2,n-1
            !        this%trans_mat%sub(i-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x)
            !        this%trans_mat%super(i)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x)
            !    end do
            !    this%trans_mat%sub(n-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x)
            else if (mesh%scheme.eq.3) then !> upwind
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1))
                    this%trans_mat%sub=1d0/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*1d0/mesh%Delta_x
                    this%trans_mat%super=1d0/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*1d0/mesh%Delta_x
                end if
            else
                error stop "Scheme not implemented yet"
            end if
        else if (mesh%targets_flag .eq. 0) then !> cell centred, dimensional
            if (mesh%scheme.eq.1) then !> traditional centered finite differences
                this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2) + this%tpt_props_heterog%flux_cent(2:n)&
                    /(2*mesh%Delta_x)
                this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2) - &
                this%tpt_props_heterog%flux_cent(1:n-1)/(2*mesh%Delta_x)
            else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024)
                this%trans_mat%super(1)=this%tpt_props_heterog%dispersion(2)/(mesh%Delta_x**2) - &
                    this%tpt_props_heterog%flux_int(2)/(2*mesh%Delta_x)
                do i=2,n-1
                    this%trans_mat%sub(i-1)=this%tpt_props_heterog%dispersion(i)/(mesh%Delta_x**2) + &
                        this%tpt_props_heterog%flux_int(i)/(2*mesh%Delta_x)
                    this%trans_mat%super(i)=this%tpt_props_heterog%dispersion(i+1)/(mesh%Delta_x**2) - &
                        this%tpt_props_heterog%flux_int(i+1)/(2*mesh%Delta_x)
                end do
                this%trans_mat%sub(n-1)=this%tpt_props_heterog%dispersion(n)/(mesh%Delta_x**2) + &
                    this%tpt_props_heterog%flux_int(n)/(2*mesh%Delta_x)
            else if (mesh%scheme.eq.3) then !> upwind
                if (minval(this%tpt_props_heterog%flux_cent)>=0d0 .or. maxval(this%tpt_props_heterog%flux_cent)<0d0) then
                    sign_flux=sign(1d0,this%tpt_props_heterog%flux_cent(1))
                    this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*&
                    this%tpt_props_heterog%flux_cent(2:n)/mesh%Delta_x
                    this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*&
                    this%tpt_props_heterog%flux_cent(1:n-1)/mesh%Delta_x
                end if
            else
                error stop "Scheme not implemented yet"
            end if
        !else if (this%dimless.eqv..false. .and. this%tpt_props_heterog%homog_flag.eqv..true.) &
        !    then
        !    this%trans_mat%sub=this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2) + &
        !        this%tpt_props_heterog%flux(1)/(2*mesh%Delta_x)
        !    this%trans_mat%super=this%tpt_props_heterog%dispersion(1)/(mesh%Delta_x**2) - &
        !        this%tpt_props_heterog%flux(1)/(2*mesh%Delta_x)
        end if
    type is (spatial_discr_rad_c)
        if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..true.) then
            !if (mesh%scheme.eq.1) then
            !    this%trans_mat%sub=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x)
            !    this%trans_mat%super=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x)
            !else if (mesh%scheme.eq.2) then
            !    this%trans_mat%super(1)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x)
            !    do i=2,n-1
            !        this%trans_mat%sub(i-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x)
            !        this%trans_mat%super(i)=1d0/(mesh%Delta_x**2) - 1d0/(2*mesh%Delta_x)
            !    end do
            !    this%trans_mat%sub(n-1)=1d0/(mesh%Delta_x**2) + 1d0/(2*mesh%Delta_x)
            !else if (mesh%scheme.eq.3) then !> upwind
            !    if (minval(this%tpt_props_heterog%flux)>=0d0 .or. maxval(this%tpt_props_heterog%flux)<0d0) then
            !        sign_flux=sign(1d0,this%tpt_props_heterog%flux(1))
            !        this%trans_mat%sub=1d0/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*1d0/mesh%Delta_x
            !        this%trans_mat%super=1d0/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*1d0/mesh%Delta_x
            !    end if
            !else
            !    error stop "Scheme not implemented yet"
            !end if
        else if (mesh%targets_flag.eq.0 .and. this%dimless.eqv..false.) then !> non-dimensionless, cell centred
            if (mesh%scheme.eq.1) then !> centered finite differences
                !this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2) + this%tpt_props_heterog%flux(2:n)/&
                !(2*mesh%Delta_x)
                !this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2) - &
                !this%tpt_props_heterog%flux(1:n-1)/(2*mesh%Delta_x)
            else if (mesh%scheme.eq.2) then !> Petchamé & Carrera (2024)
                this%trans_mat%super(1)=2d0*this%tpt_props_heterog%dispersion(2)/(0.5*(mesh%Delta_r(1)+mesh%Delta_r(2))*&
                    (1.5*mesh%Delta_r(1)+0.5*mesh%Delta_r(2))) - &
                    this%tpt_props_heterog%flux_int(2)/(mesh%Delta_r(1)+mesh%Delta_r(2))
                do i=2,n-1
                    a=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i-1))
                    b=0.5*(mesh%Delta_r(i)+mesh%Delta_r(i+1))
                    this%trans_mat%sub(i-1)=2d0*this%tpt_props_heterog%dispersion(i)/(a*(a+b)) + &
                        this%tpt_props_heterog%flux_int(i)/(2*a)
                    this%trans_mat%super(i)=2d0*this%tpt_props_heterog%dispersion(i+1)/(b*(a+b)) - &
                        this%tpt_props_heterog%flux_int(i+1)/(2*b)
                end do
                this%trans_mat%sub(n-1)=2d0*this%tpt_props_heterog%dispersion(n)/(0.5*(mesh%Delta_r(n-1)+mesh%Delta_r(n))*&
                    (1.5*mesh%Delta_r(n)+0.5*mesh%Delta_r(n-1))) + &
                    this%tpt_props_heterog%flux_int(n)/(mesh%Delta_r(n-1)+mesh%Delta_r(n))
            else if (mesh%scheme.eq.3) then !> upwind
                !if (minval(this%tpt_props_heterog%flux)>=0d0 .or. maxval(this%tpt_props_heterog%flux)<0d0) then
                !    sign_flux=sign(1d0,this%tpt_props_heterog%flux(1))
                !    this%trans_mat%sub=this%tpt_props_heterog%dispersion(2:n)/(mesh%Delta_x**2)+((sign_flux+1d0)/2)*&
                !    this%tpt_props_heterog%flux(2:n)/mesh%Delta_x
                !    this%trans_mat%super=this%tpt_props_heterog%dispersion(1:n-1)/(mesh%Delta_x**2)+((sign_flux-1d0)/2)*&
                !    this%tpt_props_heterog%flux(1:n-1)/mesh%Delta_x
                !end if
            else
                error stop "Scheme not implemented yet"
            end if
        else
            error stop "Radial transition matrix not fully implemented yet"
        end if
    end select
    this%trans_mat%diag=-this%tpt_props_heterog%source_term_flag*this%tpt_props_heterog%source_term
    this%trans_mat%diag(2:n-1)=this%trans_mat%diag(2:n-1)-this%trans_mat%sub(1:n-2)-this%trans_mat%super(2:n-1)
end subroutine 