!subroutine compute_E_mat_tpt(this,E_mat,k)
!!> E=Delta_t(k)*F^(-1)*T (tridiagonal, negative stoich_mat_react_zonemi-definite)
!!> rows sum = 0 if r=0
!>    use transport_transient_1_m
!>    implicit none
!>    
!>    class(transport_1D_transient_c) :: this
!>    class(array_c), pointer, intent(out) :: E_mat
!>    integer(kind=4), intent(in), optional :: k
!>    
!>    real(kind=8) :: sign_vel !> sign of velocity parameter
!>    integer(kind=4) :: i,n
!>    type(tridiag_matrix_c), target :: E_mat_tpt
!>    
!>    n=this%spatial_discr%Num_targets
!>    
!>    select type (trans_mat=>this%trans_mat)
!>    type is (tridiag_matrix_c)
!>        select type (stab_params=>this%stab_params)
!>        type is (stab_params_tpt_homog_c)
!>            if (this%spatial_discr%scheme.eqv.1 .or. this%spatial_discr%scheme.eqv.2) then
!>                E_mat_tpt%sub=stab_params%beta+stab_params%courant/2
!>                E_mat_tpt%super=stab_params%beta-stab_params%courant/2
!>            else if (this%spatial_discr%scheme.eqv.3) then
!>                select type (props=>this%props)
!>                type is (tpt_props_homog_c)
!>                    sign_vel=sign(1d0,props%velocity)
!>                end select
!>                E_mat_tpt%sub=((sign_vel+1)/2d0)*stab_params%courant+stab_params%beta
!>                E_mat_tpt%super=stab_params%beta+((sign_vel-1)/2d0)*stab_params%courant
!>            else
!>                error stop "Scheme not implemented yet"
!>            end if
!>        !type is (stab_params_heterog_c)
!>        !>    select type (props=>this%props)
!>        !>    type is (props_heterog_c)
!>        !>        if (this%spatial_discr%scheme.eqv.1) then
!>        !>            diag=-2d0*stab_params%beta(:,k)
!>        !>            sub=stab_params%beta(2:n,k)+stab_params%alpha(2:n,k)*props%flux(2:n)
!>        !>            super=stab_params%beta(1:n-1,k)-stab_params%alpha(1:n-1,k)*props%flux(1:n-1)
!>        !>        else if (this%spatial_discr%scheme.eqv.2) then
!>        !>            sub=stab_params%beta(2:n,k)+stab_params%alpha(2:n,k)*props%flux(2:n)
!>        !>            super=stab_params%beta(1:n-1,k)-stab_params%alpha(1:n-1,k)*props%flux(2:n)
!>        !>        else
!>        !>            !upwind
!>        !>        end if
!>        !>    end select
!>        end select
!>        E_mat_tpt%diag(1)=-E_mat_tpt%super(1)
!>        E_mat_tpt%diag(2:n-1)=-E_mat_tpt%sub(1:n-2)-E_mat_tpt%super(2:n-1)
!>        E_mat_tpt%diag(n)=-E_mat_tpt%sub(n-1)
!>        select type (time_discr=>this%time_discr)
!>        type is (time_discr_homog_c)
!>            E_mat_tpt%diag=E_mat_tpt%diag-this%props%source_term_flag*time_discr%Delta_t*this%props%source_term
!>        type is (time_discr_heterog_c)
!>            E_mat_tpt%diag=E_mat_tpt%diag-this%props%source_term_flag*time_discr%Delta_t(k)*this%props%source_term
!>        end select
!>        E_mat_tpt%diag=E_mat_tpt%diag/this%F_mat%diag
!>    end select
!>    E_mat=>E_mat_tpt
!end subroutine 