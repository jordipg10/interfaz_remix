!> Computes critical time step for RT model
subroutine compute_Delta_t_crit_RT(this)
    use RT_1D_m, only: RT_1D_c, RT_1D_transient_c, tridiag_matrix_c, norm_mat_inf
    !use chem_system_m, only: chem_system_c

    implicit none
    class(RT_1D_c) :: this
            
    type(tridiag_matrix_c) :: X_mat
    real(kind=8) :: B_norm_inf,B_norm_1,B_T_norm,lambda_sum,lipschitz,lipschitz_max,Delta_t_crit
    integer(kind=4) :: flag,i,j
    !class(chem_system_c), pointer :: chem_syst=>null()
    real(kind=8), allocatable :: Sk(:,:)
            
    procedure(norm_mat_inf), pointer :: p_norm=>null()
            
    select type (this)
    type is (RT_1D_transient_c)
        !select type (transport=>this%transport)
        !type is (transport_1D_transient_c)
            !call this%transport%compute_X_mat(0d0,X_mat)
        if (this%chemistry%chem_syst%num_kin_reacts>0) then
            do i=1,this%chemistry%num_reactive_zones
                !if (this%chemistry%reactive_zones(i)%num_kin_reacts>0) then
                    !call this%chemistry%reactive_zones(i)%compute_Delta_t_crit_reactive_zone(this%transport%X_mat,this%transport%F_mat,Delta_t_crit)
                    if (Delta_t_crit<this%Delta_t_crit) then
                        this%Delta_t_crit=Delta_t_crit
                    end if
                !end if
            end do
        end if
        !end select
    end select
end subroutine