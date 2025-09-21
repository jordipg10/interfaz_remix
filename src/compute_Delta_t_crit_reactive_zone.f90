!!> Computes critical time step for RT in a reactive zone
!subroutine compute_Delta_t_crit_reactive_zone(this,B_mat,F_mat,Delta_t_crit)
!    use reactive_zone_Lagr_m
!    use redox_kin_reaction_m
!    implicit none
!    
!    class(reactive_zone_c) :: this
!    class(tridiag_matrix_c), intent(in) :: B_mat !> explicit mixing ratios
!    class(diag_matrix_c), intent(in) :: F_mat !> storage matrix
!    real(kind=8), intent(out) :: Delta_t_crit !> critical time step
!!> Variables
!    real(kind=8) :: B_norm_inf,B_norm_1,B_T_norm,lambda_sum,lipschitz,lipschitz_max
!    integer(kind=4) :: flag,i,j
!    real(kind=8), allocatable :: Sk(:,:)
!            
!    procedure(norm_mat_inf), pointer :: p_norm=>null()
!    
!    B_norm_1=B_mat%compute_norm_1()
!    B_norm_inf=B_mat%compute_norm_inf()
!    if (B_norm_1<=1d0) then
!        p_norm=>norm_mat_inf
!        B_T_norm=B_norm_1
!        flag=0
!    else if (B_norm_inf<=1d0) then
!        p_norm=>norm_mat_1
!        B_T_norm=B_norm_inf
!        flag=1
!    end if
!    
!    lipschitz_max=0d0
!    
!    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
!        lipschitz=this%solid_chemistry%reactive_zone%chem_syst%lin_kin_reacts(i)%lambda
!        if (lipschitz>lipschitz_max) then
!            lipschitz_max=lipschitz
!        end if
!    end do
!    !do i=1,this%num_kin_reactions
!    !>    select type (kin_react=>this%kin_reactions(i)%kin_reaction)
!    !>    type is (redox_kin_c)
!    !>        lipschitz=0d0 !> Lipschitz constant of Jacobian of $r_k$
!    !>            if (kin_react%params%rate_cst>lipschitz) then
!    !>                lipschitz=kin_react%params%rate_cst
!    !>            end if
!    !>    end select
!    !>    if (lipschitz>lipschitz_max) then
!    !>        lipschitz_max=lipschitz
!    !>    end if
!    !end do
!    Delta_t_crit=(1d0-B_T_norm)/(p_norm(transpose(this%solid_chemistry%reactive_zone%chem_syst%Sk))*(1d0/maxval(F_mat%diag))*lipschitz_max)
!end subroutine