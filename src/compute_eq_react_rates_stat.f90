subroutine compute_eq_react_rates_stat(this)
    use RT_m, only: RT_c
    use metodos_sist_lin_m, only: LU_lin_syst
    use vectors_m, only: inf_norm_vec_real
    implicit none
    class(RT_c) :: this
    !> Local variables
    integer :: i
    !> Local arrays
!> Process
    ! select type(this)
    ! class is (RT_1D_transient_c)
    !     ! do i=1,this%num_waters
    !     !     call this%waters(i)%compute_eq_react_rates_stat()
    !     ! end do
    ! end select
    !> Variables for non-lumped formulation (commented out - CHAPUZA -
    !> Extract current state variables and stoichiometric matrices from aqueous chemistry object
    !c2v=this%get_c2v() !> Retrieve current concentrations of secondary non-constant activity species at time k+1 [M/L³]
    ! rk=this%get_rk() !> Commented out (CHAPUZA - workaround): retrieval of current kinetic reaction rates (not used in non-lumped formulation) [M/(L³·T)]
    ! rk_old=this%get_rk_old() !> Commented out (CHAPUZA - workaround): retrieval of old kinetic reaction rates (not used in non-lumped formulation) [M/(L³·T)]
    !Se_nc=this%solid_chemistry%reactive_zone%get_Se_nc_react_zone() !> Retrieve stoichiometric matrix for kinetic reactions (columns for non-constant activity species) - retrieved but not used in current implementation [-]
    !> Construct normal equation matrix A = lambda_r * S_nc · S_nc^T
    !A=lambda_r*matmul(Se_nc,transpose(& !> Transpose the extracted sub-matrix S_nc to get S_nc^T
    !    Se_nc)) !> Result: A = lambda_r * S_nc · S_nc^T (symmetric matrix) [-]
    !!> Construct right-hand side vector b = S_nc · (c2v - c2v_mix-r2nc_mix)
    !b=matmul(Se_nc,c2v-c2v_mix-Delta_t*r2nc_mix) !> Compute S_nc · (c2v - c2v_mix) = change in secondary species concentrations due to equilibrium reactions [M/L³]
        ! matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,& !> Commented out: additional term that would include kinetic reaction contributions
        ! this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Commented out: extraction of S_nc for kinetic reaction term
        ! this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species),Delta_t*(rk_tilde+matmul(transpose(& !> Commented out: multiply time step by sum of mixed kinetic rates and weighted current rates
        ! Sk_nc(:,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:& !> Commented out: extract non-constant activity species columns from Sk_nc
        ! this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)),& !> Commented out: transpose of Sk_nc sub-matrix
        ! theta*rk+(1d0-theta)*rk_old))) !> Commented out: time-weighted kinetic reaction rates θ·rk_new + (1-θ)·rk_old; overall comment indicates b would include kinetic contributions in lumped version
    !> Solve linear system A · Re = b using LU decomposition (with check for near-zero right-hand side)
    !if (inf_norm_vec_real(b)<this%solid_chemistry%reactive_zone%CV_params%zero) then !> Check if infinity norm of b is below zero tolerance (no significant concentration change)
    !    Re=0d0 !> If b ≈ 0, set all equilibrium reaction amounts to zero (no equilibrium reactions needed) [M/L³]
    !else !> Right-hand side is non-zero, solve the linear system
    !    call LU_lin_syst(A,b,this%solid_chemistry%reactive_zone%CV_params%zero,Re) !> Solve A · Re = b using LU decomposition with specified zero tolerance, result stored in Re [M/L³]
    !end if
    !do i=1,this%num_waters
    !    allocate(r_eq(this%waters(i)%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)) !> Allocate equilibrium reaction amounts array with size = total number of equilibrium reactions [-]
    !    
    !end do
end subroutine compute_eq_react_rates_stat
