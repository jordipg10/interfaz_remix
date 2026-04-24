subroutine compute_f_vec(this,Delta_t)                                               !< Compute source term vector f for PDE right-hand side
!> Mathematical formulation: f=Delta_t(k)*inv(F)*g                         !< f = time_step * inverse(F_matrix) * source_term
    use PDE_transient_m, only: PDE_transient_c                          !< Import transient PDE class
    use time_discr_m, only: time_discr_homog_c, time_discr_heterog_c       !< Import time discretization classes
    implicit none                                                           !< Require explicit variable declaration
    class(PDE_transient_c) :: this                                      !< PDE object containing source terms and matrices
    real(kind=8), intent(in) :: Delta_t                             !< Optional time step index (for adaptive time stepping)
    
    this%f_vec=this%source_term_PDE*Delta_t/this%F_mat%diag                                        !< Initialize f_vector with source term (external forcing, reactions, etc.)
    ! select type (time=>this%time_discr)                                    !< Polymorphic selection based on time discretization type
    ! type is (time_discr_homog_c)                                           !< Homogeneous time stepping (constant time step)
    !     this%f_vec=this%f_vec*time%Delta_t/this%F_mat%diag                 !< Scale by constant time step / F_matrix diagonal (mass matrix inversion)
    ! type is (time_discr_heterog_c)                                         !< Heterogeneous time stepping (adaptive time step)
    !     this%f_vec=this%f_vec*time%Delta_t(k)/this%F_mat%diag              !< Scale by variable time step Delta_t(k) / F_matrix diagonal
    ! end select                                                              !< End time discretization type selection
end subroutine                                                             !< End source vector computation subroutine