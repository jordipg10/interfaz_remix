!> \file compute_dc2nc_dc1_ideal.f90
!> \brief Computes the Jacobian of secondary variable activity species concentrations with respect to primary species concentrations for ideal solutions
!> \details Given the mass-action law for ideal solutions (\f$ \gamma = 1 \f$):
!> \f[
!>   c_{2,v,i} = K_i \prod_{j=1}^{n_p} c_{1,j}^{S^*_{e,v,1}(i,j)}
!> \f]
!> the analytical Jacobian entries are:
!> \f[
!>   \frac{\partial c_{2,v,i}}{\partial c_{1,j}} = \frac{c_{2,v,i} \cdot S^*_{e,v,1}(i,j)}{c_{1,j}}
!> \f]
!> where:
!> - \f$ c_{1,j} \f$ = primary species concentrations
!> - \f$ c_{2,v,i} \f$ = secondary variable activity species concentrations (from mass-action law)
!> - \f$ S^*_{e,v,1} \f$ = stoichiometric matrix of equilibrium reactions with respect to primary species
!> - \f$ K_i \f$ = equilibrium constant of reaction \f$ i \f$
!>
!> \param[in] this Aqueous chemistry object containing speciation algebra and stoichiometric matrices
!> \param[in] c1 Primary species concentrations (n_p)
!> \param[in] c2v Secondary variable activity species concentrations (n_e)
!> \param[out] dc2v_dc1 Jacobian matrix \f$ \partial \mathbf{c}_{2,v}/\partial \mathbf{c}_1 \f$ (n_e x n_p)

subroutine compute_dc2v_dc1_ideal(this,c1,c2v,dc2v_dc1)
    use aqueous_chemistry_m, only: aqueous_chemistry_c  !< Import aqueous chemistry class for speciation algebra access
    implicit none                                       !< Enforce explicit variable declarations
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this     !< Aqueous chemistry object containing speciation algebra and stoichiometric matrices [-]
    real(kind=8), intent(in) :: c1(:)                  !< Primary species concentrations \f$ \mathbf{c}_1 \f$ (n_p)
    real(kind=8), intent(in) :: c2v(:)                 !< Secondary variable activity species concentrations \f$ \mathbf{c}_{2,v} \f$ (n_e)
    real(kind=8), intent(out) :: dc2v_dc1(:,:)         !< Output Jacobian \f$ \partial \mathbf{c}_{2,v}/\partial \mathbf{c}_1 \f$ (n_e x n_p)
!> Local variables
    integer(kind=4) :: i                               !< Row loop index over equilibrium reactions (1..n_e) [-]
    integer(kind=4) :: j                               !< Column loop index over primary species (1..n_p) [-]
!> Process: compute each Jacobian entry analytically from mass-action law
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions      !< Loop over each equilibrium reaction \f$ i = 1, \ldots, n_e \f$
        do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species  !< Loop over each primary species \f$ j = 1, \ldots, n_p \f$
            dc2v_dc1(i,j)=c2v(i)*this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star(i,j)/max(c1(j),1d-30)  !< \f$ \frac{\partial c_{2,v,i}}{\partial c_{1,j}} = \frac{c_{2,v,i} \cdot S^*_{e,v,1}(i,j)}{c_{1,j}} \f$
        end do                                          !< End primary species loop
    end do                                              !< End equilibrium reactions loop
end subroutine                                          !< End compute_dc2v_dc1_ideal