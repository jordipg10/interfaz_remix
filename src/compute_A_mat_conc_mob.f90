!> @file compute_A_mat_conc_mob.f90
!> @brief Computes coefficient matrix A for implicit linear system in MRMT mobile zone
!> @details This subroutine constructs the coefficient matrix A for the theta-weighted implicit
!> linear system governing mobile zone concentrations in Multi-Rate Mass Transfer (MRMT) models:
!>
!> A * c_mob^(k+1) = b
!>
!> The matrix A combines transport operators (advection, dispersion) in the mobile zone with
!> mass transfer terms accounting for exchange with multiple immobile zones. Each immobile zone
!> contributes to the diagonal based on its exchange rate, porosity, and probability.
!>
!> @par Mathematical Formulation:
!> The matrix components are:
!>   - A_sub = -(θ*Δt/φ_m) * T_sub (subdiagonal: upstream coupling)
!>   - A_super = -(θ*Δt/φ_m) * T_super (superdiagonal: downstream coupling)
!>   - A_diag = 1 - (θ*Δt/φ_m) * T_diag + (Δt/φ_m) * Σ_j [φ_j * α_j * (θ - α_j*Δt*θ²/(1+α_j*Δt*θ))]
!>
!> where:
!>   - φ_m = mobile zone porosity [-]
!>   - φ_j = immobile zone j porosity [-]
!>   - α_j = exchange rate for immobile zone j [1/T]
!>   - θ = time weighting factor (0 ≤ θ ≤ 1) [-]
!>   - Δt = time step [T]
!>   - T = transport matrix (tridiagonal) [1/T]
!>
!> @par Physical Interpretation:
!> The diagonal term includes three contributions:
!>   1. Identity (accumulation term)
!>   2. Transport operator (negative diagonal for diffusion/dispersion stability)
!>   3. Mass transfer correction accounting for implicit exchange with all immobile zones
!>
!> @see compute_b_vec_conc_mob For right-hand side vector computation
!> @see compute_conc_imm_MRMT For immobile zone concentration updates
!> @see MRMT_m Multi-Rate Mass Transfer module
!> @author Generated documentation
!> @date November 2025

!> @brief Construct coefficient matrix A for mobile zone implicit linear system in MRMT
!> @param[in] this MRMT object containing mobile zone, immobile zones, and transport properties
!> @param[in] theta Time weighting factor θ ∈ [0,1] for temporal discretization [-]
!> @param[in] Delta_t Time step size Δt [T]
!> @param[out] A_mat Tridiagonal coefficient matrix A for linear system A*c_mob^(k+1)=b [dimensionless]
subroutine compute_A_mat_conc_mob(this,theta,Delta_t,A_mat)
    use MRMT_m, only: MRMT_1D_trans_c !< Import Multi-Rate Mass Transfer class definition
    use arrays_m, only: tridiag_matrix_c !< Import tridiagonal matrix class for efficient storage
    implicit none !< Enforce explicit variable declarations for type safety
    
    class(MRMT_1D_trans_c), intent(in) :: this !< MRMT object with mobile zone porosity, immobile zones, PDE transport matrix [-]
    real(kind=8), intent(in) :: theta !< Time weighting factor θ: 0 (explicit) ≤ θ ≤ 1 (implicit) [-]
    real(kind=8), intent(in) :: Delta_t !< Time step size Δt [T]
    class(tridiag_matrix_c), intent(out) :: A_mat !< Output: tridiagonal coefficient matrix A (sub, diag, super bands) [-]
    
    integer(kind=8) :: i !< Loop counter for iterating over immobile zones (1 to n_imm) [-]
    integer(kind=8) :: j !< Unused loop counter (legacy variable, can be removed) [-]
    real(kind=8) :: sum !< Accumulator for summing immobile zone contributions to diagonal [1/T]
    
    !> Validate time weighting factor is within physically meaningful range
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1" !< Ensure 0 ≤ θ ≤ 1 (stability requirement)
    
    !> Compute subdiagonal entries of A matrix from transport operator
    !> Accounts for implicit treatment of advection/dispersion coupling with upstream cell
    !> Scaling: -(θ*Δt/φ_m) converts transport rate [1/T] to dimensionless coefficient
    A_mat%sub=(-theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%sub !< [-] = [-]*[T]*[-]*[1/T] (subdiagonal: upstream cell coupling)
    
    !> Compute superdiagonal entries of A matrix from transport operator
    !> Accounts for implicit treatment of advection/dispersion coupling with downstream cell
    !> Scaling: -(θ*Δt/φ_m) converts transport rate [1/T] to dimensionless coefficient
    A_mat%super=(-theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%super !< [-] = [-]*[T]*[-]*[1/T] (superdiagonal: downstream cell coupling)
    
    sum=0d0 !< Initialize accumulator for immobile zone contributions to zero [1/T]
    
    !> Loop over all immobile zones to compute their collective contribution to diagonal
    do i=1,this%n_imm !< Iterate through all immobile zones (each with distinct α_i, φ_i, P_i)
        !> Accumulate contribution from immobile zone i to diagonal term
        !> Each zone contributes: φ_i * α_i * [θ - α_i*Δt*θ²/(1+α_i*Δt*θ)]
        !> The bracketed term represents the effective exchange coefficient accounting for
        !> implicit treatment at time k+1. The division by (1+α_i*Δt*θ) arises from
        !> solving the immobile zone equation for c_imm^(k+1) and substituting into mobile zone equation.
        sum=sum+this%imm_zones(i)%imm_por*this%imm_zones(i)%exch_rate*(theta-(this%imm_zones(i)%exch_rate*&
            Delta_t*theta**2)/(1+this%imm_zones(i)%exch_rate*Delta_t*theta)) !< [1/T] += [-]*[-]*[1/T]*[-] (accumulate exchange contributions)
    end do !< End loop over immobile zones
    
    !> Compute diagonal entries of A matrix combining three terms:
    !> 1. Identity term: 1.0 (accumulation/storage)
    !> 2. Transport operator diagonal: -(θ*Δt/φ_m)*T_diag (diffusion/dispersion + advection self-coupling)
    !> 3. Immobile zone exchange: +(Δt/φ_m)*sum (mass transfer sink/source from all immobile zones)
    !> The positive sign on sum reflects that exchange acts as a sink for mobile zone in implicit formulation
    A_mat%diag=1d0-(theta*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%diag+(Delta_t/this%mob_zone%mob_por)*sum !< [-] = [-] - [-]*[T]*[-]*[1/T] + [T]*[-]*[1/T] (diagonal: self-coupling + exchange)
end subroutine !< End of compute_A_mat_conc_mob subroutine