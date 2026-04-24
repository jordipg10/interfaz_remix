!> @file compute_b_vec_conc_mob.f90
!> @brief Computes right-hand side vector b for implicit linear system in MRMT mobile zone
!> @details This subroutine constructs the right-hand side vector b for the theta-weighted implicit
!> linear system governing mobile zone concentrations in Multi-Rate Mass Transfer (MRMT) models:
!>
!> A * c_mob^(k+1) = b
!>
!> The vector b combines explicit transport of old mobile concentrations with mass transfer
!> contributions from all immobile zones. The theta-weighting allows flexible temporal
!> discretization from fully explicit (θ=0) to fully implicit (θ=1) schemes.
!>
!> @par Mathematical Formulation:
!> b = X * c_mob^k - f
!>
!> where:
!>   - X_sub = ((1-θ)*Δt/φ_m) * T_sub (subdiagonal of explicit transport matrix)
!>   - X_super = ((1-θ)*Δt/φ_m) * T_super (superdiagonal of explicit transport matrix)
!>   - X_diag = 1 + ((1-θ)*Δt/φ_m) * T_diag - (Δt/φ_m) * Σ_j [φ_j * P_j * α_j * (1-θ - θ*α_j*Δt*(1-θ)/(1+α_j*Δt*θ))]
!>   - f = (Δt/φ_m) * Σ_j [φ_j * P_j * α_j * (θ-1 - θ*(1-α_j*Δt*(1-θ))/(1+α_j*Δt*θ)) * c_imm_j^k]
!>
!> @par Physical Interpretation:
!> The RHS vector accounts for:
!>   1. Explicit part of transport operator acting on old mobile concentrations (X * c_mob^k)
!>   2. Mass transfer from immobile zones weighted by their old concentrations (-f term)
!>   3. Theta-weighting provides stability control and accuracy trade-off
!>
!> @see compute_A_mat_conc_mob For coefficient matrix A computation
!> @see compute_conc_imm_MRMT For immobile zone concentration updates
!> @see MRMT_m Multi-Rate Mass Transfer module
!> @author Generated documentation
!> @date November 2025

!> @brief Construct right-hand side vector b for mobile zone implicit linear system in MRMT
!> @param[in] this MRMT object containing mobile zone, immobile zones, and transport properties
!> @param[in] theta Time weighting factor θ ∈ [0,1]: 0=explicit, 1=implicit [-]
!> @param[in] Delta_t Time step size Δt [T]
!> @param[in] conc_mob_old Mobile zone concentrations at previous time step c_mob^k [M/L³]
!> @param[in] conc_imm_old Immobile zone concentrations at previous time step c_imm^k (all zones) [M/L³]
!> @param[out] b_vec Right-hand side vector b for linear system A*c_mob^(k+1)=b [M/L³]
subroutine compute_b_vec_conc_mob(this,theta,Delta_t,conc_mob_old,conc_imm_old,b_vec)
    use MRMT_m, only: MRMT_1D_trans_c !< Import Multi-Rate Mass Transfer class definition
    use arrays_m, only: prod_mat_vec,tridiag_matrix_c !< Import matrix-vector product and tridiagonal matrix class
    implicit none !< Enforce explicit variable declarations for type safety
    class(MRMT_1D_trans_c), intent(in) :: this !< MRMT object with mobile zone porosity, immobile zones (n_imm, properties), PDE transport matrix [-]
    real(kind=8), intent(in) :: theta !< Time weighting factor θ: 0 (explicit) ≤ θ ≤ 1 (implicit) for temporal discretization [-]
    real(kind=8), intent(in) :: Delta_t !< Time step size Δt for temporal integration [T]
    real(kind=8), intent(in) :: conc_mob_old(:) !< Mobile zone concentrations at time k: c_mob^k (spatial distribution) [M/L³]
    real(kind=8), intent(in) :: conc_imm_old(:) !< Immobile zones concentrations at time k: c_imm^k (one value per immobile zone) [M/L³]
    real(kind=8), intent(out) :: b_vec(:) !< Output: right-hand side vector b = X*c_mob^k - f (spatial distribution) [M/L³]
    
    integer(kind=8) :: i,j,n_imm !< Loop counters: i for iterating over immobile zones, j unused (can be removed), n_imm unused (can be removed) [-]
    real(kind=8) :: sum1 !< Accumulator for diagonal correction from all immobile zones (mass transfer contribution to X_diag) [1/T]
    type(tridiag_matrix_c) :: X_mat !< Explicit transport matrix X with bands (sub, diag, super) for computing explicit part [dimensionless]
    real(kind=8), allocatable :: f(:),sum2(:) !< f: mass transfer vector from immobile zones [M/L³], sum2: intermediate accumulator for f computation [M/L³]
    
    !> Validate time weighting factor is within physically meaningful range
    if (theta<0d0 .or. theta>1d0) error stop "Theta must be between 0 and 1" !< Ensure 0 ≤ θ ≤ 1 (required for stability and accuracy of temporal discretization)
    
    !> Compute subdiagonal entries of explicit transport matrix X
    !> Scaling: ((1-θ)*Δt/φ_m) weights the explicit part of transport operator
    !> For θ=0 (fully explicit), coefficient is Δt/φ_m; for θ=1 (fully implicit), coefficient is 0
    X_mat%sub=((1d0-theta)*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%sub !< [-] = [-]*[T]*[-]*[1/T] (subdiagonal: explicit upstream cell coupling)
    
    !> Compute superdiagonal entries of explicit transport matrix X
    !> Scaling: ((1-θ)*Δt/φ_m) weights the explicit part of transport operator
    !> For θ=0 (fully explicit), coefficient is Δt/φ_m; for θ=1 (fully implicit), coefficient is 0
    X_mat%super=((1d0-theta)*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%super !< [-] = [-]*[T]*[-]*[1/T] (superdiagonal: explicit downstream cell coupling)
    
    sum1=0d0 !< Initialize accumulator for immobile zone contributions to diagonal to zero [1/T]
    
    !> First loop: Compute sum of immobile zone contributions to diagonal term of X matrix
    do i=1,this%n_imm !< Iterate through all n_imm immobile zones
        !> Accumulate contribution from immobile zone i to diagonal correction term
        !> Term structure: φ_i * P_i * α_i * [1-θ - θ*α_i*Δt*(1-θ)/(1+α_i*Δt*θ)]
        !> This represents the explicit part of mass transfer exchange weighted by immobile zone properties
        !> The bracketed term accounts for the partitioning between explicit and implicit treatments
        sum1=sum1+this%imm_zones(i)%imm_por*this%imm_zones(i)%prob*this%imm_zones(i)%exch_rate*(1d0-theta-(&
            theta*this%imm_zones(i)%exch_rate*Delta_t*(1d0-theta))/(1d0+this%imm_zones(i)%exch_rate*Delta_t*theta)) !< [1/T] += [-]*[-]*[1/T]*[-] (accumulate explicit exchange contributions)
    end do !< End loop over immobile zones for diagonal correction
    
    !> Compute diagonal entries of explicit transport matrix X combining three terms:
    !> 1. Identity term: 1.0 (accumulation/storage)
    !> 2. Explicit transport operator diagonal: +((1-θ)*Δt/φ_m)*T_diag (note: positive sign for explicit part)
    !> 3. Explicit immobile zone exchange correction: -(Δt/φ_m)*sum1 (mass transfer sink/source)
    !> Note: Sign convention differs from implicit A matrix due to explicit vs implicit formulation
    X_mat%diag=1d0+((1d0-theta)*Delta_t/this%mob_zone%mob_por)*this%PDE%trans_mat%diag-(Delta_t/this%mob_zone%mob_por)*sum1 !< [-] = [-] + [-]*[T]*[-]*[1/T] - [T]*[-]*[1/T] (diagonal: self-coupling + explicit exchange)
    
    !> Allocate sum2 vector with same size as immobile zone concentrations array
    allocate(sum2(size(conc_imm_old))) !< Allocate temporary array for accumulating immobile zone concentration contributions [M/L³]
    
    sum2=0d0 !< Initialize sum2 vector to zero for accumulation [M/L³]
    
    !> Second loop: Compute weighted sum of immobile zone concentrations for mass transfer term
    do i=1,this%n_imm !< Iterate through all n_imm immobile zones to build mass transfer source term
        !> Accumulate weighted contribution from immobile zone i
        !> Term structure: φ_i * P_i * α_i * [θ-1 - θ*(1-α_i*Δt*(1-θ))/(1+α_i*Δt*θ)] * c_imm_i^k
        !> The first bracketed term is always negative (θ-1 ≤ 0), representing mass leaving immobile zones
        !> This contribution depends on old immobile zone concentration c_imm_i^k
        !> The weighting accounts for theta-temporal discretization of the exchange term
        sum2=sum2+this%imm_zones(i)%imm_por*this%imm_zones(i)%prob*this%imm_zones(i)%exch_rate*(theta-1-theta*(&
            1-this%imm_zones(i)%exch_rate*Delta_t*(1-theta))/(1+this%imm_zones(i)%exch_rate*Delta_t*theta))*conc_imm_old(i) !< [M/L³] += [-]*[-]*[1/T]*[-]*[M/L³] (accumulate mass from immobile zones)
    end do !< End loop over immobile zones for mass transfer source
    
    !> Compute mass transfer correction vector f by scaling sum2 with time step and mobile porosity
    !> f represents the total mass contribution from all immobile zones to the RHS
    f=(Delta_t/this%mob_zone%mob_por)*sum2 !< [M/L³] = [T]*[-]*[M/L³] (scale accumulated mass transfer by Δt/φ_m)
    
    !> Compute final right-hand side vector b = X * c_mob^k - f
    !> This combines the explicit transport of old mobile concentrations with mass transfer from immobile zones
    !> prod_mat_vec performs tridiagonal matrix-vector product: X * c_mob^k
    !> Subtracting f adds the mass exchange contribution (note: f has negative components, so -f adds mass)
    b_vec=prod_mat_vec(X_mat,conc_mob_old)-f !< [M/L³] = [M/L³] - [M/L³] (explicit transport minus immobile zone mass transfer)
end subroutine !< End of compute_b_vec_conc_mob subroutine