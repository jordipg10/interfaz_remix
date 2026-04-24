!> \file speciation_algebra_m.f90
!> \brief Algebraic framework for chemical speciation calculations
!> \details
!>   This module provides the algebraic structure for solving chemical equilibrium
!>   speciation problems using matrix-based methods.
!>
!>   Theoretical foundation:
!>   - Component matrix approach (De Simoni et al., 2005; Saaltink et al., 1998)
!>   - Stoichiometric matrix decomposition
!>   - Primary-secondary species partitioning
!>   - Constant vs. variable activity species
!>
!>   Key matrices and transformations:
!>   
!>   Stoichiometric matrix S partitioned as:
!>   \f[
!>     S = [S_1 | S_2]
!>   \f]
!>   where S₁ = primary species, S₂ = secondary species
!>
!>   Component matrix U (De Simoni et al., 2005):
!>   \f[
!>     U = [I | -S_2^{-1} S_1]^T
!>   \f]
!>
!>   Mass action law in matrix form:
!>   \f[
!>     \log c_2 = \log K + S_2^{-1} S_1 \log c_1
!>   \f]
!>
!>   Applications:
!>   - Reactive transport speciation
!>   - Equilibrium chemistry solvers
!>   - Multi-component systems
!>   - Cation exchange reactions
!>   - Surface complexation
!>
!>   The module handles:
!>   - Matrix inversions with singularity checks
!>   - Species swapping for numerical stability
!>   - Constant activity species elimination
!>   - Aqueous vs. solid phase partitioning
!>
!> \author Jordi Petchamé-Guerrero
!> \date 2025
module speciation_algebra_m
    use metodos_sist_lin_m, only: inv_matrix
    use arrays_m, only: compute_det, id_matrix
    implicit none
    save
    private
    !> \brief Algebraic structure for chemical speciation
    !> \details
    !>   Contains all matrix transformations and dimensions needed for
    !>   equilibrium speciation calculations.
    !>
    !>   Species classification:
    !>   - Primary species: Independent components (num_prim_species)
    !>   - Secondary species: Dependent via equilibrium (num_sec_species)
    !>   - Constant activity: Fixed activity (e.g., H₂O, pure minerals)
    !>   - Variable activity: Computed from mass action (most aqueous species)
    !>
    !>   Matrix structure:
    !>   - inv_Se_2: Inverse of S₂ for all secondary species
    !>   - inv_Se_nc_2: Inverse of S₂ for variable activity secondary species
    !>   - Se_1_star = -inv_Se_2 · S₁
    !>   - Se_nc_1_star = -inv_Se_nc_2 · S₁
    !>   - logK_star = inv_Se_nc_2 · log K (variable activity)
    !>   - logK_tilde = inv_Se_2 · log K (all species)
    !>   - comp_mat: Component matrix without constant activity species
    !>   - comp_mat_cst_act: Component matrix with constant activity species
    !>   - comp_mat_aq: Aqueous component matrix
    !>
    !>   Flags:
    !>   - flag_comp: TRUE if no constant activity species in component matrix
    !>   - flag_cat_exch: TRUE if cation exchange reactions present
    type, public :: speciation_algebra_s
        logical :: flag_comp !> TRUE if component matrix has no constant activity species (De Simoni et al, 2005)
        logical :: flag_cat_exch !> TRUE if there are cation exchange reactions in reactive zone, FALSE otherwise
        integer(kind=4) :: num_species !> total number of chemical species
        integer(kind=4) :: num_eq_reactions !> number of equilibrium reactions
        integer(kind=4) :: num_cst_act_species !> number of constant activity species (e.g., H₂O, pure solids)
        integer(kind=4) :: num_var_act_species !> number of variable activity species
        integer(kind=4) :: num_sec_aq_species !> number of secondary aqueous species
        integer(kind=4) :: num_aq_sec_var_act_species !> number of aqueous secondary variable activity species
        integer(kind=4) :: num_aq_var_act_species !> number of aqueous variable activity species
        integer(kind=4) :: num_prim_species !> number of primary species (independent components)
        integer(kind=4) :: num_aq_prim_species !> number of aqueous primary species
        integer(kind=4) :: num_sec_species !> number of secondary species (dependent via equilibrium)
        real(kind=8), allocatable :: inv_Se_2(:,:) !> inverse of S₂ (secondary species stoichiometric submatrix)
        real(kind=8), allocatable :: Se_1_star(:,:) !> = -inv_Se_2 · S₁ (relates primary to secondary)
        real(kind=8), allocatable :: inv_Se_nc_2(:,:) !> inverse of S₂ for variable activity secondary species
        real(kind=8), allocatable :: Se_nc_1_star(:,:) !> = -inv_Se_nc_2 · S₁ (variable activity transformation)
        real(kind=8), allocatable :: logK_star(:) !> = inv_Se_nc_2 · log K (transformed equilibrium constants)
        real(kind=8), allocatable :: logK_tilde(:) !> = inv_Se_2 · log K (transformed equilibrium constants, all species)
        real(kind=8), allocatable :: comp_mat_cst_act(:,:) !> component matrix with constant activity species (Saaltink et al, 1998)
        real(kind=8), allocatable :: comp_mat(:,:) !> component matrix without constant activity species (De Simoni et al, 2005)
        real(kind=8), allocatable :: comp_mat_aq(:,:) !> aqueous component matrix without constant activity species
    contains
    !> Set flags and dimensions
        procedure :: set_flag_comp                          !< Set component matrix flag
        procedure :: set_flag_cat_exch                      !< Set cation exchange flag
        procedure :: set_dimensions                         !< Set all dimensions and initialize arrays
        procedure :: set_num_eq_reactions                   !< Set number of equilibrium reactions
    !> Compute dimensions
        procedure :: compute_num_prim_species               !< Compute number of primary species
        procedure :: compute_num_sec_species                !< Compute number of secondary species
        procedure :: compute_num_sec_aq_species             !< Compute number of secondary aqueous species
        procedure :: compute_num_aq_var_act_species         !< Compute number of aqueous variable activity species
        procedure :: compute_num_aq_sec_var_act_species     !< Compute number of aqueous secondary variable activity species
    !> Compute matrices and vectors
        procedure :: compute_arrays                         !< Compute all arrays with stability checks
        procedure, private :: compute_Se_1_star                     !< Compute Se_1_star = -inv_Se_2 · S₁
        procedure, private :: compute_Se_nc_1_star                  !< Compute Se_nc_1_star = -inv_Se_nc_2 · S₁
        procedure, private :: compute_comp_mat                      !< Compute component matrix
        procedure, private :: compute_comp_mat_aq                   !< Compute aqueous component matrix
        procedure, private :: compute_comp_mat_cst_act              !< Compute component matrix with constant activity species
        procedure :: compute_logK_star                      !< Compute logK_star = inv_Se_nc_2 · log K
        procedure :: compute_logK_tilde                     !< Compute logK_tilde = inv_Se_2 · log K
        procedure :: compute_inv_Se_nc_2                    !< Compute inverse of variable activity S₂
        procedure :: compute_inv_Se_2                       !< Compute inverse of S₂
    !> Eliminate species
        procedure :: elim_cst_act_species                   !< Remove constant activity species from algebra
    end type
    
    contains
        !> \brief Set dimensions and initialize speciation algebra arrays
        !> \param[in,out] this        Speciation algebra object
        !> \param[in]     n_sp        Number of species
        !> \param[in]     n_e         Number of equilibrium reactions
        !> \param[in]     n_c         Number of constant activity species
        !> \param[in]     n_aq        Number of aqueous species
        !> \param[in]     n_v_aq     Number of aqueous variable activity species
        !> \param[in]     num_min_kin Optional number of mineral kinetic reactions
        !> \param[in]     num_gas_kin Optional number of gases not in equilibrium
        !> \details
        !>   Sets all dimension counters and initializes the algebraic structure.
        !>   
        !>   Algorithm:
        !>   1. Set basic species counts
        !>   2. Compute num_var_act_species = num_species - num_cst_act_species
        !>   3. Call compute_num_prim_species (handles kinetics)
        !>   4. Call compute_num_sec_species
        !>   5. Call compute_num_sec_aq_species
        !>   6. Call compute_num_aq_sec_var_act_species
        !>   
        !>   After this call, all matrix sizes are ready for allocation.
        subroutine set_dimensions(this,n_sp,n_e,n_c,n_aq,n_v_aq,num_min_kin,num_gas_kin)
            implicit none
            !> Set the dimensions and initialize speciation algebra arrays.
            !> This routine is essential for setting up the algebraic structure for chemical speciation calculations.
            !> It ensures all internal counters and arrays are consistent with the chemical system definition.
            !> The optional arguments allow for flexible handling of kinetic mineral and gas species, which affect the primary species count and matrix structure.
            !> Example usage:
            !>   call spec_algebra%set_dimensions(n_sp, n_e, n_c, n_aq, n_v_aq, num_min_kin, num_gas_kin)
            !>   ! After this, all matrix sizes and counters are ready for equilibrium calculations.
            class(speciation_algebra_s) :: this
            integer(kind=4), intent(in) :: n_sp !> number of species
            integer(kind=4), intent(in) :: n_e !> number of equilibrium reactions
            integer(kind=4), intent(in) :: n_c !> number of constant activity species
            integer(kind=4), intent(in) :: n_aq !> number of aqueous species
            integer(kind=4), intent(in) :: n_v_aq !> number of aqueous variable activity species
            integer(kind=4), intent(in), optional :: num_min_kin !> number of mineral kinetic reactions
            integer(kind=4), intent(in), optional :: num_gas_kin !> number of gases not in equilibrium
            ! Set counters and call supporting routines
            this%num_species=n_sp
            this%num_eq_reactions=n_e
            this%num_cst_act_species=n_c
            this%num_var_act_species=this%num_species-this%num_cst_act_species
            this%num_aq_var_act_species=n_v_aq
            call this%compute_num_prim_species(num_min_kin,num_gas_kin)
            call this%compute_num_sec_species()
            call this%compute_num_sec_aq_species(n_aq)            
            call this%compute_num_aq_sec_var_act_species()
        end subroutine

        !> \brief Compute number of secondary aqueous species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     n_aq Number of aqueous species
        !> \details
        !>   Computes:
        !>   \f[
        !>     n_{\text{sec,aq}} = n_{\text{aq}} - n_{\text{prim,aq}}
        !>   \f]
        !>   
        !>   Used for sizing arrays for secondary aqueous species.
        subroutine compute_num_sec_aq_species(this,n_aq)
            implicit none
            !> Compute number of secondary aqueous species.
            !> This is used to size arrays for secondary aqueous species and is critical for matrix assembly in speciation calculations.
            !> Example:
            !>   call spec_algebra%compute_num_sec_aq_species(n_aq)
            !>   ! Now spec_algebra%num_sec_aq_species is set for array allocation.
            class(speciation_algebra_s) :: this
            integer(kind=4), intent(in) :: n_aq !> number of aqueous species
            this%num_sec_aq_species=n_aq-this%num_aq_prim_species
        end subroutine

        !> \brief Compute number of secondary species
        !> \param[in,out] this Speciation algebra object
        !> \details
        !>   Determines number of secondary species based on flag_comp:
        !>   
        !>   If flag_comp = TRUE (no constant activity in component matrix):
        !>   \f[
        !>     n_{\text{sec}} = n_{\text{eq}} + n_{\text{cst}}
        !>   \f]
        !>   
        !>   If flag_comp = FALSE (constant activity included):
        !>   \f[
        !>     n_{\text{sec}} = n_{\text{eq}}
        !>   \f]
        subroutine compute_num_sec_species(this)
            implicit none
            !> Compute number of secondary species.
            !> Determines the number of secondary species, considering whether the component matrix includes constant activity species.
            !> If constant activity species are present, the secondary species count includes them; otherwise, only equilibrium reactions are counted.
            !> Example:
            !>   call spec_algebra%compute_num_sec_species()
            !>   ! Now spec_algebra%num_sec_species is set for matrix sizing.
            class(speciation_algebra_s) :: this
            if (this%flag_comp.eqv..true.) then
                this%num_sec_species=this%num_eq_reactions+this%num_cst_act_species
            else
                this%num_sec_species=this%num_eq_reactions
            end if
        end subroutine
        
        !> \brief Compute component matrix without constant activity species
        !> \param[in,out] this Speciation algebra object
        !> \details
        !>   Constructs the component matrix U (De Simoni et al., 2005):
        !>   \f[
        !>     U = \begin{bmatrix} I \\ (-S_{nc,2}^{-1} S_1)^T \end{bmatrix}
        !>   \f]
        !>   
        !>   Structure:
        !>   - Left block: Identity matrix (num_prim_species × num_prim_species)
        !>   - Right block: Se_nc_1_star transposed
        !>   
        !>   Size: num_prim_species × num_var_act_species
        subroutine compute_comp_mat(this) !> see De Simoni et al (2005) for defintion of component matrix
            implicit none
            !> Compute component matrix without constant activity species.
            !> The component matrix is central to the algebraic solution of chemical speciation, representing the relationship between primary and variable activity species.
            !> Example:
            !>   call spec_algebra%compute_comp_mat()
            !>   ! spec_algebra%comp_mat is now allocated and filled for equilibrium calculations.
            class(speciation_algebra_s) :: this
            if (allocated(this%comp_mat)) then
                deallocate(this%comp_mat)
            end if
            allocate(this%comp_mat(this%num_prim_species,this%num_var_act_species))
            this%comp_mat(:,1:this%num_prim_species)=id_matrix(this%num_prim_species)
            this%comp_mat(:,this%num_prim_species+1:this%num_var_act_species)=transpose(this%Se_nc_1_star)
        end subroutine
        
        !> \brief Compute aqueous component matrix without constant activity species
        !> \param[in,out] this Speciation algebra object
        !> \details
        !>   Constructs aqueous-only component matrix by extracting aqueous
        !>   portion from full component matrix.
        !>   
        !>   Size: num_aq_prim_species × num_aq_var_act_species
        !>   
        !>   Used for aqueous phase equilibrium calculations.
        subroutine compute_comp_mat_aq(this) !> see De Simoni et al (2005) for defintion of component matrix
            implicit none
            !> Compute aqueous component matrix without constant activity species.
            !> This matrix is used for aqueous phase speciation, ensuring correct mapping between aqueous primary and variable activity species.
            !> Example:
            !>   call spec_algebra%compute_comp_mat_aq()
            !>   ! spec_algebra%comp_mat_aq is now allocated and filled for aqueous equilibrium calculations.
            class(speciation_algebra_s) :: this
            if (allocated(this%comp_mat_aq)) then
                deallocate(this%comp_mat_aq)
            end if
            allocate(this%comp_mat_aq(this%num_aq_prim_species,this%num_aq_var_act_species))
            this%comp_mat_aq(:,1:this%num_aq_prim_species)=id_matrix(this%num_aq_prim_species)
            this%comp_mat_aq(:,this%num_aq_prim_species+1:this%num_aq_var_act_species)=this%comp_mat(:,this%num_prim_species+1:&
                this%num_prim_species+this%num_aq_sec_var_act_species)
        end subroutine
        
        !> \brief Compute component matrix with constant activity species
        !> \param[in,out] this Speciation algebra object
        !> \details
        !>   Constructs component matrix including constant activity species
        !>   (Saaltink et al., 1998).
        !>   
        !>   Structure similar to comp_mat but includes constant activity columns.
        !>   
        !>   Size: num_prim_species × num_species
        subroutine compute_comp_mat_cst_act(this) !> computes component matrix with constant activity species (Saaltink et al, 1998)
            implicit none
            !> Compute component matrix with constant activity species.
            !> This matrix is used when constant activity species (e.g., water, inert) are present, affecting the algebraic structure and solution of the chemical system.
            !> Example:
            !>   call spec_algebra%compute_comp_mat_cst_act()
            !>   ! spec_algebra%comp_mat_cst_act is now allocated and filled for systems with constant activity species.
            class(speciation_algebra_s) :: this
            if (allocated(this%comp_mat_cst_act)) then
                deallocate(this%comp_mat_cst_act)
            end if
            allocate(this%comp_mat_cst_act(this%num_prim_species,this%num_species))
            this%comp_mat_cst_act(:,1:this%num_prim_species)=id_matrix(this%num_prim_species)
            this%comp_mat_cst_act(:,this%num_prim_species+1:this%num_species)=transpose(this%Se_1_star)
        end subroutine

        !> \brief Compute inverse of secondary equilibrium stoichiometric submatrix
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     Se_2 Secondary equilibrium stoichiometric submatrix
        !> \param[in]     tol  Tolerance for computing inverse matrix
        !> \details
        !>   Computes \f$ S_2^{-1} \f$ using matrix inversion for equilibrium reactions.
        !>   
        !>   Algorithm:
        !>   1. Check matrix is square (num_eq_reactions × num_eq_reactions)
        !>   2. Handle special case: 1×1 matrix (direct inversion)
        !>   3. General case: Call inv_matrix for Gauss-Jordan elimination
        !>   
        !>   Critical for component matrix computation: \f$ S_1^* = -S_2^{-1} S_1 \f$
        subroutine compute_inv_Se_2(this,Se_2,tol)
            implicit none
            !> Compute inverse of secondary equilibrium stoichiometric submatrix.
            !> This is a key step in matrix algebra for speciation, enabling the calculation of dependent species concentrations from equilibrium relationships.
            !> Example:
            !>   call spec_algebra%compute_inv_Se_2(Se_2, tol)
            !>   ! spec_algebra%inv_Se_2 now contains the inverse for equilibrium calculations.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(in) :: Se_2(:,:) !> secondary equilibrium stoichiometic sbmatrix
            real(kind=8), intent(in) :: tol !> for computing inverse matrix
            if (size(Se_2,1)/=size(Se_2,2)) then
                error stop "Matrix must be square (compute_inv_Se_2)"
            end if
            if (allocated(this%inv_Se_2)) then
                deallocate(this%inv_Se_2)
            end if
            allocate(this%inv_Se_2(this%num_eq_reactions,this%num_eq_reactions))
            if (size(Se_2,1)==1) then
                this%inv_Se_2(1,1)=1d0/Se_2(1,1)
            else
                call inv_matrix(Se_2,tol,this%inv_Se_2)
            end if
        end subroutine
        
        !> \brief Compute inverse for variable activity secondary species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     Se_nc_2 Secondary variable activity equilibrium stoichiometric submatrix
        !> \param[in]     tol     Tolerance for computing inverse matrix
        !> \details
        !>   Computes \f$ S_{nc,2}^{-1} \f$ for variable activity equilibrium reactions.
        !>   
        !>   Used when constant activity species are excluded from equilibrium matrix.
        !>   
        !>   Size: num_eq_reactions × num_eq_reactions (excluding constant activity)
        subroutine compute_inv_Se_nc_2(this,Se_nc_2,tol)
            implicit none
            !> Compute inverse of secondary variable activity equilibrium stoichiometric submatrix.
            !> This is used for systems where some species have variable activity, affecting the algebraic solution and matrix structure.
            !> Example:
            !>   call spec_algebra%compute_inv_Se_nc_2(Se_nc_2, tol)
            !>   ! spec_algebra%inv_Se_nc_2 now contains the inverse for variable activity equilibrium calculations.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(in) :: Se_nc_2(:,:) !> secondary variable activity equilibrium stoichiometic sbmatrix
            real(kind=8), intent(in) :: tol !> for computing inverse matrix
            if (allocated(this%inv_Se_nc_2)) then
                deallocate(this%inv_Se_nc_2)
            end if
            allocate(this%inv_Se_nc_2(this%num_eq_reactions,this%num_eq_reactions))
            call inv_matrix(Se_nc_2,tol,this%inv_Se_nc_2)
        end subroutine
        
        !> \brief Compute all speciation algebra arrays and matrices
        !> \param[in,out] this Speciation algebra object
        !> \param[in,out] Se   Equilibrium stoichiometric matrix (may be modified by swapping)
        !> \param[in,out] K    Equilibrium constants (may be swapped)
        !> \param[in]     tol  Tolerance for computing inverse matrix and determinant
        !> \param[out]    flag TRUE if stoichiometric matrix has been modified, FALSE otherwise
        !> \param[out]    swap Indices of swapped species or reactions (2 elements)
        !> \details
        !>   Main routine for speciation algebra setup. Ensures matrices are invertible
        !>   by swapping species/reactions as needed for numerical stability.
        !>   
        !>   Algorithm (flag_comp = TRUE, variable activity case):
        !>   1. Extract \f$ S_{nc,2} \f$ (secondary variable activity submatrix)
        !>   2. Compute determinant of \f$ S_{nc,2} \f$
        !>   3. If singular (det < tol or error):
        !>      - Swap a primary species with first secondary variable activity species
        !>      - Repeat until invertible
        !>   4. Compute all derived matrices: inv_Se_nc_2, Se_nc_1_star, logK_star, comp_mat
        !>   
        !>   Algorithm (flag_comp = FALSE, constant activity case):
        !>   1. Extract \f$ S_2 \f$ (all secondary species submatrix)
        !>   2. Compute determinant of \f$ S_2 \f$
        !>   3. If singular:
        !>      - Swap equilibrium reactions
        !>      - Repeat until invertible
        !>   4. Compute all derived matrices: inv_Se_2, Se_1_star, logK_tilde, comp_mat_cst_act
        !>   
        !>   Swapping strategy avoids numerical issues from near-singular matrices.
        subroutine compute_arrays(this,Se,K,tol,flag,swap)
            implicit none
            !> Compute all speciation algebra arrays and matrices.
            !> This routine is the core of the algebraic setup, ensuring all matrices are invertible and correctly assembled.
            !> It performs error checking, swaps species/reactions if needed, and updates all dependent arrays for subsequent equilibrium calculations.
            !> Example:
            !>   call spec_algebra%compute_arrays(Se, K, tol, flag, swap)
            !>   ! All matrices and vectors are now ready for equilibrium calculations.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(inout) :: Se(:,:) !> equilibrium stoichiometric matrix
            real(kind=8), intent(inout) :: K(:) !> equilibrium constants
            real(kind=8), intent(in) :: tol !> for computing inverse matrix
            logical, intent(out) :: flag !> TRUE if stoichiometric matrix has been modified, FALSE otherwise
            integer(kind=4), intent(out) :: swap(:) !> indices of swapped species or reactions in reactive zone
            
            integer(kind=4) :: i
            real(kind=8) :: det
            logical :: error
            real(kind=8), allocatable :: Se_nc_2(:,:),Se_2(:,:),aux_Se(:,:),aux_K(:)
            
            flag=.false. !> by default
            !aux_Se=Se
            
            !allocate(swap(2)) !> chapuza
            !swap=0 !> no columns are swapped by default
            
            i=0
            !> esto de abajo deberia ir dentro del if
            aux_Se=Se
            aux_K=K
            if (this%flag_comp .eqv. .true.) then
                allocate(Se_nc_2(this%num_eq_reactions,this%num_eq_reactions))
                do
                    Se_nc_2=Se(:,this%num_prim_species+1:this%num_var_act_species) !> secondary variable activity equilibrium stoichiometric submatrix
                    call compute_det(Se_nc_2,tol,det,error)
                    if ((error .eqv. .true.) .or. (abs(det)<tol)) then
                        flag=.true.
                        i=i+1
                        if (i > this%num_prim_species) then
                            print *, "ERROR in compute_arrays: Se_nc_2 submatrix is singular for all column swaps"
                            print *, "  num_prim_species =", this%num_prim_species
                            print *, "  num_eq_reactions =", this%num_eq_reactions
                            print *, "  num_var_act_species =", this%num_var_act_species
                            print *, "  num_cst_act_species =", this%num_cst_act_species
                            print *, "  Se shape:", shape(Se)
                            error stop "compute_arrays: Se_nc_2 submatrix is singular for all possible column swaps"
                        end if
                        Se=aux_Se
                        Se(:,this%num_prim_species+1)=aux_Se(:,this%num_prim_species-i+1)
                        Se(:,this%num_prim_species-i+1)=aux_Se(:,this%num_prim_species+1)
                    else
                        call this%compute_inv_Se_nc_2(Se_nc_2,tol)
                        swap(1)=this%num_prim_species-i+1 !> index of primary species to swap
                        swap(2)=this%num_prim_species+1 !> index of secondary variable activity species to swap
                        call this%compute_Se_nc_1_star(Se(:,1:this%num_prim_species))
                        call this%compute_logK_star(K)
                        call this%compute_comp_mat()
                        call this%compute_comp_mat_aq()
                        exit
                    end if
                end do
                deallocate(Se_nc_2)
            else
                allocate(Se_2(this%num_eq_reactions,this%num_eq_reactions))
                do
                    Se_2=Se(:,this%num_prim_species+1:this%num_species) !> secondary equilibrium stoichiometric submatrix
                    call compute_det(Se_2,tol,det,error)
                    if ((error .eqv. .true.) .or. (abs(det)<tol)) then
                        flag=.true.
                        i=i+1
                        if (i >= this%num_eq_reactions) then
                            print *, "ERROR in compute_arrays: Se_2 submatrix is singular for all row swaps"
                            print *, "  num_prim_species =", this%num_prim_species
                            print *, "  num_eq_reactions =", this%num_eq_reactions
                            print *, "  num_species =", this%num_species
                            print *, "  Se shape:", shape(Se)
                            error stop "compute_arrays: Se_2 submatrix is singular for all possible row swaps"
                        end if
                        Se=aux_Se
                        !> We swap equilibrium reactions
                        Se(this%num_eq_reactions,:)=aux_Se(this%num_eq_reactions-i,:) !> index of first equilibrium reaction to swap
                        Se(this%num_eq_reactions-i,:)=aux_Se(this%num_eq_reactions,:) !> index of second equilibrium reaction to swap
                        !Se(:,this%num_species)=aux_Se(:,this%num_species-i)
                        !Se(:,this%num_species-i)=aux_Se(:,this%num_species)
                    else
                        call this%compute_inv_Se_2(Se_2,tol)
                        swap(1)=this%num_eq_reactions-i
                        swap(2)=this%num_eq_reactions
                        !swap(1)=this%num_species-i
                        !swap(2)=this%num_species
                        call this%compute_Se_1_star(Se(:,1:this%num_prim_species))
                        K(swap(1))=aux_K(swap(2))
                        K(swap(2))=aux_K(swap(1))
                        call this%compute_logK_tilde(K)
                        call this%compute_comp_mat_cst_act()
                        exit
                    end if
                end do
                deallocate(Se_2)
            end if
        end subroutine
        
        
        !> \brief Compute Se_nc_1_star matrix for variable activity species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     Se_1 Primary equilibrium stoichiometric submatrix
        !> \details
        !>   Computes the matrix:
        !>   \f[
        !>     S_{nc,1}^* = -S_{nc,2}^{-1} \cdot S_1
        !>   \f]
        !>   
        !>   Used to express secondary variable activity species in terms of primary species.
        !>   
        !>   Size: num_eq_reactions × num_prim_species
        subroutine compute_Se_nc_1_star(this,Se_1)
            implicit none
            !> Compute Se_nc_1_star matrix.
            !> This matrix is used to relate primary species to secondary variable activity species in the algebraic solution.
            !> Example:
            !>   call spec_algebra%compute_Se_nc_1_star(Se_1)
            !>   ! spec_algebra%Se_nc_1_star is now allocated and filled.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(in) :: Se_1(:,:) !> primary equilibrium stoichiometric submatrix
            if (allocated(this%Se_nc_1_star)) then
                deallocate(this%Se_nc_1_star)
            end if
            allocate(this%Se_nc_1_star(this%num_eq_reactions,this%num_prim_species))
            this%Se_nc_1_star=-matmul(this%inv_Se_nc_2,Se_1)
        end subroutine
        
        !> \brief Compute Se_1_star matrix for all secondary species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     Se_1 Primary equilibrium stoichiometric submatrix
        !> \details
        !>   Computes the matrix:
        !>   \f[
        !>     S_1^* = -S_2^{-1} \cdot S_1
        !>   \f]
        !>   
        !>   Used to express all secondary species (including constant activity) in terms of primary species.
        !>   
        !>   Size: num_eq_reactions × num_prim_species
        subroutine compute_Se_1_star(this,Se_1)
            implicit none
            !> Compute Se_1_star matrix.
            !> This matrix is used to relate primary species to secondary species in the algebraic solution for systems without constant activity species.
            !> Example:
            !>   call spec_algebra%compute_Se_1_star(Se_1)
            !>   ! spec_algebra%Se_1_star is now allocated and filled.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(in) :: Se_1(:,:) !> primary equilibrium stoichiometric submatrix
            if (allocated(this%Se_1_star)) then
                deallocate(this%Se_1_star)
            end if
            allocate(this%Se_1_star(this%num_eq_reactions,this%num_prim_species))
            this%Se_1_star=-matmul(this%inv_Se_2,Se_1)
        end subroutine
        
        !> \brief Compute transformed equilibrium constants for variable activity species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     K    Equilibrium constants
        !> \details
        !>   Computes the transformed equilibrium constants:
        !>   \f[
        !>     \log K^* = S_{nc,2}^{-1} \cdot \log K
        !>   \f]
        !>   
        !>   Used in mass action law calculations for variable activity systems.
        subroutine compute_logK_star(this,K)
            implicit none
            !> Compute logK_star vector.
            !> This vector is used in equilibrium calculations for systems with variable activity species, representing transformed equilibrium constants.
            !> Example:
            !>   call spec_algebra%compute_logK_star(K)
            !>   ! spec_algebra%logK_star is now ready for equilibrium calculations.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(in) :: K(:) !> equilibrium constants
            this%logK_star=matmul(this%inv_Se_nc_2,log10(K))
        end subroutine
        
        !> \brief Compute transformed equilibrium constants for all species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     K    Equilibrium constants
        !> \details
        !>   Computes the transformed equilibrium constants:
        !>   \f[
        !>     \log \tilde{K} = S_2^{-1} \cdot \log K
        !>   \f]
        !>   
        !>   Used in mass action law calculations including constant activity species.
        subroutine compute_logK_tilde(this,K)
            implicit none
            !> Compute logK_tilde vector.
            !> This vector is used in equilibrium calculations for systems without constant activity species, representing transformed equilibrium constants.
            !> Example:
            !>   call spec_algebra%compute_logK_tilde(K)
            !>   ! spec_algebra%logK_tilde is now ready for equilibrium calculations.
            class(speciation_algebra_s) :: this
            real(kind=8), intent(in) :: K(:) !> equilibrium constants
            !print *, log10(K)
            this%logK_tilde=matmul(this%inv_Se_2,log10(K))
        end subroutine
        
        !> \brief Compute number of primary species
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     num_min_kin Number of kinetic minerals (optional)
        !> \param[in]     num_gas_kin Number of kinetic gases (optional)
        !> \details
        !>   Determines number of primary species based on flag_comp:
        !>   
        !>   If flag_comp = TRUE (no constant activity in component matrix):
        !>   \f[
        !>     n_{\text{prim}} = n_{\text{tot}} - n_{\text{eq}} - n_{\text{cst}}
        !>   \f]
        !>   
        !>   If flag_comp = FALSE (constant activity included):
        !>   \f[
        !>     n_{\text{prim}} = n_{\text{tot}} - n_{\text{eq}}
        !>   \f]
        !>   
        !>   Then computes num_aq_prim_species by subtracting:
        !>   - 1 if cation exchange present (flag_cat_exch = TRUE)
        !>   - num_min_kin if kinetic minerals present and flag_comp = FALSE
        !>   - num_gas_kin if kinetic gases present and flag_comp = FALSE
        subroutine compute_num_prim_species(this,num_min_kin,num_gas_kin)
            implicit none
            !> Compute number of primary species and aqueous primary species.
            !> This routine ensures the correct partitioning of species for matrix assembly and equilibrium calculations.
            !> It handles special cases for cation exchange, mineral kinetics, and gas kinetics, which affect the algebraic structure.
            !> Example:
            !>   call spec_algebra%compute_num_prim_species(num_min_kin, num_gas_kin)
            !>   ! Now spec_algebra%num_prim_species and %num_aq_prim_species are set.
            class(speciation_algebra_s) :: this
            integer(kind=4), intent(in), optional :: num_min_kin
            integer(kind=4), intent(in), optional :: num_gas_kin
            if (this%flag_comp.eqv..true.) then
                this%num_prim_species=this%num_species-this%num_eq_reactions-this%num_cst_act_species
            else
                this%num_prim_species=this%num_species-this%num_eq_reactions
            end if
            this%num_aq_prim_species=this%num_prim_species
            if (this%flag_cat_exch.eqv..true.) then
                this%num_aq_prim_species=this%num_aq_prim_species-1
            end if
            if (present(num_min_kin) .and. present(num_gas_kin)) then
                if (num_min_kin>0 .and. this%flag_comp.eqv..false.) then
                    this%num_aq_prim_species=this%num_aq_prim_species-num_min_kin
                end if
                if (num_gas_kin>0 .and. this%flag_comp.eqv..false.) then
                    this%num_aq_prim_species=this%num_aq_prim_species-num_gas_kin
                end if
            end if
        end subroutine
        
        !> \brief Compute number of aqueous variable activity species
        !> \param[in,out] this Speciation algebra object
        !> \details
        !>   For variable activity component matrix (flag_comp = TRUE):
        !>   \f[
        !>     n_{\text{aq,var}} = n_{\text{aq,prim}} + n_{\text{aq,sec,var}}
        !>   \f]
        subroutine compute_num_aq_var_act_species(this)
            implicit none
            !> Compute number of aqueous variable activity species.
            !> This is used for matrix sizing and equilibrium calculations in systems with variable activity aqueous species.
            !> Example:
            !>   call spec_algebra%compute_num_aq_var_act_species()
            !>   ! Now spec_algebra%num_aq_var_act_species is set.
            class(speciation_algebra_s) :: this
            if (this%flag_comp.eqv..true.) then
                this%num_aq_var_act_species=this%num_aq_prim_species+this%num_aq_sec_var_act_species
            end if
        end subroutine
        
        !> \brief Compute number of aqueous secondary variable activity species
        !> \param[in,out] this Speciation algebra object
        !> \details
        !>   Computes number of aqueous secondary variable activity species based on flag_comp:
        !>   
        !>   If flag_comp = TRUE:
        !>   \f[
        !>     n_{\text{aq,sec,var}} = n_{\text{aq,var}} - n_{\text{aq,prim}}
        !>   \f]
        !>   
        !>   If flag_comp = FALSE:
        !>   \f[
        !>     n_{\text{aq,sec,var}} = n_{\text{aq,var}} - n_{\text{aq,prim}} + 1
        !>   \f]
        subroutine compute_num_aq_sec_var_act_species(this)
            implicit none
            !> Compute number of aqueous secondary variable activity species.
            !> This is used for matrix sizing and equilibrium calculations in systems with secondary variable activity aqueous species.
            !> Example:
            !>   call spec_algebra%compute_num_aq_sec_var_act_species()
            !>   ! Now spec_algebra%num_aq_sec_var_act_species is set.
            class(speciation_algebra_s) :: this
            !integer(kind=4), intent(in) :: wat_flag !> cst act water flag
            if (this%flag_comp.eqv..true.) then !> chapuza
                this%num_aq_sec_var_act_species=this%num_aq_var_act_species-this%num_aq_prim_species
            else
                this%num_aq_sec_var_act_species=this%num_aq_var_act_species-this%num_aq_prim_species+1
            end if
        end subroutine
        
        !> \brief Set component matrix flag
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     flag_comp Logical flag (TRUE = no constant activity in component matrix, FALSE = include constant activity)
        !> \details
        !>   Sets flag_comp which determines:
        !>   - Matrix sizing
        !>   - Which equations to use (Se_nc_* vs Se_*)
        !>   - Component matrix structure
        subroutine set_flag_comp(this,flag_comp)
            implicit none
            !> Set flag for component matrix type.
            !> This affects the algebraic structure and solution of the chemical system.
            !> Example:
            !>   call spec_algebra%set_flag_comp(.true.)
            !>   ! Now spec_algebra%flag_comp is set for matrix assembly.
            class(speciation_algebra_s) :: this
            logical, intent(in) :: flag_comp
            this%flag_comp=flag_comp
        end subroutine
        
        !> \brief Set cation exchange flag
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     flag_cat_exch Logical flag (TRUE = cation exchange present, FALSE = no cation exchange)
        !> \details
        !>   Sets flag_cat_exch which affects:
        !>   - num_aq_prim_species calculation (subtract 1 if TRUE)
        !>   - Exchange site component in mass balance
        subroutine set_flag_cat_exch(this,flag_cat_exch)
            implicit none
            !> Set flag for cation exchange reactions.
            !> This affects the algebraic structure and partitioning of aqueous primary species.
            !> Example:
            !>   call spec_algebra%set_flag_cat_exch(.true.)
            !>   ! Now spec_algebra%flag_cat_exch is set for matrix assembly.
            class(speciation_algebra_s) :: this
            logical, intent(in) :: flag_cat_exch
            this%flag_cat_exch=flag_cat_exch
        end subroutine

    
        !> \brief Eliminate constant activity species from speciation algebra
        !> \param[in,out] this Speciation algebra object
        !> \param[in]     num_surf_compl Number of surface complexes
        !> \param[in]     n_aq           Number of aqueous species
        !> \param[in]     num_min_kin    Number of kinetic minerals (optional)
        !> \param[in]     num_gas_kin    Number of kinetic gases (optional)
        !> \details
        !>   Comprehensive routine to set up speciation algebra excluding constant activity species.
        !>   
        !>   Algorithm:
        !>   1. Set flag_comp = TRUE (no constant activity in component matrix)
        !>   2. Set flag_cat_exch based on surface complexes (TRUE if num_surf_compl > 0)
        !>   3. Compute num_prim_species considering kinetics
        !>   4. Compute num_sec_species
        !>   5. Compute num_sec_aq_species
        !>   6. Compute num_aq_sec_var_act_species
        !>   
        !>   Simplifies algebraic structure for systems where constant activity species
        !>   (e.g., H₂O, buffers) do not participate in component matrix.
        subroutine elim_cst_act_species(this,num_surf_compl,n_aq,num_min_kin,num_gas_kin)
            implicit none
            !> Eliminate constant activity species from speciation algebra.
            !> This routine is used to simplify the algebraic structure for systems where constant activity species are not relevant, improving computational efficiency and clarity.
            !> Example:
            !>   call spec_algebra%elim_cst_act_species(num_surf_compl, n_aq, num_min_kin, num_gas_kin)
            !>   ! Now constant activity species are removed from the algebraic structure.
            class(speciation_algebra_s) :: this
            integer(kind=4), intent(in) :: num_surf_compl !> number of surface complexes
            integer(kind=4), intent(in) :: n_aq !> number of aqueous species
            integer(kind=4), intent(in), optional :: num_min_kin !> number of mineral kinetics
            integer(kind=4), intent(in), optional :: num_gas_kin !> number of gas kinetics
            
            
            logical :: flag_surf

            call this%set_flag_comp(.true.)
            if (num_surf_compl>0) then
                flag_surf=.true.
            else
                flag_surf=.false.
            end if
            call this%set_flag_cat_exch(flag_surf)
            call this%compute_num_prim_species(num_min_kin,num_gas_kin)
            call this%compute_num_sec_species()
            call this%compute_num_sec_aq_species(n_aq)
            call this%compute_num_aq_sec_var_act_species()
        end subroutine
        
        subroutine set_num_eq_reactions(this,num_eq_reactions)
            implicit none
            !> Set number of equilibrium reactions.
            !> This is used for matrix sizing and equilibrium calculations.
            !> Example:
            !>   call spec_algebra%set_num_eq_reactions(num_eq_reactions)
            !>   ! Now spec_algebra%num_eq_reactions is set.
            class(speciation_algebra_s) :: this
            integer(kind=4), intent(in) :: num_eq_reactions
            if (num_eq_reactions<0) then
                error stop "Number of equilibrium reactions must be non-negative (set_num_eq_reactions)"
            end if
            this%num_eq_reactions=num_eq_reactions
        end subroutine

end module