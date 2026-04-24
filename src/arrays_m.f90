!> \file arrays_m.f90
!> \brief Matrix operations module providing comprehensive linear algebra functionality
!> \details
!>   - Defines abstract and concrete matrix classes for various matrix types
!>   - Implements matrix-vector and matrix-matrix operations
!>   - Provides eigenvalue/eigenvector computation for special matrix types
!>   - Supports diagonal, tridiagonal, symmetric, and Toeplitz matrices
!>   - Includes matrix norms, decompositions (LU, Cholesky), and other linear algebra operations
!> \author Jordi Petchamé-Guerrero
!> \date 2024

module arrays_m
    use vectors_m                                                              !< Import vector classes and operations for matrix-vector interactions
    implicit none                                                              !< Require explicit declaration of all variables (Fortran safety)
    save                                                                       !< Preserve module variable values between procedure calls
    private                                                                    !< Default visibility is private; only specified entities are public
    public :: copy_real_array, id_matrix, norm_mat_inf, forward_substitution, &
        backward_substitution, LU, compute_det, prod_tridiag_diag_mat, prod_diag_tridiag_mat,&
        prod_mat_vec, prod_tridiag_mat_mat, inv_SPD_mat !< Public procedures for real array operations
    !> \brief Abstract base class for all array/matrix types
    !> \details
    !>   - Provides polymorphism for different matrix implementations
    !>   - Defines common interface for matrix operations
    !>   - Cannot be instantiated directly (abstract type)
    type, public, abstract :: array_c
        integer(kind=4) :: num_cols                                            !< Number of columns in the array/matrix structure
    contains
        procedure :: allocate_array                                   !< Allocate memory for array structure
        procedure :: allocate_columns                                 !< Allocate memory for individual columns
        procedure :: prod_mat_vec                                     !< Matrix-vector multiplication: y = A*x
        !procedure :: prod_mat_mat                                    !< [DISABLED] Matrix-matrix multiplication
        procedure :: get_diag                                         !< Extract main diagonal elements
        procedure :: get_sub                                          !< Extract sub-diagonal elements (below main diagonal)
        procedure :: get_super                                        !< Extract super-diagonal elements (above main diagonal)
        procedure :: compute_norm_inf                                 !< Compute infinity norm (maximum row sum)
        procedure :: compute_norm_1                                   !< Compute 1-norm (maximum column sum)
    end type
    
    !> \brief Non-square matrix class extending array_c
    !> \details
    !>   - Represents matrices where m≠n (rows ≠ columns)
    !>   - Placeholder for future expansion
    type, public, extends(array_c) :: non_sq_matrix_c
    end type
    
    !> \brief Real-valued array/matrix implementation
    !> \details
    !>   - Stores matrix as array of real vector columns
    !>   - Column-wise storage for efficient access patterns
    !>   - Supports full matrix operations
    type, public, extends(non_sq_matrix_c) :: real_array_c
        type(vector_real_c), allocatable :: cols(:)                           !< Array of real vector columns (column-wise storage for efficient access)
    contains
        procedure :: is_matrix_full                                   !< Check if matrix is fully populated (no empty elements)
        procedure :: get_can_vec                                      !< Extract canonical vector information (basis vectors)
        procedure :: copy_real_array                                   !< Copy real array information (basis vectors)
        procedure :: compute_mix_ratios_R_opt4                             !< Compute mixing ratios for reaction amounts (optimized version 4)
    end type
    
    !> \brief Integer-valued array/matrix implementation
    !> \details
    !>   - Stores matrix as array of integer vector columns
    !>   - Used for index matrices and integer computations
    type, public, extends(non_sq_matrix_c) :: int_array_c
        type(vector_int_c), allocatable :: cols(:)                            !< Array of integer vector columns for integer matrix operations
    contains
        procedure :: get_vector_int                                   !< Extract integer vector from array structure
    end type
    
    !> \brief Square matrix class with eigenvalue/eigenvector support
    !> \details
    !>   - Base class for all square matrices
    !>   - Provides eigenvalue and eigenvector computation
    type, public, extends(array_c) :: sq_matrix_c
        !integer(kind=4) :: dim                                               !< [DISABLED] Dimension of square matrix
        real(kind=8), allocatable :: eigenvalues(:)                           !< Eigenvalues of the matrix
        real(kind=8), allocatable :: eigenvectors(:,:)                        !< Eigenvectors of the matrix (column-wise)
    contains
        procedure :: compute_eigenvalues                              !< Compute all eigenvalues
        procedure :: compute_eigenvectors                             !< Compute all eigenvectors
    end type
!****************************************************************************************************************************************************
    !> \brief Symmetric tridiagonal Toeplitz matrix class
    !> \details
    !>   - Toeplitz matrix: constant along diagonals
    !>   - Tridiagonal: only main, sub, and super diagonals are non-zero
    !>   - Symmetric: sub-diagonal equals super-diagonal
    !>   - Fully defined by two scalars (sub/super and diag)
    type, public, extends(sq_matrix_c) :: tridiag_sym_Toeplitz_matrix_c
        real(kind=8) :: sub                                                   !< Sub-diagonal element (constant for all positions)
        real(kind=8) :: diag                                                  !< Main diagonal element (constant for all positions)
    end type
    
    !> \brief General tridiagonal Toeplitz matrix class
    !> \details
    !>   - Toeplitz matrix: constant along diagonals
    !>   - Tridiagonal structure with potentially different sub/super diagonals
    !>   - Fully defined by three scalars (sub, diag, super)
    type, public, extends(tridiag_sym_Toeplitz_matrix_c) :: tridiag_Toeplitz_matrix_c
        real(kind=8) :: super                                                 !< Super-diagonal element (constant for all positions)
    contains
        procedure :: set_tridiag_Toeplitz_matrix                     !< Set all three diagonal elements
    end type
!****************************************************************************************************************************************************
    !> \brief Diagonal matrix class
    !> \details
    !>   - Only main diagonal is non-zero
    !>   - Efficient storage: O(n) instead of O(n²)
    !>   - Fast matrix-matrix and matrix-vector products
    type, public, extends(sq_matrix_c) :: diag_matrix_c
        real(kind=8), allocatable :: diag(:)                                  !< Main diagonal elements (size n for n×n matrix)
    contains
        procedure :: set_diag_matrix                                 !< Set diagonal elements from array
        procedure :: prod_mat_diag_mat                               !< Multiply general matrix by diagonal matrix
    end type
    
    !> \brief Symmetric tridiagonal matrix class
    !> \details
    !>   - Three non-zero diagonals: main, sub, super
    !>   - Symmetric: sub(i) = super(i) for all i
    !>   - Efficient storage: O(2n-1) instead of O(n²)
    type, public, extends(diag_matrix_c) :: tridiag_sym_matrix_c
        real(kind=8), allocatable :: sub(:)                                   !< Sub-diagonal elements (size n-1)
    contains
        procedure :: check_eigenvectors_tridiag_sym_matrix           !< Verify eigenvector correctness
    end type

    !> \brief General tridiagonal matrix class
    !> \details
    !>   - Three non-zero diagonals: main, sub, super
    !>   - Sub-diagonal may differ from super-diagonal (not symmetric)
    !>   - Most general tridiagonal representation
    !>   - Efficient storage: O(3n-2) instead of O(n²)
    type, public, extends(tridiag_sym_matrix_c) :: tridiag_matrix_c
        real(kind=8), allocatable :: super(:)                                 !< Super-diagonal elements (size n-1)
    contains
        procedure :: set_tridiag_matrix                              !< Set all three diagonals from arrays
        procedure :: compute_transpose_tridiag_matrix                !< Compute transpose by swapping sub/super
        ! procedure :: compute_inverse_tridiag_matrix                !< [DISABLED] Compute matrix inverse
        procedure :: prod_tridiag_mat_mat                            !< Tridiagonal matrix times general matrix
        procedure :: prod_tridiag_mat_vec                            !< Tridiagonal matrix times general vector
    end type
    !> \brief Pentadiagonal matrix class
    !> \details
    type, public, extends(tridiag_matrix_c) :: pentadiag_matrix_c
        real(kind=8), allocatable :: sub2(:)                         !< Radial positions associated with tridiagonal entries
        real(kind=8), allocatable :: super2(:)                       !< Radial positions associated with tridiagonal entries
    contains
        procedure :: set_pentadiag_matrix                           !< Set all five diagonals from arrays
        procedure :: fill_pentadiag_mat                         !< Pentadiagonal matrix times general vector
    end type
    !> \brief [DISABLED] Tridiagonal matrix with associated vector
    !> \details
    !>   - Would extend tridiagonal matrix with additional vector storage
    !>   - Currently commented out - placeholder for future development
    !type, public, extends(tridiag_matrix_c) :: tridiag_matrix_vec_c
    !>    real(kind=8), allocatable :: vector(:)                              !< Associated vector
    !end type
!****************************************************************************************************************************************************

!****************************************************************************************************************************************************
    contains
        !> \brief Allocate memory for array/matrix structure
        !> \details
        !>   - Allocates appropriate storage based on concrete type
        !>   - For diagonal matrices: allocates diag array
        !>   - For tridiagonal matrices: allocates diag, sub, and super arrays
        !>   - For real/int arrays: allocates column vector array
        !>   - Uses Fortran SELECT TYPE for polymorphic dispatch
        !> \param[inout] this Array/matrix object to allocate
        !> \param[in] n Optional number of columns (if not present, uses existing value)
        subroutine allocate_array(this,n)
            implicit none
            class(array_c) :: this                                            !< Array/matrix object
            integer(kind=4), intent(in), optional :: n                        !< Optional column count
            if (present(n)) then                                              !< Check if n was provided
                this%num_cols=n                                               !< Set number of columns
            end if
            select type (this)                                                !< Polymorphic type selection
            class is (diag_matrix_c)                                          !< Diagonal matrix branch
                allocate(this%diag(this%num_cols))                            !< Allocate diagonal array (size n)
                this%diag=0.0d0                                                   !< Initialize diagonal elements to zero
                select type (this)                                            !< Nested type selection for subclasses
                class is (tridiag_sym_matrix_c)                               !< Symmetric tridiagonal branch
                    allocate(this%sub(this%num_cols-1))                       !< Allocate sub-diagonal (size n-1)
                    this%sub=0.0d0                                                   !< Initialize sub-diagonal elements to zero
                    select type (this)                                        !< Further nested selection
                    class is (tridiag_matrix_c)                                !< General tridiagonal branch
                        allocate(this%super(this%num_cols-1))                 !< Allocate super-diagonal (size n-1)
                        this%super=0.0d0                                                   !< Initialize super-diagonal elements to zero
                        select type (this)                                        !< [DISABLED] Further nesting
                        class is (pentadiag_matrix_c)                       !< [DISABLED] Tridiag with vector
                            allocate(this%sub2(this%num_cols-1))                        !< [DISABLED] Allocate associated vector
                            allocate(this%super2(this%num_cols-1))                      !< [DISABLED] Allocate associated vector
                            this%sub2=0.0d0                                                   !< [DISABLED] Initialize associated vector to zero
                            this%super2=0.0d0                                                   !< [DISABLED] Initialize
                        end select                                           !< [DISABLED] End innermost select
                    end select                                                !< End tridiag_matrix_c select
                end select                                                    !< End tridiag_sym_matrix_c select
            type is (real_array_c)                                            !< Real array branch
                allocate(this%cols(this%num_cols))                            !< Allocate array of real column vectors
                !do i=1,this%num_cols                                         !< [DISABLED] Loop to allocate each column
                !    call this%cols(i)%allocate_vector()                      !< [DISABLED] Allocate individual vectors
            type is (int_array_c)                                             !< Integer array branch
                if (allocated(this%cols)) then                                !< Check if already allocated
                    deallocate(this%cols)                                     !< Deallocate existing columns
                end if
                allocate(this%cols(this%num_cols))                            !< Allocate array of integer column vectors
            end select                                                        !< End type selection
        end subroutine
        
        !> \brief Allocate memory for individual columns in array structure
        !> \details
        !>   - Iterates through all columns and allocates each vector
        !>   - Only applies to real_array_c and int_array_c types
        !>   - Must be called after allocate_array
        !> \param[inout] this Array object with columns to allocate
        subroutine allocate_columns(this)
            implicit none
            class(array_c) :: this                                            !< Array object
            integer(kind=4) :: i                                              !< Loop index for columns
            select type (this)                                                !< Polymorphic type selection
            type is (real_array_c)                                            !< Real array branch
                do i=1,this%num_cols                                          !< Loop over all columns
                    call this%cols(i)%allocate_vector()                       !< Allocate each real column vector
                end do
            type is (int_array_c)                                             !< Integer array branch
                do i=1,this%num_cols                                          !< Loop over all columns
                    call this%cols(i)%allocate_vector()                       !< Allocate each integer column vector
                end do
            end select                                                        !< End type selection
        end subroutine
        
        !> \brief Set diagonal matrix elements
        !> \details
        !>   - Copies diagonal elements from input array
        !>   - Automatically sets num_cols from array size
        !> \param[inout] this Diagonal matrix object
        !> \param[in] diag Array of diagonal elements
        subroutine set_diag_matrix(this,diag)
            implicit none
            class(diag_matrix_c) :: this                                      !< Diagonal matrix object
            real(kind=8), intent(in) :: diag(:)                               !< Diagonal elements to set
            this%diag=diag                                                    !< Copy diagonal elements
            this%num_cols=size(diag)                                          !< Set matrix dimension from array size
        end subroutine
        
        !> \brief Set tridiagonal Toeplitz matrix elements
        !> \details
        !>   - Toeplitz: constant values along diagonals
        !>   - Sets three scalar values for sub, diag, super
        !> \param[inout] this Tridiagonal Toeplitz matrix object
        !> \param[in] sub Sub-diagonal constant value
        !> \param[in] diag Main diagonal constant value
        !> \param[in] super Super-diagonal constant value
        subroutine set_tridiag_Toeplitz_matrix(this,sub,diag,super)
            implicit none
            class(tridiag_Toeplitz_matrix_c) :: this                          !< Tridiag Toeplitz matrix object
            real(kind=8), intent(in) :: sub,diag,super                        !< Diagonal constant values
            this%sub=sub                                                      !< Set sub-diagonal constant
            this%diag=diag                                                    !< Set main diagonal constant
            this%super=super                                                  !< Set super-diagonal constant
        end subroutine
        
        !> \brief Set general tridiagonal matrix elements
        !> \details
        !>   - Sets all three diagonals from arrays
        !>   - Validates dimension compatibility
        !>   - Requires: size(sub) = size(super) = size(diag) - 1
        !> \param[inout] this Tridiagonal matrix object
        !> \param[in] sub Sub-diagonal elements (size n-1)
        !> \param[in] diag Main diagonal elements (size n)
        !> \param[in] super Super-diagonal elements (size n-1)
        subroutine set_tridiag_matrix(this,sub,diag,super)
            implicit none
            class(tridiag_matrix_c) :: this                                   !< Tridiagonal matrix object
            real(kind=8), intent(in) :: sub(:),diag(:),super(:)               !< Diagonal arrays
            if (size(sub)/=size(super) .or. size(diag)/=(size(sub)+1)) error stop "Dimension error in set_tridiag_matrix"  !< Validate dimensions
            this%sub=sub                                                      !< Copy sub-diagonal elements
            this%diag=diag                                                    !< Copy main diagonal elements
            this%super=super                                                  !< Copy super-diagonal elements
        end subroutine

        subroutine set_pentadiag_matrix(this,sub2,sub,diag,super,super2)
            implicit none
            class(pentadiag_matrix_c) :: this                                   !< Pentadiagonal matrix object
            real(kind=8), intent(in) :: sub2(:),sub(:),diag(:),super(:),super2(:)               !< Diagonal arrays
            if (size(sub)/=size(super) .or. size(diag)/=(size(sub)+1)) then
                error stop "Dimension error in set_pentadiag_matrix"  !< Validate dimensions
            else if (size(sub2)/=size(super2) .or. size(sub)/=(size(sub2)+1)) then
                error stop "Dimension error in set_pentadiag_matrix"
            end if
            this%sub2=sub2                                                      !< Copy second sub-diagonal elements
            this%sub=sub                                                      !< Copy sub-diagonal elements
            this%diag=diag                                                    !< Copy main diagonal elements
            this%super=super                                                  !< Copy super-diagonal elements
            this%super2=super2                                                  !< Copy second super-diagonal elements
        end subroutine
        
        !> \brief Extract main diagonal elements from matrix
        !> \details
        !>   - Only works for diagonal and tridiagonal matrix types
        !>   - Returns allocatable array
        !> \param[in] this Matrix object
        !> \return diag Allocatable array of diagonal elements
        function get_diag(this) result(diag)
            implicit none
            class(array_c), intent(in) :: this                                !< Matrix object
            real(kind=8), allocatable :: diag(:)                              !< Diagonal elements
            select type (this)                                                !< Polymorphic type selection
            class is (diag_matrix_c)                                          !< Diagonal matrix branch
                diag=this%diag                                                !< Copy diagonal elements
            end select                                                        !< End type selection
        end function
        
        !> \brief Extract sub-diagonal elements from matrix
        !> \details
        !>   - Only works for symmetric tridiagonal and tridiagonal matrices
        !>   - Returns allocatable array of size n-1
        !> \param[in] this Matrix object
        !> \return sub Allocatable array of sub-diagonal elements
        function get_sub(this) result(sub)
            implicit none
            class(array_c), intent(in) :: this                                !< Matrix object
            real(kind=8), allocatable :: sub(:)                               !< Sub-diagonal elements
            select type (this)                                                !< Polymorphic type selection
            class is (tridiag_sym_matrix_c)                                   !< Symmetric tridiag branch
                sub=this%sub                                                  !< Copy sub-diagonal elements
            end select                                                        !< End type selection
        end function
        
        !> \brief Extract super-diagonal elements from matrix
        !> \details
        !>   - Only works for general tridiagonal matrices
        !>   - Returns allocatable array of size n-1
        !> \param[in] this Matrix object
        !> \return super Allocatable array of super-diagonal elements
        function get_super(this) result(super)
            implicit none
            class(array_c), intent(in) :: this                                !< Matrix object
            real(kind=8), allocatable :: super(:)                             !< Super-diagonal elements
            select type (this)                                                !< Polymorphic type selection
            class is (tridiag_matrix_c)                                       !< Tridiagonal matrix branch
                super=this%super                                              !< Copy super-diagonal elements
            end select                                                        !< End type selection
        end function
        
        !> \brief Compute infinity norm (maximum absolute row sum) of matrix
        !> \details
        !>   - ||A||_∞ = max_i Σ_j |a_ij|
        !>   - For diagonal matrix: max of absolute diagonal values
        !>   - For tridiagonal: max of row sums (at most 3 elements per row)
        !>   - Efficient O(n) computation for sparse matrices
        !> \param[in] this Matrix object
        !> \return norm Infinity norm value
        function compute_norm_inf(this) result(norm)
            implicit none
            class(array_c), intent(in) :: this                                !< Matrix object
            real(kind=8) :: norm                                              !< Infinity norm result
            
            integer(kind=4) :: i,n                                            !< Loop index and matrix dimension
            real(kind=8) :: norm_i                                            !< Row sum for current row
            
            select type (this)                                                !< Polymorphic type selection
            type is (diag_matrix_c)                                           !< Diagonal matrix branch
                norm=maxval(abs(this%diag))                                   !< Maximum absolute diagonal value
            type is (tridiag_matrix_c)                                        !< Tridiagonal matrix branch
                n=size(this%diag)                                             !< Get matrix dimension
                norm=abs(this%diag(1))+abs(this%super(1))                     !< First row (diag + super only)
                do i=2,n-1                                                    !< Loop over interior rows
                    norm_i=abs(this%sub(i-1))+abs(this%diag(i))+abs(this%super(i))  !< Row sum (sub + diag + super)
                    if (norm_i>norm) then                                     !< Check if larger than current max
                        norm=norm_i                                           !< Update maximum
                    end if
                end do
                norm_i=abs(this%sub(n-1))+abs(this%diag(n))                   !< Last row (sub + diag only)
                if (norm_i>norm) then                                         !< Check if larger than current max
                    norm=norm_i                                               !< Update maximum
                end if
            end select                                                        !< End type selection
        end function
        
        !> \brief Compute 1-norm (maximum absolute column sum) of matrix
        !> \details
        !>   - ||A||_1 = max_j Σ_i |a_ij|
        !>   - For tridiagonal: computes via transpose and infinity norm
        !>   - Uses property: ||A||_1 = ||A^T||_∞
        !> \param[in] this Matrix object
        !> \return norm 1-norm value
        function compute_norm_1(this) result(norm)
            implicit none
            class(array_c), intent(in) :: this                                !< Matrix object
            real(kind=8) :: norm                                              !< 1-norm result
            
            type(tridiag_matrix_c) :: transpose                               !< Transpose matrix
            
            select type (this)                                                !< Polymorphic type selection
            type is (tridiag_matrix_c)                                        !< Tridiagonal matrix branch
                call this%compute_transpose_tridiag_matrix(transpose)         !< Compute transpose
                norm=transpose%compute_norm_inf()                             !< Compute infinity norm of transpose
            end select                                                        !< End type selection
        end function
        
        !> \brief Compute transpose of tridiagonal matrix
        !> \details
        !>   - For tridiagonal: transpose swaps sub and super diagonals
        !>   - Diagonal remains unchanged
        !>   - Efficient O(1) operation (just pointer/reference swap)
        !> \param[in] this Original tridiagonal matrix
        !> \param[out] transpose Transposed tridiagonal matrix
        subroutine compute_transpose_tridiag_matrix(this,transpose)
            implicit none
            class(tridiag_matrix_c), intent(in) :: this                       !< Original matrix
            class(tridiag_matrix_c), intent(out) :: transpose                 !< Transpose matrix
            !> Note: Would be better to allocate transpose here               !< TODO: Allocate transpose
            transpose%sub=this%super                                          !< Swap: sub ← super
            transpose%diag=this%diag                                          !< Diagonal unchanged
            transpose%super=this%sub                                          !< Swap: super ← sub
        end subroutine
        
        !> \brief Extract submatrix from general matrix
        !> \details
        !>   - Extracts specified rows and columns
        !>   - Output submatrix must be pre-allocated
        !>   - Useful for block matrix operations
        !> \param[in] mat Input matrix
        !> \param[in] rows Row indices to extract
        !> \param[in] cols Column indices to extract
        !> \param[out] submat Extracted submatrix (must be allocated)
        subroutine get_submatrix(mat,rows,cols,submat)
            implicit none
            real(kind=8), intent(in) :: mat(:,:)                              !< Input matrix
            integer(kind=4), intent(in) :: rows(:),cols(:)                    !< Row and column indices
            real(kind=8), intent(out) :: submat(:,:)                          !< Output submatrix (must be allocated)
            
            integer(kind=4) :: i,j,n_r,n_c,i_r,i_c                            !< Loop indices and dimensions
            n_r=size(rows)                                                    !< Number of rows to extract
            n_c=size(cols)                                                    !< Number of columns to extract
            do i=1,n_r                                                        !< Loop over target rows
                do j=1,n_c                                                    !< Loop over target columns
                    submat(i,j)=mat(rows(i),cols(j))                          !< Copy element
                end do
            end do
        end subroutine
        
        !> \brief Multiply general matrix by diagonal matrix: C = A * D
        !> \details
        !>   - A is general m×n matrix, D is diagonal n×n
        !>   - Result is m×n matrix where each column of A is scaled
        !>   - Efficient: O(mn) instead of O(mn²) for general matrix product
        !> \param[in] this Diagonal matrix object
        !> \param[in] mat General matrix to multiply
        !> \return prod Result of matrix-diagonal product
        function prod_mat_diag_mat(this,mat) result(prod)
            implicit none
            class(diag_matrix_c), intent(in) :: this                          !< Diagonal matrix
            real(kind=8), intent(in) :: mat(:,:)                              !< General matrix
            real(kind=8), allocatable :: prod(:,:)                            !< Product matrix
            integer(kind=4) :: j                                              !< Column loop index
            if (size(mat,2)/=this%num_cols) then                              !< Validate dimensions
                error stop "Dimension error in prod_mat_diag_mat"
            end if
            allocate(prod(size(mat,1),this%num_cols))                         !< Allocate result
            do j=1,this%num_cols                                              !< Loop over columns
                prod(:,j)=mat(:,j)*this%diag(j)                               !< Scale column j by diagonal element j
            end do
        end function
        
        function get_vector_int(this) result(vector_int)
            implicit none
            class(int_array_c), intent(in) :: this
            integer(kind=4), allocatable :: vector_int(:)
            
            integer(kind=4) :: i,j,dim_vec
            
            dim_vec=0
            do i=1,this%num_cols
                dim_vec=dim_vec+this%cols(i)%dim
            end do
            allocate(vector_int(dim_vec))
            
            j=1
            do i=1,this%num_cols
                vector_int(j:j+this%cols(i)%dim-1)=this%cols(i)%col_1
                j=j+this%cols(i)%dim
            end do
        end function
        
        subroutine is_matrix_full(this,flag)
            implicit none
            class(real_array_c), intent(in) :: this
            logical, intent(out) :: flag
            
            integer(kind=4) :: i
            
            flag=.true.
            do i=1,this%num_cols
                if (this%cols(i)%dim<this%num_cols) then
                    flag=.false.
                    exit
                end if
            end do
        end subroutine

        !> \brief Compute spectral radius of matrix from eigenvalues
        !> \details
        !>   - Spectral radius = maximum absolute eigenvalue
        !>   - ρ(A) = max|λᵢ|
        !>   - Important for iterative method convergence analysis
        !> \param[in] lambda Vector of eigenvalues (real or complex parts)
        !> \return rho Spectral radius (largest absolute eigenvalue)
function radio_espectral(lambda) result(rho)
    implicit none
    real(kind=8), intent(in) :: lambda(:)                                     !< Eigenvalue vector
    real(kind=8) :: rho                                                       !< Spectral radius
    rho=max(maxval(lambda),maxval(-lambda))                                   !< Maximum of positive and negative extrema
end function 

        !> \brief Multiply tridiagonal matrix by general matrix: C = T * X
        !> \details
        !>   - T is tridiagonal n×n, X is general n×m
        !>   - Result C is general n×m matrix
        !>   - Efficient O(nm) algorithm using tridiagonal structure
        !>   - Each column of C computed independently
        !>   - Uses vectorized operations where possible
        !> \param[in] this Tridiagonal matrix T
        !> \param[in] X_mat General matrix X (n×m)
        !> \return C_mat Product matrix C = TX
function prod_tridiag_mat_mat(this,X_mat) result(C_mat)
    implicit none
    class(tridiag_matrix_c), intent(in) :: this                               !< Tridiagonal matrix n×n
    real(kind=8), intent(in) :: X_mat(:,:)                                    !< General matrix n×m
    real(kind=8), allocatable :: C_mat(:,:)                                   !< Product matrix n×m
    
    integer(kind=4) :: i,j,n,m                                                !< Loop indices and dimensions
    
    m=size(X_mat,2)                                                           !< Number of columns in X
    !select type (this)
    !type is (tridiag_matrix_c)
        n=size(this%diag)                                                     !< Matrix dimension
        if (n/=size(X_mat,1)) error stop "Dimension error in prod_mat_mat"   !< Validate compatibility
        if (.not. allocated(C_mat)) then                                      !< Allocate result if needed
            allocate(C_mat(n,m))
        end if
        !> First row: diagonal and superdiagonal only
        C_mat(1,:)=this%diag(1)*X_mat(1,:)+this%super(1)*X_mat(2,:)
        !> Interior rows: subdiagonal, diagonal, superdiagonal (column-wise for cache efficiency)
        do i=2,n-1
            do j=1,m
                C_mat(i,j)=this%sub(i-1)*X_mat(i-1,j)+this%diag(i)*X_mat(i,j)+this%super(i)*X_mat(i+1,j)
            end do
        end do
        !> Last row: subdiagonal and diagonal only
        C_mat(n,:)=this%sub(n-1)*X_mat(n-1,:)+this%diag(n)*X_mat(n,:)
    !end select
end function 

        !> \brief Multiply tridiagonal matrix by diagonal matrix: C = T * D
        !> \details
        !>   - T is tridiagonal n×n, D is diagonal n×n
        !>   - Result C is tridiagonal n×n
        !>   - Efficient: preserves tridiagonal structure
        !>   - Each band multiplied by corresponding diagonal elements
        !> \param[in] A Tridiagonal matrix
        !> \param[in] B Diagonal matrix
        !> \return C Tridiagonal product matrix C = A*B
function prod_tridiag_diag_mat(A,B) result(C)
    class(tridiag_matrix_c), intent(in) :: A                                  !< Tridiagonal matrix
    class(diag_matrix_c), intent(in) :: B                                     !< Diagonal matrix
    type(tridiag_matrix_c) :: C                                               !< Tridiagonal result

    integer(kind=4) :: i,n                                                    !< Loop index and dimension
    
    n=size(B%diag)                                                            !< Matrix dimension
    
    call C%allocate_array(n)                                                  !< Allocate tridiagonal structure
        
    C%diag(1)=A%diag(1)*B%diag(1)                                             !< First diagonal element
    
    if (n>1) then                                                             !< Boundary elements
        C%super(1)=A%super(1)*B%diag(2)                                       !< First superdiagonal
        C%sub(n-1)=A%sub(n-1)*B%diag(n-1)                                     !< Last subdiagonal
        C%diag(n)=A%diag(n)*B%diag(n)                                         !< Last diagonal
    end if
    if (n>2) then                                                             !< Interior elements
        do i=2,n-1
            C%diag(i)=A%diag(i)*B%diag(i)                                     !< Interior diagonal
            C%super(i)=A%super(i)*B%diag(i+1)                                 !< Interior superdiagonal
            C%sub(i-1)=A%sub(i-1)*B%diag(i-1)                                 !< Interior subdiagonal
        end do
    end if
end function

        !> \brief Compute analytical solution to ODE system: dy/dt = -Ay + b
        !> \details
        !>   - Solves first-order linear ODE system with constant coefficients
        !>   - Uses eigenvalue decomposition: y(t) = exp(-At)y₀ + integral
        !>   - For each eigenvalue λ: contribution is exp(-λt) and (1-exp(-λt))/λ
        !>   - Handles zero eigenvalues with limit case
        !>   - Requires precomputed eigenvalues and eigenvectors
        !> \param[in] A Coefficient matrix (must have eigenvalues/eigenvectors)
        !> \param[in] y0 Initial condition vector
        !> \param[in] b Source/forcing term vector
        !> \param[in] time Time at which to evaluate solution
        !> \return y Solution vector at specified time
function prod_total_sym_mat(A,y0,b,time) result(y)
        implicit none
        class(sq_matrix_c), intent(in) :: A                                   !< Coefficient matrix
        real(kind=8), intent(in) :: y0(:)                                     !< Initial condition
        real(kind=8), intent(in) :: b(:)                                      !< Forcing term
        real(kind=8), intent(in) :: time                                      !< Evaluation time
        real(kind=8), allocatable :: y(:)                                     !< Solution at time t
        
        real(kind=8), allocatable :: Q_lambda(:,:),x(:),Pt_x(:)               !< Intermediate computation arrays
        integer(kind=4) :: i,n,time_step                                      !< Loop index, dimension, time step
        real(kind=8), parameter :: epsilon=1d-6                               !< Zero eigenvalue threshold
        
        n=size(A%eigenvalues)                                                 !< System dimension
        allocate(Q_lambda(2,n),x(2*n),Pt_x(2*n))                              !< Allocate work arrays
        !> Compute exponential and integral factors for each eigenvalue
        do i=1,n
            if (abs(A%eigenvalues(i))<epsilon) then                       !< Handle zero eigenvalue limit
                Q_lambda(2,i)=0d0                                         !< lim (1-exp(-λt))/λ as λ→0 = 0
            else
                Q_lambda(2,i)=(1d0-exp(-A%eigenvalues(i)*time))/A%eigenvalues(i)  !< Integral factor
            end if
        end do
        Q_lambda(1,:)=exp(-A%eigenvalues*time)                            !< Exponential decay factors
        x(1:n)=y0                                                             !< Pack initial condition
        x(n+1:2*n)=b                                                          !< Pack forcing term
        Pt_x(1:n)=matmul(transpose(A%eigenvectors),x(1:n))                    !< Project y0 onto eigenbasis
        Pt_x(n+1:2*n)=matmul(transpose(A%eigenvectors),x(n+1:2*n))            !< Project b onto eigenbasis
        !> Combine homogeneous and particular solutions
        y=matmul(A%eigenvectors,Q_lambda(1,:)*Pt_x(1:n))+matmul(A%eigenvectors,Q_lambda(2,:)*Pt_x(n+1:2*n))
    end function

        !> \brief Multiply general matrix by vector: x = A * b
        !> \details
        !>   - Polymorphic: handles different matrix types
        !>   - Commented code shows tridiag_matrix_vec_c implementation
        !>   - Currently uses SELECT TYPE for type-specific operations
        !> \param[in] A Matrix (various types supported)
        !> \param[in] b Input vector
        !> \return x Result vector x = Ab
    function prod_mat_vec(A,b) result(x)
        implicit none
        class(array_c), intent(in) :: A                                       !< Matrix (polymorphic)
        real(kind=8), intent(in) :: b(:)                                      !< Input vector
        real(kind=8), allocatable :: x(:)                                     !< Result vector
        
        integer(kind=4) :: i,j,n,m                                            !< Loop indices and dimensions
        real(kind=8) :: sum                                                   !< Accumulator for dot product
        
        select type (A)
        !type is (tridiag_matrix_vec_c)
        !>    n=size(A%diag)
        !>    if (size(b)/=2*n) error stop "Dimension error in b"
        !>    if (n.eqv.1) then
        !>        x=A%diag(1)*b(1)+A%vector(1)*b(n+1)
        !>    else if (n.eqv.2) then
        !>        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)+A%vector(1)*b(n+1)
        !>        x(2)=A%sub(1)*b(1)+A%diag(2)*b(2)+A%vector(2)*b(n+2)
        !>    else
        !>        x(1)=A%diag(1)*b(1)+A%super(1)*b(2)+A%vector(1)*b(n+1)
        !>        do i=2,n-1
        !>            x(i)=A%sub(i-1)*b(i-1)+A%diag(i)*b(i)+A%super(i)*b(i+1)
        !>            x(i)=x(i)+A%vector(i)*b(n+i)
        !>        end do
        !>        x(n)=A%sub(n-1)*b(n-1)+A%diag(n)*b(n)+A%vector(n)*b(2*n)
        !>    end if
        type is (tridiag_matrix_c)
            n=size(A%diag)
            allocate(x(n))
            if (size(b)/=n) error stop "Dimension error in b"
            x(1)=A%diag(1)*b(1)+A%super(1)*b(2)
            do i=2,n-1
                x(i)=A%sub(i-1)*b(i-1)+A%diag(i)*b(i)+A%super(i)*b(i+1)
            end do
            x(n)=A%sub(n-1)*b(n-1)+A%diag(n)*b(n)                             !< Last element
        end select
    end function 

        !> \brief Multiply general matrix by matrix: C = A * X (polymorphic)
        !> \details
        !>   - Polymorphic implementation for different matrix types
        !>   - Tridiagonal case: efficient O(nm) algorithm
        !>   - General case would use BLAS (commented out)
        !> \param[in] this Matrix A (polymorphic - tridiag_matrix_c shown)
        !> \param[in] X_mat Matrix X (n×m)
        !> \return C_mat Product matrix C = AX
    function prod_mat_mat(this,X_mat) result(C_mat)
        implicit none
        class(array_c), intent(in) :: this                                    !< Matrix A (n×n, polymorphic)
        real(kind=8), intent(in) :: X_mat(:,:)                                !< Matrix X (n×m)
        real(kind=8), allocatable :: C_mat(:,:)                               !< Product C (n×m)
        
        integer(kind=4) :: i,j,n,m                                            !< Loop indices and dimensions
        
        m=size(X_mat,2)                                                       !< Number of columns in X
        select type (this)
        type is (tridiag_matrix_c)                                            !< Tridiagonal case
            n=size(this%diag)                                                 !< Matrix dimension
            if (n/=size(X_mat,1)) error stop "Dimension error in prod_mat_mat"  !< Validate compatibility
            if (.not. allocated(C_mat)) then                                  !< Allocate if needed
                allocate(C_mat(n,m))
            end if
            !> First row: diagonal and superdiagonal
            C_mat(1,:)=this%diag(1)*X_mat(1,:)+this%super(1)*X_mat(2,:)
            !> Interior rows: all three bands
            do i=2,n-1
                do j=1,m
                    C_mat(i,j)=this%sub(i-1)*X_mat(i-1,j)+this%diag(i)*X_mat(i,j)+this%super(i)*X_mat(i+1,j)
                end do
            end do
            !> Last row: subdiagonal and diagonal
            C_mat(n,:)=this%sub(n-1)*X_mat(n-1,:)+this%diag(n)*X_mat(n,:)
        end select
    end function 

        !> \brief Multiply diagonal matrix by tridiagonal matrix: C = D * T
        !> \details
        !>   - D is diagonal n×n, T is tridiagonal n×n
        !>   - Result C is tridiagonal (preserves structure)
        !>   - Each row i of T scaled by diagonal element D(i,i)
        !>   - Efficient: O(n) operations
        !> \param[in] A Diagonal matrix D
        !> \param[in] B Tridiagonal matrix T
        !> \return C Tridiagonal product C = DT
    function prod_diag_tridiag_mat(A,B) result(C)
        class(diag_matrix_c), intent(in) :: A                                 !< Diagonal matrix
        class(tridiag_matrix_c), intent(in) :: B                              !< Tridiagonal matrix
        type(tridiag_matrix_c) :: C                                           !< Tridiagonal result

        integer(kind=4) :: i,n                                                !< Loop index and dimension
        
        n=size(A%diag)                                                        !< Matrix dimension
        
        call C%allocate_array(n)                                              !< Allocate tridiagonal structure
        
        C%diag(1)=A%diag(1)*B%diag(1)                                             !< First diagonal
        C%super(1)=A%diag(1)*B%super(1)                                       !< First superdiagonal (row 1 scaled by D(1,1))
        
        C%sub(n-1)=A%diag(n)*B%sub(n-1)                                       !< Last subdiagonal (row n scaled by D(n,n))
        C%diag(n)=A%diag(n)*B%diag(n)                                         !< Last diagonal
        if (n>2) then                                                         !< Interior elements
            do i=2,n-1
                C%diag(i)=A%diag(i)*B%diag(i)                                 !< Interior diagonal
                C%super(i)=A%diag(i)*B%super(i)                               !< Interior superdiagonal
                C%sub(i-1)=A%diag(i)*B%sub(i-1)                               !< Interior subdiagonal
            end do
        end if
        
    end function

        !> \brief Compute Gram matrix: AᵀA
        !> \details
        !>   - Forms symmetric positive semi-definite matrix
        !>   - Result is m×m where A is n×m
        !>   - Used in least squares, normal equations
        !>   - Entry (i,j) = dot product of columns i and j of A
        !> \param[in] A Input matrix (n×m)
        !> \return AT_A Gram matrix AᵀA (m×m, symmetric)
    function prod_AT_A(A) result(AT_A)
        implicit none
        real(kind=8), intent(in) :: A(:,:)                                    !< Input matrix n×m
        real(kind=8), allocatable :: AT_A(:,:)                                !< Gram matrix m×m
        
        integer(kind=4) :: n,m,i,j                                            !< Dimensions and loop indices
        
        n=size(A,1)                                                           !< Number of rows in A
        m=size(A,2)                                                           !< Number of columns in A
        
        allocate(AT_A(m,m))                                                   !< Allocate result
    
        do j=1,m                                                              !< Loop over columns of A
            do i=1,j                                                          !< Exploit symmetry: compute upper triangle
                AT_A(i,j)=dot_product(A(:,i),A(:,j))                          !< (AᵀA)ᵢⱼ = column i · column j
                AT_A(j,i)=AT_A(i,j)                                           !< Mirror to lower triangle
            end do
        end do
    end function

        !> \brief Power iteration method for computing spectral radius
        !> \details
        !>   - Iterative method: zₖ₊₁ = Azₖ, ρₖ = ||Azₖ||/||zₖ||
        !>   - Converges to largest eigenvalue (spectral radius)
        !>   - Updates initial vector z0 in place
        !>   - Stops when |ρₖ - ρₖ₋₁| < tolerance or max iterations reached
        !> \param[in] A Matrix to analyze
        !> \param[inout] z0 Initial/final eigenvector estimate
        !> \param[in] tolerance Convergence tolerance
        !> \param[out] rho Computed spectral radius
        !> \param[out] niter Number of iterations performed
    subroutine potencia(A,z0,tolerance,rho,niter)
            implicit none
            real(kind=8), intent(in) :: A(:,:)                                !< Input matrix
            real(kind=8), intent(inout) :: z0(:)                              !< Initial/final eigenvector
            real(kind=8), intent(in) :: tolerance                             !< Convergence tolerance
            real(kind=8), intent(out) :: rho                                  !< Spectral radius (dominant eigenvalue)
            integer(kind=4), intent(out) :: niter                             !< Number of iterations
            
            integer(kind=4) :: n,k                                            !< Dimension and iteration counter
            integer(kind=4), parameter :: iter_max=100                        !< Maximum iterations
            real(kind=8) :: rho_old,rho_new                                   !< Successive estimates
            real(kind=8), allocatable :: zk(:)                                !< Work vector
            n=size(z0)                                                        !< Vector dimension
            allocate(zk(n))                                                   !< Allocate work space
            niter=0                                                           !< Initialize iteration counter
            rho_old=0d0                                                       !< Initialize previous estimate
            
            do                                                                !< Power iteration loop
                niter=niter+1                                                 !< Increment iteration counter
                zk=matmul(A,z0)                                               !< Compute zₖ = Azₖ₋₁
                rho_new=p_norm_vec(zk,2)/p_norm_vec(z0,2)                    !< Estimate ρ = ||Azₖ||/||zₖ₋₁||
                if (abs(rho_new-rho_old)<tolerance .or. niter.eq.iter_max) exit  !< Check convergence or max iters
                z0=zk                                                         !< Update vector for next iteration
                rho_old=rho_new                                               !< Update previous estimate
            end do
            rho=rho_new                                                       !< Return final estimate
        end subroutine

        !> \brief Compute infinity norm (max row sum) of matrix
        !> \details
        !>   - ||A||∞ = maxᵢ Σⱼ |aᵢⱼ|
        !>   - Maximum absolute row sum
        !>   - Induced matrix norm for p=∞ vector norm
        !> \param[in] A Input matrix (m×n)
        !> \return norm Infinity norm of A
function norm_mat_inf(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:)                                        !< Input matrix m×n
    real(kind=8) :: norm                                                      !< Infinity norm
            
    integer(kind=4) :: i,j,n,m                                                !< Loop indices and dimensions
    real(kind=8) :: sum                                                       !< Row sum accumulator
    
    m=size(A,1)
    n=size(A,2)
    sum=0d0
    do i=1,n
        sum=sum+abs(A(1,i))
    end do
    norm=sum
    do i=2,m
        sum=0d0
        do j=1,n
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function

!> Computes $L_1$ norm of a matrix
function norm_mat_1(A) result(norm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> m x n matrix
    real(kind=8) :: norm
            
    integer(kind=4) :: i,j,n,m
    real(kind=8) :: sum
    m=size(A,1)
    n=size(A,2)
    sum=0d0
    do i=1,m
        sum=sum+abs(A(i,1))
    end do
    norm=sum
    do j=2,n
        sum=0d0
        do i=1,m
            sum=sum+abs(A(i,j))
        end do
        if (sum>norm) then
            norm=sum
        end if
    end do
end function

!> LU decomposition
subroutine LU(A,tol,L,U,error,num_swaps,perm)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: tol !> tolerance
    real(kind=8), intent(out) :: L(:,:), U(:,:)
    logical, intent(out) :: error
    integer(kind=4), intent(out), optional :: num_swaps !> number of row swaps
    integer(kind=4), intent(out), optional :: perm(:) !> row permutation vector
    
    integer(kind=4) :: n,i,j,pivot_row,nswaps,tmp_int
    real(kind=8) :: factor, max_val
    real(kind=8), allocatable :: tmp_row(:)
    integer(kind=4), allocatable :: p(:)
    n=size(A,1)
    U=A
    L=0d0
    error=.false.
    nswaps=0
    allocate(tmp_row(n),p(n))
    do i=1,n
        p(i)=i
    end do
    do j=1,n
        !> Partial pivoting: find row with largest absolute value in column j
        max_val=abs(U(j,j))
        pivot_row=j
        do i=j+1,n
            if (abs(U(i,j))>max_val) then
                max_val=abs(U(i,j))
                pivot_row=i
            end if
        end do
        if (pivot_row/=j) then
            nswaps=nswaps+1
            !> Swap rows in U
            tmp_row(j:n)=U(j,j:n)
            U(j,j:n)=U(pivot_row,j:n)
            U(pivot_row,j:n)=tmp_row(j:n)
            !> Swap rows in L (only columns 1:j-1)
            if (j>1) then
                tmp_row(1:j-1)=L(j,1:j-1)
                L(j,1:j-1)=L(pivot_row,1:j-1)
                L(pivot_row,1:j-1)=tmp_row(1:j-1)
            end if
            !> Track permutation
            tmp_int=p(j)
            p(j)=p(pivot_row)
            p(pivot_row)=tmp_int
        end if
        if (abs(U(j,j))<tol) then
            error=.true.
        end if
        L(j,j)=1d0
        do i=j+1,n
            factor=U(i,j)/U(j,j)
            L(i,j)=factor
            U(i,j:n)=U(i,j:n)-factor*U(j,j:n)
        end do
    end do
    if (present(num_swaps)) num_swaps=nswaps
    if (present(perm)) perm=p
    deallocate(tmp_row,p)
end subroutine LU







function id_matrix(n)
    !> Identity nxn matrix
        implicit none
        integer(kind=4), intent(in) :: n
        integer(kind=4) :: i
        real(kind=8) :: id_matrix(n,n)
        id_matrix=0d0
        do i=1,n
            id_matrix(i,i)=1d0
        end do
    end function

    function Householder(x)
        implicit none
        real(kind=8), intent(in) :: x(:)
        real(kind=8), allocatable :: Householder(:,:)
        integer(kind=4) :: n
        real(kind=8) :: beta, norm_x
        real(kind=8), parameter :: epsilon=1d-6
        real(kind=8), allocatable :: u(:)
        n=size(x)
        allocate(Householder(n,n),u(n))
        norm_x=p_norm_vec(x,2)
        beta=1d0/(norm_x*(abs(x(1))+norm_x))
        u(1)=sign(1d0,x(1))*(abs(x(1))+norm_x)
        u(2:n)=x(2:n)
        Householder=id_matrix(n)-beta*outer_prod_vec(u,u)
        !if (abs(det(Householder))+epsilon<1d0) error stop "Not orthogonal"
    end function
    
    subroutine Gram_Schmidt_mat(v,u)
        implicit none
        real(kind=8), intent(in) :: v(:,:)
        real(kind=8), intent(out) :: u(:,:) !> base ortogonal
        integer(kind=4) :: k,j,n
        real(kind=8), allocatable :: sum_proy_ortog(:)
        n=size(v,1)
        allocate(sum_proy_ortog(n))
        u(1:n,1)=v(1:n,1)
        do k=2,n
            sum_proy_ortog=0d0
            do j=1,k-1
                sum_proy_ortog=sum_proy_ortog+proy_ortog(v(1:n,k),u(1:n,j))
            end do
            u(1:n,k)=v(1:n,k)-sum_proy_ortog
        end do
    end subroutine
    
    subroutine forward_substitution(L,b,x)
        implicit none
        real(kind=8), intent(in) :: L(:,:) !> lower triangular matrix
        real(kind=8), intent(in) :: b(:) !> vector
        real(kind=8), intent(out) :: x(:) !> solution of linear system
        integer(kind=4) :: j,k,n
        real(kind=8) :: sum
        n=size(b)
        x(1)=b(1)/L(1,1)
        do j=2,n
            sum=0d0
            do k=1,j-1
                sum=sum+L(j,k)*x(k)
            end do
            x(j)=(b(j)-sum)/L(j,j)
        end do
    end subroutine 


    !> Determinant of square matrix using LU decomposition
subroutine compute_det(A,tol,det,error)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(in) :: tol !> tolerance for determinant
    real(kind=8), intent(out) :: det !> determinant
    logical, intent(out) :: error
    
    real(kind=8) :: det_U
    real(kind=8), allocatable :: L(:,:), U(:,:)
    integer(kind=4) :: n,i,nswaps
    
    error=.false.
    n=size(A,1)
    if (n/=size(A,2)) then
        error=.true.
        error stop "Matrix must be square (det)"
    end if
    if (n==1) then
        det=A(1,1)
    else if (n==2) then
        det=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    else
        allocate(L(n,n),U(n,n))
        call LU(A,tol,L,U,error,num_swaps=nswaps)
        det_U=1d0
        do i=1,n
            det_U=det_U*U(i,i)
        end do
        !> L has unit diagonal so det_L=1; sign from row swaps
        det=((-1d0)**nswaps)*det_U
    end if
end subroutine



!> Cholesky decomposition
subroutine Cholesky(A,L)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> positive definite, symmetric matrix
    real(kind=8), intent(out) :: L(:,:) !> lower triangular matrix with positive diagonal terms
    
    real(kind=8), allocatable :: prod_L_Lt(:,:) !> L*L^T
    integer(kind=4) :: i,j,k,n
    real(kind=8) :: sum1, sum2
    real(kind=8), parameter :: tol=1d-8
    real(kind=8) :: rel_err
    
    n=size(A,1)
    if (n/=size(A,2)) error stop "Cholesky: matrix must be square"
    !> Check symmetry
    do i=1,n
        do j=i+1,n
            if (abs(A(i,j)-A(j,i))>tol*(abs(A(i,j))+abs(A(j,i))+1d-300)) then
                print *, "Cholesky: A is not symmetric. A(", i, ",", j, ")=", A(i,j), " /= A(", j, ",", i, ")=", A(j,i)
                error stop "Cholesky: matrix is not symmetric"
            end if
        end do
    end do
    !> Check positive definiteness (all diagonal elements must be positive)
    do i=1,n
        if (A(i,i)<=0d0) then
            print *, "Cholesky: A(", i, ",", i, ")=", A(i,i), " is not positive"
            error stop "Cholesky: matrix is not positive definite (non-positive diagonal)"
        end if
    end do
    L=0d0
    do j=1,n
        if (j == 1) then
            L(j,j)=sqrt(A(j,j))
        else
        sum1=0d0
        do k=1,j-1
            sum1=sum1+L(j,k)**2
        end do
        if (A(j,j)-sum1<=0d0) then
            print *, "Cholesky: pivot A(", j, ",", j, ")-sum=", A(j,j)-sum1, " is not positive"
            error stop "Cholesky: matrix is not positive definite"
        end if
        L(j,j)=sqrt(A(j,j)-sum1)
        end if
        do i=j+1,n
            sum2=0d0
            do k=1,j-1
                sum2=sum2+L(i,k)*L(j,k)
            end do
            L(i,j)=(A(i,j)-sum2)/L(j,j)
        end do
    end do
    prod_L_Lt=matmul(L,transpose(L))
    rel_err=norm_mat_inf(A-prod_L_Lt)/(norm_mat_inf(A)+1d-300)
    if (rel_err > tol) then
        print *, "Relative Cholesky error:", rel_err
        error stop "Wrong Cholesky decomposition"
    end if
end subroutine

    !> This subroutine performs backward substitution in LU decomposition of matrix
subroutine backward_substitution(U,b,x)
    implicit none
    real(kind=8), intent(in) :: U(:,:) !> upper triangular matrix
    real(kind=8), intent(in) :: b(:) !> independent term
    real(kind=8), intent(out) :: x(:) !> solution of linear system
    
    integer(kind=4) :: j,k,n
    real(kind=8) :: sum
    
    n=size(b)
    x(n)=b(n)/U(n,n)
    do j=1,n-1
        sum=0d0
        do k=n-j+1,n
            sum=sum+U(n-j,k)*x(k)
        end do
        x(n-j)=(b(n-j)-sum)/U(n-j,n-j)
    end do
end subroutine

subroutine eigenvectors_tridiag_toeplitz_matrix(A)
    class(tridiag_Toeplitz_matrix_c) :: A
    
    integer(kind=4) :: i,j
    real(kind=8) :: L2_norm_vj
    real(kind=8), parameter :: pi=4d0*atan(1d0)
    
    if (A%sub*A%super<=0d0) error stop "a*c must be positive"
    allocate(A%eigenvectors(A%num_cols,A%num_cols))
    do j=1,A%num_cols
        do i=1,A%num_cols
            A%eigenvectors(i,j)=(A%sub/A%super)**((i-1d0)/2d0)*sin((A%num_cols-j+1)*pi*i/(A%num_cols+1d0))
        end do
        L2_norm_vj=p_norm_vec(A%eigenvectors(:,j),2)
        A%eigenvectors(:,j)=A%eigenvectors(:,j)/L2_norm_vj
    end do
end subroutine

subroutine eigenvectors_tridiag_sym_matrix(a,b,lambda,v)
    real(kind=8), intent(in) :: a(:) !> diagonal elements
    real(kind=8), intent(in) :: b(:) !> non-diagonal elements
    real(kind=8), intent(in) :: lambda(:) !> eigenvalues
    real(kind=8), intent(out) :: v(:,:) !> eigenvectors
    
    real(kind=8), parameter :: epsilon=1d-12
    real(kind=8) :: L2_norm_vj
    
    integer(kind=4) :: i,j,k,n,m
    integer(kind=4) :: l(2)
    
    n=size(lambda)
    
    do j=1,n
        if (abs(a(1)-lambda(j))>=epsilon .and. abs(a(n)-lambda(j))>=epsilon) then
            v(1,j)=1d0
            v(2,j)=v(1,j)*(lambda(j)-a(1))/b(1)
            do k=3,n-1
                v(k,j)=v(k-1,j)*(lambda(j)-a(k-1))/b(k-1) - v(k-2,j)*b(k-2)/b(k-1)
            end do
            v(n,j)=v(n-1,j)*b(n-1)/(lambda(j)-a(n))
        else if (abs(a(1)-lambda(j))<epsilon) then
            v(2,j)=0d0
            v(1,j)=1d0
            do k=3,n-1
                v(k,j)=v(k-1,j)*(lambda(j)-a(k-1))/b(k-1) - v(k-2,j)*b(k-2)/b(k-1)
            end do
            v(n,j)=v(n-1,j)*b(n-1)/(lambda(j)-a(n))
        else
            v(n-1,j)=0d0
            v(n,j)=1d0
            do k=2,n-2
                v(n-k,j)=v(n-k+1,j)*(lambda(j)-a(n-k+1))/b(n-k) - v(n-k+2,j)*b(n-k+1)/b(n-k)
            end do
            v(1,j)=v(2,j)*b(1)/(lambda(j)-a(1))
        end if
    L2_norm_vj=p_norm_vec(v(:,j),2)
    v(:,j)=v(:,j)/L2_norm_vj
    end do
    do i=1,n
        do j=1,n
            if (abs(lambda(i)-lambda(j))<epsilon .and. i/=j) then
                l=[i,j]
            end if
        end do
    end do
    if (mod(n,2).eq.0) then
        m=n
    else
        m=n-1
    end if
end subroutine

subroutine check_eigenvectors(A,lambda,v,tolerance)
    real(kind=8), intent(in) :: A(:,:)
    real(kind=8), intent(in) :: lambda(:)
    real(kind=8), intent(in) :: v(:,:)
    real(kind=8), intent(in) :: tolerance
    
    real(kind=8), allocatable :: A_lambda(:,:,:),A_lambda_P(:),prod(:,:)
    integer(kind=4) :: i,j,n

    n=size(lambda)
    allocate(A_lambda(n,n,n))
    do i=1,n
        A_lambda(:,:,i)=A
        do j=1,n
            A_lambda(j,j,i)=A_lambda(j,j,i)-lambda(i)
        end do
    end do
    do i=1,n
        A_lambda_P=matmul(A_lambda(:,:,i),v(:,i))
        if (inf_norm_vec_real(A_lambda_P)>=tolerance) then
            print *, "Error in eigenvector", i, inf_norm_vec_real(A_lambda_P)
            !error stop
        end if
    end do
    prod=matmul(v,transpose(v))
    do i=1,n
       if (i.eq.1) then
            if (abs(prod(i,i)-1d0)>=tolerance .or. inf_norm_vec_real(prod(i+1:n,i))>=tolerance) print *, "Error in first eigenvector"
        else if (i>1 .and. i<n) then
            if (abs(prod(i,i)-1d0)>=tolerance .or. inf_norm_vec_real(prod(1:i-1,i))>=tolerance .or. &
                inf_norm_vec_real(prod(i+1:n,i))>=tolerance) then
                print *, "Error in eigenvector", i
                error stop
            end if
        else
            if (abs(prod(i,i)-1d0)>=tolerance .or. inf_norm_vec_real(prod(1:i-1,i))>=tolerance) print *, "Error in last eigenvector"
        end if
    end do
end subroutine


subroutine check_eigenvectors_tridiag_sym_matrix(this,tolerance)
class(tridiag_sym_matrix_c), intent(in) :: this !> tridiagonal symmetric matrix object
real(kind=8), intent(in) :: tolerance

real(kind=8), allocatable :: A_lambda(:,:,:),A_lambda_P(:),prod(:,:)
integer(kind=4) :: i,j,n

    n=size(this%eigenvalues)
    allocate(A_lambda(n,n,n))
    do i=1,n
        A_lambda(:,:,i)=0d0
        A_lambda(1,1:2,i)=[this%diag(1)-this%eigenvalues(i),this%sub(1)]
        do j=2,n-1
            A_lambda(j,j-1,i)=this%sub(j-1)
            A_lambda(j,j+1,i)=this%sub(j)
            A_lambda(j,j,i)=this%diag(j)-this%eigenvalues(i)
        end do
        A_lambda(n,n-1:n,i)=[this%sub(n-1),this%diag(n)-this%eigenvalues(i)]
    end do
    do i=1,n
        A_lambda_P=matmul(A_lambda(:,:,i),this%eigenvectors(:,i))
        if (inf_norm_vec_real(A_lambda_P)>=tolerance) then
            print *, "Error in eigenvector", i, inf_norm_vec_real(A_lambda_P)
            error stop
        end if
    end do
    prod=matmul(this%eigenvectors,transpose(this%eigenvectors))
    !if ((abs(det(this%eigenvectors))-1d0)>=tolerance) then
    !    print *, "abs(det(P)) is not 1"
    !    error stop
    !end if
end subroutine

subroutine compute_eigenvalues(this)
    class(sq_matrix_c) :: this
    
    real(kind=8), allocatable :: d(:),e(:)
    integer(kind=4) :: info
!> Lapack
    select type (this)
    class is (tridiag_sym_matrix_c)
        d=this%diag
        e=this%sub
        !call dsterf(size(d),d,e,info)
        this%eigenvalues=d
    end select
end subroutine

!> USAR LAPACK
subroutine compute_eigenvectors(this)
    class(sq_matrix_c) :: this
end subroutine

subroutine get_can_vec(this,tol,num_can_vec,ind_can_vec,num_non_can_vec,ind_non_can_vec)
!> Computes the number and indices of canonical and non-canonical vectors
implicit none
class(real_array_c), intent(in) :: this !> array object
real(kind=8), intent(in) :: tol !> tolerance
integer(kind=4), intent(out) :: num_can_vec !> number of canonical vectors
integer(kind=4), intent(out), allocatable :: ind_can_vec(:) !> indices of canonical vectors
integer(kind=4), intent(out) :: num_non_can_vec !> number of non-canonical vectors
integer(kind=4), intent(out), allocatable :: ind_non_can_vec(:) !> indices of non-canonical vectors
integer(kind=4) :: i,j
num_can_vec=0 !> we initialize number of canonical vectors
num_non_can_vec=0 !> we initialize number of non-canonical vectors
!> First loop: we count the number of canonical and non-canonical vectors
do i=1,this%num_cols
    if (abs(this%cols(i)%col_1(1)-1d0)<tol .and. &
        inf_norm_vec_real(this%cols(i)%col_1(2:this%cols(i)%dim))<tol) then
        num_can_vec=num_can_vec+1 !> we update number of canonical vectors
    else
        num_non_can_vec=num_non_can_vec+1 !> we update number of non-canonical vectors
    end if
end do
allocate(ind_can_vec(num_can_vec)) !> we allocate the array of indices of canonical vectors
allocate(ind_non_can_vec(num_non_can_vec)) !> we allocate the array of indices of canonical vectors
!> Second loop: we store the indices of canonical and non-canonical vectors
num_can_vec=0 !> we initialize number of canonical vectors
num_non_can_vec=0 !> we initialize number of non-canonical vectors
do i=1,this%num_cols
    if (abs(this%cols(i)%col_1(1)-1d0)<tol .and. &
        inf_norm_vec_real(this%cols(i)%col_1(2:this%cols(i)%dim))<tol) then
        num_can_vec=num_can_vec+1 !> we update number of canonical vectors
        ind_can_vec(num_can_vec)=i !> we store the index of the canonical vector
    else
        num_non_can_vec=num_non_can_vec+1 !> we update number of non-canonical vectors
        ind_non_can_vec(num_non_can_vec)=i !> we store the index of the non-canonical vector
    end if
end do
end subroutine

subroutine prod_tridiag_mat_vec(this,b,x)
    implicit none
    class(tridiag_matrix_c), intent(in) :: this                               !< Tridiagonal matrix
    real(kind=8), intent(in) :: b(:)                                       !< Input vector
    real(kind=8), intent(out) :: x(:)                                      !< Result vector
    integer(kind=4) :: i,n                                                !< Loop indices and dimensions
    n=size(this%diag)                                                        !< Matrix dimension
    if (size(b)/=n) error stop "Dimension error in b"
    x(1)=this%diag(1)*b(1)+this%super(1)*b(2)
    do i=2,n-1
        x(i)=this%sub(i-1)*b(i-1)+this%diag(i)*b(i)+this%super(i)*b(i+1)
    end do
    x(n)=this%sub(n-1)*b(n-1)+this%diag(n)*b(n)                                 !< Last element
end subroutine

subroutine copy_real_array(src,dest)
    implicit none
    class(real_array_c), intent(in) :: src  !< Source array
    class(real_array_c), intent(out) :: dest !< Destination array
    
    integer(kind=4) :: i

    !if (size(src%col_1) /= size(dest%col_1)) error stop "Dimension error in copy_real_array"
    !dest%num_cols = src%num_cols
    call dest%allocate_array(src%num_cols)
    !allocate(dest%cols(src%num_cols))
    do i=1,dest%num_cols
        call dest%cols(i)%allocate_vector(src%cols(i)%dim)
        dest%cols(i)%col_1 = src%cols(i)%col_1
    end do
end subroutine

subroutine QR_Householder(A,Q,R)
    implicit none
    real(kind=8), intent(in) :: A(:,:) !> square matrix
    real(kind=8), intent(out) :: Q(:,:), R(:,:)
    !real(kind=8), external :: p_norm_vec, outer_prod_vec
    real(kind=8), allocatable :: Ak(:,:), Ak_tilde(:,:), Hk(:,:), Hk_tilde(:,:), Q_aux(:,:), a1(:)
    real(kind=8) :: beta
    integer(kind=4) :: n,i
    n=size(A,1)
    allocate(Ak(n,n),Hk(n,n),Q_aux(n,n),a1(n))
    Ak=A !A0
    a1=Ak(1:n,1)
    Hk=Householder(a1) !H1
    Q=id_matrix(n)
    do i=1,n-1
        R=matmul(Hk,Ak)
        Ak=R !A1
        Ak_tilde=Ak((i+1):n,(i+1):n) !Atilde_1
        a1=Ak_tilde(1:(n-i),1)
        Hk_tilde=Householder(a1) !Htilde_2
        Hk=0d0
        Hk(1:i,1:i)=id_matrix(i)
        Hk((i+1):n,(i+1):n)=Hk_tilde !H2
        Q_aux=matmul(Q,Hk)
        Q=Q_aux
    end do
end subroutine QR_Householder

        function compute_mix_ratios_R_opt4(this) result(mix_ratios_R)
            !> \brief Computes the mixing ratios for reaction amounts in next time step at computation targets (option 4 estimation downstream waters)
            !> \details
            !>   - Computes mixing ratios for reaction amounts in next time step at computation targets
            !>   - Option 4: sums contributions from computation target and downstream waters
            !>   - Assumes mixing ratios already sorted in current, upstream, downstream order
            !>   - Efficient computation using pre-stored mixing ratios structure
            !> \param[in] this Mixing ratios array object
            implicit none
            class(real_array_c), intent(in) :: this !> mixing ratios array object
            !integer(kind=4), intent(in) :: i !> theta parameter
            !class(int_array_c), intent(in) :: mix_target_waters_indices !> mixing target waters indices
            real(kind=8), allocatable :: mix_ratios_R(:) !> new mixing ratios
            
            integer(kind=4) :: j,i,num_up,num_down,n_targets,dim
            real(kind=8) :: sum

            allocate(mix_ratios_R(this%num_cols))

            ! dim=size(mix_target_waters_indices) !> dimension of mixing waters indices

            ! num_down=mix_target_waters_indices(dim) !> number of downstream waters
            ! num_up=mix_target_waters_indices(dim-1) !> number of upstream waters
        !     !mix_ratio_R=0d0
            
        !     !allocate(mix_ratios_R(this%num_cols))
        !     !n_targets=this%spatial_discr%Num_targets
        !     !> Mixing ratios for reaction amounts in next time step at computation targets
            do i=1,this%num_cols
                sum=this%cols(i)%col_1(1) !> computation target contribution
                !> loop over downstream waters
                ! do j=1,mix_target_waters_indices%cols(i)%col_1(mix_target_waters_indices%cols(i)%dim) !> number of downstream waters
                !     num_up=mix_target_waters_indices%cols(i)%col_1(&
                !         mix_target_waters_indices%cols(i)%dim-1) !> number of upstream waters
                !     sum=sum+this%cols(i)%col_1(1+num_up+j) !> we update sum with downstream waters
                ! end do
                do j=1,this%cols(i)%dim-i
                    sum=sum+this%cols(i)%col_1(i+j) !> we update sum with downstream waters
                end do
                mix_ratios_R(i)=sum !> final mixing ratio for target i
            end do
            !mix_ratios_R=mix_ratios_R*theta !> multiply by theta parameter
        end function compute_mix_ratios_R_opt4
        
        function fill_pentadiag_mat(this, Nx) result(M)
            implicit none
            class(pentadiag_matrix_c) :: this !> pentadiagonal matrix object
            integer(kind=4), intent(in) :: Nx !> number of cells in x direction (matrix dimension)
            real(kind=8), allocatable :: M(:,:) !> output matrix to fill pentadiagonal structure
            
            integer(kind=4) :: i,j,Ny
                        
            Ny=this%num_cols/Nx !> number of cells in y direction
            allocate(M(this%num_cols,this%num_cols)) !> Allocate output matrix
            M=0d0 !> Initialize matrix to zero
            !this%diag(1)=this%diag(1) !> First diagonal element
            !this%super(1)=this%super(1) !> First superdiagonal element
            !this%super2(1)=this%super2(1) !> First second superdiagonal element
            
            do i=1,this%num_cols
                M(i,i)=this%diag(i) !> Diagonal elements
                !this%super(i)=this%super(i) !> Superdiagonal elements
                !this%sub(i-1)=this%sub(i-1) !> Subdiagonal elements
                !if (i<n-1) then
                !    this%super2(i)=this%super2(i) !> Second superdiagonal elements
                !    this%sub2(i-1)=this%sub2(i-1) !> Second subdiagonal elements
                !end if
            end do
            do j=1,Ny
                do i=1,Nx-1
                    M((j-1)*Nx+i,(j-1)*Nx+i+1)=this%super((j-1)*(Nx-1)+i) !> Superdiagonal elements
                    M((j-1)*Nx+i+1,(j-1)*Nx+i)=this%sub((j-1)*(Nx-1)+i) !> Subdiagonal elements
                end do
            end do
            do i=1,Nx
                do j=1,Ny-1
                    M(j*Nx+i,(j-1)*Nx+i)=this%sub2((i-1)*(Ny-1)+j) !> Second subdiagonal elements
                    M((j-1)*Nx+i,j*Nx+i)=this%super2((i-1)*(Ny-1)+j) !> Second superdiagonal elements
                end do
            end do
            !this%diag(n)=this%diag(n) !> Last diagonal element
            !this%sub(n-1)=this%sub(n-1) !> Last subdiagonal element
        end function fill_pentadiag_mat
        
        subroutine inv_SPD_mat(A,inv_A) !> Invert symmetric positive definite matrix using Cholesky decomposition
            implicit none
            real(kind=8), intent(in) :: A(:,:) !> symmetric positive definite matrix
            real(kind=8), intent(out) :: inv_A(:,:) !> inverse of A (must be pre-allocated)
            
            real(kind=8), allocatable :: L(:,:), y(:), x(:), prod(:,:)
            integer(kind=4) :: n,i
            real(kind=8) :: rel_err
            real(kind=8), parameter :: tol=1d-8
            
            if (size(A,1) /= size(A,2)) error stop "Matrix must be square for inversion"
            
            n=size(A,1)
            allocate(L(n,n),y(n),x(n))
            
            call Cholesky(A,L) !> Compute Cholesky decomposition A = LL^T
            
            do i=1,n
                y=0d0
                y(i)=1d0 !> Unit vector for forward substitution
                call forward_substitution(L,y,x) !> Solve Ly = e_i for y
                call backward_substitution(transpose(L),x,inv_A(:,i)) !> Solve L^T x = y for column i of inv_A
            end do
            !> Verify inverse: ||A * inv_A - I||_inf / ||A||_inf < tol
            prod=matmul(A,inv_A)
            do i=1,n
                prod(i,i)=prod(i,i)-1d0
            end do
            rel_err=norm_mat_inf(prod)/(norm_mat_inf(A)+1d-300)
            if (rel_err > tol) then
                print *, "Relative inverse error:", rel_err
                error stop "Wrong SPD matrix inversion"
            end if
        end subroutine inv_SPD_mat

end module