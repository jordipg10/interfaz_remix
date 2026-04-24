!> \file vectors_m.f90
!> \brief Vector operations module for mathematical vector computations
!> \details 
!>   Provides abstract and concrete vector classes with common operations.
!>   Supports both real-valued and integer-valued vectors with polymorphism.
!>
!>   Key Features:
!>   - Abstract base class (vector_c) for polymorphic vector operations
!>   - Concrete implementations: vector_real_c (double precision), vector_int_c (integer)
!>   - Vector norms: p-norm, infinity norm
!>   - Vector products: outer product, dot product, orthogonal projection
!>   - Dynamic memory allocation for runtime-sized vectors
!>
!>   Vector Norms:
!>   1. **p-norm**: \f$ ||x||_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p} \f$
!>   2. **Euclidean norm** (p=2): \f$ ||x||_2 = \sqrt{\sum_{i=1}^n x_i^2} \f$
!>   3. **Infinity norm**: \f$ ||x||_\infty = \max_i |x_i| \f$
!>
!>   Vector Products:
!>   - **Outer product**: \f$ (u \otimes v)_{ij} = u_i v_j \f$ (produces matrix)
!>   - **Dot product**: \f$ u \cdot v = \sum_{i=1}^n u_i v_i \f$ (produces scalar)
!>   - **Orthogonal projection**: \f$ \text{proj}_v(u) = \frac{u \cdot v}{v \cdot v} v \f$
!>
!>   Class Hierarchy:
!>   - vector_c (abstract base)
!>     ├─ vector_real_c (real*8)
!>     └─ vector_int_c (integer*4)
!>
!>   \see arrays_m, array_ops_m
module vectors_m
    implicit none
    save
    
    !> \class vector_c
    !> \brief Abstract base class for all vector types
    !> \details 
    !>   Provides common interface for polymorphic vector operations.
    !>   Cannot be instantiated directly - use concrete derived types.
    type, public, abstract :: vector_c
        !> \var dim
        !> Vector dimension (number of elements) [-]
        integer(kind=4) :: dim
    contains
        !> Allocate memory for vector elements
        procedure :: allocate_vector
        !> Set the vector dimension
        procedure :: set_dim
    end type
    
    !> \class vector_real_c
    !> \brief Real-valued vector implementation
    !> \extends vector_c
    !> \details 
    !>   Concrete vector class for double-precision real numbers.
    !>   Provides dynamic allocation for runtime sizing.
    type, public, extends(vector_c) :: vector_real_c
        !> \var col_1
        !> Array of double-precision real values
        real(kind=8), allocatable :: col_1(:)
    end type
    
    !> \class vector_int_c
    !> \brief Integer-valued vector implementation
    !> \extends vector_c
    !> \details 
    !>   Concrete vector class for 4-byte integers.
    !>   Provides dynamic allocation for runtime sizing.
    type, public, extends(vector_c) :: vector_int_c
        !> \var col_1
        !> Array of 4-byte integer values
        integer(kind=4), allocatable :: col_1(:)
    end type
    
    contains
        !> \brief Compute outer product of two vectors
        !> \param[in] x First vector (size n)
        !> \param[in] y Second vector (size m)
        !> \return outer_prod_vec Matrix of size n×m where result(i,j) = x(i)·y(j)
        !> \details 
        !>   Computes tensor product (outer product) of two vectors.
        !>   
        !>   Mathematical Definition:
        !>   \f[
        !>     (x \otimes y)_{ij} = x_i \cdot y_j
        !>   \f]
        !>   
        !>   Properties:
        !>   - Result is always a rank-2 tensor (matrix)
        !>   - Dimension: n×m where n=size(x), m=size(y)
        !>   - Result is NOT commutative: x⊗y ≠ y⊗x (different dimensions)
        !>   
        !>   Example:
        !>   \code{.f90}
        !>   x = [1, 2]
        !>   y = [3, 4, 5]
        !>   result = [[3,  4,  5],
        !>             [6,  8, 10]]
        !>   \endcode
        function outer_prod_vec(x,y)
            implicit none
            real(kind=8), intent(in) :: x(:), y(:)                             !< Input vectors: assumed-shape arrays (size determined at runtime)
            real(kind=8), allocatable :: outer_prod_vec(:,:)                   !< Return value: allocatable 2D matrix for outer product result
            integer(kind=4) :: n,m,i,j                                         !< Local variables: vector sizes (n,m) and loop indices (i,j)
            n=size(x)                                                          !< Get size of first vector using intrinsic size() function
            m=size(y)                                                          !< Get size of second vector using intrinsic size() function
            allocate(outer_prod_vec(n,m))                                      !< Allocate result matrix with dimensions n×m (n rows, m columns)
            do i=1,n                                                           !< Outer loop: iterate through rows (elements of vector x)
                do j=1,m                                                       !< Inner loop: iterate through columns (elements of vector y)
                    outer_prod_vec(i,j)=x(i)*y(j)                              !< Compute outer product element: matrix(i,j) = x_i * y_j
                end do                                                         !< End inner loop over columns
            end do                                                             !< End outer loop over rows
        end function outer_prod_vec
        
        !> \brief Compute p-norm of a vector
        !> \param[in] x Vector for which to compute norm
        !> \param[in] p Norm parameter (must be ≥ 1)
        !> \return p_norm_vec The p-norm value ||x||_p
        !> \details 
        !>   Computes the p-norm (Hölder norm) of a vector.
        !>   
        !>   Mathematical Definition:
        !>   \f[
        !>     ||x||_p = \left(\sum_{i=1}^n |x_i|^p\right)^{1/p}
        !>   \f]
        !>   
        !>   Special Cases:
        !>   - p=1: Manhattan/Taxicab norm \f$ ||x||_1 = \sum |x_i| \f$
        !>   - p=2: Euclidean norm \f$ ||x||_2 = \sqrt{\sum x_i^2} \f$
        !>   - p→∞: Maximum norm \f$ ||x||_\infty = \max |x_i| \f$
        !>   
        !>   Implementation Notes:
        !>   - p=2 uses optimized sqrt() for efficiency
        !>   - General case uses power operation (1/p)
        !>   - Validates p ≥ 1 (mathematical requirement)
        !>   
        !>   \warning Halts execution if p < 1
        function p_norm_vec(x,p)
            implicit none
            real(kind=8), intent(in) :: x(:)                                   !< Input vector: assumed-shape array for flexible sizing
            integer(kind=4), intent(in) :: p                                   !< Norm parameter: p=1 for L1, p=2 for L2 (Euclidean), etc.
            real(kind=8) :: p_norm_vec                                          !< Function return value: computed p-norm scalar
            
            integer(kind=4) :: i,n                                             !< Local variables: loop index (i) and vector size (n)
            real(kind=8) :: sum                                                 !< Accumulator variable for sum of |x_i|^p terms
            real(kind=8), parameter :: epsilon=1d-9                            !< Small tolerance constant (10^-9) for numerical precision checks
            if (p<1) error stop "p must be >= 1"                               !< Input validation: p-norm only defined for p ≥ 1, halt execution if violated
            n=size(x)                                                          !< Get vector size using intrinsic size() function
            sum=0d0                                                            !< Initialize sum accumulator to zero (0d0 = double precision zero)
            do i=1,n                                                           !< Loop through all vector elements
                sum=sum+abs(x(i))**p                                           !< Accumulate: sum += |x_i|^p (absolute value raised to power p)
            end do                                                             !< End summation loop
            if (p==2) then                                                     !< Special case: p=2 is Euclidean norm (most common)
                p_norm_vec=sqrt(sum)                                           !< For p=2: ||x||_2 = sqrt(∑|x_i|^2) using optimized sqrt function
            else                                                               !< General case: any p ≠ 2
                p_norm_vec=sum**(1d0/p)                                        !< General formula: ||x||_p = (∑|x_i|^p)^(1/p)
            end if                                                             !< End conditional for norm computation
        end function p_norm_vec
        
        !> \brief Compute infinity norm (maximum absolute value) for real vectors
        !> \param[in] x Real-valued vector
        !> \return inf_norm_vec_real The infinity norm ||x||_∞
        !> \details 
        !>   Computes the infinity norm (maximum norm) of a real vector.
        !>   
        !>   Mathematical Definition:
        !>   \f[
        !>     ||x||_\infty = \max_{1 \leq i \leq n} |x_i|
        !>   \f]
        !>   
        !>   Properties:
        !>   - Limit case of p-norm as p→∞
        !>   - Computationally efficient (single pass through vector)
        !>   - Always well-defined for finite vectors
        !>   
        !>   Implementation: Uses intrinsic maxval() and abs() functions
        function inf_norm_vec_real(x)
            real(kind=8), intent(in) :: x(:)                                   !< Input: real vector with assumed-shape for flexibility
            real(kind=8) :: inf_norm_vec_real                                   !< Return value: infinity norm (scalar)
            inf_norm_vec_real=maxval(abs(x))                                   !< Compute: ||x||_∞ = max_i |x_i| using intrinsic maxval and abs functions
        end function
        
        !> \brief Compute infinity norm for integer vectors
        !> \param[in] x Integer-valued vector
        !> \return inf_norm_vec_int The infinity norm ||x||_∞
        !> \details 
        !>   Computes the infinity norm (maximum norm) of an integer vector.
        !>   
        !>   Mathematical Definition:
        !>   \f[
        !>     ||x||_\infty = \max_{1 \leq i \leq n} |x_i|
        !>   \f]
        !>   
        !>   Implementation: Uses intrinsic maxval() and abs() functions
        function inf_norm_vec_int(x)
            integer(kind=4), intent(in) :: x(:)                                !< Input: integer vector with assumed-shape
            integer(kind=4) :: inf_norm_vec_int                                 !< Return value: infinity norm (integer)
            inf_norm_vec_int=maxval(abs(x))                                    !< Compute: ||x||_∞ = max_i |x_i| for integer elements
        end function
        
        !> \brief Compute orthogonal projection of vector u onto vector v
        !> \param[in] u Vector to be projected
        !> \param[in] v Vector defining the projection direction
        !> \return proy_ortog Orthogonal projection of u onto v
        !> \details 
        !>   Computes the orthogonal (scalar) projection of u onto v.
        !>   
        !>   Mathematical Definition:
        !>   \f[
        !>     \text{proj}_v(u) = \frac{u \cdot v}{v \cdot v} v
        !>   \f]
        !>   
        !>   Geometric Interpretation:
        !>   - Result is the component of u in the direction of v
        !>   - Result is parallel to v
        !>   - Minimizes ||u - proj_v(u)||
        !>   
        !>   Properties:
        !>   - If u ⊥ v (orthogonal), then proj_v(u) = 0
        !>   - If u || v (parallel), then proj_v(u) = ±u
        !>   - Length: ||proj_v(u)|| = |u·v| / ||v||
        !>   
        !>   \warning Vectors must have same dimension; halts if sizes differ
        !>   \warning Undefined if v = 0 (division by zero)
        function proy_ortog(u,v)
            real(kind=8), intent(in) :: u(:), v(:)                             !< Input vectors: u (to project) and v (projection direction)
            real(kind=8), allocatable :: proy_ortog(:)                         !< Return: projection vector (same size as u and v)
            integer(kind=4) :: n                                               !< Vector dimension
            if (size(u)/=size(v)) error stop "u and v must have same dimension"  !< Validate: vectors must have equal dimensions
            n=size(u)                                                          !< Get common vector dimension
            allocate(proy_ortog(n))                                            !< Allocate result vector
            proy_ortog=(dot_product(u,v)/dot_product(v,v))*v                   !< Compute: proj_v(u) = [(u·v)/(v·v)]·v using intrinsic dot_product
        end function
        
        !> \brief Compute sum of squared elements
        !> \param[in] x Input vector
        !> \return sum_squares Sum of x_i² for all elements
        !> \details 
        !>   Computes the sum of squares of all vector elements.
        !>   Equivalent to the squared Euclidean norm.
        !>   
        !>   Mathematical Definition:
        !>   \f[
        !>     \text{sum\_squares}(x) = \sum_{i=1}^n x_i^2
        !>   \f]
        !>   
        !>   Relationship to Norms:
        !>   - sum_squares(x) = ||x||_2²
        !>   - ||x||_2 = sqrt(sum_squares(x))
        !>   
        !>   Used for:
        !>   - Variance calculations
        !>   - Residual error computations
        !>   - Avoiding sqrt() when only squared norm needed
        function sum_squares(x)
            implicit none
            real(kind=8), intent(in) :: x(:)                                   !< Input vector
            real(kind=8) :: sum_squares                                        !< Return: sum of squared elements
            integer(kind=4) :: i                                               !< Loop index
            sum_squares=0d0                                                    !< Initialize accumulator to zero
            do i=1,size(x)                                                     !< Loop over all vector elements
                sum_squares=sum_squares+x(i)**2                                !< Accumulate: sum += x_i²
            end do                                                             !< End summation loop
        end function
        
        !> \brief Allocate memory for vector elements
        !> \param[in,out] this   Vector object to allocate
        !> \param[in]     dim    Optional dimension (if omitted, uses this%dim)
        !> \details 
        !>   Allocates storage for vector elements based on concrete type.
        !>   Uses polymorphic dispatch to allocate appropriate array type.
        !>   
        !>   Behavior:
        !>   - If dim provided: sets this%dim and allocates with that size
        !>   - If dim omitted: uses existing this%dim value
        !>   
        !>   Polymorphic Allocation:
        !>   - vector_real_c: allocates real(8) array
        !>   - vector_int_c: allocates integer(4) array
        !>   
        !>   \note Values are uninitialized after allocation
        !>   \warning Overwrites existing allocation (potential memory leak if not deallocated first)
        subroutine allocate_vector(this,dim)
            implicit none
            class(vector_c) :: this                                            !< Vector object (polymorphic)
            integer(kind=4), intent(in), optional :: dim                       !< Optional: dimension to allocate
            integer(kind=4) :: alloc_stat                                      !< Allocation status code
            if (present(dim)) then                                             !< Check if dimension argument provided
                this%dim=dim                                                   !< Set vector dimension from argument
            end if                                                             !< Otherwise use existing this%dim
            if (this%dim < 0) then
                write(*,*) 'ERROR in allocate_vector: invalid dimension =', this%dim
                error stop "allocate_vector: negative dimension"
            end if
            select type (this)                                                 !< Polymorphic dispatch based on actual type
            class is (vector_real_c)                                           !< Case: real vector
                if (allocated(this%col_1)) deallocate(this%col_1)              !< Prevent memory leak from double allocation
                allocate(this%col_1(this%dim), stat=alloc_stat)                !< Allocate real array of size dim
                if (alloc_stat /= 0) then
                    write(*,*) 'ERROR in allocate_vector (real): failed to allocate', this%dim, 'elements'&
                        &, ' (', int(this%dim,8)*8_8, ' bytes)'
                    error stop "allocate_vector: allocation failed for vector_real_c"
                end if
            class is (vector_int_c)                                            !< Case: integer vector
                if (allocated(this%col_1)) deallocate(this%col_1)              !< Prevent memory leak from double allocation
                allocate(this%col_1(this%dim), stat=alloc_stat)                !< Allocate integer array of size dim
                if (alloc_stat /= 0) then
                    write(*,*) 'ERROR in allocate_vector (int): failed to allocate', this%dim, 'elements'&
                        &, ' (', int(this%dim,8)*4_8, ' bytes)'
                    error stop "allocate_vector: allocation failed for vector_int_c"
                end if
            end select                                                         !< End type-specific allocation
        end subroutine
        
        !> \brief Set the vector dimension
        !> \param[in,out] this Vector object to modify
        !> \param[in]     dim  Dimension value to set
        !> \details 
        !>   Sets the dimension attribute of the vector.
        !>   Does NOT allocate memory - call allocate_vector() for that.
        !>   
        !>   Typical Usage:
        !>   \code{.f90}
        !>   call vec%set_dim(10)
        !>   call vec%allocate_vector()  ! Uses dim=10
        !>   \endcode
        subroutine set_dim(this,dim)
            implicit none
            class(vector_c) :: this                                            !< Vector object (polymorphic)
            integer(kind=4), intent(in) :: dim                                 !< Dimension to set
            this%dim=dim                                                       !< Store dimension value
        end subroutine
end module