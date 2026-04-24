!> \file penta_Cholesky_inv_m.f90
!> \brief Inverse of a symmetric positive definite pentadiagonal matrix via
!>        banded Cholesky decomposition (LAPACK DPBTRF / DPBTRS).
!> \details
!> Used for the 2D REMIX transition matrix \f$ T \f$ of size
!> \f$ N = N_x \times N_y \f$ with five diagonals at offsets
!> \f$ \{-N_x,\,-1,\,0,\,+1,\,+N_x\} \f$.
!>
!> Algorithm:
!>  1. Pack the pentadiagonal entries into LAPACK upper-banded storage.
!>  2. Cholesky factor \f$ T = L L^T \f$ via DPBTRF.
!>  3. Solve \f$ T \, X = I \f$ via DPBTRS to obtain \f$ T^{-1} \f$.
!>  4. Verify \f$ \|T \, T^{-1} - I\|_\infty < \varepsilon \f$.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup algebra
!> \see metodos_sist_lin_m
!===============================================================================
! Module: penta_cholesky_inv
!
! Computes the inverse of a symmetric positive definite (SPD) pentadiagonal
! matrix using Cholesky decomposition via LAPACK banded routines (DPBTRF/DPBTRS).
!
! Context (REMIX Theory, Section 16):
!   The 2D transport transition matrix T is pentadiagonal of size N x N,
!   where N = Nx * Ny. Its 5 diagonals are:
!
!     offset -Nx : lower far sub-diagonal  (y-direction coupling)
!     offset -1  : lower near sub-diagonal (x-direction coupling)
!     offset  0  : main diagonal
!     offset +1  : upper near super-diagonal (x-direction coupling)
!     offset +Nx : upper far super-diagonal  (y-direction coupling)
!
!   Three diagonals (offsets -1, 0, +1) are contiguous; the other two
!   (offsets -Nx, +Nx) are separated by a gap of Nx rows.
!
!   The near sub-/super-diagonal (offset ±1) has structural zeros at
!   y-layer boundaries (every Nx-th position), so only Ny*(Nx-1) non-zero
!   entries exist.  The far sub-/super-diagonal (offset ±Nx) couples
!   consecutive y-layers, yielding Nx*(Ny-1) entries.
!
!   For an SPD matrix, T = L * L^T  (Cholesky), where L is banded lower
!   triangular with bandwidth Nx.  The inverse is computed by solving
!   T * X = I using LAPACK's banded Cholesky routines.
!
! Subroutine signature:
!   call penta_spd_inverse(n, Nx, Ny, d0, d1, dm, Ainv)
!
! Arguments:
!   n             [in]  integer   - matrix dimension (= Nx * Ny)
!   Nx            [in]  integer   - x-direction grid size (= bandwidth gap)
!   Ny            [in]  integer   - y-direction grid size
!   d0(n)         [in]  real(dp), allocatable - main diagonal (offset 0)
!   d1(Ny*(Nx-1)) [in]  real(dp), allocatable - near sub-/super-diagonal (offset ±1)
!                                   Stored contiguously per y-layer: entries
!                                   for layer j are d1((j-1)*(Nx-1)+1 : j*(Nx-1)),
!                                   coupling nodes (i,j)-(i+1,j) for i=1..Nx-1.
!   dm(Nx*(Ny-1)) [in]  real(dp), allocatable - far sub-/super-diagonal (offset ±Nx)
!                                   Stored contiguously: dm((j-1)*Nx+1 : j*Nx)
!                                   couples layer j to layer j+1, for j=1..Ny-1.
!   Ainv(n,n)     [out] real(dp), allocatable - the full inverse matrix
!
! Note: Uses LAPACK DPBTRF (banded Cholesky factorisation) and DPBTRS
!       (banded Cholesky solve) for optimal performance.
!===============================================================================

module penta_cholesky_inv_m
  implicit none

  integer, parameter :: dp = selected_real_kind(15, 307)   ! double precision

  private
  public :: penta_spd_inverse

contains

  !-----------------------------------------------------------------------------
  ! penta_spd_inverse
  !
  ! Computes Ainv = T^{-1} for an n x n SPD pentadiagonal matrix T whose
  ! lower triangle is described by three diagonal arrays: d0, d1, dm.
  !
  ! Algorithm:
  !   1. Load the pentadiagonal entries into LAPACK banded storage AB(kd+1, n).
  !   2. Call DPBTRF to compute the banded Cholesky factorisation in-place.
  !   3. Set up Ainv = Identity, then call DPBTRS to solve T * Ainv = I.
  !-----------------------------------------------------------------------------
  subroutine penta_spd_inverse(n, Nx, Ny, d0, d1, dm, Ainv)
    implicit none
    integer(kind=4),  intent(in)  :: n, Nx, Ny
    real(kind=8), allocatable, intent(in)  :: d0(:)
    real(kind=8), allocatable, intent(in)  :: d1(:)
    real(kind=8), allocatable, intent(in)  :: dm(:)
    real(kind=8), allocatable, intent(out) :: Ainv(:,:)

    ! --- Local variables ---
    real(kind=8), allocatable :: AB(:,:)   ! Banded storage for LAPACK (kd+1, n)
    real(kind=8), allocatable :: Tv(:)     ! Work vector for sparse T * Ainv(:,j)
    real(kind=8) :: err_norm               ! ||T*Ainv - I||_inf
    real(kind=8) :: row_sum                ! Row sum accumulator
    real(kind=8), parameter :: verify_tol = 1.0d-6  ! Verification tolerance
    integer(kind=4) :: i, j, k, kd, ldab, info
    integer(kind=4) :: row, ix, jy, idx_d1

    ! ------------------------------------------------------------------
    ! Input validation
    ! ------------------------------------------------------------------
    if (n /= Nx * Ny) then
      write(*,*) 'ERROR: n /= Nx*Ny:', n, Nx, Ny
      error stop 'penta_spd_inverse: inconsistent dimensions n /= Nx*Ny'
    end if
    if (.not. allocated(d0) .or. size(d0) /= n) then
      error stop 'penta_spd_inverse: d0 must be allocated with size n'
    end if
    if (.not. allocated(d1) .or. size(d1) /= Ny*(Nx-1)) then
      write(*,*) 'ERROR: size(d1) =', size(d1), ' expected', Ny*(Nx-1)
      error stop 'penta_spd_inverse: d1 must be allocated with size Ny*(Nx-1)'
    end if
    if (.not. allocated(dm) .or. size(dm) /= Nx*(Ny-1)) then
      write(*,*) 'ERROR: size(dm) =', size(dm), ' expected', Nx*(Ny-1)
      error stop 'penta_spd_inverse: dm must be allocated with size Nx*(Ny-1)'
    end if
    do i = 1, n
      if (d0(i) <= 0.0_dp) then
        write(*,*) 'ERROR: d0(', i, ') =', d0(i), ' is not positive'
        error stop 'penta_spd_inverse: main diagonal must be positive for SPD matrix'
      end if
    end do

    kd = Nx     ! half-bandwidth (number of sub-diagonals)
    ldab = kd + 1

    ! ------------------------------------------------------------------
    ! Step 1:  Load T into LAPACK banded storage (upper triangle)
    !
    ! LAPACK DPBTRF with UPLO='U' uses upper-triangle banded storage:
    !   AB(kd+1+i-j, j) = A(i,j)  for max(1,j-kd) <= i <= j
    !
    ! Equivalently for the upper triangle:
    !   AB(kd+1, j)      = A(j,j)          (diagonal, row kd+1)
    !   AB(kd+1-(k), j)  = A(j-k, j)      for k = 1..min(kd, j-1)
    !
    ! For our symmetric matrix, A(i,j) = A(j,i), so:
    !   super-diagonal offset +1:  A(i, i+1) stored at AB(kd, i+1)
    !   super-diagonal offset +Nx: A(i, i+Nx) stored at AB(1, i+Nx)
    ! ------------------------------------------------------------------
    allocate(AB(ldab, n))
    AB = 0.0_dp

    ! --- Main diagonal (d0) -> row kd+1 ---
    do i = 1, n
      AB(kd + 1, i) = d0(i)
    end do

    ! --- Near super-diagonal (offset +1): A(i, i+1) = d1 ---
    do jy = 1, Ny
      do ix = 1, Nx - 1
        row = (jy - 1) * Nx + ix              ! global row i
        idx_d1 = (jy - 1) * (Nx - 1) + ix
        AB(kd, row + 1) = d1(idx_d1)          ! AB(kd, i+1) = A(i, i+1)
      end do
    end do

    ! --- Far super-diagonal (offset +Nx = kd): A(i, i+Nx) = dm ---
    do jy = 1, Ny - 1
      do ix = 1, Nx
        row = (jy - 1) * Nx + ix              ! global row i
        k   = (ix - 1) * (Ny - 1) + jy        ! index into dm (matches sub2/super2 storage)
        AB(1, row + Nx) = dm(k)               ! AB(1, i+Nx) = A(i, i+Nx)
      end do
    end do

    ! ------------------------------------------------------------------
    ! Step 2:  Banded Cholesky factorisation via LAPACK DPBTRF
    ! ------------------------------------------------------------------
    call dpbtrf('U', n, kd, AB, ldab, info)
    if (info /= 0) then
      write(*,*) 'ERROR: DPBTRF failed with info =', info
      error stop 'penta_spd_inverse: Cholesky factorisation failed'
    end if

    ! ------------------------------------------------------------------
    ! Step 3:  Solve T * Ainv = I via LAPACK DPBTRS
    ! ------------------------------------------------------------------
    allocate(Ainv(n, n))
    Ainv = 0.0_dp
    do i = 1, n
      Ainv(i, i) = 1.0_dp
    end do

    call dpbtrs('U', n, kd, n, AB, ldab, Ainv, n, info)
    if (info /= 0) then
      write(*,*) 'ERROR: DPBTRS failed with info =', info
      error stop 'penta_spd_inverse: banded Cholesky solve failed'
    end if

    deallocate(AB)

    ! ------------------------------------------------------------------
    ! Step 4:  Verification  ||T * Ainv - I||_inf < tol
    !
    ! Compute (T * Ainv - I) column by column using the sparse
    ! pentadiagonal structure.  Only one O(n) work vector is needed.
    ! ------------------------------------------------------------------
    allocate(Tv(n))
    err_norm = 0.0_dp

    do j = 1, n
      ! Compute Tv = T * Ainv(:,j) using the sparse pentadiagonal bands
      do i = 1, n
        Tv(i) = d0(i) * Ainv(i, j)
      end do
      ! Near sub-/super-diagonal (offset ±1)
      do jy = 1, Ny
        do ix = 1, Nx - 1
          row = (jy - 1) * Nx + ix
          idx_d1 = (jy - 1) * (Nx - 1) + ix
          Tv(row)     = Tv(row)     + d1(idx_d1) * Ainv(row + 1, j)
          Tv(row + 1) = Tv(row + 1) + d1(idx_d1) * Ainv(row, j)
        end do
      end do
      ! Far sub-/super-diagonal (offset ±Nx)
      do jy = 1, Ny - 1
        do ix = 1, Nx
          row = (jy - 1) * Nx + ix
          k   = (ix - 1) * (Ny - 1) + jy      ! matches sub2/super2 storage
          Tv(row)      = Tv(row)      + dm(k) * Ainv(row + Nx, j)
          Tv(row + Nx) = Tv(row + Nx) + dm(k) * Ainv(row, j)
        end do
      end do
      ! Subtract column j of identity: (T*Ainv - I)(:,j)
      Tv(j) = Tv(j) - 1.0_dp
      ! Update max-entry norm
      do i = 1, n
        err_norm = max(err_norm, abs(Tv(i)))
      end do
    end do

    deallocate(Tv)

    if (err_norm >= verify_tol) then
      write(*,*) 'ERROR: ||T*Ainv - I||_inf =', err_norm, ' >= tol =', verify_tol
      error stop 'penta_spd_inverse: inverse verification failed'
    end if

  end subroutine penta_spd_inverse

end module penta_cholesky_inv_m
