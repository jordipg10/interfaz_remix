!> \file properties_m.f90
!> \brief Abstract base module for physical properties in PDE-based transport models.
!> \details
!> Defines the abstract superclass `props_c` for physical properties (diffusion, flow,
!> transport) that are used in PDE discretizations. Properties include source/sink terms,
!> homogeneity flags, and stationarity flags. Concrete subclasses
!> (e.g., diff_props_heterog_m, transport_properties_heterog_m, flow_props_heterog_m)
!> implement the deferred `read_props` procedure.
!>
!> \author Jordi
!> \date Unknown
!> \ingroup transport

!> \brief Abstract properties module for PDE-based transport.
!> \details Contains the abstract type `props_c` and common property setter/allocator routines.
module properties_m
    use BCs_m, only: BCs_1D_c !< BCs_1D_c used for evaporation flag in source term
    use spatial_discr_m, only: spatial_discr_c !< Spatial discretization class
    implicit none
    save
    private
    !> \brief Abstract superclass for physical properties.
    !> \details Stores source/sink term arrays and flags that indicate whether
    !> properties are spatially homogeneous or temporally stationary. Subclasses
    !> must implement `read_props` to load property data from files.
    type, public, abstract :: props_c
        integer(kind=4) :: source_term_order                    !< [#] Polynomial degree of source term (if applicable)
        real(kind=8), allocatable :: source_term(:)             !< [1/T] Source/sink term array r(x) at each cell
        integer(kind=4), allocatable :: source_term_flag(:)     !< [-] Flag per cell: 0 if r<0 and no evaporation, 1 otherwise
        logical :: homog_flag                                   !< [-] .true. if all properties are spatially homogeneous
        logical :: stat_flag                                    !< [-] .true. if properties are temporally stationary
    contains
        procedure :: set_source_term           !< Set source/sink term array
        procedure :: set_source_term_order     !< Set polynomial degree of source term
        procedure :: set_source_term_flag      !< Set source term flags based on BCs
        procedure :: set_homog_flag            !< Set spatial homogeneity flag
        procedure :: set_stat_flag             !< Set temporal stationarity flag
        procedure(read_props), public, deferred :: read_props  !< Read properties from file (deferred)
        procedure :: allocate_props            !< Allocate property arrays for given number of cells
    end type
    
    abstract interface
    
        subroutine read_props(this,root,spatial_discr)
            import props_c
            import spatial_discr_c
            class(props_c) :: this
            character(len=*), intent(in) :: root
            class(spatial_discr_c), intent(in), optional :: spatial_discr
        end subroutine
        
        subroutine are_props_homog(this)
            import props_c
            class(props_c) :: this
        end subroutine
            
    end interface
    
    contains
        !> \brief Set source/sink term array.
        !> \param[in,out] this  Properties object
        !> \param[in]     source_term  Source/sink term values at each cell [1/T]
        subroutine set_source_term(this,source_term)
            implicit none
            class(props_c) :: this
            real(kind=8), intent(in) :: source_term(:)
            this%source_term=source_term
        end subroutine
        
        !> \brief Set polynomial degree of source term.
        !> \param[in,out] this  Properties object
        !> \param[in]     source_term_order  Polynomial degree [#]
        subroutine set_source_term_order(this,source_term_order)
            implicit none
            class(props_c) :: this
            integer(kind=4), intent(in) :: source_term_order
            this%source_term_order=source_term_order
        end subroutine
        
         !> \brief Set source term flags based on boundary conditions.
         !> \details Flags cells where the source term is active. A flag of 0 indicates
         !> discharge (negative source) without evaporation; 1 otherwise.
         !> \param[in,out] this  Properties object
         !> \param[in]     BCs   1D boundary conditions (used for evaporation flag)
         subroutine set_source_term_flag(this,BCs)
            implicit none
            class(props_c) :: this
            class(BCs_1D_c), intent(in) :: BCs
            
            integer(kind=4) :: i
            if (.not. allocated(this%source_term_flag)) then
                allocate(this%source_term_flag(size(this%source_term)))
            end if
            this%source_term_flag=1
            do i=1,size(this%source_term)
                if (this%source_term(i)<0 .and. BCs%evap.eqv..false.) then !> discharge
                    this%source_term_flag(i)=0
                end if
            end do
         end subroutine
         
         !> \brief Set spatial homogeneity flag.
         !> \param[in,out] this  Properties object
         !> \param[in]     flag  .true. if all properties are spatially homogeneous
         subroutine set_homog_flag(this,flag)
         implicit none
         class(props_c) :: this
         logical, intent(in) :: flag
         this%homog_flag=flag 
         end subroutine
         
        !> \brief Set temporal stationarity flag.
        !> \param[in,out] this  Properties object
        !> \param[in]     flag  .true. if properties do not change in time
        subroutine set_stat_flag(this,flag)
        implicit none
        class(props_c) :: this
        logical, intent(in) :: flag
        this%stat_flag=flag
        end subroutine
        
        !> \brief Allocate property arrays for a given number of mesh cells.
        !> \param[in,out] this      Properties object
        !> \param[in]     Num_cells Number of mesh cells [#]
        subroutine allocate_props(this,Num_cells)
            class(props_c) :: this
            integer(kind=4), intent(in) :: Num_cells            !< [#] Number of mesh cells
            if (.not. allocated(this%source_term)) allocate(this%source_term(Num_cells))
            if (.not. allocated(this%source_term_flag)) allocate(this%source_term_flag(Num_cells))
        end subroutine
end module