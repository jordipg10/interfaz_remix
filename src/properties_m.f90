module properties_m
    use BCs_m, only: BCs_t !> BCs_t is used to set the source term flag
    use spatial_discr_m, only: spatial_discr_c !> spatial discretization class
    implicit none
    save
    type, public, abstract :: props_c !> properties superclass
        integer(kind=4) :: source_term_order                    !> if source term is polynomial, this is the degree
        real(kind=8), allocatable :: source_term(:)             !> (r)
        integer(kind=4), allocatable :: source_term_flag(:)     !> 0 if r<0 & no evaporation
                                                                !> 1 otherwise
        logical :: homog_flag                                  !> True if all properties are homogenous, False otherwise
        logical :: stat_flag !> TRUE if properties are stationary, FALSE otherwise
    contains
        procedure, public :: set_source_term
        procedure, public :: set_source_term_order
        procedure, public :: set_source_term_flag
        procedure, public :: set_homog_flag
        procedure, public :: set_stat_flag !> set stationary flag
        procedure(read_props), public, deferred :: read_props
        procedure(are_props_homog), public, deferred :: are_props_homog
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
        subroutine set_source_term(this,source_term)
            implicit none
            class(props_c) :: this
            real(kind=8), intent(in) :: source_term(:)
            this%source_term=source_term
        end subroutine
        
        subroutine set_source_term_order(this,source_term_order)
            implicit none
            class(props_c) :: this
            integer(kind=4), intent(in) :: source_term_order
            this%source_term_order=source_term_order
        end subroutine
        
         subroutine set_source_term_flag(this,BCs)
            implicit none
            class(props_c) :: this
            class(BCs_t), intent(in) :: BCs
            
            integer(kind=4) :: i
            allocate(this%source_term_flag(size(this%source_term)))
            this%source_term_flag=1
            do i=1,size(this%source_term)
                if (this%source_term(i)<0 .and. BCs%evap.eqv..false.) then !> discharge
                    this%source_term_flag(i)=0
                end if
            end do
         end subroutine
         
         subroutine set_homog_flag(this,flag) !> Set homogeneity flag
         implicit none
         class(props_c) :: this !> properties object
         logical, intent(in) :: flag !> True if all properties are homogenous, False otherwise
         this%homog_flag=flag 
         end subroutine
         
        subroutine set_stat_flag(this,flag)
        implicit none
        class(props_c) :: this !> properties object
        logical, intent(in) :: flag !> TRUE if properties are stationary, FALSE otherwise
        this%stat_flag=flag !> set stationary flag
        end subroutine
        
end module