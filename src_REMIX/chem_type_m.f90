!!> Chemistry type module: contains the properties of different types of waters, solids and gases
!module chem_type_m
!    use aqueous_chemistry_m
!    implicit none
!    save
!    type, public, abstract :: chem_type_c !> chemistry type superclass
!        character(len=256) :: name
!    contains
!        procedure, public :: set_chem_type_name
!    end type
!    
!    type, public, extends(chem_type_c) :: water_type_c !> water type subclass
!        type(aqueous_chemistry_c) :: aq_chem !> aqueous chemistry object
!    contains
!        procedure, public :: set_aq_chem
!        procedure, public :: read_wat_type_CHEPROO
!    end type
!    
!    type, public, extends(chem_type_c) :: solid_type_c !> solid type subclass
!        type(solid_chemistry_c) :: solid_chem !> solid chemistry object
!    contains
!        procedure, public :: set_solid_chem
!    end type
!    
!    type, public, extends(chem_type_c) :: gas_type_c !> gas type subclass
!        type(gas_chemistry_c) :: gas_chem !> gas chemistry object
!    contains
!        procedure, public :: set_gas_chem
!    end type
!    
!    interface
!        subroutine read_wat_type_CHEPROO(this,n_p_aq,num_cstr,model,Jac_flag,unit,niter,CV_flag,surf_chem)
!            import water_type_c
!            import chem_system_c
!            import solid_chemistry_c
!            implicit none
!            class(water_type_c) :: this
!            integer(kind=4), intent(in) :: n_p_aq !> number of primary aqueous species
!            integer(kind=4), intent(in) :: num_cstr !> number of constrains
!            integer(kind=4), intent(in) :: model !> activity coefficients model
!            integer(kind=4), intent(in) :: Jac_flag !> 0: incremental coeffficinets, 1: analtical
!            integer(kind=4), intent(in) :: unit !> file unit
!            integer(kind=4), intent(out) :: niter !> number of iterations in Newton method
!            logical, intent(out) :: CV_flag !> TRUE if Newton method converges, FALSE otherwise
!            !real(kind=8), intent(in), optional :: c1_surf !> chapuza
!            !real(kind=8), intent(in), optional :: CEC !> chapuza
!            !real(kind=8), intent(out), optional :: conc_exch(:) !> chapuza
!            class(solid_chemistry_c), intent(inout), optional :: surf_chem
!        end subroutine
!    end interface
!    
!    contains
!      
!        
!        subroutine set_aq_chem(this,aq_chem) !> sets aqueous chemistry object in water type
!            implicit none
!            class(water_type_c) :: this
!            class(aqueous_chemistry_c), intent(in) :: aq_chem
!            this%aq_chem=aq_chem
!        end subroutine
!        
!        subroutine set_gas_chem(this,gas_chem) !> sets gas chemistry object in gas type
!            implicit none
!            class(gas_type_c) :: this
!            class(gas_chemistry_c), intent(in) :: gas_chem
!            this%gas_chem=gas_chem
!        end subroutine
!        
!        subroutine set_solid_chem(this,solid_chem) !> sets solid chemistry object in solid type
!            implicit none
!            class(solid_type_c) :: this
!            class(solid_chemistry_c), intent(in) :: solid_chem
!            this%solid_chem=solid_chem
!        end subroutine
!        
!        subroutine set_chem_type_name(this,name) !> sets chemistry type name
!            implicit none
!            class(chem_type_c) :: this
!            character(len=*), intent(in) :: name
!            this%name=name
!        end subroutine
!       
!end module 