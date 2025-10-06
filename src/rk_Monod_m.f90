!module rk_Monod_m
!    use rk_m
!    use Monod_m
!    implicit none
!    save
!    type, public, extends(rk_c) :: rk_redox_kin_c
!        !!class(glob_stoich_c), allocatable :: stoichiometry(:)
!        !integer(kind=4) :: n_paths ! numero de reacciones Monod
!        !integer(kind=4) :: n_t ! numero de terminos
!        !!integer(kind=4) :: n_M
!        !!real(kind=8), allocatable :: conc_cat(:)
!        !!type(species_c), allocatable :: TEAs(:) ! terminal electron acceptors
!        !!type(species_c), allocatable :: Monod_terms(:)
!        !type(species_c), allocatable :: catalysts(:)
!        !type(species_c) :: TEA
!        !!type(species_c) :: catalyst
!        !!type(species_c) :: inhibitor
!        !integer(kind=4) :: n_cat
!        !!real(kind=8) :: k_cat
!        !real(kind=8), allocatable :: k_cat(:)
!        !integer(kind=4) :: n_inh
!        !!type(aqueous_chemistry_c) :: inh
!        !!real(kind=8), allocatable :: conc_inh(:)
!        !type(species_c), allocatable :: inhibitors(:)
!        !!type(aqueous_chemistry_c) :: kappa_inh
!        !real(kind=8), allocatable :: k_inh(:)
!        !real(kind=8) :: rate_cst ! mu
!    contains
!        !procedure, public :: read_rk=>read_Monod
!        procedure, public :: compute_rk=>compute_Monod_bis
!        procedure, public :: compute_drk_dc=>compute_drk_dc_Monod
!        !procedure, public :: read_Monod_syst
!        !procedure, public :: set_parameters
!    end type
!        
!    interface
!        !subroutine read_Monod(this,filename)
!        !    import rk_redox_kin_c
!        !    implicit none
!        !    class(rk_redox_kin_c) :: this
!        !    character(len=*), intent(in) :: filename
!        !end subroutine
!        
!       ! subroutine read_Monod_syst(this,name)
!       !     import rk_redox_kin_c
!       !     implicit none
!       !     class(rk_redox_kin_c) :: this
!       !     character(len=*), intent(in) :: name
!       !     !type(glob_stoich_c), intent(out), allocatable :: stoich(:)
!       !end subroutine
!        
!        subroutine compute_Monod_bis(this,conc,kin_params)
!            import rk_redox_kin_c
!            implicit none
!            class(rk_redox_kin_c) :: this
!            real(kind=8), intent(in) :: conc(:)
!            class(kin_params_c), intent(in) :: kin_params
!            !integer(kind=4), intent(in), optional :: i
!        end subroutine
!        
!        subroutine compute_drk_dc_Monod(this,chem_syst,drk_dc,conc)
!            import rk_redox_kin_c
!            import chem_system_c
!            implicit none
!            class(rk_redox_kin_c) :: this
!            class(chem_system_c), intent(in) :: chem_syst
!            real(kind=8), intent(out) :: drk_dc(:)
!            real(kind=8), intent(in), optional :: conc(:)
!            !integer(kind=4), intent(in), optional :: i
!        end subroutine
!    end interface
!    
!    contains
!        !subroutine read_Monod(this,filename)
!        !    implicit none
!        !    class(rk_redox_kin_c) :: this
!        !    character(len=*), intent(in) :: filename
!        !    
!        !    integer(kind=4) :: i
!        !    
!        !    open(unit=1,file=filename,status='old',action='read')
!        !    read(1,"(/,F10.2)") this%k_cat
!        !    read(1,*) this%n_cat
!        !    allocate(this%species_cat(this%n_cat))
!        !    read(1,*) this%n_inh
!        !    allocate(this%species_inh(this%n_inh))
!        !    do i=1,this%n_cat
!        !        read(1,*) this%species_cat(i)%name
!        !    end do
!        !    do i=1,this%n_inh
!        !        read(1,*) this%species_inh(i)%name
!        !    end do
!        !    read(1,*) this%kappa_cat
!        !    read(1,*) this%kappa_inh
!        !    close(1)
!        !end subroutine
!    
!        !subroutine set_parameters(this,Monod)
!        !    implicit none
!        !    class(rk_redox_kin_c) :: this
!        !    class(redox_kin_c), intent(in) :: Monod
!        !    this%n_t=Monod%n_t
!        !    this%TEA=Monod%TEA
!        !    this%n_cat=Monod%n_cat
!        !    this%catalysts=Monod%catalysts
!        !    this%k_cat=Monod%k_cat
!        !    this%n_inh=Monod%n_inh
!        !    this%inhibitors=Monod%inhibitors
!        !    this%k_inh=Monod%k_inh
!        !    this%rate_cst=Monod%rate_cst
!        !end subroutine
!end module