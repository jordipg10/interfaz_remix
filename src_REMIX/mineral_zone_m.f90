module mineral_zone_m
    use mineral_m, only: mineral_c
    use chem_system_m, only: chem_system_c
    implicit none
    save
!**********************************************************************************************************************************!
    type, public :: mineral_zone_c !> mineral zone class
        integer(kind=4) :: id=0 !> identifier of mineral zone
        integer(kind=4) :: num_minerals=0
        integer(kind=4) :: num_minerals_eq=0 !> number of minerals in equilibrium
        integer(kind=4) :: num_minerals_kin=0 !> number of minerals NOT in equilibrium
        integer(kind=4) :: num_minerals_kin_cst_act=0 !> number of minerals with constant activity
        integer(kind=4) :: num_minerals_kin_var_act=0 !> number of minerals with variable activity
        !type(mineral_c), allocatable :: minerals(:) !> first kinetic, then equilibrium
        integer(kind=4), allocatable :: ind_min_chem_syst(:) !> indices of minerals in chemical system (first kinetic var act, then kinetic cst act, then equilibrium var act, then equilibrium cst act)
        integer(kind=4), allocatable :: ind_min_Sk(:) !> indices of minerals in kinetic stoichiometric matrix of chemical system (first kinetic var act, then kinetic cst act)
        integer(kind=4) :: num_minerals_eq_cst_act=0 !> number of minerals in equilibrium with constant activity
        integer(kind=4) :: num_minerals_eq_var_act=0 !> number of minerals in equilibrium with variable activity
        class(chem_system_c), pointer :: chem_syst => null() !> chemical system associated to mineral zone
    contains
        procedure, public :: set_num_mins_min_zone
        procedure, public :: set_num_mins_eq_min_zone
        procedure, public :: set_num_mins_kin_min_zone
        procedure, public :: set_num_mins_kin_cst_act_min_zone
        procedure, public :: set_num_mins_eq_cst_act_min_zone
        procedure, public :: set_num_mins_eq_var_act_min_zone
        procedure, public :: set_num_mins_kin_var_act_min_zone
        procedure, public :: set_chem_syst_min_zone
        procedure, public :: set_ind_min_Sk
        procedure, public :: set_id
        !procedure, public :: allocate_minerals_eq
        !procedure, public :: allocate_minerals_min_zone
        procedure, public :: allocate_ind_chem_syst_min_zone
        procedure, public :: allocate_ind_min_Sk
        !procedure, public :: set_minerals_eq
        !procedure, public :: update_mineral_zone
        !procedure, public :: compare_min_zones
        procedure, public :: is_mineral_in_min_zone
    end type
!**************************************************************************************************
    interface
        subroutine update_mineral_zone(this,old_min_ind)
            import mineral_zone_c
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in) :: old_min_ind(:)
        end subroutine
    end interface
    
    
    
    contains
        subroutine set_num_mins_min_zone(this,num_minerals)
            implicit none
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_minerals
            if (present(num_minerals)) then
                this%num_minerals=num_minerals
            else
                this%num_minerals=this%chem_syst%num_minerals
            end if
        end subroutine
        
        subroutine set_num_mins_kin_min_zone(this,num_min_kin)
            implicit none
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_min_kin
            if (present(num_min_kin)) then
                this%num_minerals_kin=num_min_kin
            else
                this%num_minerals_kin=this%num_minerals-this%num_minerals_eq
            end if
        end subroutine
        
        subroutine set_num_mins_eq_min_zone(this,num_min_eq)
            implicit none
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in) :: num_min_eq
            this%num_minerals_eq=num_min_eq
        end subroutine
        
        subroutine set_num_mins_kin_cst_act_min_zone(this,num_mins_kin_cst_act)
            implicit none
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_kin_cst_act
            this%num_minerals_kin_cst_act=num_mins_kin_cst_act
        end subroutine
        
        subroutine set_num_mins_kin_var_act_min_zone(this,num_mins_kin_var_act)
            implicit none
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_kin_var_act
            this%num_minerals_kin_var_act=num_mins_kin_var_act
        end subroutine

        subroutine set_num_mins_eq_cst_act_min_zone(this,num_mins_eq_cst_act)
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_eq_cst_act
            this%num_minerals_eq_cst_act=num_mins_eq_cst_act
        end subroutine
        
        subroutine set_chem_syst_min_zone(this,chem_syst)
            implicit none
            class(mineral_zone_c) :: this
            class(chem_system_c), intent(in), target :: chem_syst
            this%chem_syst=>chem_syst
        end subroutine
        
        subroutine allocate_ind_chem_syst_min_zone(this,num_mins)
           implicit none
           class(mineral_zone_c) :: this
           integer(kind=4), intent(in), optional :: num_mins
           if (present(num_mins)) then
               this%num_minerals=num_mins
            end if
            allocate(this%ind_min_chem_syst(this%num_minerals))
        end subroutine
        
        subroutine allocate_ind_min_Sk(this,num_min_kin)
            implicit none
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_min_kin !> number of kinetic minerals
            if (present(num_min_kin)) then
                call this%set_num_mins_kin_min_zone(num_min_kin)
            end if
            !allocate(this%ind_min_kin_chem_syst(this%num_minerals_kin))
            allocate(this%ind_min_Sk(this%num_minerals_kin))
        end subroutine
        
      
        
        !subroutine set_minerals_eq(this,minerals_eq)
        !    implicit none
        !    class(mineral_zone_c) :: this
        !    class(mineral_c), intent(in) :: minerals_eq(:)
        !    integer(kind=4) :: i,j
        !    
        !    if (size(minerals_eq)>this%num_minerals) error stop "Dimension error in minerals in equilibrium"
        !    this%minerals_eq=minerals_eq
        !end subroutine
        
        ! subroutine compare_min_zones(min_zone_1,min_zone_2,flag)
        !     implicit none
        !     class(mineral_zone_c), intent(in) :: min_zone_1
        !     class(mineral_zone_c), intent(in) :: min_zone_2
        !     logical, intent(out) :: flag !> TRUE if same minerals & same properties, FALSE otherwise
            
        !     integer(kind=4) :: i,j
            
        !     flag=.true.
        !     if (min_zone_1%num_minerals/=min_zone_2%num_minerals) then
        !         flag=.false.
        !     else if (min_zone_1%num_minerals_eq/=min_zone_2%num_minerals_eq) then
        !         flag=.false.
        !     else if (min_zone_1%num_minerals_eq_cst_act /= min_zone_2%num_minerals_eq_cst_act) then
        !         flag = .false.
        !     else
        !         j=1 !> index of mineral in min_zone_2
        !         do i=1,min_zone_1%num_minerals
        !             if (min_zone_1%ind_chem_syst(i)/=min_zone_2%ind_chem_syst(j)) then
        !                 j=j+1
        !                 if (j>min_zone_2%num_minerals) then
        !                     flag=.false.
        !                     exit
        !                 end if
        !             else
        !                 j=1
        !             end if
        !         end do
        !     end if
        ! end subroutine
        
        ! subroutine allocate_ind_min_kin(this)
        !     implicit none
        !     class(mineral_zone_c) :: this
        !     allocate(this%ind_min_kin(this%num_minerals_kin))
        ! end subroutine
        
        subroutine is_mineral_in_min_zone(this,mineral,flag,mineral_ind)
        !> This subroutine checks if a mineral belongs to the mineral zone
            class(mineral_zone_c), intent(in) :: this                   !< mineral zone
            class(mineral_c), intent(in) :: mineral                     !< mineral
            logical, intent(out) :: flag                                !< TRUE if mineral belongs to mineral zone, FALSE otherwise
            integer(kind=4), intent(out), optional :: mineral_ind       !> index of mineral in "ind_chem_syst" attribute (if not belongs: 0)
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(mineral_ind)) then
                mineral_ind=0
            end if
            do i=1,this%num_minerals
                if (mineral%name==this%chem_syst%minerals(this%ind_min_chem_syst(i))%name) then
                    flag=.true.
                    if (present(mineral_ind)) then
                        mineral_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine

        subroutine set_ind_min_Sk(this,ind_min_Sk)
            !> This subroutine sets the indices of kinetic minerals in the kinetic stoichiometric matrix of the chemical system
            class(mineral_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_min_Sk(:)
            integer(kind=4) :: i,num_min_kin,num_sp
            if (present(ind_min_Sk)) then
                this%ind_min_Sk=ind_min_Sk
            else
                num_min_kin=0 !> counter kinetic minerals
                num_sp=this%chem_syst%speciation_alg%num_aq_prim_species
                !> Variable activity kinetic minerals
                do i=1,this%num_minerals_kin_var_act
                    this%ind_min_Sk(num_min_kin+i)=num_sp+i
                end do
                num_min_kin=num_min_kin+this%num_minerals_kin_var_act
                num_sp=num_sp+this%num_minerals_kin_var_act+this%chem_syst%cat_exch%num_surf_compl+&
                    this%chem_syst%gas_phase%num_var_act_species+this%chem_syst%speciation_alg%num_sec_aq_species
                !> Constant activity kinetic minerals
                do i=1,this%num_minerals_kin_cst_act
                    this%ind_min_Sk(num_min_kin+i)=num_sp+this%ind_min_chem_syst(i)
                end do
            end if
        end subroutine
        
        subroutine set_num_mins_eq_var_act_min_zone(this,num_mins_eq_var_act)
        !> This subroutine sets the number of minerals in equilibrium with variable activity
        class(mineral_zone_c) :: this
        integer(kind=4), intent(in) :: num_mins_eq_var_act
        this%num_minerals_eq_var_act=num_mins_eq_var_act
        end subroutine
        
        subroutine set_id(this,id)
        class(mineral_zone_c) :: this
        integer(kind=4), intent(in) :: id
        this%id=id
        end subroutine
end module