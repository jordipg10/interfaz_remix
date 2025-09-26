!> Set of heterogeneous reactions in equilibrium and "non flowing" species that define these reactions
module reactive_zone_Lagr_m
    use chem_system_m, only: chem_system_c, species_c, eq_reaction_c, gas_phase_c, cat_exch_c, speciation_algebra_c, aq_phase_c, &
        mineral_c, aq_species_c, species_c, gas_c, inf_norm_vec_real, reaction_c, inf_norm_vec_real, LU_lin_syst, id_matrix
    use CV_params_m, only: CV_params_s
    implicit none
    save
    type, public :: reactive_zone_c
        integer(kind=4) :: num_non_flowing_species=0 !> number of non-flowing species
        type(species_c), allocatable :: non_flowing_species(:) !> non-flowing species
        integer(kind=4) :: num_minerals=0 !> number of minerals in equilibrium
        integer(kind=4) :: num_minerals_cst_act=0 !> number of minerals with constant activity
        integer(kind=4) :: num_minerals_var_act=0 !> number of minerals with variable activity
        !type(mineral_c), allocatable :: minerals(:) !> minerals in equilibrium
        integer(kind=4), allocatable :: ind_mins_chem_syst(:) !> indices of minerals in chemical system (first var act, then cst act)
        type(cat_exch_c) :: cat_exch_zone !> cation exchange zone (AQUI PODRIAS USAR INDICES TAMBIEN)
        integer(kind=4) :: num_solids=0 !> number of solids
        type(gas_phase_c) :: gas_phase !> gas phase (AQUI PODRIAS USAR INDICES TAMBIEN)
        real(kind=8), allocatable :: stoich_mat(:,:) !> stoichiometric matrix
        integer(kind=4), allocatable :: ind_mins_stoich_mat(:) !> indices minerals in stoichiometric matrix (same order as in "ind_mins_chem_syst")
        integer(kind=4), allocatable :: ind_gases_stoich_mat(:) !> indices gases in stoichiometric matrix
        !real(kind=8), allocatable :: stoich_mat_sol(:,:) !> solid stoichiometric matrix
        !real(kind=8), allocatable :: stoich_mat_gas(:,:) !> gas stoichiometric matrix
        !integer(kind=4) :: num_eq_reactions=0 !> number of equilibrium reactions
        !type(eq_reaction_c), allocatable :: eq_reactions(:) !> equilibrium heterogeneous reactions
        integer(kind=4), allocatable :: ind_eq_reacts(:) !> indices of equilibrium reactions in chemical system
        class(chem_system_c), pointer :: chem_syst !>  (same chemical system as chemistry class)
        type(speciation_algebra_c) :: speciation_alg !> speciation algebra object
        class(CV_params_s), pointer :: CV_params !> convergence parameters for speciation and reactive mixing computations (chapuza)
        real(kind=8), allocatable :: U_SkT_prod(:,:) !> =U*S_k,nc^T
        integer(kind=4), allocatable :: ind_aq_var_act_species(:) !> indices of aqueous variable activty species in aqueous phase of chemical system (sorted in primary and secondary species)
    contains
    !> Set
        procedure, public :: set_ind_mins_chem_syst
        procedure, public :: set_non_flowing_species
        procedure, public :: set_non_flowing_species_from_chem_syst
        procedure, public :: set_single_non_flowing_species
        procedure, public :: set_num_non_flowing_species
        procedure, public :: set_num_solids
        procedure, public :: set_chem_syst_react_zone
        !procedure, public :: set_num_eq_reactions
        procedure, public :: set_cat_exch_zone
        procedure, public :: set_ind_eq_reacts
        procedure, public :: set_stoich_mat_react_zone
        procedure, public :: set_ind_mins_stoich_mat
        procedure, public :: set_ind_gases_stoich_mat
        procedure, public :: set_num_mins_cst_act
        procedure, public :: set_num_mins_var_act
        procedure, public :: set_num_mins
        procedure, public :: set_gas_phase
        procedure, public :: set_speciation_alg_dimensions
        procedure, public :: set_CV_params
        procedure, public :: set_ind_aq_var_act_species
    !> Allocate/deallocate
        procedure, public :: allocate_non_flowing_species
        procedure, public :: allocate_ind_mins
        procedure, public :: allocate_ind_eq_reacts
        !procedure, public :: allocate_ind_mins_stoich_mat
        procedure, public :: allocate_ind_gases_stoich_mat
        procedure, public :: deallocate_react_zone
    !> Update
        procedure, public :: update_num_eq_reacts
        !procedure, public :: update_reactive_zone
        procedure, public :: update_eq_reactions
    !> Get
        procedure, public :: get_eq_csts_react_zone
    !> Compute
        !procedure, public :: compute_num_eq_reactions
        !procedure, public :: compute_num_cst_act_species_react_zone
        procedure, public :: compute_U_SkT_prod
        procedure, public :: compute_speciation_alg_arrays
    !> Is
        procedure, public :: is_nf_species_in_react_zone
        procedure, public :: is_mineral_in_react_zone
    !> ASsign
        procedure, public :: assign_react_zone
    !> Rearrange
        procedure, public :: rearrange_non_flowing_species
    end type
!**************************************************************************************************
    interface
         subroutine set_stoich_mat_react_zone(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
         end subroutine
        
         subroutine set_stoich_mat_sol_rz(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
         end subroutine
         
         subroutine set_stoich_mat_gas_rz(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
        end subroutine
        
        !subroutine compute_rk_mat(this,conc,rk_mat)
        !>    import reactive_zone_c
        !>    implicit none
        !>    class(reactive_zone_c), intent(in) :: this
        !>    real(kind=8), intent(in) :: conc(:,:)
        !>    real(kind=8), intent(out) :: rk_mat(:,:)
        !end subroutine
        !
        !subroutine compute_rk_vec(this,conc,rk_vec)
        !>    import reactive_zone_c
        !>    implicit none
        !>    class(reactive_zone_c), intent(in) :: this
        !>    real(kind=8), intent(in) :: conc(:)
        !>    real(kind=8), intent(out) :: rk_vec(:)
        !end subroutine
        
        
        

        
       
        
        !subroutine compute_Delta_t_crit_reactive_zone(this,X_mat,F_mat,Delta_t_crit)
        !    import reactive_zone_c
        !    import tridiag_matrix_c
        !    import diag_matrix_c
        !    implicit none
        !    class(reactive_zone_c) :: this
        !    class(tridiag_matrix_c), intent(in) :: X_mat
        !    class(diag_matrix_c), intent(in) :: F_mat
        !    real(kind=8), intent(out) :: Delta_t_crit
        !end subroutine
        
        subroutine read_reactive_zone_Lagr(this,filename,line)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
            character(len=*), intent(in) :: filename !> nombre del archivo de entrada
            integer(kind=4), intent(in) :: line
        end subroutine
        
        subroutine write_reactive_zone(this)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
        end subroutine
        
       
        

        
       
        
        ! subroutine update_reactive_zone(this,old_nf_ind,new_react_zone)
        !     import reactive_zone_c
        !     implicit none
        !     class(reactive_zone_c) :: this
        !     integer(kind=4), intent(in) :: old_nf_ind(:)
        !     class(reactive_zone_c), intent(out) :: new_react_zone
        ! end subroutine
        
        subroutine update_eq_reactions(this,old_eq_reacts_ind)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: old_eq_reacts_ind(:)
        end subroutine
        
        subroutine compare_react_zones(react_zone_1,react_zone_2,flag)
            import reactive_zone_c
            implicit none
            class(reactive_zone_c), intent(in) :: react_zone_1
            class(reactive_zone_c), intent(in) :: react_zone_2
            logical, intent(out) :: flag !> true if same non flowing species, false otherwise
        end subroutine
    end interface
    
    
    
    contains
        subroutine set_num_non_flowing_species(this,num_non_flowing_species)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_non_flowing_species
            if (present(num_non_flowing_species)) then
                this%num_non_flowing_species=num_non_flowing_species
            else
                this%num_non_flowing_species=this%num_minerals+this%cat_exch_zone%num_surf_compl+this%gas_phase%num_gases_eq
            end if
        end subroutine
        
        subroutine set_num_solids(this,num_solids)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_solids
            if (present(num_solids)) then
                this%num_solids=num_solids
            else
                this%num_solids=this%num_minerals+this%cat_exch_zone%num_surf_compl
            end if
        end subroutine
        
        subroutine set_num_mins_cst_act(this,num_mins_cst_act)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_cst_act
            this%num_minerals_cst_act=num_mins_cst_act
        end subroutine
        
        subroutine set_num_mins_var_act(this,num_mins_var_act)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: num_mins_var_act
            this%num_minerals_var_act=num_mins_var_act
        end subroutine
        
        subroutine set_gas_phase(this,gas_phase)
            implicit none
            class(reactive_zone_c) :: this
            class(gas_phase_c), intent(in) :: gas_phase
            this%gas_phase=gas_phase
        end subroutine
        
        subroutine check_num_solids(this)
            implicit none 
            class(reactive_zone_c) :: this
            if (this%num_solids/=this%num_minerals+this%cat_exch_zone%num_surf_compl) then
                error stop "Wrong number of solids in reactive zone object"
            end if
        end subroutine
        
        subroutine allocate_non_flowing_species(this,num_non_flowing_species)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_non_flowing_species
            if (present(num_non_flowing_species)) then
                call this%set_num_non_flowing_species(num_non_flowing_species)
            else
                call this%set_num_non_flowing_species()
            end if
            if (allocated(this%non_flowing_species)) then
                deallocate(this%non_flowing_species)
            end if
            allocate(this%non_flowing_species(this%num_non_flowing_species))
        end subroutine
        
        subroutine set_non_flowing_species_from_chem_syst(this,non_flowing_species_ind)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: non_flowing_species_ind(:)
            integer(kind=4) :: i
            call this%set_num_non_flowing_species(size(non_flowing_species_ind))
            call this%allocate_non_flowing_species()
            do i=1,this%num_non_flowing_species
                call this%non_flowing_species(i)%assign_species(this%chem_syst%species(non_flowing_species_ind(i)))
            end do
        end subroutine
        !
        
        subroutine set_non_flowing_species(this,non_flowing_species)
        !<      constant activity minerals
        !<      constant activity gases
        !<      variable activity minerals
        !<      surface complexes
        !<      variable activity gases
            implicit none
            class(reactive_zone_c) :: this
            class(species_c), intent(in), optional :: non_flowing_species(:)
            
            integer(kind=4) :: i
            
            if (present(non_flowing_species)) then
                this%non_flowing_species=non_flowing_species
            else
                ! j=0 !> counter constant non flowing species
                ! k=0 !> counter variable non flowing species
                call this%allocate_non_flowing_species(this%num_solids+this%gas_phase%num_gases_eq)
                !> First: minerals with constant activity
                do i=1,this%num_minerals_cst_act
                    call this%non_flowing_species(i)%assign_species(this%chem_syst%minerals(this%ind_mins_chem_syst(i))%mineral)
                end do
                !> Second: minerals with variable activity
                do i=1,this%num_minerals_var_act
                    call this%non_flowing_species(this%num_minerals_cst_act+i)%assign_species(this%chem_syst%minerals(&
                        this%ind_mins_chem_syst(this%num_minerals_cst_act+i))%mineral)
                end do
                !> Third: surface complexes
                do i=1,this%cat_exch_zone%num_surf_compl
                    call this%non_flowing_species(this%num_minerals+i)%assign_species(this%cat_exch_zone%surf_compl(i))
                end do
                !> Fourth: gases with constant activity
                do i=1,this%gas_phase%num_gases_eq_cst_Act
                    call this%non_flowing_species(this%num_solids+i)%assign_species(this%gas_phase%gases(i))
                end do
                !> Fifth: gases with variable activity
                do i=1,this%gas_phase%num_gases_eq_var_act
                    call this%non_flowing_species(this%num_solids+this%gas_phase%num_gases_eq_cst_act+i)%assign_species(&
                        this%gas_phase%gases(this%gas_phase%num_gases_eq_cst_act+i))
                end do
            end if
        end subroutine
        
        subroutine rearrange_non_flowing_species(this)
        !< This subroutine rearranges the "non_flowing_species" attribute in the following order:
        !<      constant activity minerals
        !<      constant activity gases
        !<      variable activity minerals
        !<      surface complexes
        !<      variable activity gases
            implicit none
            class(reactive_zone_c) :: this
            
            integer(kind=4) :: i,j,k,l,m,n
            type(species_c), allocatable :: old_non_flow_species(:)
            LOGICAL :: flag_gas,flag_surf
            
            old_non_flow_species=this%non_flowing_species
            deallocate(this%non_flowing_species)
            call this%allocate_non_flowing_species()
            
            j=0 !> counter constant activity gases
            k=0 !> counter variable activity gases
            l=0 !> counter surface complexes
            m=0 !> counter constant activity minerals
            n=0 !> counter variable activity minerals
            do i=1,this%num_non_flowing_species
                call old_non_flow_species(I)%is_gas(flag_gas)
                if (flag_gas .eqv. .true.) then
                    if (old_non_flow_species(i)%cst_act_flag.eqv..true.) then
                        j=j+1
                        call this%non_flowing_species(this%num_minerals_cst_act+j)%assign_species(old_non_flow_species(I))
                    else
                        k=k+1
                        call this%non_flowing_species(&
                        this%num_non_flowing_species-this%gas_phase%num_gases_eq_var_act+k)%assign_species(&
                        old_non_flow_species(I))
                    end if
                else
                    call old_non_flow_species(I)%is_surf_compl(flag_surf)
                    if (flag_surf .eqv. .true.) then
                        l=l+1
                        call this%non_flowing_species(this%num_minerals+this%gas_phase%num_gases_eq_cst_act+l)%assign_species(&
                        old_non_flow_species(I))
                    else if (old_non_flow_species(i)%cst_act_flag .eqv. .true.) then
                        m=m+1
                        call this%non_flowing_species(m)%assign_species(old_non_flow_species(I))
                    else
                        n=n+1
                    call this%non_flowing_species(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act+n)%assign_species(&
                        old_non_flow_species(I))
                    end if
                end if
            end do
        end subroutine
        
        subroutine set_single_non_flowing_species(this,non_flowing_species_ind,chem_syst_ind)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: non_flowing_species_ind !> index in reactive zone
            integer(kind=4), intent(in) :: chem_syst_ind !> index in chemical system
            if(non_flowing_species_ind<1 .or. non_flowing_species_ind>chem_syst_ind .or. chem_syst_ind>this%chem_syst%num_species) &
                error stop "Index out of bounds"
            call this%non_flowing_species(non_flowing_species_ind)%assign_species(this%chem_syst%species(chem_syst_ind))
        end subroutine
        
        subroutine set_chem_syst_react_zone(this,chem_syst)
            implicit none
            class(reactive_zone_c) :: this
            class(chem_system_c), intent(in), target :: chem_syst
            this%chem_syst=>chem_syst
        end subroutine
        
        function get_eq_csts_react_zone(this,cst_act_gases) result(K)
            implicit none
            class(reactive_zone_c), intent(in) :: this !> reactive zone object
            real(kind=8), intent(in), optional :: cst_act_gases(:) !> constant activity gases in gas phase
            real(kind=8), allocatable :: K(:) !> equilibrium constants
            
            integer(kind=4) :: i,ind_eq_react
          
            allocate(K(this%speciation_alg%num_eq_reactions))
            
            ind_eq_react=0 !> counter equilibrium reactions
            
            do i=1,this%num_minerals_cst_act
                K(i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(i))%eq_cst
            end do
            ind_eq_react=ind_eq_react+this%num_minerals_cst_act
            if (present(cst_act_gases)) then
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    K(this%num_minerals_cst_act+i) = this%chem_syst%eq_reacts(this%ind_eq_reacts(&
                        this%num_minerals_cst_act+i))%eq_cst/cst_act_gases(i)
                end do
            else
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    K(this%num_minerals_cst_act+i) = this%chem_syst%eq_reacts(this%ind_eq_reacts(&
                        this%num_minerals_cst_act+i))%eq_cst
                end do
            end if
            ind_eq_react=ind_eq_react+this%gas_phase%num_gases_eq_cst_act
            do i=ind_eq_react+1,this%speciation_alg%num_eq_reactions
                K(i)=this%chem_syst%eq_reacts(this%ind_eq_reacts(i))%eq_cst
            end do
           
        end function
        
        !subroutine get_minerals_eq(this)
        !>    implicit none
        !>    class(reactive_zone_c) :: this
        !>    
        !>    integer(kind=4) :: i
        !>    !do i=1,this%mineral_zone%num_minerals_eq
        !>    !>    this%non_flowing_species(i)%name=this%mineral_zone%minerals_eq(i)%name
        !>    !>    this%non_flowing_species(i)%cst_act_flag=.true.
        !>    !>    this%non_flowing_species(i)%valence=0
        !>    !end do
        !>    i=1 !> counter minerals reactive zone
        !>    j=1 !> counter equilibrium reactions chemical system
        !>    do
        !>        call this%chem_syst%eq_reacts(j)%is_species_in_react(this%non_flowing_species(i),flag,sp_ind)
        !>        !call this%chem_syst%is_eq_reaction_in_chem_syst(this%non_flowing_species(i),flag,eq_react_ind)
        !>        if (flag.eqv..true.) then
        !>            this%eq_reactions(i)=this%chem_syst%eq_reacts(j)
        !>            if (i<this%num_eq_recactions) then
        !>                i=i+1
        !>                j=1
        !>            else
        !>                exit
        !>            end if
        !>            !call append_int_1D_array(eq_react_indices,eq_react_ind)
        !>        else if (j<this%chem_syst%num_eq_reacts) then
        !>            j=j+1
        !>        else
        !>            error stop "This equilibrium reaction is not in the chemical system"
        !>        end if
        !>    end do
        !end subroutine
       
        
        
        
      
       
        
        subroutine set_ind_eq_reacts(this,ind_eq_reacts)
        !< This subroutine sets the "ind_eq_reacts" attribute in the following order:
        !< If "flag_comp" attribute in "speciation_alg" attribute of this class is TRUE:
        !<      constant activity minerals
        !<      constant activity gases
        !<      redox
        !<      aqueous complexes
        !<      variable activity minerals
        !<      cation exchange
        !<      variable activity gases
        !< Else:
        !<      redox
        !<      aqueous complexes
        !<      minerals
        !<      cation exchange
        !<      gases
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_eq_reacts(:) !> indices of equilibrium reactions in chemical system
            
            integer(kind=4) :: i,l,j,k,eq_react_ind,n_eq,sp_ind
            integer(kind=4), allocatable :: eq_react_indices(:)
            logical :: flag

            call this%allocate_ind_eq_reacts()
            if (present(ind_eq_reacts)) then
                !call this%allocate_eq_reactions(size(eq_reactions_ind))
                if (size(ind_eq_reacts)/=this%speciation_alg%num_eq_reactions) error stop
                this%ind_eq_reacts=ind_eq_reacts
            else if ((this%speciation_alg%flag_comp .eqv. .true.) .and. (this%chem_syst%speciation_alg%flag_comp .eqv. .true.)) then
                do i=1,this%chem_syst%num_aq_eq_reacts
                    this%ind_eq_reacts(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act+i)=&
                        this%chem_syst%num_minerals_eq_cst_act+this%chem_syst%gas_phase%num_gases_eq_cst_act+i
                end do
                i=1 !> counter constant activity minerals & gases reactive zone
                l=1 !> counter variable activity minerals, surface complexes & gases reactive zone
                j=1 !> counter equilibrium reactions chemical system
                k=1 !> counter non flowing species   
                if (this%num_non_flowing_species>0) then
                    do
                        if (this%non_flowing_species(k)%name=='x-') then !> chapuza
                            k=k+1
                        end if
                        call this%chem_syst%eq_reacts(j)%is_species_in_react(this%non_flowing_species(k),flag,sp_ind)
                        if (flag .eqv. .true.) then
                            if (this%non_flowing_species(k)%cst_act_flag .eqv. .true.) then
                                this%ind_eq_reacts(i)=j
                                i=i+1
                            else
                                this%ind_eq_reacts(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act + & 
                                    this%chem_syst%num_aq_eq_reacts+l)=j
                                l=l+1
                            end if
                            if (k<this%num_non_flowing_species) then
                                k=k+1
                                j=1
                            else
                                exit
                            end if
                        else if (j<this%chem_syst%num_eq_reacts) then
                            j=j+1
                        else
                            print *, this%non_flowing_species(k)%name
                            error stop "This equilibrium reaction is not in the chemical system"
                        end if
                    end do
                end if
            else if ((this%speciation_alg%flag_comp .eqv. .true.) .and. (this%chem_syst%speciation_alg%flag_comp .eqv. .false.))then
                do i=1,this%chem_syst%num_aq_eq_reacts
                    this%ind_eq_reacts(i)=i
                end do
                i=1 !> counter constant activity minerals & gases reactive zone
                l=1 !> counter variable activity minerals, surface complexes & gases reactive zone
                j=1 !> counter equilibrium reactions chemical system
                k=1 !> counter non flowing species   
                if (this%num_non_flowing_species>0) then
                    do
                        if (this%non_flowing_species(k)%name=='x-') then !> chapuza
                            k=k+1
                        end if
                        call this%chem_syst%eq_reacts(j)%is_species_in_react(this%non_flowing_species(k),flag,sp_ind)
                        if (flag .eqv. .true.) then
                            if (this%non_flowing_species(k)%cst_act_flag .eqv. .true.) then
                                this%ind_eq_reacts(this%chem_syst%num_aq_eq_reacts+i)=j
                                i=i+1
                            else
                                this%ind_eq_reacts(this%num_minerals_cst_act+this%gas_phase%num_gases_eq_cst_act + & 
                                    this%chem_syst%num_aq_eq_reacts+l)=j
                                l=l+1
                            end if
                            if (k<this%num_non_flowing_species) then
                                k=k+1
                                j=1
                            else
                                exit
                            end if
                        else if (j<this%chem_syst%num_eq_reacts) then
                            j=j+1
                        else
                            print *, this%non_flowing_species(k)%name
                            error stop "This equilibrium reaction is not in the chemical system"
                        end if
                    end do
                end if
            else
                i=0 !> counter equilibrium reactions reactive zone
                j=1 !> counter equilibrium reactions chemical system
                k=1 !> counter non flowing species
                do l=1,this%chem_syst%num_aq_eq_reacts
                    this%ind_eq_reacts(l)=l
                end do
                i=i+this%chem_syst%num_aq_eq_reacts+1
                if (this%num_non_flowing_species>0) then
                    do
                        if (this%non_flowing_species(k)%name=='x-') then !> chapuza
                            k=k+1
                        end if
                        call this%chem_syst%eq_reacts(j)%is_species_in_react(this%non_flowing_species(k),flag,sp_ind)
                        if (flag .eqv. .true.) then
                            this%ind_eq_reacts(i)=j
                            i=i+1
                            if (k<this%num_non_flowing_species) then
                                k=k+1
                                j=1
                            else
                                exit
                            end if
                        else if (j<this%chem_syst%num_eq_reacts) then
                            j=j+1
                        else
                            error stop "This equilibrium reaction is not in the chemical system"
                        end if
                    end do
                end if
            end if
        end subroutine
                
        subroutine set_ind_mins_chem_syst(this,ind_mins_chem_syst)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_mins_chem_syst(:) !> indices of minerals in chemical system
            integer(kind=4) :: i
            if (present(ind_mins_chem_syst)) then
                if (size(ind_mins_chem_syst)>this%chem_syst%num_minerals_eq) error stop
                this%num_minerals=size(ind_mins_chem_syst)
                this%ind_mins_chem_syst=ind_mins_chem_syst
            else !> default: all minerals in chemical system
                this%num_minerals=this%chem_syst%num_minerals_eq
                do i=1,this%num_minerals
                    this%ind_mins_chem_syst(i)=this%chem_syst%num_minerals_kin+i
                end do
                ! do i=1,this%num_minerals
                !     if (this%minerals(i)%mineral%cst_act_flag.eqv..true.) then
                !         this%num_minerals_cst_act=this%num_minerals_cst_act+1
                !     else
                !         this%num_minerals_var_act=this%num_minerals_var_act+1
                !     end if
                ! end do
            end if
        end subroutine
        
        subroutine allocate_ind_mins(this,num_minerals)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: num_minerals !> number of minerals in reactive zone
            if (present(num_minerals)) then
                if (num_minerals<0 .or. num_minerals>this%chem_syst%num_minerals_eq) error stop
                this%num_minerals=num_minerals
            end if
            if (allocated(this%ind_mins_chem_syst)) deallocate(this%ind_mins_chem_syst)
            if (allocated(this%ind_mins_stoich_mat)) deallocate(this%ind_mins_stoich_mat)
            allocate(this%ind_mins_chem_syst(this%num_minerals),this%ind_mins_stoich_mat(this%num_minerals))
        end subroutine
        
        !subroutine set_mineral_zone(this,mineral_zone)
        !>    implicit none
        !>    class(reactive_zone_c) :: this
        !>    class(mineral_zone_c), intent(in) :: mineral_zone
        !>    this%mineral_zone=mineral_zone
        !end subroutine
        
        !subroutine set_num_eq_reactions(this,num_eq_reacts)
        !    implicit none
        !    class(reactive_zone_c) :: this
        !    integer(kind=4), intent(in) :: num_eq_reacts
        !    this%speciation_alg%num_eq_reactions=num_eq_reacts
        !end subroutine
        
        !subroutine set_num_kin_reactions(this,num_kin_reacts)
        !>    implicit none
        !>    class(reactive_zone_c) :: this
        !>    integer(kind=4), intent(in) :: num_kin_reacts
        !>    this%num_kin_reactions=num_kin_reacts
        !end subroutine
        
        subroutine update_num_eq_reacts(this,num_old_eq_reacts)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in) :: num_old_eq_reacts
            this%speciation_alg%num_eq_reactions=this%speciation_alg%num_eq_reactions-num_old_eq_reacts
        end subroutine
        
        subroutine allocate_ind_eq_reacts(this)
            implicit none
            class(reactive_zone_c) :: this
            if (allocated(this%ind_eq_reacts)) then
                deallocate(this%ind_eq_reacts)
            end if
            allocate(this%ind_eq_reacts(this%speciation_alg%num_eq_reactions))
        end subroutine
        
        !subroutine allocate_kin_reactions(this,num_kin_reactions)
        !>    implicit none
        !>    class(reactive_zone_c) :: this
        !>    integer(kind=4), intent(in), optional :: num_kin_reactions
        !>    if (present(num_kin_reactions)) then
        !>        call this%set_num_kin_reactions(num_kin_reactions)
        !>    end if
        !>    allocate(this%kin_reactions(this%num_kin_reactions))
        !end subroutine
        
        subroutine is_nf_species_in_react_zone(this,nf_species,flag,nf_species_ind)
            implicit none
            class(reactive_zone_c), intent(in) :: this
            class(species_c), intent(in) :: nf_species
            logical, intent(out) :: flag
            integer(kind=4), intent(out) :: nf_species_ind
            
            integer(kind=4) :: i
            
            nf_species_ind=0
            flag=.false.
            do i=1,this%num_non_flowing_species
                if (nf_species%name==this%non_flowing_species(i)%name) then
                    flag=.true.
                    nf_species_ind=i
                    exit
                end if
            end do
        end subroutine
        
        !function compute_num_species_react_zone(this) result(n_sp)
        !    implicit none
        !    class(reactive_zone_c), intent(in) :: this
        !    integer(kind=4) :: n_sp
        !    n_sp=this%chem_syst%aq_phase%num_species+this%num_non_flowing_species+this%chem_syst%num_minerals_kin
        !end function
        
        !function compute_num_cst_act_species_react_zone(this) result(n_c)
        !    implicit none
        !    class(reactive_zone_c), intent(in) :: this
        !    integer(kind=4) :: n_c
        !    logical :: flag
        !    n_c=this%num_non_flowing_species    !> we assume:    1) all non-flowing species are pure minerals
        !    call this%chem_syst%aq_phase%is_water_in_aq_phase(flag)
        !    if (flag.eqv..true.) then
        !        n_c=n_c+1 
        !    end if
        !end function
        
        
        
        subroutine is_mineral_in_react_zone(this,mineral,flag,index)
            class(reactive_zone_c), intent(in) :: this !> reactive zone
            class(mineral_c), intent(in) :: mineral !> mineral
            logical, intent(out) :: flag !> true if mineral is in reactive zone, false otherwise
            integer(kind=4), intent(out), optional :: index !> index of mineral in reactive zone
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(index)) then
                index=0
            end if
            do i=1,this%num_minerals
                if (mineral%name==this%chem_syst%minerals(this%ind_mins_chem_syst(i))%name .and. mineral%mineral%cst_act_flag.eqv.&
                    this%chem_syst%minerals(this%ind_mins_chem_syst(i))%mineral%cst_act_flag) then
                    flag=.true.
                    if (present(index)) then
                        index=i
                    end if
                    exit
                end if
            end do
        end subroutine
        
        !subroutine set_prim_species_indices(this,num_prim_species)
        !>    implicit none
        !>    class(reactive_zone_c) :: this
        !>    integer(kind=4), intent(in) :: num_prim_species
        !>    
        !>    integer(kind=4) :: i,j,k
        !>    logical :: flag
        !>    
        !>    allocate(this%prim_species_indices(num_prim_species))
        !>    
        !>    k=1 !> counter primary species
        !>    do i=1,this%num_eq_reactions
        !>        do j=1,num_prim_species
        !>            call this%eq_reactions(i)%is_species_in_react(this%chem_syst%aq_phase%aq_species(j),flag)
        !>            if (flag.eqv..true.) then
        !>                this%prim_species_indices(k)=j
        !>                k=k+1
        !>            else
        !>                continue
        !>            end if
        !>        end do
        !>    end do
        !end subroutine
        
        subroutine set_cat_exch_zone(this,cat_exch_zone)
            implicit none
            class(reactive_zone_c) :: this
            class(cat_exch_c), intent(in), optional :: cat_exch_zone
            
            integer(kind=4) :: i
            logical :: flag
            
            if (present(cat_exch_zone)) then
                do i=1,cat_exch_zone%num_surf_compl
                    call this%chem_syst%cat_exch%is_surf_compl_in(cat_exch_zone%surf_compl(i),flag)
                    if (flag.eqv..false.) then
                        error stop "Surface complex not in chemical system"
                    end if
                end do
                this%cat_exch_zone=cat_exch_zone
            else
                this%cat_exch_zone=this%chem_syst%cat_exch
            end if
        end subroutine
        
        subroutine deallocate_react_zone(this)
            implicit none
            class(reactive_zone_c) :: this
            deallocate(this%ind_mins_chem_syst)
            deallocate(this%non_flowing_species)
            deallocate(this%stoich_mat)
            deallocate(this%ind_mins_stoich_mat)
            deallocate(this%ind_gases_stoich_mat)
            deallocate(this%ind_eq_reacts)
            deallocate(this%U_SkT_prod)
        end subroutine
        
        subroutine assign_react_zone(this,react_zone) !> chapuza 
            implicit none
            class(reactive_zone_c) :: this !< reactive zone
            class(reactive_zone_c), intent(in) :: react_zone !< reactive zone to be assigned
            if (associated(react_zone%chem_syst)) then
                this%chem_syst=>react_zone%chem_syst
            else
                error stop "Chemical system not associated with reactive zone"
            end if
            if (associated(react_zone%CV_params)) then
                this%CV_params=>react_zone%CV_params
            else
                error stop "CV parameters not associated with reactive zone"
            end if
            if (allocated(react_zone%ind_mins_chem_syst)) then
                this%ind_mins_chem_syst=react_zone%ind_mins_chem_syst
            end if
            if (allocated(react_zone%non_flowing_species)) then
                this%non_flowing_species=react_zone%non_flowing_species
            end if
            this%num_non_flowing_species=react_zone%num_non_flowing_species
            this%gas_phase=react_zone%gas_phase
            this%num_minerals=react_zone%num_minerals
            this%num_minerals_cst_Act=react_zone%num_minerals_cst_Act
            this%num_minerals_var_Act=react_zone%num_minerals_var_Act
            this%cat_exch_zone=react_zone%cat_exch_zone
            this%num_solids=react_zone%num_solids
            if (allocated(react_zone%stoich_mat)) then
                this%stoich_mat=react_zone%stoich_mat
            end if
            if (allocated(react_zone%ind_mins_stoich_mat)) then
                this%ind_mins_stoich_mat=react_zone%ind_mins_stoich_mat
            end if
            if (allocated(react_zone%ind_gases_stoich_mat)) then
                this%ind_gases_stoich_mat=react_zone%ind_gases_stoich_mat
            end if
            this%speciation_alg=react_zone%speciation_alg
            if (allocated(react_zone%ind_eq_reacts)) then
                call this%set_ind_eq_reacts(react_zone%ind_eq_reacts)
            end if
            if (react_zone%chem_syst%num_kin_reacts>0 .and. allocated(react_zone%U_SkT_prod)) then
                this%U_SkT_prod=react_zone%U_SkT_prod
            end if
        end subroutine
        
        subroutine set_speciation_alg(this,speciation_alg)
            implicit none
            class(reactive_zone_c) :: this
            type(speciation_algebra_c), intent(in) :: speciation_alg
            this%speciation_alg=speciation_alg
        end subroutine
        
        subroutine set_speciation_alg_dimensions(this,flag_comp)
            implicit none
            class(reactive_zone_c) :: this
            logical, intent(in), optional :: flag_comp !> TRUE if component matrix has no constant activity species (De Simoni et al, 2005), FALSE otherwise
            
            integer(kind=4) :: i,n_sp,n_c,n_eq,n_gas_kin
            logical :: flag_cat_exch
            
            n_gas_kin=0
            
            if (.not. associated(this%chem_syst)) then
                error stop "Chemical system not associated with reactive zone"
            else if (this%num_non_flowing_species>0) then
                n_sp=this%chem_syst%aq_phase%num_species+this%num_non_flowing_species!+this%chem_syst%num_minerals_kin
                n_c=this%chem_syst%aq_phase%wat_flag+this%num_minerals_cst_act+this%gas_phase%num_cst_act_species
                !do i=1,this%num_minerals
                !    if (this%minerals(I)%mineral%cst_act_flag.eqv..true.) then
                !        n_c=n_c+1
                !    end if
                !end do
                !do i=1,this%chem_syst%num_minerals_kin
                !    if (this%chem_syst%minerals(i)%mineral%cst_act_flag.eqv..true.) then
                !        n_c=n_c+1
                !    end if
                !end do
                !do i=1,this%gas_phase%num_species
                !    if (this%gas_phase%gases(i)%cst_act_flag.eqv..true.) then
                !        n_c=n_c+1
                !    end if
                !end do
                n_eq=this%num_minerals+this%gas_phase%num_gases_eq+this%cat_exch_zone%num_exch_cats+ & 
                    this%chem_syst%num_aq_eq_reacts
                if (this%cat_exch_zone%num_surf_compl>0) then
                    flag_cat_exch=.true.
                else
                    flag_cat_exch=.false.
                end if
            else !> all equilibrium reactions are aqueous
                n_sp=this%chem_syst%aq_phase%num_species
                n_c=this%chem_syst%aq_phase%num_cst_act_species
                n_eq=this%chem_syst%num_aq_eq_reacts
                if (this%chem_syst%cat_exch%num_surf_compl>0) then
                    flag_cat_exch=.true.
                else
                    flag_cat_exch=.false.
                end if
            end if
            if (present(flag_comp)) then
                call this%speciation_alg%set_flag_comp(flag_comp)
            end if
            call this%speciation_alg%set_flag_cat_exch(flag_cat_exch)
            call this%speciation_alg%set_dimensions(n_sp,n_eq,n_c,this%chem_syst%aq_phase%num_species, & 
                this%chem_syst%aq_phase%num_species-this%chem_syst%aq_phase%wat_flag)
        end subroutine
        
        subroutine compute_speciation_alg_arrays(this,flag,swap,cst_act_gases)
            implicit none
            class(reactive_zone_c) :: this
            logical, intent(out) :: flag !> TRUE if stoichiometric matrix has been modified, FALSE otherwise
            integer(kind=4), intent(out) :: swap(:) !> species or reactions to swap in stoichiometric matrix (already allocated)
            real(kind=8), intent(in), optional :: cst_act_gases(:) !> constant activity of gases in reactive zone (chapuza)

            real(kind=8), allocatable :: Se(:,:),K(:),aux_Se(:,:),aux_Sk(:,:)
            integer(kind=4) :: aux_col
            type(species_c), allocatable :: aux_species(:)
            type(eq_reaction_c), allocatable :: aux_eq_reacts(:)
            
            !logical :: flag
            !type(aq_phase_c), target :: aux_aq_phase
                        
            !call aq_phase_new%copy_attributes(this%aq_phase)

            flag=.false. !> false by default
            swap=0 !> no swap by default
            
            if (this%speciation_alg%num_eq_reactions>0) then
                !Se=this%stoich_mat
                K=this%get_eq_csts_react_zone(cst_act_gases)
                call this%speciation_alg%compute_arrays(this%stoich_mat,K,this%CV_params%zero,flag,swap)
                if (flag .eqv. .true.) then
                    !allocate(aux_Se(this%speciation_Alg%num_eq_reactions,this%speciation_alg%num_Species))
                    !aux_Se=Se
                    if (this%speciation_alg%flag_comp .eqv. .false.) then
                        !allocate(aux_ind_eq_reacts(2))
                        !aux_ind_eq_reacts=this%ind_eq_reacts(swap)
                        this%ind_eq_reacts(swap(1))=swap(2)
                        this%ind_eq_reacts(swap(2))=swap(1)
                        !this%stoich_mat(swap(1),:)=aux_Se(swap(2),:)
                        !this%stoich_mat(swap(2),:)=aux_Se(swap(1),:)
                        !deallocate(aux_ind_eq_reacts)
                    else
                        allocate(aux_Sk(this%chem_syst%num_species,this%chem_syst%num_species),aux_species(2))
                        aux_species=this%chem_syst%species(swap)
                        this%chem_syst%species(swap(1))=aux_species(2)
                        this%chem_syst%species(swap(2))=aux_species(1)
                        aux_Sk=this%chem_syst%Sk
                        this%chem_syst%Sk(:,swap(1))=aux_Sk(:,swap(2))
                        this%chem_syst%Sk(:,swap(2))=aux_Sk(:,swap(1))
                        ! this%stoich_mat(:,swap(1))=aux_Se(:,swap(2))
                        ! this%stoich_mat(:,swap(2))=aux_Se(:,swap(1))
                        deallocate(aux_species,aux_Sk)
                        ! this%ind_aq_var_act_species(swap(1))=swap(2)
                        ! this%ind_aq_var_act_species(swap(2))=swap(1)
                    end if
                    !deallocate(aux_Se)
                end if
                deallocate(K)
            else if (associated(this%chem_syst)) then
                Se=this%chem_syst%Se
                K=this%chem_syst%get_eq_csts()
                if (this%speciation_alg%flag_comp .eqv. .false.) then
                    this%speciation_alg%comp_mat_cst_act=id_matrix(this%speciation_alg%num_prim_species)
                    call this%speciation_alg%compute_inv_Se_2(Se(&
                        :,this%speciation_alg%num_prim_species+1:this%speciation_alg%num_species),&
                        this%CV_params%zero)
                    call this%speciation_alg%compute_logK_tilde(K)
                else
                    this%speciation_alg%comp_mat=id_matrix(this%speciation_alg%num_prim_species)
                    this%speciation_alg%comp_mat_aq=id_matrix(this%speciation_alg%num_aq_prim_species)
                    call this%speciation_alg%compute_inv_Se_nc_2(Se(&
                        :,this%speciation_alg%num_prim_species+1:this%speciation_alg%num_var_act_species),&
                        this%CV_params%zero)
                    !call this%speciation_alg%compute_logK_tilde(K)
                end if
                !call this%speciation_alg%compute_comp_mat()
                !call this%speciation_alg%compute_arrays(this%chem_syst%Se,K,this%CV_params%zero,flag,swap)
            else
                error stop "Chemical system not associated with reactive zone"
            end if
        end subroutine
        
        subroutine set_CV_params(this,CV_params)
            implicit none
            class(reactive_zone_c) :: this
            class(CV_params_s), intent(in), target :: CV_params
            this%CV_params=>CV_params
        end subroutine
        
        subroutine compute_U_SkT_prod(this,ind_kin)
            implicit none
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_kin(:) !> indices of kinetic reactions in chemical system
            
            if (present(ind_kin)) then
                if (this%speciation_alg%num_eq_reactions==0) then
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_var_act_species,size(ind_kin)))
                    end if
                    this%U_SkT_prod=transpose(this%chem_syst%Sk(ind_kin,1:this%speciation_alg%num_var_act_species)) !> chapuza
                else
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_prim_species,size(ind_kin)))
                    end if
                    this%U_SkT_prod = matmul(this%speciation_alg%comp_mat, &
                                    transpose(this%chem_syst%Sk(ind_kin, 1:this%speciation_alg%num_var_act_species)))
                end if
            else !> we consider all kinetic reactions in the chemical system
                if (this%speciation_alg%num_eq_reactions==0) then
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_var_act_species,this%chem_syst%num_kin_reacts))
                    end if
                    this%U_SkT_prod=transpose(this%chem_syst%Sk(:,1:this%speciation_alg%num_var_act_species)) !> chapuza
                else if (this%chem_syst%num_kin_reacts>0) then
                    if (.not. allocated(this%U_SkT_prod)) then
                        allocate(this%U_SkT_prod(this%speciation_alg%num_prim_species,this%chem_syst%num_kin_reacts))
                    end if
                    this%U_SkT_prod = matmul(this%speciation_alg%comp_mat, &
                                             transpose(this%chem_syst%Sk(:, 1:this%speciation_alg%num_var_act_species)))
                else if (.not. allocated(this%U_SkT_prod)) then
                    allocate(this%U_SkT_prod(this%speciation_alg%num_prim_species,0))
                    this%U_SkT_prod=0d0
                end if
            end if
        end subroutine
        
        !subroutine compute_num_eq_reactions(this)
        !    implicit none
        !    class(reactive_zone_c) :: this

        subroutine set_ind_aq_var_act_species(this) !> sets the indices of variable activity aqueous species in the chemical system
        !! first primary species, then secondary var act species
        !! we assume that the species in the chemical system are already ordered in variable & constant activity
            class(reactive_zone_c) :: this
            integer(kind=4) :: i
            do i=1,this%speciation_alg%num_aq_var_act_species
                this%ind_aq_var_act_species(i)=i !> by default
            end do
        end subroutine

        subroutine set_ind_mins_stoich_mat(this,ind_mins_stoich_mat)
        !> This subroutine sets the "ind_mins_stoich_mat" attribute
            !! We assume that the indices of minerals are already allocated
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_mins_stoich_mat(:) !> indices of solids in stoichiometric matrix
            integer(kind=4) :: i,num_mins,num_sp
            !call this%allocate_ind_mins()
            if (present(ind_mins_stoich_mat)) then
                if (size(ind_mins_stoich_mat)>this%num_solids) error stop
                !this%num_solids=size(ind_mins_stoich_mat)
                this%ind_mins_stoich_mat=ind_mins_stoich_mat
            else
                num_mins=0 !> counter minerals
                num_sp=this%speciation_alg%num_prim_species+this%speciation_alg%num_sec_aq_species !> counter species in stoichiometric matrix
                ! !> Free site
                ! if (this%speciation_alg%flag_cat_exch .eqv. .true.) then
                !     num_sol=num_sol+1
                !     num_sp=num_sp+1
                !     this%ind_mins_stoich_mat(num_sol)=num_sp+1
                ! end if
                !> Variable activity minerals in equilibrium
                do i=1,this%num_minerals_var_act
                    this%ind_mins_stoich_mat(num_mins+i)=num_sp+i
                end do
                num_mins=num_mins+this%num_minerals_var_act
                num_sp=num_sp+this%num_minerals_var_act+this%cat_exch_zone%num_exch_cats+this%gas_phase%num_gases_eq_var_act
                ! !> Surface complexes
                ! do i=1,this%cat_exch_zone%num_exch_cats
                !     this%ind_solids(num_mins+i)=num_sp+i
                ! end do
                ! num_mins=num_mins+this%cat_exch_zone%num_exch_cats
                ! num_sp=num_sp+this%cat_exch_zone%num_exch_cats+this%gas_phase%num_gases_eq_var_act+this%chem_syst%aq_phase%wat_flag
                !> Constant activity minerals in equilibrium
                do i=1,this%num_minerals_cst_act
                    this%ind_mins_stoich_mat(num_mins+i)=num_sp+i
                end do
            end if
        end subroutine

        subroutine set_ind_gases_stoich_mat(this,ind_gases_stoich_mat)
        !> This subroutine sets the "ind_gases_stoich_mat" attribute
            !! We assume that the number of gases in the gas phase attribute is already set
            class(reactive_zone_c) :: this
            integer(kind=4), intent(in), optional :: ind_gases_stoich_mat(:) !> indices of gases in stoichiometric matrix
            integer(kind=4) :: i,num_gas,num_sp
            call this%allocate_ind_gases_stoich_mat()
            if (present(ind_gases_stoich_mat)) then
                if (size(ind_gases_stoich_mat)>this%gas_phase%num_gases_eq) error stop
                !this%num_gases=size(ind_gases)
                this%ind_gases_stoich_mat=ind_gases_stoich_mat
            else
                !this%num_gases=this%gas_phase%num_gases_eq
                num_gas=0 !> counter gases
                num_sp=this%speciation_alg%num_var_act_species-this%gas_phase%num_gases_eq_var_act !> counter species in stoichiometric matrix
                !> Variable activity gases in equilibrium
                do i=1,this%gas_phase%num_gases_eq_var_act
                    this%ind_gases_stoich_mat(num_gas+i)=num_sp+i
                end do
                num_gas=num_gas+this%gas_phase%num_gases_eq_var_act
                num_sp=num_sp+this%speciation_alg%num_cst_act_species-this%gas_phase%num_gases_eq_cst_act
                !> Constant activity gases in equilibrium
                do i=1,this%gas_phase%num_gases_eq_cst_act
                    this%ind_gases_stoich_mat(num_gas+i)=num_sp+i
                end do
            end if
        end subroutine

        ! subroutine allocate_ind_mins_stoich_mat(this)
        ! !> This subroutine allocates the "ind_mins_stoich_mat" attribute
        !     class(reactive_zone_c) :: this
        !     if (allocated(this%ind_mins_stoich_mat)) then
        !         deallocate(this%ind_mins_stoich_mat)
        !     end if
        !     allocate(this%ind_mins_stoich_mat(this%num_solids))
        ! end subroutine

        subroutine allocate_ind_gases_stoich_mat(this)
        !> This subroutine allocates the "ind_gases_stoich_mat" attribute
            class(reactive_zone_c) :: this
            if (allocated(this%ind_gases_stoich_mat)) then
                deallocate(this%ind_gases_stoich_mat)
            end if
            allocate(this%ind_gases_stoich_mat(this%gas_phase%num_gases_eq))
        end subroutine

        subroutine set_num_mins(this,num_mins)
        class(reactive_zone_c) :: this
        integer(kind=4), intent(in) :: num_mins !> number of minerals in reactive zone
        this%num_minerals=num_mins
        end subroutine
end module