!!> Reads reactive zones (species that define equilibrium reactions)
!!> Suponemos que el archivo ya ha sido abierto
!subroutine read_reactive_zones_Lagr(this,unit)
!    use chemistry_m, only: chemistry_c
!    use aqueous_chemistry_m, only: aqueous_chemistry_c, solid_chemistry_c, mineral_c, &
!    eq_reaction_c, chem_system_c, species_c, speciation_algebra_c
!    !use redox_kin_reaction_m, only: redox_kin_reaction_c
!    !use reaction_m, only: reaction_c
!    implicit none
!    class(chemistry_c) :: this
!    integer(kind=4), intent(in) :: unit
!    
!    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:),c2nc_init_guess(:)
!    integer(kind=4) :: i,j,k,l,num_rz,num_minerals,num_mins_rz,ind_mins_rz,niwtype, species_ind, ind_nf_sp,ind_rz,kin_react_type,&
!    num_nf_sp,n_k,chem_syst_type,num_non_flowing_species,num_react_zones,nischem,solid_chem_ind,ind_min
!    integer(kind=4), allocatable :: n_tar(:),min_eq_ind(:),n_tar_sol(:)
!    real(kind=8) :: aux,conc,temp
!    character(len=256) :: label,react_zone_file,file_eq_react,name,str,str1,str2,str3,str4,str5,path_name,Monod_name,&
!    file_kin_params,chem_syst_file
!    logical :: flag
!    type(aqueous_chemistry_c), allocatable :: initial_water_types(:)
!    type(solid_chemistry_c), allocatable :: init_solid_chems(:)
!    
!    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
!    solid_species_str(:),kin_react_names(:)
!    type(species_c) :: species,water,constrain
!    type(species_c), allocatable :: aq_species(:),cst_act_species(:),prim_species(:),non_flowing_species(:),species_array(:)
!    type(species_c), allocatable, target :: aq_species_tgt(:)
!    type(mineral_c) :: mineral
!    class(mineral_c), allocatable :: minerals(:)
!    ! class(kin_params_c), pointer :: p_kin_params=>null()
!    ! type(Monod_params_c), allocatable, target :: Monod_params_array(:)
!    ! type(Monod_params_c), target :: Monod_params
!    ! class(kin_reaction_c), pointer :: p_kin_react=>null()
!    ! type(kin_reaction_poly_c) :: kin_react_ptr
!    ! class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
!    ! type(redox_kin_c), target :: Monod
!    type(eq_reaction_c) :: eq_react
!    type(eq_reaction_c), allocatable :: eq_reacts(:)
!    
!    class(chem_system_c), pointer :: p_chem_syst=>null()
!    type(chem_system_c), target :: my_chem_syst
!   
!    type(speciation_algebra_c) :: spec_alg
!!> First iteration
!    do
!        read(unit,*) label
!        if (label=='end') then
!            rewind(unit)
!        else if (label=='REACTIVE ZONES') then
!            read(unit,*) num_rz !> number of reactive zones
!            call this%allocate_reactive_zones(num_rz)
!            i=1 !> counter reactive zones
!            do 
!                read(unit,*) ind_rz !> index reactive zone
!                if (ind_rz<1) then
!                    error stop "Index of reactive zone must be positive integer"
!                else if (ind_rz>num_rz) then
!                    error stop "Index of reactive zone cannot be larger than number of reactive zones"
!                end if
!                !num_rz=num_rz+1
!                num_nf_sp=0 !> counter non flowing species reactive zone
!                num_mins_rz=0 !> counter minerals reactive zone
!                do
!                    read(unit,*) name
!                    if (name=='*') then
!                        call this%reactive_zones(ind_rz)%allocate_non_flowing_species(num_nf_sp)
!                        call this%reactive_zones(ind_rz)%allocate_ind_mins(num_mins_rz)
!                        call this%reactive_zones(ind_rz)%set_chem_syst_react_zone(this%chem_syst)
!                        exit
!                    end if
!                    call mineral%set_phase_name(name)
!                    call this%chem_syst%is_mineral_in_chem_syst(mineral,flag,ind_min)
!                    if (flag.eqv..true.) then
!                        num_mins_rz=num_mins_rz+1
!                    end if
!                    num_nf_sp=num_nf_sp+1
!                end do
!                if (i<num_rz) then
!                    i=i+1
!                else
!                    exit
!                end if
!            end do
!        else
!            continue
!        end if
!    end do
!!> second iteration
!    do
!        read(unit,*) label
!        if (label=='end') then
!            exit
!        else if (label=='REACTIVE ZONES') then
!            read(unit,*) num_rz !> number of reactive zones
!            !call this%allocate_reactive_zones(num_rz)
!            i=1 !> counter reactive zones
!            do 
!                read(unit,*) ind_rz !> index reactive zone
!                ind_nf_sp=0 !> counter non flowing species
!                ind_mins_rz=0 !> counter minerals reactive zone
!                !if (rz_ind<1) then
!                !>    error stop "Index of reactive zone must be positive integer"
!                !else if (rz_ind>num_rz) then
!                !>    error stop "Index of reactive zone cannot be larger than number of reactive zones"
!                !end if
!                !num_rz=num_rz+1
!                !num_nf_sp=0
!                do
!                    read(unit,*) species%name
!                    if (species%name=='*') then
!                        exit
!                    end if
!                    call this%chem_syst%is_species_in_chem_syst(species,flag,species_ind)
!                    if (flag.eqv..true.) then
!                        ind_nf_sp=ind_nf_sp+1
!                        call this%reactive_zones(ind_rz)%set_single_non_flowing_species(ind_nf_sp,species_ind)
!                    end if
!                    call this%chem_syst%is_mineral_in_chem_syst(mineral,flag,ind_min)
!                    if (flag.eqv..true.) then
!                        ind_mins_rz=ind_mins_rz+1
!                        !> set mineral
!                        this%reactive_zones(ind_rz)%minerals(ind_mins_rz)=this%chem_syst%minerals(ind_min)
!                    end if
!                end do
!            !> Equilibrium reactions
!                call this%reactive_zones(ind_rz)%set_eq_reactions()
!                if (i<num_rz) then
!                    i=i+1
!                else
!                    exit
!                end if
!            end do
!        else
!            continue
!        end if
!    end do
!
!    
!    
!    
!    
!    !water_flag=0
!
!    !backspace(5)
!    !>    read(5,*) str
!    !>    print *, str
!    !>    if (str=='REACTIVE_ZONES') then
!    !>        line=line+1
!    !>        read(5,*) this%num_reactive_zones
!    !>        if (this%num_reactive_zones<1) then
!    !>            error stop "There are no reactive zones"
!    !>        else
!    !>            allocate(this%reactive_zones(this%num_reactive_zones))
!    !>            read(5,*) react_zone_file
!    !>            do l=1,this%num_reactive_zones
!    !>                line=line+1
!    !>                call this%reactive_zones(l)%set_chem_syst_react_zone(this%solid_chemistry%reactive_zone%chem_syst)
!    !>                call this%reactive_zones(l)%read_reactive_zone_Lagr(react_zone_file,line+1+l)
!    !>            end do
!    !>        end if
!    !>    else
!    !>        error stop "Wrong record in read_reactive_zones"
!    !>    end if
!end subroutine!!> Reads reactive zones (species that define equilibrium reactions)
!!> Suponemos que el archivo ya ha sido abierto
!subroutine read_reactive_zones_Lagr(this,unit)
!    use chemistry_m, only: chemistry_c
!    use aqueous_chemistry_m, only: aqueous_chemistry_c, solid_chemistry_c, mineral_c, &
!    eq_reaction_c, chem_system_c, species_c, speciation_algebra_c
!    !use redox_kin_reaction_m, only: redox_kin_reaction_c
!    !use reaction_m, only: reaction_c
!    implicit none
!    class(chemistry_c) :: this
!    integer(kind=4), intent(in) :: unit
!    
!    real(kind=8), allocatable :: Sk(:,:),logK(:),gamma_1(:),gamma_2(:),c2nc_init_guess(:)
!    integer(kind=4) :: i,j,k,l,num_rz,num_minerals,num_mins_rz,ind_mins_rz,niwtype, species_ind, ind_nf_sp,ind_rz,kin_react_type,&
!    num_nf_sp,n_k,chem_syst_type,num_non_flowing_species,num_react_zones,nischem,solid_chem_ind,ind_min
!    integer(kind=4), allocatable :: n_tar(:),min_eq_ind(:),n_tar_sol(:)
!    real(kind=8) :: aux,conc,temp
!    character(len=256) :: label,react_zone_file,file_eq_react,name,str,str1,str2,str3,str4,str5,path_name,Monod_name,&
!    file_kin_params,chem_syst_file
!    logical :: flag
!    type(aqueous_chemistry_c), allocatable :: initial_water_types(:)
!    type(solid_chemistry_c), allocatable :: init_solid_chems(:)
!    
!    character(len=256), allocatable :: aq_species_str(:),prim_species_str(:),cst_act_species_str(:),minerals_str(:),&
!    solid_species_str(:),kin_react_names(:)
!    type(species_c) :: species,water,constrain
!    type(species_c), allocatable :: aq_species(:),cst_act_species(:),prim_species(:),non_flowing_species(:),species_array(:)
!    type(species_c), allocatable, target :: aq_species_tgt(:)
!    type(mineral_c) :: mineral
!    class(mineral_c), allocatable :: minerals(:)
!    ! class(kin_params_c), pointer :: p_kin_params=>null()
!    ! type(Monod_params_c), allocatable, target :: Monod_params_array(:)
!    ! type(Monod_params_c), target :: Monod_params
!    ! class(kin_reaction_c), pointer :: p_kin_react=>null()
!    ! type(kin_reaction_poly_c) :: kin_react_ptr
!    ! class(kin_reaction_poly_c), allocatable :: kin_reacts(:)
!    ! type(redox_kin_c), target :: Monod
!    type(eq_reaction_c) :: eq_react
!    type(eq_reaction_c), allocatable :: eq_reacts(:)
!    
!    class(chem_system_c), pointer :: p_chem_syst=>null()
!    type(chem_system_c), target :: my_chem_syst
!   
!    type(speciation_algebra_c) :: spec_alg
!!> First iteration
!    do
!        read(unit,*) label
!        if (label=='end') then
!            rewind(unit)
!        else if (label=='REACTIVE ZONES') then
!            read(unit,*) num_rz !> number of reactive zones
!            call this%allocate_reactive_zones(num_rz)
!            i=1 !> counter reactive zones
!            do 
!                read(unit,*) ind_rz !> index reactive zone
!                if (ind_rz<1) then
!                    error stop "Index of reactive zone must be positive integer"
!                else if (ind_rz>num_rz) then
!                    error stop "Index of reactive zone cannot be larger than number of reactive zones"
!                end if
!                !num_rz=num_rz+1
!                num_nf_sp=0 !> counter non flowing species reactive zone
!                num_mins_rz=0 !> counter minerals reactive zone
!                do
!                    read(unit,*) name
!                    if (name=='*') then
!                        call this%reactive_zones(ind_rz)%allocate_non_flowing_species(num_nf_sp)
!                        call this%reactive_zones(ind_rz)%allocate_ind_mins(num_mins_rz)
!                        call this%reactive_zones(ind_rz)%set_chem_syst_react_zone(this%chem_syst)
!                        exit
!                    end if
!                    call mineral%set_phase_name(name)
!                    call this%chem_syst%is_mineral_in_chem_syst(mineral,flag,ind_min)
!                    if (flag.eqv..true.) then
!                        num_mins_rz=num_mins_rz+1
!                    end if
!                    num_nf_sp=num_nf_sp+1
!                end do
!                if (i<num_rz) then
!                    i=i+1
!                else
!                    exit
!                end if
!            end do
!        else
!            continue
!        end if
!    end do
!!> second iteration
!    do
!        read(unit,*) label
!        if (label=='end') then
!            exit
!        else if (label=='REACTIVE ZONES') then
!            read(unit,*) num_rz !> number of reactive zones
!            !call this%allocate_reactive_zones(num_rz)
!            i=1 !> counter reactive zones
!            do 
!                read(unit,*) ind_rz !> index reactive zone
!                ind_nf_sp=0 !> counter non flowing species
!                ind_mins_rz=0 !> counter minerals reactive zone
!                !if (rz_ind<1) then
!                !>    error stop "Index of reactive zone must be positive integer"
!                !else if (rz_ind>num_rz) then
!                !>    error stop "Index of reactive zone cannot be larger than number of reactive zones"
!                !end if
!                !num_rz=num_rz+1
!                !num_nf_sp=0
!                do
!                    read(unit,*) species%name
!                    if (species%name=='*') then
!                        exit
!                    end if
!                    call this%chem_syst%is_species_in_chem_syst(species,flag,species_ind)
!                    if (flag.eqv..true.) then
!                        ind_nf_sp=ind_nf_sp+1
!                        call this%reactive_zones(ind_rz)%set_single_non_flowing_species(ind_nf_sp,species_ind)
!                    end if
!                    call this%chem_syst%is_mineral_in_chem_syst(mineral,flag,ind_min)
!                    if (flag.eqv..true.) then
!                        ind_mins_rz=ind_mins_rz+1
!                        !> set mineral
!                        this%reactive_zones(ind_rz)%minerals(ind_mins_rz)=this%chem_syst%minerals(ind_min)
!                    end if
!                end do
!            !> Equilibrium reactions
!                call this%reactive_zones(ind_rz)%set_eq_reactions()
!                if (i<num_rz) then
!                    i=i+1
!                else
!                    exit
!                end if
!            end do
!        else
!            continue
!        end if
!    end do
!
!    
!    
!    
!    
!    !water_flag=0
!
!    !backspace(5)
!    !>    read(5,*) str
!    !>    print *, str
!    !>    if (str=='REACTIVE_ZONES') then
!    !>        line=line+1
!    !>        read(5,*) this%num_reactive_zones
!    !>        if (this%num_reactive_zones<1) then
!    !>            error stop "There are no reactive zones"
!    !>        else
!    !>            allocate(this%reactive_zones(this%num_reactive_zones))
!    !>            read(5,*) react_zone_file
!    !>            do l=1,this%num_reactive_zones
!    !>                line=line+1
!    !>                call this%reactive_zones(l)%set_chem_syst_react_zone(this%solid_chemistry%reactive_zone%chem_syst)
!    !>                call this%reactive_zones(l)%read_reactive_zone_Lagr(react_zone_file,line+1+l)
!    !>            end do
!    !>        end if
!    !>    else
!    !>        error stop "Wrong record in read_reactive_zones"
!    !>    end if
!end subroutine