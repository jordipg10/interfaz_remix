!subroutine compute_conc_biomass(this,time,F_mat,conc_c) !> F*dc_c/dt=S_c^T*r_k
!>    use chemistry_m
!>    implicit none
!>    class(chemistry_c) :: this
!>    !real(kind=8), intent(in) :: conc_init(:,:)
!>    real(kind=8), intent(in) :: time
!>    class(diag_matrix_c), intent(in) :: F_mat
!>    real(kind=8), intent(out), optional :: conc_c(:,:) !> tiene que estar alocatado
!>            
!>    integer(kind=4) :: i,j,water_flag
!>    real(kind=8), allocatable :: rk_mat(:,:),req_mat(:,:)
!>    type(species_c) :: water
!>    logical :: flag
!>    
!>    !do j=1,size(this%solid_chem)
!>    !>    call this%solid_chem(j)%allocate_concentrations(this%solid_chemistry%reactive_zone%chem_syst)
!>    !end do
!>    
!>    !req_mat=this%get_req_mat()
!>    !select type (this)
!>    !type is (chemistry_transient_c)
!>        !select type (transport=>this%transport)
!>        !class is (diffusion_transient_c)
!>        water%name='h2o'
!>        call water%belongs_to(this%solid_chemistry%reactive_zone%chem_syst%cst_act_species,flag)
!>        if (flag==.true.) then
!>            water_flag=1
!>        else
!>            water_flag=0
!>        end if
!>        select type (chem_syst=>this%solid_chemistry%reactive_zone%chem_syst)
!>        type is (chem_system_eq_kin_c)
!>            !allocate(rk_mat(chem_syst%num_reactions,transport%spatial_discr%Num_targets))
!>            !rk_mat=1d0
!>            !rk_mat(:,1)=1d0
!>            if (water_flag==1) then
!>                do j=1,size(this%aqueous_chem)
!>                    !allocate(this%solid_zone%minerals_eq(i)%conc_solid(transport%spatial_discr%Num_targets))
!>                    !allocate(this%solid_zone%solid_species(i)%conc_solid(transport%spatial_discr%Num_targets))
!>                    !print *, this%solid_zone%solid_species(i)%conc_solid
!>                    this%aqueous_chem(j)%concentrations(chem_syst%num_species)=this%aqueous_chem_init(j)%concentrations(chem_syst%num_species)+time*dot_product(chem_syst%Sc(:,1),[this%r_eq_mat(:,j),this%rk_mat(:,j)])/F_mat%diag(j)
!>                end do
!>            end if
!>            do i=1,chem_syst%num_solid_species
!>                do j=1,size(this%solid_chem)
!>                    !this%aqueous_chem(j)%concentrations(chem_syst%num_species)=this%aqueous_chem_init(j)%concentrations(chem_syst%num_species)+time*dot_product(chem_syst%Sc(:,1))
!>                    !this%solid_zone%minerals_eq(i)%conc_solid(j)=this%solid_zone_init%minerals_eq(i)%conc_solid(j)+transport%time_discr%Final_time*dot_product(chem_syst%Sc(:,i),rk_mat(:,j))/transport%F_mat%diag(j)
!>                    !this%solid_zone%solid_species(i)%conc_solid(j)=this%solid_zone_init%solid_species(i)%conc_solid(j)+transport%time_discr%Final_time*dot_product(chem_syst%Sc(:,i),this%r_eq(:,j)%reaction_rate)/transport%F_mat%diag(j)
!>                    !this%solid_chem(j)%concentrations(i)=this%solid_chem_init(j)%concentrations(i)+time*dot_product(chem_syst%Sc(:,i+water_flag),this%r_eq(:,j)%reaction_rate)/F_mat%diag(j)
!>                    this%solid_chem(j)%concentrations(i)=this%solid_chem_init(j)%concentrations(i)+time*dot_product(chem_syst%Sc(:,i+water_flag),[this%r_eq_mat(:,j),this%rk_mat(:,j)])/F_mat%diag(j)
!>                end do
!>            end do
!>        end select
!>        if (present(conc_c)) then
!>            !allocate(conc_c(this%solid_chemistry%reactive_zone%chem_syst%num_cst_act_species,size(this%aqueous_chem)))
!>            if (water_flag==1) then
!>                do j=1,size(this%aqueous_chem)
!>                    conc_c(1,j)=this%aqueous_chem(j)%concentrations(this%solid_chemistry%reactive_zone%chem_syst%num_species)
!>                end do
!>            end if
!>            do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_solid_species
!>                do j=1,size(this%aqueous_chem)
!>                    conc_c(i+water_flag,j)=this%solid_chem(j)%concentrations(i)
!>                end do
!>            end do
!>        end if
!>    !end select
!end subroutine