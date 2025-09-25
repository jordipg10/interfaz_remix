!> Writes data and results of chemistry object
subroutine write_chemistry(this,unit,Lagr_flag)
    use chemistry_m, only: chemistry_c, inf_norm_vec_real
    
    implicit none
    class(chemistry_c), intent(in) :: this !> chemistry object
    integer(kind=4), intent(in) :: unit !> file unit
    logical, intent(in) :: Lagr_flag !> flag for Lagrangian particles
    
    integer(kind=4) :: i,j,l,num_dom,num_ext
    integer(kind=4), allocatable :: dom_indices(:),ext_indices(:)
    
    write(unit,"(/,2x,'Chemical input option:',I25)") this%read_opt
    if (this%lump_flag .eqv. .true. ) then
        write(unit,"(/,2x,'WMA option:',10x,'Lumped')") 
    else
        write(unit,"(/,2x,'WMA option:',10x,'Consistent')")
        write(unit,"(/,2x,'Consistent WMA option:',I25)") this%cons_opt
        write(unit,"(/,2x,'rk downstream waters option:',I25)") this%rk_down_opt
        if (this%cons_opt>1 .and. this%rk_down_opt>1) then
            write(unit,"(/,2x,'Parameter mu:',ES15.5)") this%CV_params%mu
        end if
    end if
    write(unit,"(/,2x,'rk average option:',I25)") this%rk_avg_opt
    write(unit,"(/,2x,'Number of lumpings:',I25)") this%num_lump
    write(unit,"(/,2x,'Aqueous species:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,A15)") this%chem_syst%aq_phase%aq_species(i)%name
    end do
    write(unit,"(2x,'Minerals (name + molar volume):'/)")
    do i=1,this%chem_syst%num_minerals
        write(unit,"(10x,A15,ES15.5)") this%chem_syst%minerals(i)%name, this%chem_syst%minerals(i)%mineral%mol_vol
    end do
    write(unit,"(2x,'Gases:'/)")
    do i=1,this%chem_syst%gas_phase%num_species
        write(unit,"(10x,A15)") this%chem_syst%gas_phase%gases(i)%name
    end do
    write(unit,"(2x,'Surface complexes:'/)")
    do i=1,this%chem_syst%cat_exch%num_surf_compl
        write(unit,"(10x,A15)") this%chem_syst%cat_exch%surf_compl(i)%name
    end do
    write(unit,"(2x,'Equilibrium reactions:'/)")
    do i=1,this%chem_syst%num_eq_reacts
        write(unit,"(10x,A30/)") this%chem_syst%eq_reacts(i)%name
    end do
    write(unit,"(2x,'Kinetic reactions:'/)")
    do i=1,this%chem_syst%num_kin_reacts
        call this%chem_syst%kin_reacts(i)%kin_reaction%write_reaction(unit)
    end do
    write(unit,"(2x,'Global stoichiometric matrix:',/)")
    do i=1,this%chem_syst%num_reacts
        write(unit,"(10x,*(F15.5))") (this%chem_syst%stoich_mat(i,j), j=1,this%chem_syst%num_species)
    end do
    write(unit,"(2x,'Equilibrium constants:',/)")
    do i=1,this%chem_syst%num_eq_reacts
        write(unit,"(10x,*(ES15.5))") this%chem_syst%eq_reacts(i)%eq_cst
    end do
    write(unit,"(/,2x,'Global component matrix:'/)")
    ! do i=1,this%chem_syst%speciation_alg%num_prim_species
    !     write(unit,"(10x,*(F15.5))") (this%chem_syst%speciation_alg%comp_mat(i,j), &
    !     j=1,this%chem_syst%speciation_alg%num_var_act_species)
    ! end do
    if (Lagr_flag) then
        ! write(unit,"(2x,'Initial water particles coordinates (by rows):'/)")
        ! do i=1,this%num_target_waters
        !     write(unit,"(2x,*(ES15.5))") this%target_waters_init(i)%solid_chemistry%tar%coord !> Coordinates of the targets
        ! end do
        ! write(unit,"(/,2x,'Final water particles coordinates:'/)")
        ! do i=1,this%num_target_waters
        !     write(unit,"(2x,*(ES15.5))") this%target_waters(i)%solid_chemistry%tar%coord !> Coordinates of the targets
        ! end do
        write(unit,"(2x,'Initial water particles targets (by rows):'/)")
        do i=1,this%num_target_waters_init
            write(unit,"(2x,I5)") this%target_waters_init(i)%solid_chemistry%tar%id !> ID of the targets
        end do
        write(unit,"(/,2x,'Final water particles targets:'/)")
        do i=1,this%num_target_waters
            write(unit,"(2x,I5)") this%target_waters(i)%solid_chemistry%tar%id !> ID of the targets
        end do
        write(unit,"(/,2x,'Initial water particles volumes (by rows):'/)")
        do i=1,this%num_target_waters_init
            write(unit,"(2x,ES15.5)") this%target_waters_init(i)%volume
        end do
        write(unit,"(/,2x,'Final water particles volumes:'/)")
        do i=1,this%num_target_waters
            write(unit,"(2x,ES15.5)") this%target_waters(i)%volume
        end do
    end if
    write(unit,"(/,2x,'Concentration of recharge waters (in molalities): &
        (rows -> aqueous species in the same order as before, columns -> targets)'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%ext_waters_indices(j))%concentrations(&
            this%target_waters(this%ext_waters_indices(j))%indices_aq_species(i)), j=1,this%num_ext_waters)
    end do
    write(unit,"(/,2x,'Activities of recharge waters:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%ext_waters_indices(j))%activities(&
            this%target_waters(this%ext_waters_indices(j))%indices_aq_species(i)), j=1,this%num_ext_waters)
    end do
    write(unit,"(/,2x,'Concentration of boundary waters (in molalities): &
        (rows -> aqueous species in the same order as before, columns -> targets)'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%bd_waters_indices(j))%concentrations(&
            this%target_waters(this%bd_waters_indices(j))%indices_aq_species(i)), j=1,this%num_bd_waters)
    end do
    write(unit,"(/,2x,'Activities of boundary waters:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%bd_waters_indices(j))%activities(&
            this%target_waters(this%bd_waters_indices(j))%indices_aq_species(i)), j=1,this%num_bd_waters)
    end do
    write(unit,"(/,2x,'Initial concentration of aqueous species in the domain:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%concentrations(&
        this%target_waters_init(this%dom_tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Initial activities aqueous species domain:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%activities(&
        this%target_waters_init(this%dom_tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Final concentration of aqueous species in the domain:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%concentrations(&
            this%target_waters(this%dom_tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Final activities aqueous species domain:'/)")
    do i=1,this%chem_syst%aq_phase%num_species
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%activities(&
            this%target_waters(this%dom_tar_wat_indices(j))%indices_aq_species(i)), j=1,this%num_target_waters_dom)
    end do
    !if (this%chem_syst%aq_phase%ind_wat>0) then
        ! write(unit,"(/,2x,'Initial activity water domain:'/)")
        ! write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%activities(&
        !     this%target_waters_init(this%dom_tar_wat_indices(j))%indices_aq_species(this%chem_syst%aq_phase%ind_wat)),&
        !     j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Initial salinity domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%salinity,j=1,this%num_target_waters_dom)
        ! write(unit,"(/,2x,'Activity water domain:'/)")
        ! write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%activities(&
        !     this%target_waters(this%dom_tar_wat_indices(j))%indices_aq_species(this%chem_syst%aq_phase%ind_wat)),&
        !     j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Final salinity domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%salinity, &
            j=1,this%num_target_waters_dom)
    !end if
    if (this%chem_syst%aq_phase%ind_prot>0) then
        write(unit,"(/,2x,'Initial pH domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%pH, j=1,this%num_target_waters_dom)
        write(unit,"(/,2x,'Final pH domain:'/)")
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%pH, j=1,this%num_target_waters_dom)
    end if
    write(unit,"(/,2x,'Aqueous mean equilibrium reaction amounts:'/)")
    do i=1,this%chem_syst%aq_phase%num_aq_complexes
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%Re_mean(i), j=1,this%num_target_waters_dom)
    end do
    ! write(unit,"(/,2x,'Aqueous equilibrium reaction rates:'/)")
    ! do i=1,this%chem_syst%aq_phase%num_aq_complexes
    !     write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%r_eq(i), j=1,this%num_target_waters_dom)
    ! end do
    write(unit,"(/,2x,'Initial aqueous kinetic reaction rates:'/)")
    do i=1,this%chem_syst%num_aq_kin_reacts
        write(unit,"(10x,*(ES15.5))") (this%target_waters_init(this%dom_tar_wat_indices(j))%rk(i), j=1,this%num_target_waters_dom)
    end do
    
    write(unit,"(/,2x,'Aqueous kinetic reaction rates:'/)")
    do i=1,this%chem_syst%num_aq_kin_reacts
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%rk(i), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Aqueous mean kinetic reaction amounts:'/)")
    do i=1,this%chem_syst%num_aq_kin_reacts
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%Rk_mean(i), j=1,this%num_target_waters_dom)
    end do
    write(unit,"(/,2x,'Aqueous estimated kinetic reaction amounts:'/)")
    do i=1,this%chem_syst%num_aq_kin_reacts
        write(unit,"(10x,*(ES15.5))") (this%target_waters(this%dom_tar_wat_indices(j))%Rk_est(i), j=1,this%num_target_waters_dom)
    end do
    do l=1,this%num_mineral_zones
        write(unit,"(/,2x,'Mineral zone',I5,':'/)") l
        write(unit,"(10x,'Minerals in equilibrium:'/)")
        do i=1,this%mineral_zones(l)%num_minerals_eq
            write(unit,"(20x,A15/)") this%mineral_zones(l)%chem_syst%minerals(&
                this%mineral_zones(l)%ind_min_chem_syst(this%mineral_zones(l)%num_minerals_kin+i))%name
        end do
        call this%link_target_waters_mineral_zone(l,dom_indices,ext_indices)
        !write(unit,"(/,10x,'Mineral equilibrium reaction rates:'/)")
        !do i=1,this%mineral_zones(l)%num_minerals_eq
        !    write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%r_eq(i), j=1,size(dom_indices))
        !end do
        write(unit,"(10x,'Kinetic minerals:'/)")
        do i=1,this%mineral_zones(l)%num_minerals_kin
            write(unit,"(20x,A15/)")this%mineral_zones(l)%chem_syst%minerals(this%mineral_zones(l)%ind_min_chem_syst(i))%name
        end do
        write(unit,"(/,10x,'Mineral mean kinetic reaction amounts:'/)")
        do i=1,this%mineral_zones(l)%num_minerals_kin
            !print *, this%target_waters(dom_indices(1))%solid_chemistry%Rk_mean
            write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%Rk_mean(i), j=1,SIZE(dom_indices))
        end do
        write(unit,"(/,10x,'Mineral estimated kinetic reaction amounts:'/)")
        do i=1,this%mineral_zones(l)%num_minerals_kin
            write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%Rk_est(i), j=1,SIZE(dom_indices))
        end do
        write(unit,"(/,10x,'Absolute error kinetic reaction amounts:'/)")
        do i=1,this%mineral_zones(l)%num_minerals_kin
            write(unit,"(10x,*(ES15.5))") (abs(this%target_waters(dom_indices(j))%solid_chemistry%Rk_est(i)-&
            this%target_waters(dom_indices(j))%solid_chemistry%Rk_mean(i)), j=1,SIZE(dom_indices))
        end do
    end do
    do l=1,this%num_reactive_zones
        call this%link_target_waters_reactive_zone(l,dom_indices,ext_indices)
        num_dom=size(dom_indices)
        num_ext=size(ext_indices)
        if (num_dom>0) then
            write(unit,"(/,2x,'Reactive zone',I5,':'/)") l
            write(unit,"(10x,'Non-flowing species:'/)")
            do i=1,this%reactive_zones(l)%num_non_flowing_species
                write(unit,"(20x,A15)") this%reactive_zones(l)%non_flowing_species(i)%name
            end do
            write(unit,"(10x,'Number of equilibrium reactions:',I10/)") this%reactive_zones(l)%speciation_alg%num_eq_reactions
            write(unit,"(10x,'Se:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_eq_reactions
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%stoich_mat(i,j), &
                j=1,this%reactive_zones(l)%speciation_alg%num_species)
            end do
            write(unit,"(10x,'Sk:',/)")
            do i=1,this%reactive_zones(l)%chem_syst%num_kin_reacts
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%chem_syst%Sk(i,j), &
                j=1,this%reactive_zones(l)%chem_syst%num_species)
            end do
            write(unit,"(10x,'U*Sk^T:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_aq_prim_species
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%U_SkT_prod(i,j), j=1,size(this%reactive_zones(l)%U_SkT_prod,2)) !> chapuza
            end do
            write(unit,"(/,10x,'Equilibrium constants:',/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_eq_reactions
                write(unit,"(10x,*(ES15.5))") this%reactive_zones(l)%chem_syst%eq_reacts(&
                    this%reactive_zones(l)%ind_eq_reacts(i))%eq_cst
            end do
            !write(unit,"(/,10x,'Aqueous species:',/)")
            !do i=1,this%chem_syst%aq_phase%num_species
            !    write(unit,"(20x,A15)") this%reactive_zones(l)%chem_syst%aq_phase%aq_species(this%target_waters(dom_indices(1))%indices_aq_species(i))%name
            !end do
            !write(unit,"(/,10x,'Primary species:',/)")
            !do i=1,this%reactive_zones(l)%speciation_alg%num_aq_prim_species
            !    ! write(unit,"(20x,A15)") this%reactive_zones(l)%chem_syst%species(&
            !    !     this%reactive_zones(l)%ind_aq_var_act_species(i))%name
            !end do
            if (this%reactive_zones(l)%cat_exch_zone%num_surf_compl>0) then
                write(unit,"(20x,A15)") this%reactive_zones(l)%cat_exch_zone%surf_compl(1)%name
            end if
            write(unit,"(/,10x,'Component matrix:'/)")
            !print *, l
            !print *, this%reactive_zones(l)%speciation_alg%num_prim_species
            !print *, this%reactive_zones(l)%speciation_alg%num_var_act_species
            !print *, this%reactive_zones(l)%speciation_alg%comp_mat
            do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species
                write(unit,"(10x,*(F15.5))") (this%reactive_zones(l)%speciation_alg%comp_mat(i,j), &
                    j=1,this%reactive_zones(l)%speciation_alg%num_var_act_species)
            end do           
            write(unit,"(/,10x,'Initial concentration of components:'/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species
                write(unit,"(10x,*(ES15.5))") (dot_product(this%reactive_zones(l)%speciation_alg%comp_mat(i,:), &
                    this%target_waters_init(dom_indices(j))%get_conc_nc()), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Final concentration of components:'/)")
            do i=1,this%reactive_zones(l)%speciation_alg%num_prim_species
                write(unit,"(10x,*(ES15.5))") (dot_product(this%reactive_zones(l)%speciation_alg%comp_mat(i,:), &
                    this%target_waters(dom_indices(j))%get_conc_nc()), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Final concentration of aqueous species:'/)")
            do i=1,this%chem_syst%aq_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%concentrations(this%target_waters(&
                    dom_indices(j))%indices_aq_species(i)), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Final activities aqueous species:'/)")
            do i=1,this%chem_syst%aq_phase%num_species
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%activities(this%target_waters(&
                    dom_indices(j))%indices_aq_species(i)), j=1,size(dom_indices))
            end do
             write(unit,"(/,10x,'Initial volumetric fractions of minerals:'/)")
             do i=1,this%reactive_zones(l)%num_minerals
                  write(unit,"(10x,*(ES15.5))") (this%target_solids_init(dom_indices(j))%vol_fracts(this%target_solids_init(&
                      dom_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(dom_indices))
             end do
             write(unit,"(/,10x,'Volumetric fractions of minerals:'/)")
             do i=1,this%reactive_zones(l)%num_minerals
                 write(unit,"(10x,*(ES15.5))") (this%target_solids(dom_indices(j))%vol_fracts(this%target_solids(&
                     dom_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(dom_indices))
             end do
             write(unit,"(/,10x,'Initial concentration of minerals:'/)")
             do i=1,this%reactive_zones(l)%num_minerals
                 write(unit,"(10x,*(ES15.5))") (this%target_solids_init(dom_indices(j))%concentrations(this%target_solids_init(&
                     dom_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(dom_indices))
             end do
             write(unit,"(/,10x,'Concentration of minerals:'/)")
             do i=1,this%reactive_zones(l)%num_minerals
                 write(unit,"(10x,*(ES15.5))") (this%target_solids(dom_indices(j))%concentrations(this%target_solids(&
                     dom_indices(j))%mineral_zone%num_minerals_kin+i), j=1,size(dom_indices))
             end do
            write(unit,"(/,10x,'Initial concentration of surface complexes:'/)")
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl
                write(unit,"(10x,*(ES15.5))") (this%target_solids_init(dom_indices(j))%concentrations(&
                    this%target_solids_init(dom_indices(j))%mineral_zone%num_minerals+i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Concentration of surface complexes:'/)")
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl
                write(unit,"(10x,*(ES15.5))") (this%target_solids(dom_indices(j))%concentrations(&
                    this%target_solids(dom_indices(j))%mineral_zone%num_minerals+i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial activity of surface complexes:'/)")
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl
                write(unit,"(10x,*(ES15.5))") (this%target_solids_init(dom_indices(j))%activities(&
                this%target_solids_init(dom_indices(j))%mineral_zone%num_minerals+i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Activity of surface complexes:'/)")
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_surf_compl
                write(unit,"(10x,*(ES15.5))") (this%target_solids(dom_indices(j))%activities(&
                    this%target_solids(dom_indices(j))%mineral_zone%num_minerals+i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial concentration of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(dom_indices(j))%concentrations(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Concentration of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%concentrations(i), &
                j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial partial pressures of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(dom_indices(j))%activities(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Partial pressures of gases:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%activities(i),j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Initial volume of gas:'/)")
            if (this%reactive_zones(l)%gas_phase%num_species>0) then
                write(unit,"(10x,*(ES15.5))") (this%target_gases_init(dom_indices(j))%volume, j=1,size(dom_indices))
            end if
            write(unit,"(/,10x,'Volume of gas:'/)")
            if (this%reactive_zones(l)%gas_phase%num_species>0) then
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%volume, j=1,size(dom_indices))
            end if
            write(unit,"(/,10x,'Aqueous equilibrium reaction amounts:'/)")
            do i=1,this%reactive_zones(l)%chem_syst%num_aq_eq_reacts
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%Re_mean(i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Mineral mean equilibrium reaction amounts:'/)")
            do i=1,this%reactive_zones(l)%num_minerals
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%Re_mean(i), j=1,size(dom_indices))
            end do
            !> chapuza
            !write(unit,"(/,10x,'Mineral kinetic reaction rates:'/)")
            !do i=1,this%target_waters(dom_indices(1))%solid_chemistry%num_min_kin
            !    write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%rk(i), j=1,size(dom_indices))
            !end do
            write(unit,"(/,10x,'Exchange equilibrium reaction rates:'/)")
            do i=1,this%reactive_zones(l)%cat_exch_zone%num_exch_cats
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%solid_chemistry%r_eq(&
                this%reactive_zones(l)%num_minerals+i), j=1,size(dom_indices))
            end do
            write(unit,"(/,10x,'Gas equilibrium reaction rates:'/)")
            do i=1,this%reactive_zones(l)%gas_phase%num_gases_eq
                write(unit,"(10x,*(ES15.5))") (this%target_waters(dom_indices(j))%gas_chemistry%r_eq(i), j=1,size(dom_indices))
            end do
            deallocate(dom_indices,ext_indices)
        end if
    end do
end subroutine write_chemistry