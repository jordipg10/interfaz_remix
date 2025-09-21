!> Solid chemistry subclass:This class contains the concentrations, activities, volumetric fractions and specific surfaces of solid species
!! We assume that the pointers to reactive zone and mineral zone classes are always set
module solid_chemistry_m
    use local_chemistry_m, only: local_chemistry_c
    use reactive_zone_Lagr_m, only: reactive_zone_c, inf_norm_vec_real, species_c, LU_lin_syst
    use mineral_zone_m, only: mineral_zone_c, mineral_c
    use target_m, only: target_c
    implicit none
    save
!***********************************************************************************************************************************
    type, public, extends(local_chemistry_c) :: solid_chemistry_c
        integer(kind=4) :: num_solids=0 !> number of solid species (chapucilla)
        real(kind=8), allocatable :: equivalents(:) !> equivalents of adsorbed cations
        real(kind=8), allocatable :: vol_fracts(:) !> volumetric fractions minerals (dimensionless) (first kinetic, then equilibrium)
        real(kind=8), allocatable :: react_surfaces(:) !> specific surface minerals per unit volume
        class(reactive_zone_c), pointer :: reactive_zone !> solids in equilibrium belong to a certain reactive zone
        class(mineral_zone_c), pointer :: mineral_zone !> minerals belong to a certain mineral zone
        real(kind=8) :: CEC !> Cation exchange capacity [eq/L] (we assume we have one or zero adsorption surfaces)
        class(target_c), pointer :: tar !> target associated to solid chemistry (we assume bijection)
    contains
    !> Set
        procedure, public :: set_concentrations=>set_conc_solids
        procedure, public :: set_indices_solids
        procedure, public :: set_vol_fracts
        procedure, public :: set_react_surfaces
        procedure, public :: set_reactive_zone
        procedure, public :: set_mineral_zone
        procedure, public :: set_CEC
        procedure, public :: set_conc_free_site
        procedure, public :: set_act_surf_compl
        procedure, public :: set_target
    !> Allocate
        procedure, public :: allocate_vol_fracts
        procedure, public :: allocate_react_surfaces
        procedure, public :: allocate_conc_solids
        procedure, public :: allocate_reaction_rates_solid_chem
        procedure, public :: allocate_activities
        procedure, public :: allocate_log_act_coeffs_solid_chem
        procedure, public :: allocate_equivalents
    !> Compute
        procedure, public :: compute_activities_solids
        procedure, public :: compute_conc_ads_cats
        procedure, public :: compute_equivalents
        procedure, public :: compute_conc_minerals_iter_EE
        procedure, public :: compute_mass_bal_mins
        procedure, public :: compute_conc_surf_ideal_Newton
        procedure, public :: compute_conc_surf_ideal_Picard
        procedure, public :: compute_conc_surf_ideal
        procedure, public :: compute_num_solids_solid_chem
        !procedure, public :: compute_conc_solids
    !> Update
        procedure, public :: update_conc_ads_cats
        procedure, public :: update_act_ads_cats
    !> Check
        procedure, public :: check_solid_chemistry
    !> ASsign
        procedure, public :: assign_solid_chemistry
        procedure, public :: modify_mix_ratios_Rk
    end type
    
    interface
        subroutine compute_conc_surf_ideal_Newton(this,conc_cats,act_ads_cats_ig,niter,CV_flag)
            import solid_chemistry_c
            implicit none
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_cats(:) !> chapuza (dim=n� cation exchange half reactions)
            real(kind=8), intent(in) :: act_ads_cats_ig(:) !> surface complex concentrations (must be already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_conc_surf_ideal_Picard(this,conc_cats,act_ads_cats_ig,niter,CV_flag)
            import solid_chemistry_c
            implicit none
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_cats(:) !> chapuza (dim=n� cation exchange half reactions)
            real(kind=8), intent(in) :: act_ads_cats_ig(:) !> surface complex concentrations (must be already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_conc_surf_ideal(this,conc_cats)
            import solid_chemistry_c
            implicit none
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_cats(:) !> chapuza (dim=n� cation exchange half reactions)
            !real(kind=8), intent(in) :: act_ads_cats_ig(:) !> surface complex concentrations (must be already allocated)
            !integer(kind=4), intent(out) :: niter !> number of iterations
            !logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        !subroutine read_solid_chem_init(this,filename,reactive_zones,line,num_tar)
        !>    import reactive_zone_c
        !>    import solid_chemistry_c
        !>    implicit none
        !>    class(solid_chemistry_c) :: this
        !>    character(len=*), intent(in) :: filename
        !>    class(reactive_zone_c), intent(in) :: reactive_zones(:)
        !>    integer(kind=4), intent(inout) :: line
        !>    integer(kind=4), intent(out) :: num_tar
        !end subroutine
        
        ! subroutine compute_conc_minerals_iter_EE(this,Delta_t)
        !     import solid_chemistry_c
        !     implicit none
        !     class(solid_chemistry_c) :: this
        !     !real(kind=8), intent(in) :: r_eq(:) !> reaction rates
        !     real(kind=8), intent(in) :: Delta_t !> time step
        ! end subroutine
        
        subroutine compute_mass_bal_mins(this,Delta_t)
            import solid_chemistry_c
            implicit none
            class(solid_chemistry_c) :: this
            !real(kind=8), intent(in) :: r_eq(:) !> equilibrium reaction rates
            real(kind=8), intent(in) :: Delta_t !> time step
        end subroutine
        
        ! subroutine compute_conc_solids(this,r_vec,porosity,time)
        !     import solid_chemistry_c
        !     implicit none
        !     class(solid_chemistry_c) :: this
        !     real(kind=8), intent(in) :: r_vec(:) !> reaction rates
        !     real(kind=8), intent(in) :: porosity
        !     real(kind=8), intent(in) :: time
        ! end subroutine
    end interface
    
    
    
    contains
        
       
        
        subroutine set_conc_solids(this,conc)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
            if (this%reactive_zone%chem_syst%num_solids<size(conc)) error stop "Dimension error in set_conc_solids"
            this%concentrations=conc
        end subroutine        
        
        subroutine set_conc_free_site(this)
            class(solid_chemistry_c) :: this
            !real(kind=8), intent(in) :: conc
            if (ASSOCIATED(this%mineral_zone)) then
                this%concentrations(this%mineral_zone%num_minerals+1)=1d-32 !> by default
            else if (ASSOCIATED(this%reactive_zone)) then
                this%concentrations(1)=1d-32 !> by default
            else
                error stop "Reactive zone not associated to solid chemistry"
            end if
        end subroutine
        !subroutine set_num_min_kin(this,num_min_kin)
        !    implicit none
        !    class(solid_chemistry_c) :: this
        !    integer(kind=4), intent(in) :: num_min_kin
        !    if (num_min_kin>this%reactive_zone%chem_syst%num_minerals_kin) error stop "Dimension error in set_num_min_kin"
        !    this%num_min_kin=num_min_kin
        !end subroutine
        !
        !subroutine set_num_kin(this) !> autentica chapuza
        !    implicit none
        !    class(solid_chemistry_c) :: this
        !    this%num_kin=this%num_min_kin+this%reactive_zone%chem_syst%num_lin_kin_reacts+&
        !    this%reactive_zone%chem_syst%num_redox_kin_reacts
        !end subroutine
        
        subroutine set_vol_fracts(this,vol_fracts)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: vol_fracts(:)
            if (this%mineral_zone%num_minerals<size(vol_fracts)) error stop "Dimension error in set_vol_fracts"
            this%vol_fracts=vol_fracts
        end subroutine
        
        subroutine allocate_vol_fracts(this)
            class(solid_chemistry_c) :: this
            allocate(this%vol_fracts(this%mineral_zone%num_minerals))
            this%vol_fracts=0d0
        end subroutine
        
        subroutine set_react_surfaces(this,react_surfaces)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: react_surfaces(:)
            if (this%mineral_zone%num_minerals<size(react_surfaces)) error stop "Dimension error in set_react_surfaces"
            this%react_surfaces=react_surfaces
        end subroutine
        
        subroutine allocate_react_surfaces(this)
            class(solid_chemistry_c) :: this
            allocate(this%react_surfaces(this%mineral_zone%num_minerals))
            this%react_surfaces=0d0
        end subroutine
        
        subroutine allocate_conc_solids(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%concentrations)) then
                deallocate(this%concentrations)
            end if
            !if (associated(this%mineral_zone)) then
            !    allocate(this%concentrations(this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin))
            !else
                allocate(this%concentrations(this%num_solids),this%conc_old(this%num_solids),this%conc_old_old(this%num_solids))
            !end if
            !this%concentrations=0d0
        end subroutine
        
        subroutine allocate_activities(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%activities)) then
                deallocate(this%activities)
            end if
            !if (associated(this%mineral_zone)) then
                !allocate(this%activities(this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin))
            !else
                allocate(this%activities(this%num_solids))
            !end if
            this%activities=0d0
        end subroutine
        
        subroutine allocate_log_act_coeffs_solid_chem(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%log_act_coeffs)) then
                deallocate(this%log_act_coeffs)
            end if
            !if (associated(this%mineral_zone)) then
            !    allocate(this%log_act_coeffs(this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin))
            !else
                allocate(this%log_act_coeffs(this%num_solids))
            !end if
            this%log_act_coeffs=0d0
        end subroutine
        
        subroutine allocate_equivalents(this)
            class(solid_chemistry_c) :: this
            if (allocated(this%equivalents)) then
                deallocate(this%equivalents)
            end if
            allocate(this%equivalents(this%reactive_zone%cat_exch_zone%num_exch_cats))
        end subroutine
        
        subroutine set_reactive_zone(this,reactive_zone)
            class(solid_chemistry_c) :: this
            class(reactive_zone_c), intent(in), target :: reactive_zone
            if (associated(reactive_zone%chem_syst)) then
                this%reactive_zone=>reactive_zone
            else
                error stop "Reactive zone object is not associated to a chemical system"
            end if
        end subroutine  
        
        subroutine set_mineral_zone(this,mineral_zone)
            class(solid_chemistry_c) :: this
            class(mineral_zone_c), intent(in), target :: mineral_zone
            if (associated(mineral_zone%chem_syst)) then
                this%mineral_zone=>mineral_zone
            else
                error stop "mineral zone object is not associated to a chemical system"
            end if
        end subroutine
        
        !subroutine allocate_conc_comp_solids(this,n_p_sol)
        !    implicit none
        !    class(solid_chemistry_c) :: this
        !    integer(kind=4), intent(in) :: n_p_sol
        !    if (allocated(this%conc_comp)) then
        !        deallocate(this%conc_comp)
        !    end if
        !    allocate(this%conc_comp(n_p_sol))
        !end subroutine
        
        subroutine compute_activities_solids(this)
            class(solid_chemistry_c) :: this
            this%activities=this%concentrations*(10**this%log_act_coeffs)
        end subroutine
        
        !subroutine compute_concentrations_ads_cats(this)
        !    implicit none
        !    class(solid_chemistry_c) :: this
        !    this%concentrations=this%concentrations*(10**this%log_act_coeffs)
        !end subroutine
        
        subroutine compute_equivalents(this)
            class(solid_chemistry_c) :: this
            integer(kind=4) :: i
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats
                this%equivalents(i)=this%reactive_zone%cat_exch_zone%surf_compl(1+i)%valence*this%concentrations(&
                this%reactive_zone%num_minerals+1+i)
            end do
        end subroutine
        
        subroutine set_CEC(this,CEC)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: CEC
            if (CEC<0d0) then
                error stop "CEC cannot be negative"
            else
                this%CEC=CEC
            end if
        end subroutine
        
        subroutine set_act_surf_compl(this,act_ads_cats)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(in) :: act_ads_cats(:)
            integer(kind=4) :: i
            this%activities(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats)=1d0-SUM(act_ads_cats)
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats
                this%activities(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats+i)=act_ads_cats(i)
            end do
        end subroutine
        
        subroutine set_indices_solids(this)
            class(solid_chemistry_c) :: this
            integer(kind=4) :: i,j,k 
            j=0
            k=0
            do i=1,this%mineral_zone%num_minerals
                if (this%mineral_zone%chem_syst%minerals(i)%mineral%cst_act_flag.eqv..false.) then
                    j=j+1
                    this%var_act_species_indices(j)=i
                else
                    k=k+1
                    this%cst_act_species_indices(k)=i
                end if
            end do
        end subroutine
        
    !> Updates concentration adsorbed cations in Newton method
        subroutine update_conc_ads_cats(this,conc_ads_cats,Delta_conc_ads_cats)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(inout) :: conc_ads_cats(:) !> concentration adsorbed cations
            real(kind=8), intent(inout) :: Delta_conc_ads_cats(:) !> adsorbed cation concentration difference
    
            integer(kind=4) :: i
            real(kind=8), allocatable :: conc_old(:)
    
            if (this%reactive_zone%CV_params%control_factor>1d0 .or. this%reactive_zone%CV_params%control_factor<0d0) then
                error stop "Control factor must be in (0,1)"
            end if
            conc_old=conc_ads_cats
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats
                if (conc_ads_cats(i)+Delta_conc_ads_cats(i)<=this%reactive_zone%CV_params%control_factor*conc_ads_cats(i)) then
                    conc_ads_cats(i)=this%reactive_zone%CV_params%control_factor*conc_ads_cats(i)
                else if (conc_ads_cats(i)+Delta_conc_ads_cats(i)>=conc_ads_cats(i)/this%reactive_zone%CV_params%control_factor) then
                    conc_ads_cats(i)=conc_ads_cats(i)/this%reactive_zone%CV_params%control_factor
                else
                    conc_ads_cats(i)=conc_ads_cats(i)+Delta_conc_ads_cats(i)
                end if
                Delta_conc_ads_cats(i)=conc_ads_cats(i)-conc_old(i)
            end do
        end subroutine
        
    !> Updates activity adsorbed cations in Newton method
        subroutine update_act_ads_cats(this,act_ads_cats,Delta_act_ads_cats)
            class(solid_chemistry_c) :: this
            real(kind=8), intent(inout) :: act_ads_cats(:) !> concentration adsorbed cations
            real(kind=8), intent(inout) :: Delta_act_ads_cats(:) !> adsorbed cation concentration difference
    
            integer(kind=4) :: i
            real(kind=8), allocatable :: act_old(:)
    
            if (this%reactive_zone%CV_params%control_factor>1d0 .or. this%reactive_zone%CV_params%control_factor<0d0) then
                error stop "Control factor must be in (0,1)"
            end if
            act_old=act_ads_cats
            do
                do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats
                    if (act_ads_cats(i)+Delta_act_ads_cats(i)<=this%reactive_zone%CV_params%control_factor*act_ads_cats(i)) then
                        act_ads_cats(i)=this%reactive_zone%CV_params%control_factor*act_ads_cats(i)
                    else if (act_ads_cats(i)+Delta_act_ads_cats(i)>=act_ads_cats(i)/this%reactive_zone%CV_params%control_factor)then
                        act_ads_cats(i)=act_ads_cats(i)/this%reactive_zone%CV_params%control_factor
                    else
                        act_ads_cats(i)=act_ads_cats(i)+Delta_act_ads_cats(i)
                    end if
                    Delta_act_ads_cats(i)=act_ads_cats(i)-act_old(i)
                end do
                if (SUM(act_ads_cats)>=1d0) then
                    act_ads_cats=act_old
                    Delta_act_ads_cats=Delta_act_ads_cats/2d0
                    if (inf_norm_vec_real(Delta_act_ads_cats)<this%reactive_zone%CV_params%abs_tol) then
                        error stop "Delta_act_ads_cats is too small"
                    end if
                else
                    exit
                end if
            end do
        end subroutine
        
        subroutine compute_conc_ads_cats(this)
            class(solid_chemistry_c) :: this
            integer(kind=4) :: i
            do i=1,this%reactive_zone%cat_exch_zone%num_exch_cats
                this%concentrations(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats+i)=&
                this%activities(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats+i)/(&
                10**this%log_act_coeffs(this%num_solids-this%reactive_zone%cat_exch_zone%num_exch_cats+i))
            end do
        end subroutine

        subroutine compute_num_solids_solid_chem(this)
            class(solid_chemistry_c) :: this
            this%num_solids=this%reactive_zone%num_solids+this%mineral_zone%num_minerals_kin
        end subroutine

        !> Computes concentration of solids at a given time
!!> phi_j*dc_s/dt=S_s^T*r_j
! subroutine compute_conc_solids(this,r_vec,time) 
!     implicit none
!     class(solid_chemistry_c) :: this
!     real(kind=8), intent(in) :: r_vec(:) !> reaction rates
!     real(kind=8), intent(in) :: time

!     integer(kind=4) :: i
    
!     do i=1,this%reactive_zone%num_solids
!         this%concentrations(i)=this%concentrations(i)+time*dot_product(this%reactive_zone%stoich_mat(:,&
!             this%reactive_zone%ind_solids(i)),r_vec)/this%vol_fracts(i)
!     end do
           
! end subroutine

!> Updates concentration solids in Newton method
subroutine update_conc_solids(this,Delta_c_s,control_factor)
    implicit none
    class(solid_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c_s(:) !> solid concentration difference
    real(kind=8), intent(in) :: control_factor !> must \f$\in (0,1)\f$
    
    integer(kind=4) :: i
    real(kind=8), allocatable :: conc_old(:)
    
    if (control_factor>1d0 .or. control_factor<0d0) error stop "Control factor must be in (0,1)"
    conc_old=this%concentrations
    do i=1,size(this%concentrations)
        if (this%concentrations(i)+Delta_c_s(i)<=control_factor*this%concentrations(i)) then
            this%concentrations(i)=control_factor*this%concentrations(i)
        else if (this%concentrations(i)+Delta_c_s(i)>=this%concentrations(i)/control_factor) then
            this%concentrations(i)=this%concentrations(i)/control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c_s(i)
        end if
        Delta_c_s(i)=this%concentrations(i)-conc_old(i)
    end do
end subroutine

!> This subroutine checks for zero concentrations in a solid chemistry object
subroutine check_solid_chemistry(this,flag,indices)
    class(solid_chemistry_c) :: this
    !real(kind=8), intent(in) :: tolerance !> tolerance for concentrations of solids
    integer(kind=4), intent(out) :: flag !> 1 if no zero concentrations, 0 otherwise
    integer(kind=4), intent(out), allocatable :: indices(:) !> indices of zero concentrations
    
    integer(kind=4) :: num_new_non_flowing_species,i,j,k
    integer(kind=4), allocatable :: old_nf_ind(:),old_solid_ind(:)
    type(species_c), allocatable :: new_non_flowing_species(:)
    type(solid_chemistry_c), allocatable :: new_solid_chems(:)
    real(kind=8), parameter :: epsilon=1d-9 !> arbitrario
    
    flag=1
    do i=1,this%reactive_zone%num_solids
        if (this%concentrations(i)<this%reactive_zone%CV_params%abs_tol) then
            flag=0
            indices=[indices,j] !> indices of zero concentration solids
        else
            continue
        end if
    end do
end subroutine

!> Computes concentration of minerals after a given time step
!! We assume minerals have constant activity
subroutine compute_conc_minerals_iter_EE(this,Delta_t)
!> Arguments
    class(solid_chemistry_c) :: this !> solid chemistry object
    real(kind=8), intent(in) :: Delta_t !> time step
!> Variables
    integer(kind=4) :: i !> counter minerals
    real(kind=8), parameter :: eps=1d-16 !> chapuza
!> Process
    do i=1,this%reactive_zone%num_minerals
        if (abs(this%vol_fracts(this%mineral_zone%num_minerals_kin+i))<eps) then
            continue
        else
            this%concentrations(this%mineral_zone%num_minerals_kin+i)=this%concentrations(this%mineral_zone%num_minerals_kin+i)+&
                Delta_t*dot_product(this%reactive_zone%stoich_mat(:,this%reactive_zone%ind_mins_stoich_mat(i)),&
                this%r_eq(1:this%reactive_zone%num_minerals))/this%vol_fracts(this%mineral_zone%num_minerals_kin+i)
        end if
    end do
    do i=1,this%mineral_zone%num_minerals_kin
        if (abs(this%vol_fracts(i))<eps) then
            continue
        else
            this%concentrations(i)=this%concentrations(i)+Delta_t*dot_product(this%mineral_zone%chem_syst%stoich_mat(:,&
            this%mineral_zone%ind_min_Sk(i)),this%rk(1:this%mineral_zone%num_minerals_kin))/&
            this%vol_fracts(i)
        end if
    end do
end subroutine

subroutine modify_mix_ratios_Rk(this,mix_ratio_Rk_init,c_tilde,Delta_t,rk_tilde,mix_ratio_Rk_new,num_lump)
    !> This subroutine modifies the mixing ratios of the kinetic reaction rates for a target
    !> AQUI DEBERIAS GUARDAR LOS NUEVOS LAMBDAS
        class(solid_chemistry_c), intent(in) :: this !> solid chemistry object
        real(kind=8), intent(in) :: mix_ratio_Rk_init !> initial mixing ratio of the kinetic reaction rates in this target
        real(kind=8), intent(in) :: c_tilde(:) !> concentration vector
        real(kind=8), intent(in) :: Delta_t !> time step
        real(kind=8), intent(inout) :: rk_tilde(:) !> kinetic reaction rate contributions
        real(kind=8), intent(out) :: mix_ratio_Rk_new !> new mixing ratio of the kinetic reaction rates
        integer(kind=4), intent(out) :: num_lump !> number of lumpings

        integer(kind=4) :: i !> loop index
        integer(kind=4) :: n_nc !> number of variable activity species
        real(kind=8), parameter :: alpha=1.05 !> mixing ratio factor
        real(kind=8), allocatable :: R_tilde(:) !> reaction amount contributions after mixing

        R_tilde=Delta_t*rk_tilde !> we compute reaction part of concentrations after mixing
        n_nc=this%reactive_zone%speciation_alg%num_var_act_species
        mix_ratio_Rk_new=mix_ratio_Rk_init !> we initialise new mixing ratio (chapuza)
        !flag=.false. !> we initialise lumping flag
        num_lump=0 !> we initialise number of lumpings
        i=1
        !> We check negative concentrations after mixing
        do
            !if (c_tilde(i)+R_tilde_up(i)<=0d0) then
                !error stop "Negative concentration after mixing upstream waters in subroutine water_mixing_iter_EI_kin_aq_anal_ideal_opt2"
            !if (c_tilde(i)+R_tilde_up(i)+R_tilde_down(i)<=0d0 .and. R_tilde_down(i)<0d0) then
            if (c_tilde(i)+R_tilde(i)<=0d0) then
                !flag=.true. !> we set lumping flag
                num_lump=num_lump+1 !> we increase number of lumpings
                !print *, "Warning: negative concentration after mixing"
                !print *, "Species index: ", i
                !R_tilde_down(i)=R_tilde_down(i)/2d0 !> autentica chapuza
                !R_tilde_down=R_tilde_down*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new) !> we reduce contributions of reaction amounts in this target (chapuza)
                !R_tilde_up=R_tilde_up*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new) !> we reduce contributions of reaction amounts in this target (chapuza)
                R_tilde=R_tilde*(1d0-Delta_t*alpha*mix_ratio_Rk_new)/(1d0-Delta_t*mix_ratio_Rk_new) !> we reduce contributions of reaction amounts in this target (chapuza)
                mix_ratio_Rk_new=mix_ratio_Rk_new*alpha !> we increase mixing ratio in this target (chapuza)
                if (Delta_t*mix_ratio_Rk_new>1d0) then
                    !> chapuza
                    mix_ratio_Rk_new=1d0/Delta_t
                    !R_tilde_down=0d0
                    !R_tilde_up=0d0
                    rk_tilde=0d0
                    exit
                end if
            else if (i<n_nc) then
                i=i+1
            else
                rk_tilde=R_tilde/Delta_t !> we update kinetic reaction rate contributions
                exit
            end if
        end do
end subroutine

subroutine assign_solid_chemistry(this,solid_chemistry)
class(solid_chemistry_c) :: this !> solid chemistry object
type(solid_chemistry_c), intent(in) :: solid_chemistry !> solid chemistry object to assign
if (associated(solid_chemistry%reactive_zone)) then
    this%reactive_zone=>solid_chemistry%reactive_zone
else
    error stop "Reactive zone not associated with solid chemistry"
end if
if (associated(solid_chemistry%mineral_zone)) then
    this%mineral_zone=>solid_chemistry%mineral_zone
else
    error stop "Mineral zone not associated with solid chemistry"
end if
this%num_solids=solid_chemistry%num_solids
this%concentrations=solid_chemistry%concentrations
this%activities=solid_chemistry%activities
this%log_act_coeffs=solid_chemistry%log_act_coeffs
this%name=solid_chemistry%name
if (allocated(solid_chemistry%r_eq)) then
    this%r_eq=solid_chemistry%r_eq
end if
if (allocated(solid_chemistry%rk)) then
    this%rk=solid_chemistry%rk
end if
if (allocated(solid_chemistry%vol_fracts)) then
    this%vol_fracts=solid_chemistry%vol_fracts
end if
if (allocated(solid_chemistry%react_surfaces)) then
    this%react_surfaces=solid_chemistry%react_surfaces
end if
this%equivalents=solid_chemistry%equivalents
this%CEC=solid_chemistry%CEC
end subroutine

subroutine set_target(this,target_obj)
    class(solid_chemistry_c) :: this
    type(target_c), intent(in), target :: target_obj
    this%tar=>target_obj
end subroutine

subroutine allocate_reaction_rates_solid_chem(this)
    class(solid_chemistry_c) :: this !> solid chemistry object
    if (.not. allocated(this%r_eq)) then
        if (.not. associated(this%reactive_zone)) then
            error stop "Reactive zone not associated with solid chemistry"
        else
            allocate(this%r_eq(this%reactive_zone%num_minerals+&
                this%reactive_zone%cat_exch_zone%num_exch_cats))
            allocate(this%Re_mean(this%reactive_zone%num_minerals+&
                this%reactive_zone%cat_exch_zone%num_exch_cats))
            this%r_eq=0d0 !> by default
            this%Re_mean=0d0 !> by default
        end if
    end if
    !allocate(this%solid_chemistry%r_eq(this%solid_chemistry%reactive_zone%num_minerals+&
    !this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
    !this%solid_chemistry%r_eq=0d0 !> by default
    if (.not. allocated(this%rk)) then
        if (.not. associated(this%mineral_zone)) then
            error stop "Mineral zone not associated with solid chemistry"
        else
            allocate(this%rk(this%mineral_zone%num_minerals_kin))
            allocate(this%rk_old(this%mineral_zone%num_minerals_kin))
            allocate(this%rk_old_old(this%mineral_zone%num_minerals_kin))
            allocate(this%rk_old_old_old(this%mineral_zone%num_minerals_kin))
            allocate(this%Rk_mean(this%mineral_zone%num_minerals_kin))
            allocate(this%Rk_est(this%mineral_zone%num_minerals_kin))
            this%rk=0d0 !> by default
            this%Rk_est=0d0 !> by default
            this%rk_old=0d0 !> by default
            this%rk_old_old=0d0 !> by default
            this%rk_old_old_old=0d0 !> by default
            this%Rk_mean=0d0 !> by default
        end if
    end if
end subroutine
end module