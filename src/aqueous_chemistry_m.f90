!> Aqueous chemistry subclass:
!!>   contains local information of aqueous solution
!!>   points to solid chemistry, gas chemistry and aqueous phase classes
module aqueous_chemistry_m
    use local_chemistry_m, only: local_chemistry_c
    use solid_chemistry_m, only: solid_chemistry_c, mineral_c, mineral_zone_c
    use gas_chemistry_m, only: gas_chemistry_c
    use chem_system_m, only: chem_system_c, kin_reaction_c
    use reactive_zone_m, only: reactive_zone_c, aq_phase_c, reaction_c, eq_reaction_c, & 
        cat_exch_c, speciation_algebra_c, aq_species_c, species_c, gas_phase_c, gas_c
    use metodos_sist_lin_m, only: diag_matrix_c, real_array_c, int_array_c, tridiag_matrix_c, id_matrix, inf_norm_vec_real, &
        norm_mat_inf, LU_lin_syst, outer_prod_vec
    use params_aq_sol_m, only: params_aq_sol_t
    use array_ops_m, only: get_ind_int, is_int_in_1D_array
    implicit none
    save
    type, public, extends(local_chemistry_c) :: aqueous_chemistry_c
        real(kind=8), allocatable :: pos(:) !> position of the object in space
        type(params_aq_sol_t) :: params_aq_sol !> parameters of aqueous solution
        real(kind=8) :: ionic_act !> ionic activity I
        real(kind=8) :: pH
        real(kind=8) :: pe
        real(kind=8) :: salinity
        real(kind=8) :: alkalinity !> m_OH- + m_HCO3- + 2*m_CO3-2
        class(solid_chemistry_c), pointer :: solid_chemistry=>null()
        class(gas_chemistry_c), pointer :: gas_chemistry=>null()
        class(aq_phase_c), pointer :: aq_phase=>null()
        integer(kind=4), allocatable :: indices_aq_species(:) !> indices of aqueous species in state variables
        integer(kind=4), allocatable :: indices_aq_phase(:) !> indices of aqueous species in aqueous phase
        !integer(kind=4) :: ind_wat=0 !> index of water in aqueous phase pointer (if not present = 0) (NO NECESARIO)
        !integer(kind=4), allocatable :: ind_diss_solids(:) !> indices of dissolved solids in indices_aq_species attribute (NO NECESARIO)
        integer(kind=4), allocatable :: ind_prim_species(:) !> indices of aqueous primary species in indices_aq_species attribute
        integer(kind=4), allocatable :: ind_var_act_species(:) !> indices of variable activity species in indices_aq_species attribute
        integer(kind=4), allocatable :: ind_sec_species(:) !> indices of aqueous secondary species in indices_aq_species attribute (sorted in variable & constant activity)
        !integer(kind=4), allocatable :: ind_cst_act_species(:) !> indices of constant activity species in indices_aq_species attribute
        ! integer(kind=4) :: ind_prot=0 !> index 'H+' in indices_aq_species attribute
        ! integer(kind=4) :: ind_oh=0 !> index 'OH-' in indices_aq_species attribute
        ! integer(kind=4) :: ind_bicarb=0 !> index 'HCO3-' in indices_aq_species attribute
        ! integer(kind=4) :: ind_carb=0 !> index 'CO3-2' in indices_aq_species attribute
        type(int_array_c) :: indices_rk !> indices of species, in chemical system, relevant for kinetic reaction rates
        !type(vector_int_c) :: ind_up_down !> indices of upstream & downstream objects of this class that mix with object of this class
    contains
        procedure, public :: read_wat_type_CHEPROO
    !> Set
        procedure, public :: set_concentrations=>set_conc_aq_species
        procedure, public :: set_pos
        procedure, public :: set_act_aq_species
        procedure, public :: set_act_diss_solids
        procedure, public :: set_conc_aq_prim_species
        procedure, public :: set_conc_prim_species
        procedure, public :: set_conc_var_act_species
        procedure, public :: set_conc_sec_var_act_species
        procedure, public :: set_conc_sec_aq_species
        procedure, public :: set_conc_sec_aq_var_act_species
        procedure, public :: set_conc_sec_species
        procedure, public :: set_indices_aq_species_aq_chem
        !procedure, public :: set_indices_aq_phase_aq_chem
        procedure, public :: set_spec_alg_aq_chem
        procedure, public :: set_ind_species
        !procedure, public :: set_ind_wat_aq_chem
        !procedure, public :: set_ind_prot_aq_chem
        procedure, public :: set_pH
        procedure, public :: set_pe
        procedure, public :: set_ionic_act
        procedure, public :: set_conc_single_species
        procedure, public :: set_log_act_coeffs
        procedure, public :: set_solid_chemistry
        procedure, public :: set_gas_chemistry
        procedure, public :: set_aq_phase
        procedure, public :: set_indices_rk
        procedure, public :: set_Rk
        procedure, public :: set_Rk_est
        procedure, public :: set_rk_mean
        !procedure, public :: set_conc_old
        !procedure, public :: set_conc_old_old
    !> Allocate
        procedure, public :: allocate_reaction_rates
        procedure, public :: allocate_conc_aq_species
        procedure, public :: allocate_activities_aq_species
        procedure, public :: allocate_log_act_coeffs_aq_chem
        procedure, public :: allocate_indices_aq_species_aq_chem
        procedure, public :: allocate_indices_aq_phase_aq_chem
        !procedure, public :: allocate_ind_diss_solids_aq_chem
    !> Compute
        procedure, public :: compute_activities_diss_solids
        procedure, public :: compute_activities_aq
        procedure, public :: compute_activities
        procedure, public :: compute_activity
        procedure, public :: compute_act_from_MAL
        procedure, public :: compute_act_water
        procedure, public :: compute_molarities
        procedure, public :: compute_molalities
        procedure, public :: compute_salinity
        procedure, public :: compute_alkalinity
        procedure, public :: compute_pH
        procedure, public :: compute_conc_comp
        procedure, public :: compute_conc_comp_cst_act
        !procedure, public :: compute_re_mean
        procedure, public :: compute_Re_lump
        procedure, public :: compute_Re_rk_lump
        procedure, public :: compute_Re
        procedure, public :: compute_Re_rk
        procedure, public :: compute_rk_new
        procedure, public :: compute_Jacobian_rk_anal
        procedure, public :: compute_rk_Jac_rk_anal
        procedure, public :: compute_rk_Jac_rk_incr_coeff
        procedure, public :: compute_log_K_aq_chem
        procedure, public :: compute_ionic_act
        procedure, public :: compute_d_log_gamma_d_I_aq_chem
        procedure, public :: compute_res_Jac_res_anal
        procedure, public :: compute_res_Jac_res_anal_ideal
        procedure, public :: compute_res_Jac_res_anal_exch
        procedure, public :: compute_res_Jac_res_incr_coeff
        procedure, public :: compute_res_init
        procedure, public :: compute_log_act_coeff_wat
        procedure, public :: compute_saturation_kin_min
        procedure, public :: compute_Rk
    !> Get
        procedure, public :: get_indices_reaction
        procedure, public :: get_rk_new !> gets new kinetic reaction rates
        procedure, public :: get_rk_old !> gets old kinetic reaction rates
        procedure, public :: get_Sk_nc !> gets kinetic stoichiometric matrix var act species
        procedure, public :: get_c1
        procedure, public :: get_a1
        procedure, public :: get_conc_nc_aq_perm
        procedure, public :: get_c1_aq
        procedure, public :: get_c2
        procedure, public :: get_c2nc
        procedure, public :: get_conc_nc
        procedure, public :: get_conc_nc_old
        procedure, public :: get_conc
        procedure, public :: get_log_gamma2nc
        procedure, public :: get_log_gamma2
        procedure, public :: get_u_aq
    !> Speciation
        procedure, public :: compute_c2_from_c1_aq_ideal
        procedure, public :: compute_c2_from_a1_ideal
        procedure, public :: compute_c2nc_from_c1_aq_ideal
        procedure, public :: compute_c2nc_from_c1_ideal
        procedure, public :: compute_c2nc_from_c1_aq_Picard
        procedure, public :: compute_c2nc_from_c1_Picard
        procedure, public :: compute_c2_from_c1_aq_Picard
        procedure, public :: compute_c2_from_c1_Picard
        procedure, public :: compute_c_nc_from_u_Newton
        procedure, public :: compute_c_nc_from_u_aq_Newton
        procedure, public :: compute_c_nc_from_u_aq_Newton_ideal
        procedure, public :: compute_c_nc_from_u_Newton_ideal
        procedure, public :: compute_residual
        procedure, public :: compute_residual_cst_act
        procedure, public :: compute_dc2nc_dc1_gamma_cst
        procedure, public :: compute_dc2_dc1_ideal
        procedure, public :: compute_dc2nc_dc1_ideal
        procedure, public :: compute_dc2nc_dc1_aq_ideal
        procedure, public :: compute_dc2nc_dc1
        procedure, public :: compute_dc2nc_dc1_aq
        procedure, public :: compute_dc2_dc1
        procedure, public :: compute_c_from_u_bin_syst
    !> Mixing
        !procedure, public :: interfaz_comps_vars
        procedure, public :: mixing_iter_comp
        procedure, public :: mixing_iter_comp_lump
        procedure, public :: mixing_iter_comp_ideal
        procedure, public :: mixing_iter_comp_ideal_lump
        procedure, public :: mixing_iter_species
        procedure, public :: mixing_iter_comp_exch
        procedure, public :: mixing_iter_comp_exch_lump
        procedure, public :: mixing_iter_comp_exch_ideal
        procedure, public :: water_mixing_iter_EE_eq_kin
        procedure, public :: water_mixing_iter_EE_eq_kin_ideal
        procedure, public :: water_mixing_iter_EE_eq_kin_ideal_lump
        procedure, public :: water_mixing_iter_EE_kin_ideal
        procedure, public :: water_mixing_iter_EE_kin_ideal_lump
        procedure, public :: water_mixing_iter_EfI_eq_kin_anal_lump
        procedure, public :: water_mixing_iter_EfI_eq_kin_anal_ideal
        procedure, public :: water_mixing_iter_EI_eq_kin_anal_ideal_opt2
        procedure, public :: water_mixing_iter_EI_eq_kin_anal_ideal_opt2_lump
        !procedure, public :: water_mixing_iter_EfI_eq_kin_anal_ideal_opt2
        procedure, public :: water_mixing_iter_EfI_eq_kin_anal_ideal_lump
        procedure, public :: water_mixing_iter_EfI_kin_anal
        procedure, public :: water_mixing_iter_EfI_kin_anal_ideal
        procedure, public :: water_mixing_iter_EI_kin_anal_ideal_opt2
        procedure, public :: water_mixing_iter_EI_kin_anal_ideal_opt2_lump
        procedure, public :: water_mixing_iter_EfI_kin_anal_ideal_lump
        !procedure, public :: compute_c_tilde
        procedure, public :: compute_u_tilde
        procedure, public :: reaction_iteration_EE_eq_kin
        procedure, public :: reaction_iteration_EE_eq_kin_lump
        procedure, public :: reaction_iteration_EE_kin_lump
        procedure, public :: reaction_iteration_EE_kin
        procedure, public :: compute_dfk_dc_aq_EfI
        procedure, public :: compute_dfk_dc_aq_EfI_lump
        procedure, public :: compute_dfk_dc1_aq_EfI_lump
        procedure, public :: compute_dfk_dc1_aq_EfI_ideal
        procedure, public :: compute_dfk_dc1_aq_EfI_ideal_lump
        procedure, public :: Newton_EfI_rk_eq_kin_aq_anal_lump
        procedure, public :: Newton_EfI_rk_eq_kin_aq_anal_ideal
        procedure, public :: Newton_EI_rk_eq_kin_aq_anal_ideal_opt2
        procedure, public :: Newton_EI_rk_eq_kin_aq_anal_ideal_opt2_lump
        procedure, public :: Newton_EI_rk_kin_aq_anal_ideal_opt2
        procedure, public :: Newton_EI_rk_kin_aq_anal_ideal_opt2_lump
        procedure, public :: Newton_EfI_rk_eq_kin_aq_anal_ideal_lump
        !procedure, public :: Newton_EfI_rk_kin_aq_anal
        procedure, public :: Newton_EfI_rk_kin_aq_anal_ideal
        procedure, public :: Newton_EfI_rk_kin_aq_anal_ideal_lump
        procedure, public :: solve_reactive_mixing_iter !> 
    !> Update
        procedure, public :: update_conc_aq_prim_species
        procedure, public :: update_conc_aq_prim_species_bis
        procedure, public :: update_conc_prim_species
        procedure, public :: update_conc_aq_species
        !procedure, public :: update_conc_sec_species
        !procedure, public :: update_conc_sec_aq_species
        !procedure, public :: update_conc_sec_var_act_species
        procedure, public :: change_spec_alg_aq_chem
        procedure, public :: update_conc_nc
    !> Check
        procedure, public :: check_conc_aq_var_act_species
        procedure, public :: check_conc_var_act_species
        procedure, public :: check_act_aq_species
        procedure, public :: check_dc2nc_dc1_aq
        procedure, public :: check_dc2nc_dc1
        procedure, public :: check_dc2_dc1
    !> Initialise
        procedure, public :: initialise_conc_incr_coeff
        procedure, public :: initialise_conc_anal
        procedure, public :: initialise_conc_anal_ideal
        procedure, public :: initialise_conc_anal_exch
        procedure, public :: initialise_conc_nc_iterative_method_ideal
        procedure, public :: initialise_c1_aq_iterative_method
    !> Rearrange
        procedure, public :: rearrange_state_vars
    end type
        
    interface
        subroutine compute_c2nc_from_c1_aq_ideal(this,c2nc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(out) :: c2nc(:)
        end subroutine
        
        subroutine compute_c2nc_from_c1_ideal(this,c1,c2nc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1(:)
            real(kind=8), intent(out) :: c2nc(:)
        end subroutine
    end interface
    
    interface
        
        subroutine compute_c2_from_c1_aq_ideal(this,c2)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(out) :: c2(:) 
        end subroutine 
        
        subroutine compute_c2_from_a1_ideal(this,a1,log_gamma2,c2)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: a1(:)
            real(kind=8), intent(in) :: log_gamma2(:)
            real(kind=8), intent(out) :: c2(:) 
        end subroutine
    end interface
    
    interface
        
        subroutine compute_c2nc_from_c1_Picard(this,c1,c2nc_ig,c2nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1(:) !> primary concentrations
            real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(out) :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
    end interface
    
    interface
        subroutine solve_reactive_mixing_iter(this,c1_old,mixing_ratios_conc,conc_old,porosity,Delta_t,solver)
            import aqueous_chemistry_c
            import real_array_c
            import int_array_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:) !> concentrations primary species previous two time steps
            real(kind=8), intent(in) :: mixing_ratios_conc(:) !> mixing ratios array
            real(kind=8), intent(in) :: conc_old(:,:) !> concentrations previous time step
            real(kind=8), intent(in) :: porosity !>
            real(kind=8), intent(in) :: Delta_t !> time step
            external :: solver
        end subroutine
        
        subroutine compute_c2nc_from_c1_aq_Picard(this,c2nc_ig,c2nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(out) :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine 
        
        subroutine compute_c2_from_c1_aq_Picard(this,c2_ig,c2,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2_ig(:) !>  (dim=n_eq)
            real(kind=8), intent(out) :: c2(:) !>  (dim=n_eq)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_c2_from_c1_Picard(this,c1,c2_ig,c2,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1(:) !>  (dim=n_p)
            real(kind=8), intent(in) :: c2_ig(:) !>  (dim=n_eq)
            real(kind=8), intent(out) :: c2(:) !>  (dim=n_eq)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_c_nc_from_u_Newton(this,c1_ig,c2nc_ig,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_ig(:) !> initial guess primary concentrations
            real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: conc_comp(:) !> component concentrations
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_c_nc_from_u_aq_Newton(this,c2nc_ig,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: conc_comp(:) !> component concentrations
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_c_nc_from_u_aq_Newton_ideal(this,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_comp(:) !> component concentrations
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
        
        subroutine compute_c_nc_from_u_Newton_ideal(this,c1_ig,conc_comp,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_ig(:) !> initial guess primary concentrations
            real(kind=8), intent(in) :: conc_comp(:) !> component concentrations
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine
                
        subroutine initialise_conc_incr_coeff(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: icon(:)
            integer(kind=4), intent(in) :: n_icon(:)
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine 
        
        subroutine initialise_conc_anal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag,path_py)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: icon(:)
            integer(kind=4), intent(in) :: n_icon(:)
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
            character(len=*), intent(in), optional :: path_py !> path output Python
        end subroutine
        
        subroutine initialise_conc_anal_ideal(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag,path_py)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: icon(:)
            integer(kind=4), intent(in) :: n_icon(:)
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
            character(len=*), intent(in), optional :: path_py !> path output Python
        end subroutine
        
        subroutine initialise_conc_anal_exch(this,icon,n_icon,indices_constrains,ctot,niter,CV_flag)
            import aqueous_chemistry_c
            import solid_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: icon(:)
            integer(kind=4), intent(in) :: n_icon(:)
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            !class(solid_chemistry_c), intent(inout) :: surf_chem !> surface chemistry
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine 
        
        ! subroutine compute_residual(this,conc_comp,c_nc,residual)
        !     import aqueous_chemistry_c
        !     implicit none
        !     class(aqueous_chemistry_c), intent(in) :: this
        !     real(kind=8), intent(in) :: conc_comp(:)
        !     real(kind=8), intent(in) :: c_nc(:)
        !     real(kind=8), intent(out) :: residual(:)
        ! end subroutine
        

        
        ! subroutine compute_residual_cst_act(this,conc_comp,conc,residual)
        !     import aqueous_chemistry_c
        !     implicit none
        !     class(aqueous_chemistry_c), intent(in) :: this
        !     real(kind=8), intent(in) :: conc_comp(:) !> component concentrations
        !     real(kind=8), intent(in) :: conc(:) !> species concentrations
        !     real(kind=8), intent(out) :: residual(:)
        ! end subroutine
        
        subroutine compute_dc2nc_dc1_gamma_cst(this,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine
        
        subroutine compute_dc2_dc1_ideal(this,c1,c2,dc2_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1(:) !>  (dim=n_p)
            real(kind=8), intent(in) :: c2(:) !>  (dim=n_eq)
            real(kind=8), intent(out) :: dc2_dc1(:,:)
        end subroutine
        
        subroutine compute_dc2nc_dc1_ideal(this,c1,c2nc,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1(:) !>  (dim=n_p)
            real(kind=8), intent(in) :: c2nc(:) !>  (dim=n_eq)
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine
        
        subroutine compute_dc2nc_dc1_aq_ideal(this,c2nc,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc(:) !>  (dim=n_eq)
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine
        
        subroutine compute_dc2nc_dc1(this,c1,c2nc,out_prod,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1(:) !> 
            real(kind=8), intent(in) :: c2nc(:) !> 
            real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_nc_d_I and z_nc^2
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine
        
        subroutine compute_dc2nc_dc1_aq(this,c2nc,out_prod,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc(:)
            real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_nc_d_I and z_nc^2
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine
                
        subroutine compute_dc2_dc1(this,out_prod,c1,c2,dc2_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_d_I and z^2
            real(kind=8), intent(in) :: c1(:) !>  (dim=n_p)
            real(kind=8), intent(in) :: c2(:) !>  (dim=n_eq)
            real(kind=8), intent(out) :: dc2_dc1(:,:)
        end subroutine
        
        
        subroutine compute_c2nc_from_c1_expl_homog(this)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
        end subroutine 
        
        subroutine compute_c2nc_from_c1_Picard_homog(this,tolerance,niter_max,niter)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: tolerance
            integer(kind=4), intent(in) :: niter_max
            integer(kind=4), intent(out) :: niter !> number of iterations
        end subroutine 
        
        subroutine compute_c_nc_from_u_Newton_homog(this,tolerance,rel_tolerance,control_factor,niter_max,niter,CV_flag)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: tolerance
            real(kind=8), intent(in) :: rel_tolerance
            real(kind=8), intent(in) :: control_factor
            integer(kind=4), intent(in) :: niter_max
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> TRUE if converges, FALSE otherwise
        end subroutine 
        
        subroutine compute_residual_homog(this,residual)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: residual(:)
        end subroutine
        
        subroutine compute_dc2nc_dc1_gamma_cst_homog(this,tol,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: tol
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine
        
        subroutine compute_dc2nc_dc1_homog(this,tol,dc2nc_dc1)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: tol
            real(kind=8), intent(out) :: dc2nc_dc1(:,:)
        end subroutine 
        
    !     subroutine update_conc_aq_species(this,Delta_c)
    !         import aqueous_chemistry_c
    !         implicit none
    !         class(aqueous_chemistry_c) :: this
    !         real(kind=8), intent(inout) :: Delta_c(:)
    !         !real(kind=8), intent(in) :: control_factor !> must $\in (0,1)$
    !     end subroutine
        
    ! !> Updates concentration aqueous primary species in iterative method
    !     subroutine update_conc_aq_prim_species(this,Delta_c1)
    !         import aqueous_chemistry_c
    !         implicit none
    !         class(aqueous_chemistry_c) :: this
    !         real(kind=8), intent(inout) :: Delta_c1(:)
    !     end subroutine
        
    !> Updates concentration aqueous and solid primary species in iterative method
        ! subroutine update_conc_prim_species(this,c1,Delta_c1)
        !     import aqueous_chemistry_c
        !     implicit none
        !     class(aqueous_chemistry_c) :: this
        !     real(kind=8), intent(inout) :: c1(:)
        !     real(kind=8), intent(inout) :: Delta_c1(:)
        ! end subroutine
        
        ! subroutine initialise_iterative_method(conc_old_old,conc_old,param,initial_guess)
        !     implicit none
        !     real(kind=8), intent(in) :: conc_old_old(:)
        !     real(kind=8), intent(in) :: conc_old(:)
        !     real(kind=8), intent(in) :: param
        !     real(kind=8), intent(out) :: initial_guess(:)
        ! end subroutine
        

        
        subroutine compute_log_act_coeffs(this)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
        end subroutine
        
        subroutine compute_log_Jacobian_act_coeffs(this,d_log_gamma_d_log_c)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(out) :: d_log_gamma_d_log_c(:,:)
        end subroutine
        
        subroutine mixing_iter_comp(this,c1_old,c2nc_ig,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            !real(kind=8), intent(in), optional :: porosity
        end subroutine

        subroutine mixing_iter_comp_lump(this,c1_old,c2nc_ig,c_tilde,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            !real(kind=8), intent(in), optional :: porosity
        end subroutine
        
        subroutine mixing_iter_comp_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(inout) :: mix_ratio_Rk
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            !real(kind=8), intent(in), optional :: porosity
            !real(kind=8), intent(in), optional :: Delta_t !> time step
        end subroutine

        subroutine mixing_iter_comp_ideal_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: rk_tilde(:)
            !real(kind=8), intent(in) :: mix_ratio_Rk
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            !real(kind=8), intent(in), optional :: porosity
            !real(kind=8), intent(in), optional :: Delta_t !> time step
        end subroutine
        
        subroutine mixing_iter_comp_exch(this,c1_old,c2nc_ig,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            !real(kind=8), intent(in), optional :: porosity
        end subroutine

        subroutine mixing_iter_comp_exch_lump(this,c1_old,c2nc_ig,c_tilde,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            !real(kind=8), intent(in), optional :: porosity
        end subroutine
        
        subroutine mixing_iter_comp_exch_ideal(this,c1_old,c_tilde,conc_nc,porosity,Delta_t)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            real(kind=8), intent(in), optional :: porosity
            real(kind=8), intent(in), optional :: Delta_t !> time step
        end subroutine
        
        subroutine mixing_iter_species(this,c1_old,c2nc_ig,c_tilde,conc_nc,conc_comp,porosity,Delta_t)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            !class(aqueous_chemistry_c), intent(in) :: this_old !> aqueous chemistry object at previous time step
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(out) :: conc_nc(:)
            real(kind=8), intent(out) :: conc_comp(:) !> concentration components
            real(kind=8), intent(in), optional :: porosity !> storage matrix
            real(kind=8), intent(in), optional :: Delta_t !> time step
        end subroutine
            
        
        
        subroutine mixing_iter_EI_aq_chem(this,c_tilde,rk_mat,porosity,Delta_t)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_tilde(:,:)
            real(kind=8), intent(in), optional :: rk_mat(:,:) !> rk matrix
            real(kind=8), intent(in), optional :: porosity !> storage matrix
            real(kind=8), intent(in), optional :: Delta_t !> time step
        end subroutine
        
        subroutine mixing_iter_EfI_aq_chem(this,porosity,Delta_t,tolerance,rel_tolerance,control_factor,niter_max)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            !class(aqueous_chemistry_c), intent(in) :: mixing_waters(:) !> all of them have same aqueous species and components
            !real(kind=8), intent(in) :: mixing_ratios_conc(:) !> A^(-T)_j
            real(kind=8), intent(in), optional :: porosity !> storage matrix
            real(kind=8), intent(in), optional :: Delta_t !> time step
            real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EE_eq_kin(this,c1_old,c2nc_ig,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
        end subroutine
        
        subroutine water_mixing_iter_EE_eq_kin_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity !> storage matrix
            
        end subroutine
        
        subroutine water_mixing_iter_EE_eq_kin_ideal_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            !real(kind=8), intent(in) :: rk_tilde(:)
            !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity !> storage matrix
            
        end subroutine
                
        subroutine water_mixing_iter_EE_kin_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity !> storage matrix
            
        end subroutine
        
        subroutine water_mixing_iter_EE_kin_ideal_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity !> storage matrix
        end subroutine
        
        subroutine water_mixing_iter_EfI_eq_kin_anal_lump(this,c1_old,c2nc_ig,c_tilde,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !class(aqueous_chemistry_c), intent(in) :: this_old !> aqueous chemistry object at previous time step
            !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            !real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            !real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            !real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            !integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EfI_eq_kin_anal_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !integer(kind=4), intent(in) :: ind         !> index of object "this" in mixing ratios array
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            
            !real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            !real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            !real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            !integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EI_eq_kin_anal_ideal_opt2(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(inout) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            
            !real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            !real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            !real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            !integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EI_eq_kin_anal_ideal_opt2_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: rk_tilde(:)
            !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            
            !real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            !real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            !real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            !integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EfI_eq_kin_anal_ideal_opt2(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,conc_nc,porosity,Delta_t)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            real(kind=8), intent(in), optional :: porosity
            real(kind=8), intent(in), optional :: Delta_t !> time step
            !real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            !real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            !real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            !integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EfI_eq_kin_anal_ideal_lump(this,c1_old,c_tilde,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !class(aqueous_chemistry_c), intent(in) :: this_old !> aqueous chemistry object at previous time step
            !integer(kind=4), intent(in) :: ind             !> index of object "this" in mixing ratios array
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            
            !real(kind=8), intent(in), optional :: tolerance !> for Newton CV
            !real(kind=8), intent(in), optional :: rel_tolerance !> for Newton/Picard CV
            !real(kind=8), intent(in), optional :: control_factor !> controls Delta_c1_j in Newton
            !integer(kind=4), intent(in), optional :: niter_max !> number maximum of iterations
        end subroutine
        
        subroutine water_mixing_iter_EfI_kin_anal(this,c1_old,c2nc_ig,c_tilde,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in) :: porosity
            
        end subroutine
        
        subroutine water_mixing_iter_EfI_kin_anal_ideal(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
            
        end subroutine

        subroutine water_mixing_iter_EI_kin_anal_ideal_opt2(this,c1_old,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: rk_tilde(:)
            !real(kind=8), intent(in) :: rk_tilde_up(:)
            !real(kind=8), intent(in) :: rk_tilde_down(:)
            real(kind=8), intent(inout) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
        end subroutine
        
        subroutine water_mixing_iter_EI_kin_anal_ideal_opt2_lump(this,c1_old,c_tilde,Delta_t,theta,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            !import time_discr_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: rk_tilde(:)
            !real(kind=8), intent(in) :: rk_tilde_up(:)
            !real(kind=8), intent(in) :: rk_tilde_down(:)
            !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
        end subroutine
        
        subroutine water_mixing_iter_EfI_kin_anal_ideal_lump(this,c1_old,c_tilde,Delta_t,conc_nc)
            import aqueous_chemistry_c
            import diag_matrix_c 
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c1_old(:)
            !real(kind=8), intent(in) :: c2nc_ig(:)
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:)
            !real(kind=8), intent(out) :: conc_comp(:)
            !real(kind=8), intent(in), optional :: porosity
        end subroutine
        
        function compute_c2nc_tilde_aq_chem(this,mixing_ratios_conc,mixing_waters) result(c2nc_tilde)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
            real(kind=8), intent(in) :: mixing_ratios_conc(:)
            class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
            real(kind=8), allocatable :: c2nc_tilde(:)
        end function
        
        ! function compute_c_tilde(this,mixing_ratios,conc_old) result(c_tilde)
        !     import aqueous_chemistry_c
        !     implicit none
        !     class(aqueous_chemistry_c), intent(in) :: this
        !     !integer(kind=4), intent(in) :: ind
        !     real(kind=8), intent(in) :: mixing_ratios(:)
        !     real(kind=8), intent(in) :: conc_old(:,:) !> concentration of variable activity species before mixing
        !     !class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
        !     real(kind=8), allocatable :: c_tilde(:)
        ! end function
        
        ! function compute_u_tilde(this,c_tilde) result(u_tilde)
        !     import aqueous_chemistry_c
        !     implicit none
        !     class(aqueous_chemistry_c), intent(in) :: this
        !     real(kind=8), intent(in) :: c_tilde(:) !> concentration of "mobile" species after mixing
        !     real(kind=8), allocatable :: u_tilde(:)
        ! end function
            
        subroutine reaction_iteration_EE_eq_kin(this,Delta_t,rk_tilde,conc_comp_react)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
            real(kind=8), intent(out) :: conc_comp_react(:) !> must be allocated
        end subroutine
        
        subroutine reaction_iteration_EE_eq_kin_lump(this,Delta_t,conc_comp_react)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            !real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
            real(kind=8), intent(out) :: conc_comp_react(:) !> must be allocated
        end subroutine
        
        subroutine reaction_iteration_EI_tpt_EE_rk_eq_kin_aq_chem(this,mixing_ratios,rk_mat,porosity,Delta_t,conc_comp_react)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: mixing_ratios(:) !> implicit transport
            real(kind=8), intent(in) :: rk_mat(:,:) !> rk matrix
            real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(out) :: conc_comp_react(:) !> must be allocated
        end subroutine
        
        subroutine reaction_iteration_EE_kin_lump(this,Delta_t,conc_react)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(out) :: conc_react(:) !> must be allocated
        end subroutine

        subroutine reaction_iteration_EE_kin(this,Delta_t,rk_tilde,conc_react)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
            real(kind=8), intent(out) :: conc_react(:) !> must be already allocated
        end subroutine
        
        subroutine compute_dfk_dc1_aq_EfI_lump(this,c2nc,drk_dc,Delta_t,dfk_dc1)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c2nc(:)
            real(kind=8), intent(in) :: drk_dc(:,:)
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(out) :: dfk_dc1(:,:)
        end subroutine
        
        subroutine compute_dfk_dc1_aq_EfI_ideal(this,c2nc,drk_dc,Delta_t,theta,mix_ratio_Rk,dfk_dc1)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c2nc(:)
            real(kind=8), intent(in) :: drk_dc(:,:)
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(in) :: mix_ratio_Rk
            real(kind=8), intent(out) :: dfk_dc1(:,:)
        end subroutine
        
        subroutine compute_dfk_dc1_aq_EfI_ideal_lump(this,c2nc,drk_dc,Delta_t,theta,dfk_dc1)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c2nc(:)
            real(kind=8), intent(in) :: drk_dc(:,:)
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t
            real(kind=8), intent(in) :: theta
            real(kind=8), intent(out) :: dfk_dc1(:,:)
        end subroutine
        
        ! subroutine compute_dfk_dc_aq_EfI(this,drk_dc,porosity,Delta_t,dfk_dc)
        !     import aqueous_chemistry_c
        !     import diag_matrix_c
        !     implicit none
        !     class(aqueous_chemistry_c), intent(in) :: this
        !     real(kind=8), intent(in) :: drk_dc(:,:)
        !     real(kind=8), intent(in) :: porosity
        !     real(kind=8), intent(in) :: Delta_t
        !     real(kind=8), intent(out) :: dfk_dc(:,:)
        ! end subroutine
        
        subroutine Newton_EfI_rk_eq_kin_aq_anal_lump(this,c2nc_ig,c_tilde,Delta_t,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        subroutine Newton_EfI_rk_eq_kin_aq_anal_ideal(this,u_tilde,u_rk_tilde,mix_ratio_Rk,Delta_t,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_tilde(:)
            !real(kind=8), intent(in) :: porosity
            !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing
            real(kind=8), intent(in) :: u_rk_tilde(:) !> reaction part of components concentrations after mixing 
            !real(kind=8), intent(in) :: rk_tilde_new(:) !> new kinetic reaction rates after mixing
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        subroutine Newton_EI_rk_eq_kin_aq_anal_ideal_opt2(this,u_tilde,u_rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_tilde(:) !> component concentrations after mixing
            !real(kind=8), intent(in) :: porosity
            !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing
            real(kind=8), intent(in) :: u_rk_tilde(:) !> reaction part of components concentrations after mixing 
            !real(kind=8), intent(in) :: rk_tilde_new(:) !> new kinetic reaction rates after mixing
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        subroutine Newton_EI_rk_eq_kin_aq_anal_ideal_opt2_lump(this,u_tilde,Delta_t,theta,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_tilde(:) !> component concentrations after mixing
            !real(kind=8), intent(in) :: porosity
            !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing
            !real(kind=8), intent(in) :: u_rk_tilde(:) !> reaction part of components concentrations after mixing 
            !real(kind=8), intent(in) :: rk_tilde_new(:) !> new kinetic reaction rates after mixing
            !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine

        subroutine Newton_EI_rk_kin_aq_anal_ideal_opt2(this,c_tilde,rk_old,Rk_tilde,mix_ratio_Rk,Delta_t,theta,conc_nc,niter,&
            CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: c_tilde(:) !> component concentrations after mixing
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: rk_old(:) !> kinetic reaction rates at previous time step
            real(kind=8), intent(in) :: Rk_tilde(:) !> upstream waters contributions of reaction amounts after mixing
            !real(kind=8), intent(in) :: Rk_tilde_down(:) !> estimation of downstream waters contributions of reaction amounts after mixing
            !real(kind=8), intent(in) :: Rk_tilde_down(:) !> estimation of downstream waters contributions of reaction amounts after mixing
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
            end subroutine
            
            subroutine Newton_EI_rk_kin_aq_anal_ideal_opt2_lump(this,c_tilde,rk_old,Delta_t,theta,conc_nc,niter,&
                CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: c_tilde(:) !> component concentrations after mixing
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: rk_old(:) !> kinetic reaction rates at previous time step
            !real(kind=8), intent(in) :: Rk_tilde(:) !> upstream waters contributions of reaction amounts after mixing
            !real(kind=8), intent(in) :: Rk_tilde_down(:) !> estimation of downstream waters contributions of reaction amounts after mixing
            !real(kind=8), intent(in) :: Rk_tilde_down(:) !> estimation of downstream waters contributions of reaction amounts after mixing
            !real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            !real(kind=8), intent(in) :: y !> sum of upstream mixing ratios
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor
            real(kind=8), intent(inout) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        subroutine Newton_EfI_rk_eq_kin_aq_anal_ideal_opt2(this,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: porosity
            !real(kind=8), intent(in) :: rk_tilde_old(:) !> old kinetic reaction rates after mixing
            real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
            !real(kind=8), intent(in) :: rk_tilde_new(:) !> new kinetic reaction rates after mixing
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine

        
        subroutine Newton_EfI_rk_eq_kin_aq_anal_ideal_lump(this,u_tilde,Delta_t,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: c2nc_ig(:) !> initial guess secondary variable activity concentrations
            real(kind=8), intent(in) :: u_tilde(:)
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(out) :: conc_nc(:) !> variable activity concentrations (already allocated)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        subroutine Newton_EI_tpt_EfI_rk_eq_kin_aq(this,c_tilde,rk_mat,porosity,Delta_t,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
            real(kind=8), intent(in) :: c_tilde(:,:)
            real(kind=8), intent(in) :: rk_mat(:,:) !> rk matrix
            real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t !> time step
            !class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        !  subroutine Newton_EfI_rk_kin_aq_anal(this,c_tilde,porosity,Delta_t,niter,CV_flag)
        !     import aqueous_chemistry_c
        !     import diag_matrix_c
        !     implicit none
        !     class(aqueous_chemistry_c) :: this
        !     real(kind=8), intent(in) :: c_tilde(:)
        !     real(kind=8), intent(in) :: porosity
        !     real(kind=8), intent(in) :: Delta_t !> time step
        !     !class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
        !     integer(kind=4), intent(out) :: niter !> number of iterations
        !     logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        !  end subroutine
         
         subroutine Newton_EfI_rk_kin_aq_anal_ideal(this,c_tilde,rk_tilde,mix_ratio_Rk,Delta_t,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_tilde(:)
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: rk_tilde(:) !> kinetic reaction rates after mixing
            real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio of kinetic reaction rate in this target
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(inout) :: conc_nc(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
         end subroutine
         
         subroutine Newton_EfI_rk_kin_aq_anal_ideal_lump(this,c_tilde,porosity,Delta_t,conc_nc,niter,CV_flag)
            import aqueous_chemistry_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_tilde(:)
            real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(inout) :: conc_nc(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        
        subroutine Picard_EfI_eq_kin_aq_chem(this,mixing_ratios,porosity,Delta_t,mixing_waters,niter,CV_flag)
            import aqueous_chemistry_c
            !import tridiag_matrix_c
            import diag_matrix_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: mixing_ratios(:)
            real(kind=8), intent(in) :: porosity
            real(kind=8), intent(in) :: Delta_t !> time step
            class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
            integer(kind=4), intent(out) :: niter !> number of iterations
            logical, intent(out) :: CV_flag !> FALSE: no CV, TRUE: CV
        end subroutine
        

        
        subroutine compute_re_mean_tilde_aq_chem(this,this_old,mixing_ratios,mixing_waters,Delta_t,re_mean_tilde)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            !integer(kind=4), intent(in) :: ind                              !> index of object "this" in mixing ratios array
            class(aqueous_chemistry_c), intent(in) :: this_old              !> nombre muy cutre
            real(kind=8), intent(in) :: mixing_ratios(:)
            class(aqueous_chemistry_c), intent(in) :: mixing_waters(:)
            real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(in) :: porosity
            real(kind=8), intent(out) :: re_mean_tilde(:)
        end subroutine
        
        subroutine compute_re_mean(this,c2nc_tilde,Delta_t,theta,rk_tilde)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_tilde(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            real(kind=8), intent(in) :: rk_tilde(:)
        end subroutine

        subroutine compute_Re(this,c2nc_tilde,Delta_t,theta)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_tilde(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            !real(kind=8), intent(in) :: rk_tilde(:)
        end subroutine

        subroutine compute_Re_rk(this,c2nc_tilde,Delta_t,theta,rk_tilde)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_tilde(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            real(kind=8), intent(in) :: rk_tilde(:)
        end subroutine
        
        subroutine compute_Re_lump(this,c2nc_tilde,Delta_t,theta)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_tilde(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            !real(kind=8), intent(in) :: porosity
        end subroutine

        subroutine compute_Re_rk_lump(this,c2nc_tilde,Delta_t,theta)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc_tilde(:) !> concentrations secondary variable activity species after mixing at time step k
            real(kind=8), intent(in) :: Delta_t !> time step
            real(kind=8), intent(in) :: theta !> time weighting factor for kinetic reactions
            !real(kind=8), intent(in) :: porosity
        end subroutine
        
        subroutine compute_rk_eq_heterog(this)!,rk)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
        end subroutine
        
        subroutine compute_rk_new(this,rk_new)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(out) :: rk_new(:)
        end subroutine
        
        subroutine compute_Jacobian_rk_anal(this,drk_dc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dc(:,:)
        end subroutine
        
        subroutine compute_rk_Jac_rk_anal(this,rk_new,drk_dc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: rk_new(:) !>  (must be allocated)
            real(kind=8), intent(out) :: drk_dc(:,:)
        end subroutine
        
        subroutine compute_rk_Jac_rk_incr_coeff(this,drk_dc)
            import aqueous_chemistry_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(out) :: drk_dc(:,:)
        end subroutine
        
        function get_indices_reaction(this,reaction) result(indices)
            import aqueous_chemistry_c
            import reaction_c
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            class(reaction_c), intent(in) :: reaction
            integer(kind=4), allocatable :: indices(:)
        end function
        
        subroutine compute_d_log_gamma_d_I_aq_chem(this,d_log_gamma_d_I)
            import aqueous_chemistry_c
            import params_aq_sol_t
            implicit none
            class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: ionic_act
            !class(params_aq_sol_t), intent(in) :: params_aq_sol
            real(kind=8), intent(inout) :: d_log_gamma_d_I(:) !> must be allocated
        end subroutine
        
        subroutine compute_res_Jac_res_anal_exch(this,conc,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,&
            log_Jacobian_act_coeffs,CEC,res,Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            real(kind=8), intent(in) :: CEC
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_Jac_res_anal(this,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,log_Jacobian_act_coeffs,res,&
            Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_Jac_res_anal_ideal(this,indices_icon,n_icon,indices_constrains,ctot,dc2_dc1,res,Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobian secondary-primary
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_Jac_res_incr_coeff(this,c2_init,indices_icon,n_icon,indices_constrains,ctot,res,Jac_res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2_init(:)
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
            real(kind=8), intent(out) :: Jac_res(:,:) !> Jacobian of residual in Newton-Raphson
        end subroutine
        
        subroutine compute_res_init(this,indices_icon,n_icon,indices_constrains,ctot,res)
            import aqueous_chemistry_c
            import int_array_c
            implicit none
            !> Pre-process
            class(aqueous_chemistry_c) :: this
            class(int_array_c), intent(in) :: indices_icon
            integer(kind=4), intent(in) :: n_icon(:) !> number of each icon
            integer(kind=4), intent(in) :: indices_constrains(:,:)
            real(kind=8), intent(in) :: ctot(:)
            real(kind=8), intent(out) :: res(:) !> residual in Newton-Raphson
        end subroutine
        
        subroutine read_wat_type_CHEPROO(this,n_p_aq,num_cstr,num_gas_zones,model,Jac_opt,unit,niter,CV_flag,gas_chem)
            import aqueous_chemistry_c
            import chem_system_c
            import solid_chemistry_c
            import gas_chemistry_c
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: n_p_aq !> number of primary aqueous species
            integer(kind=4), intent(in) :: num_cstr !> number of constrains
            integer(kind=4), intent(in) :: num_gas_zones !> number of gas zones
            integer(kind=4), intent(in) :: model !> activity coefficients model
            integer(kind=4), intent(in) :: Jac_opt !> 0: incremental coeffficinets, 1: analtical
            integer(kind=4), intent(in) :: unit !> file unit
            integer(kind=4), intent(out) :: niter !> number of iterations in Newton method
            logical, intent(out) :: CV_flag !> TRUE if Newton method converges, FALSE otherwise
            !type(solid_chemistry_c), intent(in), optional :: surf_chem
            type(gas_chemistry_c), intent(in), optional :: gas_chem
        end subroutine
        
            !subroutine interfaz_comps_vars(this,u_tilde,Delta_t,u_new)
            !import aqueous_chemistry_c
            !class(aqueous_chemistry_c) :: this
            !real(kind=8), intent(in) :: u_tilde(:) !> concentrations after solving conservative transport
            !real(kind=8), intent(in) :: Delta_t !> time step
            !real(kind=8), intent(out) :: u_new(:) !> concentrations after solving reactive mixing
            !end subroutine
    end interface
    
    contains        
    
        ! subroutine set_ind_wat_aq_chem(this,ind_wat)
        !     implicit none
        !     class(aqueous_chemistry_c) :: this
        !     integer(kind=4), intent(in) :: ind_wat
        !     this%ind_wat=ind_wat
        ! end subroutine
        
        subroutine allocate_conc_aq_species(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            allocate(this%concentrations(this%aq_phase%num_species),this%conc_old(this%aq_phase%num_species),&
                this%conc_old_old(this%aq_phase%num_species))
        end subroutine
        
        subroutine allocate_activities_aq_species(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            allocate(this%activities(this%aq_phase%num_species))
        end subroutine
        
        subroutine compute_activities_diss_solids(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            this%activities(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))*&
                (10**(this%log_act_coeffs(this%indices_aq_species(this%aq_phase%ind_diss_solids))))
        end subroutine
        
        subroutine compute_activities_aq(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            call this%compute_activities_diss_solids()
            call this%compute_act_water()
        end subroutine
        
        subroutine compute_activities(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            call this%compute_activities_diss_solids()
            call this%compute_act_water()
            call this%solid_chemistry%compute_activities_solids()
        end subroutine
        
        subroutine compute_activity(this,ind)
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: ind
            this%activities(ind)=this%concentrations(ind)*(10**this%log_act_coeffs(ind))
        end subroutine
        
        subroutine compute_act_water(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (this%aq_phase%ind_wat>0) then
                this%activities(this%indices_aq_species(this%aq_phase%ind_wat))=1d0-0.018*sum(this%concentrations(&
                    this%indices_aq_species(this%aq_phase%ind_diss_solids))) !>1st order approximation
            end if
        end subroutine
                
        subroutine allocate_log_act_coeffs_aq_chem(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            allocate(this%log_act_coeffs(this%aq_phase%num_species))
            allocate(this%log_Jacobian_act_coeffs(this%aq_phase%num_species,this%aq_phase%num_species))
            this%log_act_coeffs=0d0
            this%log_Jacobian_act_coeffs=0d0
        end subroutine
        
        subroutine set_conc_aq_species(this,conc)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:)
            if (size(conc)/=this%aq_phase%num_species) then
                error stop "Dimension error in concentration of aqueous species"
            else
                this%concentrations=conc
            end if
        end subroutine
        
        subroutine set_log_act_coeffs(this,log_act_coeffs) !> we set or initialise the logarithm of activity coefficientes of  & 
            !! aqueous species
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: log_act_coeffs(:)
            integer(kind=4) :: i
            if (present(log_act_coeffs)) then
                this%log_act_coeffs=log_act_coeffs
            else if (.not. allocated(this%log_act_coeffs)) then
                error stop "Activity coefficients not allocated"
            else
                do i=1,this%aq_phase%num_species
                    this%log_act_coeffs(i)=0d0 !> 
                end do
            end if
        end subroutine
        
        subroutine set_conc_aq_prim_species(this,conc_aq_prim_species)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_aq_prim_species(:)
            
            integer(kind=4) :: i
            if (size(conc_aq_prim_species)/=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
                error stop "Dimension error in concentration of primary species"
            else
                do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                    this%concentrations(this%ind_var_act_species(i))=conc_aq_prim_species(i)
                end do
            end if
        end subroutine
        
        subroutine set_conc_prim_species(this,conc_prim_species)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_prim_species(:)
            if (size(conc_prim_species)/=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) then
                error stop "Dimension error in concentration of primary species"
            end if
            this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1)=conc_prim_species(&
                this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            this%concentrations(this%ind_var_act_species(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))&
                =conc_prim_species(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
        end subroutine
        
        subroutine set_act_diss_solids(this,act_diss_solids)
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: act_diss_solids(:)
        if (size(act_diss_solids)/=this%aq_phase%num_species-this%aq_phase%wat_flag) then
            error stop "Dimension error in activity of dissolved solids"
        end if
        this%activities(this%ind_var_act_species)=act_diss_solids !> chapuza
        end subroutine
        
        subroutine set_conc_var_act_species(this,conc_var_act_species)
        !> variable activity species ordered in:
        !! aqueous primary species
        !! solid primary species
        !! aqueous secondary variable activity species
        !! eq minerals var act
        !! surface complexes
        !! eq gases var act
        
        !! FALTA COMPROBAR CONCENTRACIONES NEGATIVAS
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_var_act_species(:)
            if (size(conc_var_act_species)/=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species) then
                error stop "Dimension error in concentration of variable activity species"
            else
                this%concentrations(this%ind_var_act_species)=&
                    conc_var_act_species(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)
                this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                    this%solid_chemistry%num_solids)=&
                    conc_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl) 
                if (associated(this%gas_chemistry)) then
                    this%gas_chemistry%concentrations(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=&
                        conc_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act)
                    !this%gas_chemistry%concentrations(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+&
                    !    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_cst_act+1:&
                    !    this%solid_chemistry%reactive_zone%gas_phase%num_species)=&
                    !    conc_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                    !    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+1:&
                    !    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)
                end if
            end if
        end subroutine
        
        subroutine set_conc_sec_var_act_species(this,c2nc) !> sets concentration of secondary variable activity species
        !> c2nc must be ordered in:
        !! 1. aqueous secondary variable activity species
        !! 2. surface complexes
        !! 3. variable activity gases in equilibrium
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2nc(:)
            if (size(c2nc)/=this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions) then
                error stop "Dimension error in concentration of secondary variable activity species"
            else if (size(c2nc)==this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species) then
                this%concentrations(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species&
                    +1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))=c2nc
            else
                this%concentrations(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species&
                    +1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))=&
                    c2nc(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species)
                this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                    this%solid_chemistry%num_solids)=&
                c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)
                if (associated(this%gas_chemistry)) then
                    this%gas_chemistry%concentrations(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=&
                        c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)
                end if
            end if
        end subroutine
        
        subroutine set_conc_sec_aq_species(this,c2aq) !> sets concentration of secondary aqueous species
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c2aq(:) !> concentration of secondary aqueous species
            !if (present(c2aq)) then
                if (size(c2aq)/=this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)then
                    error stop "Dimension error in concentration of secondary aqueous species"
                else
                    this%concentrations(this%ind_sec_species)=c2aq
                end if
            !else
            !    this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
            !    this%aq_phase%num_species)=1d-16 !> by default
            !end if
        end subroutine
        
        !subroutine update_conc_sec_aq_species(this,c2aq) !> updates concentration of secondary aqueous species
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    real(kind=8), intent(in) :: c2aq(:) !> concentration of secondary aqueous species
        !    this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:this%aq_phase%num_species)=&
        !    c2aq
        !end subroutine
        
        !subroutine update_conc_sec_var_act_species(this,c2nc) !> updates concentration of secondary variable activity species & 
        !    !! assuming there are aqueous and solid species
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    real(kind=8), intent(in) :: c2nc(:) !> concentration of secondary variable activity species
        !    this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
        !    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)=c2nc(&
        !    1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)
        !    this%solid_chemistry%concentrations(this%solid_chemistry%num_solids-&
        !    this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats+1:this%solid_chemistry%num_solids)=&
        !    c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
        !    this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)
        !end subroutine
        
        !subroutine update_conc_sec_aq_var_act_species(this,c2nc_aq)!> updates concentration of secondary aqueous variable activity &
        !    !! species
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    real(kind=8), intent(in) :: c2nc_aq(:) !> concentration of secondary aqueous species
        !    this%concentrations(this%indices_aq_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
        !        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))=c2nc_aq
        !end subroutine
        
       subroutine set_conc_sec_species(this,c2)
        !> c2 must be ordered in:
        !! 1. aqueous complexes
        !! 2. var act minerals in equilibrium
        !! 3. cst act minerals in equilibrium
        !! 4. surface complexes
        !! 5. var act gases in equilibrium
        !! 6. cst act gases in equilibrium
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            real(kind=8), intent(in) :: c2(:) !> concentration of secondary species

            if (size(c2)/=this%solid_chemistry%reactive_zone%speciation_alg%num_species-&
                this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species) then
                error stop "Dimension error in 'set_conc_sec_species'"
            else if (this%solid_chemistry%reactive_zone%speciation_alg%flag_comp.eqv..false.) then
                this%concentrations(this%ind_sec_species)=c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)
                this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals_kin+1:&
                    this%solid_chemistry%mineral_zone%num_minerals)=c2(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                    this%solid_chemistry%reactive_zone%num_minerals)
                if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                    this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                        this%solid_chemistry%num_solids)=c2(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                        this%solid_chemistry%reactive_zone%num_minerals+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                        this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)
                end if
                if (associated(this%gas_chemistry)) then
                    this%gas_chemistry%concentrations(1:this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act)=c2(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)
                    this%gas_chemistry%concentrations(this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)=c2(&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act)
                end if
            end if
       end subroutine
       
       !subroutine update_conc_sec_species(this,c2)
       !     implicit none
       !     class(aqueous_chemistry_c) :: this
       !     real(kind=8), intent(in) :: c2(:)
       !
       !         this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
       !         this%aq_phase%num_species)=c2(1:this%aq_phase%num_aq_complexes)
       !         this%solid_chemistry%concentrations(1:this%solid_chemistry%reactive_zone%num_minerals)=&
       !         c2(this%aq_phase%num_aq_complexes+1:this%aq_phase%num_aq_complexes+this%solid_chemistry%reactive_zone%num_minerals)
       !         if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
       !             this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+2:&
       !             this%solid_chemistry%reactive_zone%num_solids)=c2(this%aq_phase%num_aq_complexes+&
       !             this%solid_chemistry%reactive_zone%num_minerals+1:this%aq_phase%num_aq_complexes+&
       !             this%solid_chemistry%reactive_zone%num_minerals+this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)
       !         end if
       !         if (associated(this%gas_chemistry)) then
       !             this%gas_chemistry%concentrations(1:this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq)=&
       !             c2(this%gas_chemistry%reactive_zone%speciation_alg%num_eq_reactions-&
       !             this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1:&
       !             this%gas_chemistry%reactive_zone%speciation_alg%num_eq_reactions)
       !         end if
       !end subroutine

       
       function compute_conc_comp(this,c_nc) result(conc_comp)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: c_nc(:) !> variable activity species
            real(kind=8), allocatable :: conc_comp(:) !> component concentrations 
            conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_nc)
        end function
        
        !subroutine compute_conc_comp_aq(this)
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !end subroutine
        
        function compute_conc_comp_cst_act(this,conc) result(conc_comp)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc(:) !> concentration of all species
            real(kind=8), allocatable :: conc_comp(:)
            conc_comp=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act,conc)
        end function
        
      
        
        subroutine set_ionic_act(this,ionic_act)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: ionic_act
            this%ionic_act=ionic_act
        end subroutine
        
        subroutine compute_ionic_act(this) !> computes ionic activity of solution (concentrations are molalities)
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4) :: i
            this%ionic_act=0d0
            do i=1,this%aq_phase%num_species
                this%ionic_act=this%ionic_act+this%concentrations(i)*this%aq_phase%aq_species(i)%valence**2
            end do
            this%ionic_act=0.5*this%ionic_act
        end subroutine
        
        subroutine compute_log_K_aq_chem(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4) :: i
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts
                call this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(i)%compute_logK_dep_T(this%temp)
            end do
            if (associated(this%solid_chemistry)) then
                do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions
                    call this%solid_chemistry%reactive_zone%chem_syst%eq_reacts(&
                        this%solid_chemistry%reactive_zone%ind_eq_reacts(i))%compute_logK_dep_T(this%temp)
                end do
            end if
        end subroutine
       
        subroutine set_solid_chemistry(this,solid_chemistry)
            implicit none
            class(aqueous_chemistry_c) :: this
            type(solid_chemistry_c), intent(in), target :: solid_chemistry
            if (associated(solid_chemistry%reactive_zone)) then
                this%solid_chemistry=>solid_chemistry
            else
                error stop "Solid chemistry object is not associated to a reactive zone"
            end if
        end subroutine
        
        subroutine set_gas_chemistry(this,gas_chemistry)
            implicit none
            class(aqueous_chemistry_c) :: this
            class(gas_chemistry_c), intent(in), target :: gas_chemistry
            if (associated(gas_chemistry%reactive_zone)) then
                this%gas_chemistry=>gas_chemistry
            else
                error stop "gas chemistry object is not associated to a reactive zone"
            end if
        end subroutine
        
        ! subroutine set_indices_aq_species_aq_chem(this,indices_aq_species)
        !     implicit none
        !     class(aqueous_chemistry_c) :: this
        !     integer(kind=4), intent(in), optional :: indices_aq_species(:)
            
        !     integer(kind=4) :: i
            
        !     if (present(indices_aq_species)) then
        !         if (size(indices_aq_species)/=this%aq_phase%num_species) then
        !             error stop "Dimension error in indices of aqueous species"
        !         else
        !             this%indices_aq_species=indices_aq_species
        !         end if
        !     else
        !         call this%allocate_indices_aq_species()
        !         do i=1,this%aq_phase%num_species
        !             this%indices_aq_species(i)=i !> by default
        !         end do
        !     end if
        ! end subroutine
        
        subroutine set_pH(this,pH)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: pH
            if (present(pH)) then
                this%pH=pH
            else
                this%pH=7d0 !> default
            end if
        end subroutine
        
        subroutine compute_pH(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (this%aq_phase%ind_prot/=0) then !> chapuza
                this%pH=-log10(this%activities(this%indices_aq_species(this%aq_phase%ind_prot)))
            end if
        end subroutine
        
        subroutine set_pe(this,pe)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: pe
            if (present(pe)) then
                this%pe=pe
            else
                this%pe=4d0 !> default
            end if
        end subroutine
        
        subroutine set_conc_single_species(this,conc_sp,sp_ind)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_sp !> concentration of species
            integer(kind=4), intent(in) :: sp_ind !> index of species
            if (sp_ind<1 .or. sp_ind>this%aq_phase%num_species) then
                error stop "Error in aqueous species index"
            else if (.not. allocated(this%concentrations)) then
                error stop "Aqueous species concentrations must be allocated"
            else if (conc_sp<0d0) then
                error stop "Concentrations cannot be negative"
            end if
            this%concentrations(sp_ind)=conc_sp
        end subroutine
        
        ! subroutine allocate_ind_diss_solids_aq_chem(this)
        !     class(aqueous_chemistry_c) :: this
        !     if (allocated(this%indices_aq_species(this%aq_phase%ind_diss_solids))) then
        !         deallocate(this%indices_aq_species(this%aq_phase%ind_diss_solids))
        !     end if
        !     allocate(this%indices_aq_species(this%aq_phase%ind_diss_solids)(this%aq_phase%num_species-this%aq_phase%wat_flag))
        ! end subroutine
        
        subroutine allocate_reaction_rates(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (associated(this%solid_chemistry%reactive_zone%chem_syst)) then
                allocate(this%Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%rk_old_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%rk_old_old_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%rk_new(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                allocate(this%Rk_accum(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                    this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts))
                this%Rk=0d0 !> by default
                this%Rk_est=0d0 !> by default
                this%rk_old=0d0 !> by default
                this%rk_new=0d0 !> by default
                this%rk_old_old=0d0 !> by default
                this%rk_old_old_old=0d0 !> by default
                this%rk_mean=0d0 !> by default
                this%Rk_accum=0d0 !> by default
            else
                error stop "Chemical system not associated to reactive zone"
            end if
            allocate(this%re_mean(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            allocate(this%Re(this%solid_chemistry%reactive_zone%chem_syst%num_aq_eq_reacts))
            this%re_mean=0d0 !> by default
            this%Re=0d0 !> by default
            if (associated(this%solid_chemistry)) then
                if (.not. allocated(this%solid_chemistry%re_mean)) then
                    allocate(this%solid_chemistry%re_mean(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    allocate(this%solid_chemistry%Re(this%solid_chemistry%reactive_zone%num_minerals+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                    this%solid_chemistry%re_mean=0d0 !> by default
                    this%solid_chemistry%Re=0d0 !> by default
                end if
                !allocate(this%solid_chemistry%re_mean(this%solid_chemistry%reactive_zone%num_minerals+&
                !this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats))
                !this%solid_chemistry%re_mean=0d0 !> by default
                if (.not. allocated(this%solid_chemistry%Rk)) then
                    allocate(this%solid_chemistry%Rk(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%rk_old(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%rk_old_old(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%rk_old_old_old(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%rk_mean(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%rk_new(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%Rk_est(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    allocate(this%solid_chemistry%Rk_accum(this%solid_chemistry%mineral_zone%num_minerals_kin))
                    this%solid_chemistry%rk_new=0d0 !> by default
                    this%solid_chemistry%Rk_est=0d0 !> by default
                    this%solid_chemistry%rk_old=0d0 !> by default
                    this%solid_chemistry%rk_old_old=0d0 !> by default
                    this%solid_chemistry%rk_old_old_old=0d0 !> by default
                    this%solid_chemistry%rk_mean=0d0 !> by default
                    this%solid_chemistry%Rk=0d0 !> by default
                    this%solid_chemistry%Rk_accum=0d0 !> by default
                end if
            else
                error stop "Solid chemistry not associated to aqueous chemistry"
            end if
            if (associated(this%gas_chemistry)) then
                if (.not. allocated(this%gas_chemistry%re_mean)) then
                    allocate(this%gas_chemistry%re_mean(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    allocate(this%gas_chemistry%Re(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                    this%gas_chemistry%re_mean=0d0 !> by default
                    this%gas_chemistry%Re=0d0 !> by default
                end if
                !allocate(this%gas_chemistry%re_mean(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                !allocate(this%gas_chemistry%Re(this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq))
                !this%gas_chemistry%re_mean=0d0 !> by default
                !this%gas_chemistry%Re=0d0 !> by default
            end if
        end subroutine
        
        subroutine set_aq_phase(this,aq_phase)
            implicit none
            class(aqueous_chemistry_c) :: this
            class(aq_phase_c), intent(in), target :: aq_phase
            this%aq_phase=>aq_phase
        end subroutine
        
        ! subroutine set_ind_diss_solids_aq_chem(this)
        !     implicit none
        !     class(aqueous_chemistry_c) :: this
            
        !     integer(kind=4) :: i,j
            
        !     call this%allocate_ind_diss_solids_aq_chem()
        !     j=1
        !     do i=1,this%aq_phase%num_species
        !         if (this%aq_phase%aq_species(i)%name/='h2o') then
        !             this%ind_diss_solids(j)=i
        !             if (j<this%aq_phase%num_species-this%aq_phase%wat_flag) then
        !                 j=j+1
        !             else
        !                 exit
        !             end if
        !         end if
        !     end do
        ! end subroutine
        
        !subroutine set_speciation_alg_dimensions(this,flag_comp)
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    logical, intent(in) :: flag_comp !> TRUE if component matrix has no constant activity species (De Simoni et al, 2005), FALSE otherwise
        !    
        !    integer(kind=4) :: i,n_sp,n_c,n_eq,n_gas_kin
        !    logical :: flag_cat_exch
        !    
        !    n_gas_kin=0
        !    
        !    if (.not. associated(this%solid_chemistry%reactive_zone%chem_syst)) then
        !        error stop
        !    else if (associated(this%solid_chemistry)) then !> aqueous chemistry is associated to a solid zone
        !        n_sp=this%aq_phase%num_species+this%solid_chemistry%reactive_zone%num_non_flowing_species+this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin
        !        n_c=this%aq_phase%wat_flag
        !        do i=1,this%solid_chemistry%reactive_zone%num_minerals
        !            if (this%solid_chemistry%reactive_zone%minerals(I)%mineral%cst_act_flag==.true.) then
        !                n_c=n_c+1
        !            end if
        !        end do
        !        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin
        !            if (this%solid_chemistry%reactive_zone%chem_syst%minerals(i)%mineral%cst_act_flag==.true.) then
        !                n_c=n_c+1
        !            end if
        !        end do
        !        n_gas_kin=this%solid_chemistry%reactive_zone%gas_phase%num_species-this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq
        !        do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_species
        !            if (this%solid_chemistry%reactive_zone%gas_phase%gases(i)%cst_act_flag==.true.) then
        !                n_c=n_c+1
        !            end if
        !        end do
        !        n_eq=this%solid_chemistry%reactive_zone%num_eq_reactions
        !        if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        !            flag_cat_exch=.true.
        !        else
        !            flag_cat_exch=.false.
        !        end if
        !    else if (associated(this%gas_chemistry)) then !> aqueous chemistry is associated to a gas zone
        !        n_gas_kin=this%gas_chemistry%reactive_zone%gas_phase%num_gases_kin
        !        n_sp=this%aq_phase%num_species+this%gas_chemistry%reactive_zone%num_non_flowing_species+this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin+n_gas_kin
        !        n_c=this%aq_phase%wat_flag
        !        n_eq=this%gas_chemistry%reactive_zone%num_eq_reactions
        !        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin
        !            if (this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(I)%mineral%mineral%cst_act_flag==.true.) then
        !                n_c=n_c+1
        !            end if
        !        end do
        !        do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_species
        !            if (this%gas_chemistry%reactive_zone%gas_phase%gases(i)%cst_act_flag==.true.) then
        !                n_c=n_c+1
        !            end if
        !        end do
        !        if (this%gas_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        !            flag_cat_exch=.true.
        !        else
        !            flag_cat_exch=.false.
        !        end if
        !    else !> all equilibrium reactions are homogeneous
        !        n_sp=this%solid_chemistry%reactive_zone%chem_syst%num_species
        !        n_c=this%solid_chemistry%reactive_zone%chem_syst%num_cst_act_species
        !        n_eq=this%solid_chemistry%reactive_zone%chem_syst%num_eq_reacts
        !        if (this%solid_chemistry%reactive_zone%chem_syst%cat_exch%num_surf_compl>0) then
        !            flag_cat_exch=.true.
        !        else
        !            flag_cat_exch=.false.
        !        end if
        !    end if
        !    call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(flag_comp)
        !    call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_cat_exch)
        !    call this%solid_chemistry%reactive_zone%speciation_alg%set_dimensions(n_sp,n_eq,n_c,this%aq_phase%num_species,this%aq_phase%num_species-this%aq_phase%wat_flag,this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin,n_gas_kin)
        !end subroutine
       
        !subroutine compute_speciation_alg_arrays(this,flag,cols)
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    logical, intent(out) :: flag !> TRUE if 
        !    !class(aq_phase_c), intent(out) :: aq_phase_new
        !    integer(kind=4), intent(out) :: cols(:)
        !    
        !    real(kind=8), allocatable :: Se(:,:),K(:),aux_Se(:,:),aux_Sk(:,:)
        !    integer(kind=4) :: aux_col
        !    !logical :: flag
        !    !type(aq_phase_c), target :: aux_aq_phase
        !                
        !    !call aq_phase_new%copy_attributes(this%aq_phase)
        !    
        !    if (associated(this%solid_chemistry)) then
        !        Se=this%solid_chemistry%reactive_zone%stoich_mat
        !        aux_Se=Se
        !        K=this%solid_chemistry%reactive_zone%get_eq_csts_react_zone()
        !        call this%solid_chemistry%reactive_zone%speciation_alg%compute_arrays(Se,K,this%solid_chemistry%reactive_zone%CV_params%zero,flag,cols)
        !        if (flag==.true.) then
        !            this%solid_chemistry%reactive_zone%stoich_mat(:,cols(1))=aux_Se(:,cols(2))
        !            this%solid_chemistry%reactive_zone%stoich_mat(:,cols(2))=aux_Se(:,cols(1))
        !            aux_Sk=this%solid_chemistry%reactive_zone%chem_syst%Sk
        !            this%solid_chemistry%reactive_zone%chem_syst%Sk(:,cols(1))=aux_Sk(:,cols(2))
        !            this%solid_chemistry%reactive_zone%chem_syst%Sk(:,cols(2))=aux_Sk(:,cols(1))
        !            !aux_col=cols(2)
        !            !this%indices_aq_species(cols(2))=cols(1)
        !            !this%indices_aq_species(cols(1))=aux_col
        !            !call this%aq_phase%aq_species(cols(1))%assign_species(aux_aq_phase%aq_species(cols(2)))
        !            !call this%aq_phase%aq_species(cols(2))%assign_species(aux_aq_phase%aq_species(cols(1)))
        !            !call aq_phase_new%allocate_aq_species(this%aq_phase%num_species)
        !            !call aq_phase_new%aq_species(cols(1))%assign_species(this%aq_phase%aq_species(cols(2)))
        !            !call aq_phase_new%aq_species(cols(2))%assign_species(this%aq_phase%aq_species(cols(1)))
        !            !nullify(this%aq_phase)
        !            !call this%set_aq_phase(aux_aq_phase)
        !        end if
        !    else if (associated(this%solid_chemistry%reactive_zone%chem_syst)) then
        !        Se=this%solid_chemistry%reactive_zone%chem_syst%Se
        !        K=this%solid_chemistry%reactive_zone%chem_syst%get_eq_csts()
        !    else
        !        error stop
        !    end if
            
                
        !end subroutine
        
        
        
        subroutine compute_act_from_MAL(this,ind_act_aq_phase,eq_react,act_non_aq_sp) !> computes activity aqueous species from mass action law
            implicit none
            class(aqueous_chemistry_c) :: this
            integer(kind=4), intent(in) :: ind_act_aq_phase !> index species unknown activity in aqueous phase
            class(eq_reaction_c), intent(in) :: eq_react !> equilibrium reaction
            real(kind=8), intent(in), optional :: act_non_aq_sp !> activity of non-aqueous species (partial pressure if gas)
            
            real(kind=8) :: log_a,dot_prod
            integer(kind=4) :: ind_act_react,i
            integer(kind=8), allocatable :: react_indices(:)
            
            react_indices=this%get_indices_reaction(eq_react)
            
            dot_prod=0d0
            do i=1,size(react_indices)
                if (ind_act_aq_phase==react_indices(i)) then
                    ind_act_react=i
                else
                    dot_prod=dot_prod+eq_react%stoichiometry(i)*log10(this%activities(react_indices(i)))
                end if
            end do
            if (present(act_non_aq_sp)) then
                if (size(react_indices)==eq_react%num_species-1) then
                    dot_prod=dot_prod+log10(act_non_aq_sp)
                else
                    error stop "Dimension error in compute_act_from_MAL"
                end if
            end if
            log_a=(log10(eq_react%eq_cst)-dot_prod)/eq_react%stoichiometry(ind_act_react)
            this%activities(ind_act_aq_phase)=10**log_a
        end subroutine
        
        subroutine compute_molarities(this) !> from molalities
            implicit none
            class(aqueous_chemistry_c) :: this
            
            ! integer(kind=4) :: i
            
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))*this%density*(1d0-this%salinity)
            
            if (this%aq_phase%wat_flag==1) then
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=this%concentrations(&
                    this%indices_aq_species(this%aq_phase%ind_wat))*this%density
            end if

        end subroutine
        
        subroutine compute_salinity(this) !> concentrations are molalities
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8) :: TDS, fluid_mass
            integer(kind=4) :: i
            TDS=0d0 !> [kg_soluto/kg_w]
            do i=1,this%aq_phase%num_species-this%aq_phase%wat_flag
                TDS=TDS+this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids(i))) * &
                    this%aq_phase%aq_species(this%aq_phase%ind_diss_solids(i))%molecular_weight
            end do
            fluid_mass=1d0+TDS
            this%salinity=TDS/fluid_mass
        end subroutine
        
        subroutine compute_alkalinity(this) !> concentrations are molalities
            implicit none
            class(aqueous_chemistry_c), intent(inout) :: this
            real(kind=8) :: alkalinity
            integer(kind=4) :: i
            alkalinity = 0d0
            do i=1, this%aq_phase%num_species
                alkalinity = alkalinity + this%concentrations(this%indices_aq_species(i)) * this%aq_phase%aq_species(i)%alk_contrib
            end do
            this%alkalinity = alkalinity
        end subroutine
        
        subroutine compute_molalities(this) !> from molarities
            implicit none
            class(aqueous_chemistry_c) :: this
            
            ! real(kind=8) :: log_a,dot_prod
            ! integer(kind=4) :: ind_act_react,i
            
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))/(this%density*(1d0-this%salinity))
            this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=&
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))/this%density
        end subroutine
        
        function get_a1(this) result(a1) !> gets primary activities
        implicit none
        class(aqueous_chemistry_c), intent(in) :: this
        real(kind=8), allocatable :: a1(:)

        integer(kind=4) :: i

        allocate(a1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
        do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            a1(i)=this%activities(this%ind_prim_species(i))
        end do
        !> faltan los solidos y gases
        
        !if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        !    a1(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=this%activities(&
        !        this%solid_chemistry%reactive_zone%num_minerals+1)
        !end if
        !if  (associated(this%gas_chemistry)) then
        !    a1(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
        !        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=this%gas_chemistry%activities(&
        !        this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1:this%gas_chemistry%reactive_zone%gas_phase%num_species)
        !end if
        end function
        
        function get_c1(this) result(c1) !> gets primary concentrations
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: c1(:)
            
            integer(kind=4) :: i
            
            allocate(c1(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            
            c1=1d0 !> chapuza
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                c1(i)=this%concentrations(this%ind_var_act_species(i))
            end do
            if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
                c1(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals+1)
            end if
            if  (associated(this%gas_chemistry)) then
                c1(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=this%gas_chemistry%concentrations(&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+1:&
                    this%gas_chemistry%reactive_zone%gas_phase%num_species)/this%volume
            end if
        end function
        
        function get_c1_aq(this) result(c1_aq) !> gets aqueous primary concentrations
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: c1_aq(:)
            
            integer(kind=4) :: i
            
            allocate(c1_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                c1_aq(i)=this%concentrations(this%ind_prim_species(i))
            end do
        end function
        
        function get_c2(this) result(c2) !> gets secondary concentrations
        !! ordered in:
        !! aqueous secondary var act species
        !! aqueous secondary cst act species
        !! var act minerals in equilibrium
        !! cst act minerals in equilibrium
        !! surface complexes
        !! var act gases in eq
        !! cst act gases in equilibrium
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: c2(:)
            integer(kind=4) :: i
            
            allocate(c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species))
        !> Aqueous secondary species (sorted in var act and cst act)
            c2(1:this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species)=this%concentrations(this%ind_sec_species)
            !do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
            !    c2(i)=this%concentrations(this%indices_aq_species(&
            !    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            !end do
        !> Minerals in equilibrium (sorted in var act and cst act)
            do i=1,this%solid_chemistry%reactive_zone%num_minerals
                c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
        !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
                c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                    this%solid_chemistry%reactive_zone%num_minerals+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i)
            end do
        !> Gases in equilibrium
            if (associated(this%gas_chemistry)) then
                !> Gases with var act
                do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                    c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq+i)=this%gas_chemistry%concentrations(&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
                end do
                !> Gases with cst act
                do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act
                    c2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_species-&
                    this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)=this%gas_chemistry%concentrations(i)
                end do
            end if
        end function
        
        function get_c2nc(this) result(c2nc) !> gets secondary variable act concentrations
        !! ordered in:
        !! sec var act aq species
        !! var act minerals in eq
        !! surf complexes
        !! var act gases in eq
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: c2nc(:)
            
            integer(kind=4) :: i
            
            allocate(c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                c2nc(i)=this%concentrations(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
            c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                this%solid_chemistry%reactive_zone%num_minerals_var_act)=this%solid_chemistry%concentrations(&
                this%solid_chemistry%mineral_zone%num_minerals_kin+1:this%solid_chemistry%mineral_zone%num_minerals_kin+&
                this%solid_chemistry%reactive_zone%num_minerals_var_act)
            if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                    this%solid_chemistry%reactive_zone%num_minerals_var_act+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=&
                    this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+2:&
                    this%solid_chemistry%num_solids)
            end if
            if (associated(this%gas_chemistry)) then
                !if (associated(this%solid_chemistry)) then
                    c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%num_solids:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)=this%gas_chemistry%concentrations(&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)
                !else
                !    c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)=this%gas_chemistry%concentrations
                !end if
            end if
        end function
        
        function get_conc_nc(this) result(conc_nc) !> gets conc var act species per unit volume of water
        !> we assume all minerals have cst act
        !> Ordered in:
        !! aqueous primary species
        !! solid primary species
        !! gas primary species
        !! aqueous secondary var act species
        !! surface complexes
        !! gases eq var act
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: conc_nc(:)
            
            integer(kind=4) :: i
            
            allocate(conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            !do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
            !    conc_nc(i)=this%concentrations(this%ind_var_act_species(i))
            !end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                conc_nc(i)=this%concentrations(this%ind_var_act_species(i))
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                !conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=&
                !! this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=this%concentrations(&
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
                if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                        this%solid_chemistry%concentrations(this%solid_chemistry%mineral_zone%num_minerals+1)
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=this%solid_chemistry%concentrations(&
                        this%solid_chemistry%mineral_zone%num_minerals+2:this%solid_chemistry%num_solids)
                end if
            if (associated(this%gas_chemistry)) then
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act)=this%gas_chemistry%concentrations(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)
                !conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+1:&
                !    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)=this%gas_chemistry%concentrations(&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+1:&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_species)
            end if
        end function

        function get_conc_nc_old(this) result(conc_nc) !> gets conc var act species per unit volume of water from previous time step
        !> we assume all minerals have cst act
        !> Ordered in:
        !! aqueous primary species
        !! solid primary species
        !! gas primary species
        !! aqueous secondary var act species
        !! surface complexes
        !! gases eq var act
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: conc_nc(:)
            
            integer(kind=4) :: i
            
            allocate(conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            !do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
            !    conc_nc(i)=this%concentrations(this%ind_var_act_species(i))
            !end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                conc_nc(i)=this%conc_old(this%ind_var_act_species(i))
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                !conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=&
                !! this%conc_old(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=this%conc_old(&
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
                if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
                        this%solid_chemistry%conc_old(this%solid_chemistry%mineral_zone%num_minerals+1)
                    conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+1:&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats)=this%solid_chemistry%conc_old(&
                        this%solid_chemistry%mineral_zone%num_minerals+2:this%solid_chemistry%num_solids)
                end if
            if (associated(this%gas_chemistry)) then
                conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                    this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+1:&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act)=this%gas_chemistry%conc_old(&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+1:&
                    this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq)
                !conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species-&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+1:&
                !    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)=this%gas_chemistry%conc_old(&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_species-&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_gases_kin_var_act+1:&
                !    this%solid_chemistry%reactive_zone%gas_phase%num_species)
            end if
        end function
        
        function get_log_gamma2nc(this) result(log_gamma2nc)!> gets log10 act coeffs sec var act species
        !> We assume minerals have cst act
        !> ordered in:
        !! aqueous sec var act species
        !! surface complexes
        !! gases var act in eq
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma2nc(:)
            integer(kind=4) :: i
            allocate(log_gamma2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                log_gamma2nc(i)=this%log_act_coeffs(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
            !if (associated(this%solid_chemistry)) then
                do i=1,this%solid_chemistry%reactive_zone%cat_Exch_zone%num_exch_cats
                    log_gamma2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+i)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals+1+i)
                end do
            !end if
            if (associated(this%gas_chemistry)) then
                !if (associated(this%solid_chemistry)) then
                !    do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_species
                !        log_gamma2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+this%solid_chemistry%reactive_zone%cat_Exch_zone%num_exch_cats+i)=this%gas_chemistry%log_act_coeffs(i)
                !    end do
                !else
                    do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_gases_eq_var_act
                        log_gamma2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+&
                        this%solid_chemistry%reactive_zone%cat_Exch_zone%num_exch_cats+i)=this%gas_chemistry%log_act_coeffs(&
                        this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq_cst_act+i)
                    end do
                !end if
            end if
        end function
        
        function get_log_gamma2(this) result(log_gamma2) !> gets log_10 secondary activity coefficients
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: log_gamma2(:)
            integer(kind=4) :: i
            allocate(log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                log_gamma2(i)=this%log_act_coeffs(this%ind_sec_species(i))
            end do
            do i=1,this%solid_chemistry%reactive_zone%num_minerals
                log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+i)=&
                    this%solid_chemistry%log_act_coeffs(this%solid_chemistry%mineral_zone%num_minerals_kin+i)
            end do
            do i=1,this%solid_chemistry%reactive_zone%cat_Exch_zone%num_exch_cats
                log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                    this%solid_chemistry%reactive_zone%num_minerals+i)=this%solid_chemistry%log_act_coeffs(&
                    this%solid_chemistry%mineral_zone%num_minerals+1+i)
            end do
            !end if
            if (associated(this%gas_chemistry)) then
                !if (associated(this%solid_chemistry)) then
                    !do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_species
                    !    log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                    !    this%solid_chemistry%reactive_zone%cat_Exch_zone%num_exch_cats+i)=this%gas_chemistry%log_act_coeffs(i)
                    !end do
                !else
                    do i=1,this%solid_chemistry%reactive_zone%gas_phase%num_gases_eq
                        log_gamma2(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species+&
                            this%solid_chemistry%reactive_zone%num_minerals+i)=this%gas_chemistry%log_act_coeffs(i)
                    end do
                !end if
                    
            end if
        end function
        
        function get_conc(this) result(conc)
        !> gets all concentrations of species
        !! We assume minerals have cst act
        !! concentrations are ordered in:
        !! aqueous primary
        !! aqueous sec var act
        !! aqueous cst act
        !! minerals kin
        !! minerals eq
        !! surf complexes
        !! gases eq
        !! gases kin
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), allocatable :: conc(:)
            
            integer(kind=4) :: i
            allocate(conc(this%solid_chemistry%reactive_zone%speciation_alg%num_species))
        !> Aqueous primary species
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                conc(i)=this%concentrations(this%ind_var_act_species(i))
            end do
        !> Aqueous secondary variable activity
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
                conc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+i)=&
                this%concentrations(this%ind_var_act_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))
            end do
        !> Solids
            !> Minerals
            do i=1,this%solid_chemistry%mineral_zone%num_minerals
                conc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+this%aq_phase%wat_flag+i)=&
                    this%solid_chemistry%concentrations(i)
            end do
            !> Surface complexes
            do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl
                conc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+this%aq_phase%wat_flag+&
                    this%solid_chemistry%mineral_zone%num_minerals+i)=this%solid_chemistry%concentrations(&
                    this%solid_chemistry%mineral_zone%num_minerals+i)
            end do
            ! !if (associated(this%solid_chemistry)) then
            ! !> Solid variable activity species
            !     if (this%solid_chemistry%reactive_zone%speciation_alg%flag_cat_exch.eqv..true.) then
            !         conc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)=&
            !         this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+1)
            !     !> Solid secondary variable activity
            !         do i=1,this%solid_chemistry%reactive_zone%cat_exch_zone%num_exch_cats
            !             conc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+&
            !             this%solid_chemistry%reactive_zone%speciation_alg%num_aq_Sec_var_act_species+i)=&
            !             this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+1+i)
            !         end do
            !     end if
            
            !end if
        !> Gases
            if (associated(this%gas_chemistry)) then
                ! if (associated(this%solid_chemistry)) then
                !     do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_species
                !         conc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+&
                !         this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl+i)=this%gas_chemistry%concentrations(i)
                !     end do
                ! else
                    do i=1,this%gas_chemistry%reactive_zone%gas_phase%num_species
                        conc(this%solid_chemistry%reactive_zone%speciation_alg%num_species-&
                            this%gas_chemistry%reactive_zone%gas_phase%num_species+i)=&
                            this%gas_chemistry%concentrations(i)
                    end do
                !end if
            end if
        end function
        
        subroutine check_zero_conc(this,zero_flag,ind_zero) !> checks zero concentrations
        implicit none
        class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry
        logical :: zero_flag !> TRUE if zero concentration, FALSE otherwise
        integer(kind=4), intent(out) :: ind_zero !> index of zero concentration
        
        integer(kind=4) :: i !> counter
        
        zero_flag=.false. !> initialize flag
        ind_zero=0 !> initialize index
        do i=1,this%aq_phase%num_species
            if (this%concentrations(i)<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then
                zero_flag=.true.
                ind_zero=i
            end if
        end do
        end subroutine
        
        subroutine check_conc_aq_var_act_species(this,conc_comp) !> checks concentration aqueous variable activity species from &
            !! components
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_comp(:)
            
            real(kind=8), allocatable :: res(:) !> residual
            
            res=conc_comp-matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,&
            this%concentrations(this%ind_var_act_species(1:&
                this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)))
            if (inf_norm_vec_real(res)>this%solid_chemistry%reactive_zone%CV_params%abs_tol) then
                error stop "Error in aqueous variable activity concentrations"
            end if
        end subroutine
        
        subroutine check_conc_var_act_species(this,conc_nc,conc_comp) !> checks concentration variable activity species from & 
            !! components
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: conc_nc(:) !> concentration variable activity species
            real(kind=8), intent(in) :: conc_comp(:) !> concentration components
            
            real(kind=8), allocatable :: res(:) !> residual

            
            res=conc_comp-matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,conc_nc)
            if (inf_norm_vec_real(res)>this%solid_chemistry%reactive_zone%CV_params%abs_tol) then
                error stop "Error in variable activity concentrations"
            end if
        end subroutine
        
        subroutine check_act_aq_species(this) !> checks activity aqueous species
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            
            real(kind=8), allocatable :: log_res(:) !> log_10(residual)
            
            log_res=log10(this%solid_chemistry%reactive_zone%get_eq_csts_react_zone())-&
            matmul(this%solid_chemistry%reactive_zone%stoich_mat(:,1:this%aq_phase%num_species),log10(this%activities))
            if (inf_norm_vec_real(log_res)>this%solid_chemistry%reactive_zone%CV_params%log_abs_tol) then
                print *, inf_norm_vec_real(log_res)
                error stop "Error in aqueous activities"
            end if
        end subroutine
        
        function get_conc_nc_aq_perm(this,perm) result(c_nc_aq) !> gets concentrations of aqueous species after a permutation
           implicit none
           class(aqueous_chemistry_c), intent(in) :: this
           integer(kind=4), intent(in) :: perm(:)
           real(kind=8), allocatable :: c_nc_aq(:)
           allocate(c_nc_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
           c_nc_aq=this%concentrations(perm)
        end function
        
        !subroutine compute_d_log_gamma_nc_d_log_c2nc_aq(this,d_log_gamma_d_I)
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    real(kind=8), intent(in) :: d_log_gamma_d_I(:)  
        !    !real(kind=8), intent(in) :: c2_aq(:)  
        !    !real(kind=8), intent(out) :: d_log_gamma_d_c2_aq(:,:) 
        !    
        !    integer(kind=4) :: j
        !                            
        !    do j=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
        !        this%log_Jacobian_act_coeffs(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+j)=5d-1*this%aq_phase%z2(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+j)*this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+j)*log(1d1)*d_log_gamma_d_I(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species)
        !    end do
        !end subroutine
        
        subroutine rearrange_state_vars(this,old_aq_phase)
            implicit none
            class(aqueous_chemistry_c) :: this
            class(aq_phase_c), intent(in) :: old_aq_phase
            
            integer(kind=4) :: i,j
            real(kind=8), allocatable :: aux(:,:)
            
            if (this%aq_phase%num_species/=old_aq_phase%num_species) error stop
            !> falta comparar especies acuosas
            allocate(aux(this%aq_phase%num_species,3)) !> chapuza
            aux(:,1)=this%concentrations
            aux(:,2)=this%activities
            aux(:,3)=this%log_act_coeffs
            i=1
            j=1
            do 
                if (this%aq_phase%aq_species(i)%name==old_aq_phase%aq_species(j)%name) then
                    this%concentrations(i)=aux(j,1)
                    this%activities(i)=aux(j,2)
                    this%log_act_coeffs(i)=aux(j,3)
                    if (i<this%aq_phase%num_species) then
                        i=i+1
                        j=1
                    else 
                        exit
                    end if
                else if (j<old_aq_phase%num_species) then
                    j=j+1
                else if (i<this%aq_phase%num_species) then
                    i=i+1
                    j=1
                else 
                    exit
                end if
            end do
        end subroutine
        
        subroutine check_dc2nc_dc1_aq(this,c2nc,dc2nc_dc1,log_Jacobian_act_coeffs)
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c2nc(:)
            real(kind=8), intent(in) :: dc2nc_dc1(:,:)
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            
            real(kind=8), allocatable :: res(:,:),lhs(:,:),rhs(:,:),lhs_1(:,:),lhs_2(:,:)
            type(diag_matrix_c) :: c1_diag,c2nc_inv_diag
            
            call c1_diag%set_diag_matrix(this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            call c2nc_inv_diag%set_diag_matrix(1d0/c2nc)
            lhs_1=c2nc_inv_diag%prod_mat_diag_mat(id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))+log_Jacobian_act_coeffs(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            lhs_2=c1_diag%prod_mat_diag_mat(dc2nc_dc1)
            lhs=matmul(lhs_1,lhs_2)
            rhs=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(&
            1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+&
            id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))-&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            res=lhs-rhs
            if (norm_mat_inf(res)>this%solid_chemistry%reactive_zone%CV_params%log_abs_tol) then
                print *, norm_mat_inf(res)
                error stop "Error in dc2nc_dc1"
            end if
        end subroutine
        
        subroutine check_dc2nc_dc1(this,c1,c2nc,dc2nc_dc1,log_Jacobian_act_coeffs)
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c1(:) !> chapuza
            real(kind=8), intent(in) :: c2nc(:) !> chapuza
            real(kind=8), intent(in) :: dc2nc_dc1(:,:) !> Jacobiano
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            
            real(kind=8), allocatable :: res(:,:),lhs(:,:),rhs(:,:),lhs_1(:,:),lhs_2(:,:)
            type(diag_matrix_c) :: c1_diag,c2_inv_diag
            
            call c1_diag%set_diag_matrix(c1)
            call c2_inv_diag%set_diag_matrix(1d0/c2nc)
            lhs_1=c2_inv_diag%prod_mat_diag_mat(id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))+&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
            lhs_2=c1_diag%prod_mat_diag_mat(dc2nc_dc1)
            lhs=matmul(lhs_1,lhs_2)
            rhs=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_nc_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+id_matrix(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))-log_Jacobian_act_coeffs(&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            res=lhs-rhs
            if (norm_mat_inf(res)>this%solid_chemistry%reactive_zone%CV_params%log_abs_tol) then
                print *, norm_mat_inf(res)
                error stop "Error in dc2nc_dc1"
            end if
        end subroutine
        
        subroutine check_dc2_dc1(this,c1,c2,dc2_dc1,log_Jacobian_act_coeffs)
            implicit none
            class(aqueous_chemistry_c), intent(in) :: this
            real(kind=8), intent(in) :: c1(:) !> chapuza
            real(kind=8), intent(in) :: c2(:) !> chapuza
            real(kind=8), intent(in) :: dc2_dc1(:,:) !> Jacobiano
            real(kind=8), intent(in) :: log_Jacobian_act_coeffs(:,:)
            
            real(kind=8), allocatable :: res(:,:),lhs(:,:),rhs(:,:),lhs_1(:,:),lhs_2(:,:)
            type(diag_matrix_c) :: c1_diag,c2_inv_diag
            
            call c1_diag%set_diag_matrix(c1(1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))
            call c2_inv_diag%set_diag_matrix(1d0/c2)
            lhs_1=c2_inv_diag%prod_mat_diag_mat(id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions)-&
            matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_Jacobian_act_coeffs(1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species))+&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species,&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species))
            lhs_2=c1_diag%prod_mat_diag_mat(dc2_dc1)
            lhs=matmul(lhs_1,lhs_2)
            rhs=matmul(this%solid_chemistry%reactive_zone%speciation_alg%Se_1_star,log_Jacobian_act_coeffs(&
            1:this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)+&
            id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species))-&
            log_Jacobian_act_coeffs(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_species,1:&
            this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
            res=lhs-rhs
            if (norm_mat_inf(res)>this%solid_chemistry%reactive_zone%CV_params%log_abs_tol) then
                print *, norm_mat_inf(res)
                error stop "Error in dc2_dc1"
            end if
        end subroutine
        
        subroutine compute_log_act_coeff_wat(this)
            implicit none
            class(aqueous_chemistry_c) :: this
            if (this%aq_phase%ind_wat>0) then
                this%log_act_coeffs(this%indices_aq_species(this%aq_phase%ind_wat))=&
                    log10(this%activities(this%indices_aq_species(this%aq_phase%ind_wat)))-&
                    log10(this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat)))
            end if
        end subroutine
        
        function compute_saturation_kin_min(this,ind_min) result(saturation)
        !> This function is meant to be used for minerals NOT in equilibrium
            implicit none
            class(aqueous_chemistry_c) :: this !> aqueous chemistry object
            integer(kind=4) :: ind_min !> index of kinetic mineral in chemical system
            real(kind=8) :: saturation
            
            real(kind=8) :: IAP
            integer(kind=4) :: i,ind_kin_react

            IAP=1d0
            ind_kin_react=this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+ind_min !> index of kinetic reaction in chemical system
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(ind_kin_react)%kin_reaction%num_species-1
                !print *, this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(ind_kin_react)%kin_reaction%indices_aq_species(i)
                IAP=IAP*this%activities(this%indices_aq_species(this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(&
                    ind_kin_react)%kin_reaction%indices_aq_phase(i)))**(-&
                    this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(ind_kin_react)%kin_reaction%stoichiometry(i))
            end do
            saturation=IAP*this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(ind_kin_react)%kin_reaction%eq_cst
            !saturation=this%solid_chemistry%reactive_zone%chem_syst%kin_reacts(ind_kin_react)%kin_reaction%eq_cst/IAP
        end function
        
        !subroutine update_conc_nc(this,c_nc)
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    real(kind=8), intent(in) :: c_nc(:)
        !    this%concentrations(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)=&
        !    c_nc(1:this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species)
        !    this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+1)=&
        !    c_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species)
        !    this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        !    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+1)=&
        !    c_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
        !    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+1)
        !    this%solid_chemistry%concentrations(this%solid_chemistry%reactive_zone%num_minerals+2:&
        !    this%solid_chemistry%reactive_zone%num_solids)=&
        !    c_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species+2:&
        !    this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)
        !end subroutine
        
        subroutine compute_c_from_u_bin_syst(this,u,K)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: u
            real(kind=8), intent(in) :: K
            this%concentrations(1)=(u+sqrt(u**2+4*K))/2d0
            this%concentrations(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1)=(-u+sqrt(u**2+4*K))/2d0
        end subroutine
        
        subroutine initialise_conc_nc_iterative_method_bis(this,conc_old_old,conc_old,param)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: conc_old_old(:)
            real(kind=8), intent(in) :: conc_old(:)
            real(kind=8), intent(in) :: param
            
            integer(kind=4) :: i
            
            if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                if ((1d0+param)*conc_old(i)-param*conc_old_old(i)<0d0) then
                    this%concentrations(this%ind_var_act_species(i))=conc_old(i)
                else
                    this%concentrations(this%ind_var_act_species(i))=(1d0+param)*conc_old(i)-param*conc_old_old(i)
                end if
            end do
        end subroutine
        
        !subroutine initialise_c1_iterative_method(this,conc_old_old,conc_old,param)
        !    implicit none
        !    class(aqueous_chemistry_c) :: this
        !    real(kind=8), intent(in) :: conc_old_old(:)
        !    real(kind=8), intent(in) :: conc_old(:)
        !    real(kind=8), intent(in) :: param
        !    
        !    integer(kind=4) :: i
        !    
        !    if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
        !    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
        !        if ((1d0+param)*conc_old(i)-param*conc_old_old(i)<0d0) then
        !            this%concentrations(this%indices_aq_species(i))=conc_old(i)
        !        else
        !            this%concentrations(this%indices_aq_species(i))=(1d0+param)*conc_old(i)-param*conc_old_old(i)
        !        end if
        !    end do
        !end subroutine
        
        subroutine set_act_aq_species(this,act)
            implicit none
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in), optional :: act(:)
            integer(kind=4) :: i
            if (present(act)) then
                do i=1,this%aq_phase%num_species
                    this%activities(i)=act(i)
                end do
            else !> ideal activities by default
                this%activities(this%indices_aq_species(this%aq_phase%ind_diss_solids))=&
                    this%concentrations(this%indices_aq_species(this%aq_phase%ind_diss_solids))
                this%activities(this%indices_aq_species(this%aq_phase%ind_wat))=1d0 !> activity of ideal water
                this%concentrations(this%indices_aq_species(this%aq_phase%ind_wat))=1d0/18d-3 !> concentration of ideal water in mol/m3
            end if
        end subroutine
        
        subroutine allocate_indices_aq_species_aq_chem(this)
        class(aqueous_chemistry_c) :: this
        if (allocated(this%indices_aq_species)) deallocate(this%indices_aq_species)
        allocate(this%indices_aq_species(this%aq_phase%num_species))
        end subroutine
        
        subroutine allocate_indices_aq_phase_aq_chem(this)
        class(aqueous_chemistry_c) :: this
        if (allocated(this%indices_aq_phase)) deallocate(this%indices_aq_phase)
        allocate(this%indices_aq_phase(this%aq_phase%num_species))
        end subroutine
        
        subroutine set_indices_aq_species_aq_chem(this)
            
            class(aqueous_chemistry_c) :: this
            
            integer(kind=4) :: i
            
            call this%allocate_indices_aq_species_aq_chem()
            call this%allocate_indices_aq_phase_aq_chem()
            do i=1,this%aq_phase%num_species
                this%indices_aq_species(i)=i !> by default
                this%indices_aq_phase(i)=i !> by default
                ! if (this%aq_phase%aq_species(i)%name=='h2o') then
                !     this%ind_wat=i
                ! else if (this%aq_phase%aq_species(i)%name=='h+') then
                !     this%ind_prot=i
                ! else if (this%aq_phase%aq_species(i)%name=='oh-') then
                !     this%ind_oh=i
                ! else if (this%aq_phase%aq_species(i)%name=='hco3-') then
                !     this%ind_bicarb=i
                ! else if (index(this%aq_phase%aq_species(i)%name,'co3-2')/=0) then
                ! !else if (this%aq_phase%aq_species(i)%name=='co3-2') then
                !     this%ind_carb=i
                ! else
                !     continue
                ! end if
            end do
        end subroutine
    
        ! subroutine set_ind_prot_aq_chem(this,ind_prot)
        !     implicit none
        !     class(aqueous_chemistry_c) :: this
        !     integer(kind=4), intent(in) :: ind_prot
        !     this%ind_prot=ind_prot
        ! end subroutine

        !> Computes Jacobian of ionic activity with respect to primary cocnentrationa
        subroutine compute_dI_dc1(this,dc2aq_dc1,dI_dc1)
            class(aqueous_chemistry_c) :: this
            real(kind=8), intent(in) :: dc2aq_dc1(:,:)
            real(kind=8), intent(out) :: dI_dc1(:) !> Jacobian of ionic activity with respect to primary cocnentrations & 
            !! (must be already allocated)
            
            integer(kind=4) :: i
            
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
                dI_dc1(i)=5d-1*this%aq_phase%aq_species(i)%valence**2 + &
                dot_product(this%solid_chemistry%reactive_zone%chem_syst%z2(&
                this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:this%aq_phase%num_species),dc2aq_dc1(:,i))
            end do

        end subroutine

!> Computes Jacobian of Newton residual with respect to variable activity species concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
!> We assume the chemical system has only kinetic reactions
!> We apply the lumping technique
subroutine compute_dfk_dc_aq_EfI_lump(this,drk_dc,Delta_t,theta,dfk_dc)
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates with respect to variable activity species concentrations
    !real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(out) :: dfk_dc(:,:) !> Jacobian Newton residual - var act concentrations (must be already allocated)
!> Process
    !> We compute Jacobian Newton residual
    dfk_dc=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)-theta*Delta_t*&
        matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,drk_dc)
end subroutine

!> Computes Jacobian of Newton residual with respect to variable activity species concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
!> We assume the chemical system has only kinetic reactions
subroutine compute_dfk_dc_aq_EfI(this,drk_dc,Delta_t,theta,mix_ratio_Rk,dfk_dc)
    !> Arguments
        class(aqueous_chemistry_c), intent(in) :: this !> aqueous chemistry object
        real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates with respect to variable activity species concentrations
        !real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
        real(kind=8), intent(in) :: Delta_t !> time step
        real(kind=8), intent(in) :: theta !> time weighting factor
        real(kind=8), intent(in) :: mix_ratio_Rk !> mixing ratio for kinetic reaction rates
        real(kind=8), intent(out) :: dfk_dc(:,:) !> Jacobian Newton residual - var act concentrations (must be already allocated)
    !> Process
        !> We compute Jacobian Newton residual
        dfk_dc=id_matrix(this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)-theta*Delta_t*mix_ratio_Rk*&
            matmul(this%solid_chemistry%reactive_zone%U_SkT_prod,drk_dc)
    end subroutine

!> This function computes aqueous component concentrations after mixing (u_tilde)
function compute_u_tilde(this,c_tilde) result(u_tilde)
    implicit none
!> Arguments
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: c_tilde(:) !> concentration of "mobile" species after mixing
    real(kind=8), allocatable :: u_tilde(:) !> component concentrations after mixing
!> Variables
    integer(kind=4) :: i !> index aqueous components
    integer(kind=4) :: j !> index mixing waters
!> Process
    u_tilde=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_tilde)
end function

!> Computes residual of component concentrations
!!> The component matrix INCLUDES constant activity species
!!> residual: \f$ res = U*c - u \f$
subroutine compute_residual_cst_act(this,conc_comp,conc,residual)
    implicit none
    !> Variables
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: conc_comp(:)
    real(kind=8), intent(in) :: conc(:)
    real(kind=8), intent(out) :: residual(:) !> residual of component concentrations
    
    !> Process
    residual=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_cst_act,conc) - conc_comp
end subroutine

!> Computes residual of component concentrations
!!> The component matrix does not include constant activity species
!!> residual: \f$ res = U*c_{nc} - u \f$
subroutine compute_residual(this,conc_comp,c_nc,residual)
    implicit none
    !> Variables
    class(aqueous_chemistry_c), intent(in) :: this
    real(kind=8), intent(in) :: conc_comp(:)
    real(kind=8), intent(in) :: c_nc(:)
    real(kind=8), intent(out) :: residual(:) !> residual of component concentrations
    
    !> Process
    residual=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat,c_nc) - conc_comp
end subroutine

!> Updates concentration aqueous primary species in iterative method
subroutine update_conc_aq_prim_species(this,Delta_c1,zero_flag)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    logical, intent(out) :: zero_flag !> flag to indicate if concentration of primary species is zero
    
    integer(kind=4) :: i,n_p_aq
    real(kind=8), allocatable :: c1_old(:)
    
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    zero_flag=.false. !> default value
    
    if (n_p_aq/=size(Delta_c1) .and. this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species/=size(Delta_c1)) then
        error stop "Dimension error in update_conc_aq_prim_species"
    end if
    c1_old=this%concentrations(this%ind_var_act_species(1:n_p_aq))
    do i=1,n_p_aq
        if (this%concentrations(this%ind_var_act_species(i))+Delta_c1(i)<=&
        this%solid_chemistry%reactive_zone%CV_params%control_factor*&
        this%concentrations(this%ind_var_act_species(i))) then
            this%concentrations(this%ind_var_act_species(i))=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            this%concentrations(this%ind_var_act_species(i))
        else if (this%concentrations(this%ind_var_act_species(i))+Delta_c1(i)>=this%concentrations(this%ind_var_act_species(i))/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(this%ind_var_act_species(i))=this%concentrations(this%ind_var_act_species(i))/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(this%ind_var_act_species(i))=this%concentrations(this%ind_var_act_species(i))+Delta_c1(i)
        end if
        if (this%concentrations(this%ind_var_act_species(i))<this%solid_chemistry%reactive_zone%CV_params%abs_tol**3) then
            this%concentrations(this%ind_var_act_species(i))=this%solid_chemistry%reactive_zone%CV_params%abs_tol**3 !> chapuza
            print *, "Warning: concentration of primary species is zero"
            print *, "Primary species: ", this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(i)))%name
            !error stop
            zero_flag=.true. !> flag to indicate that concentration of primary species is zero
        end if
        Delta_c1(i)=this%concentrations(this%ind_var_act_species(i))-c1_old(i)
    end do
end subroutine

!> Updates concentration variable activity species in iterative method
subroutine update_conc_nc(this,conc_nc,Delta_c,zero_flag,sp_ind)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: conc_nc(:) !> concentration of variable activity species
    real(kind=8), intent(inout) :: Delta_c(:) !> difference in concentration of variable activity species
    logical, intent(out) :: zero_flag !> flag to indicate if concentration of variable activity species is zero
    integer(kind=4), intent(out) :: sp_ind !> index of species with zero concentration
    
    integer(kind=4) :: i,n_nc,n_p_aq
    real(kind=8), allocatable :: c_old(:)
    
    n_nc=this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species !> number of variable activity species
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species !> number of aqueus primary species
    zero_flag=.false. !> default value
    
    if (n_nc/=size(Delta_c) .or. n_nc/=size(conc_nc)) then
        error stop "Dimension error in update_conc_nc"
    end if
    c_old=conc_nc
    do i=1,n_nc
        if (conc_nc(i)+Delta_c(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*conc_nc(i)) then
            conc_nc(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*conc_nc(i)
        else if (conc_nc(i)+Delta_c(i)>=conc_nc(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            conc_nc(i)=conc_nc(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            conc_nc(i)=conc_nc(i)+Delta_c(i)
        end if
        if (i<=n_p_aq .and. conc_nc(i)<this%solid_chemistry%reactive_zone%CV_params%abs_tol**2) then
            conc_nc(i)=this%solid_chemistry%reactive_zone%CV_params%abs_tol**2 !> chapuza
            print *, "Warning: concentration of primary species is zero"
            print *, "Primary species: ", this%aq_phase%aq_species(this%indices_aq_phase(this%ind_var_act_species(i)))%name
            zero_flag=.true. !> flag to indicate that concentration of variable activity species is zero
            sp_ind=i !> index of species with zero concentration
        end if
        Delta_c(i)=conc_nc(i)-c_old(i)
    end do
    call this%set_conc_var_act_species(conc_nc)
end subroutine

subroutine update_conc_aq_prim_species_bis(this,Delta_c1)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    
    integer(kind=4) :: i,n_p_aq
    real(kind=8), allocatable :: c1_old(:)
    
    n_p_aq=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    
    if (n_p_aq/=size(Delta_c1) .and. this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species/=size(Delta_c1)) then
        error stop "Dimension error in update_conc_aq_prim_species"
    end if
    c1_old=this%concentrations(1:n_p_aq)
    do i=1,n_p_aq
        if (this%concentrations(i)+Delta_c1(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            this%concentrations(i)) then
            this%concentrations(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
            this%concentrations(i)
        else if (this%concentrations(i)+Delta_c1(i)>=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(i)=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c1(i)
        end if
        Delta_c1(i)=this%concentrations(i)-c1_old(i)
    end do
end subroutine


!> Updates concentration aqueous species in Newton method
subroutine update_conc_aq_species(this,Delta_c_aq)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: Delta_c_aq(:)
    
    integer(kind=4) :: i
    real(kind=8), allocatable :: conc_old(:)
    
    if (this%solid_chemistry%reactive_zone%CV_params%control_factor>1d0 .or. &
    this%solid_chemistry%reactive_zone%CV_params%control_factor<0d0) error stop "Control factor must be in (0,1)"
    conc_old=this%concentrations
    do i=1,this%aq_phase%num_species
        if (this%concentrations(i)+Delta_c_aq(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*&
        this%concentrations(i)) then
            this%concentrations(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*this%concentrations(i)
        else if (this%concentrations(i)+Delta_c_aq(i)>=this%concentrations(i)/&
            this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            this%concentrations(i)=this%concentrations(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            this%concentrations(i)=this%concentrations(i)+Delta_c_aq(i)
        end if
        Delta_c_aq(i)=this%concentrations(i)-conc_old(i)
    end do
end subroutine update_conc_aq_species

!> Updates concentration aqueous and solid primary species in iterative method
subroutine update_conc_prim_species(this,c1,Delta_c1)
    class(aqueous_chemistry_c) :: this
    real(kind=8), intent(inout) :: c1(:)
    real(kind=8), intent(inout) :: Delta_c1(:) !> must be already allocated
    
    integer(kind=4) :: i,n_p_aq,n_p
    real(kind=8), allocatable :: c1_old(:)
    
    n_p=this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species
    
    if (n_p/=size(c1)) error stop "Dimension error in update_conc_prim_species"
        
    c1_old=c1 !> old primary concentrations
    do i=1,n_p
        if (c1_old(i)+Delta_c1(i)<=this%solid_chemistry%reactive_zone%CV_params%control_factor*c1_old(i)) then
            c1(i)=this%solid_chemistry%reactive_zone%CV_params%control_factor*c1_old(i)
        else if (c1_old(i)+Delta_c1(i)>=c1_old(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor) then
            c1(i)=c1_old(i)/this%solid_chemistry%reactive_zone%CV_params%control_factor
        else
            c1(i)=c1_old(i)+Delta_c1(i)
        end if
        Delta_c1(i)=c1(i)-c1_old(i)
    end do
    call this%set_conc_prim_species(c1)
end subroutine

! !> Computes Jacobian of Newton residual with respect to aqueous species concentrations in reactive mixing iteration using Euler fully implicit in chemical reactions
! !> We assume the chemical system has only kinetic reactions
! !> We assume all species are aqueous
! subroutine compute_dfk_dc_aq_EfI(this,drk_dc,porosity,Delta_t,dfk_dc)
! !> Arguments
!     class(aqueous_chemistry_c), intent(in) :: this
!     real(kind=8), intent(in) :: drk_dc(:,:) !> Jacobian of kinetic reaction rates
!     real(kind=8), intent(in) :: porosity !> in solid chemistry associated to this aqueous chemistry
!     real(kind=8), intent(in) :: Delta_t !> time step
!     real(kind=8), intent(out) :: dfk_dc(:,:) !> Jacobian Newton residual - aqueous concentrations (must be already allocated)
! !> Process
!     !> We compute Jacobian Newton residual
!         dfk_dc=id_matrix(this%aq_phase%num_species)-(Delta_t/porosity)*&
!             matmul(transpose(this%solid_chemistry%reactive_zone%chem_syst%Sk),drk_dc)
! end subroutine

!> This subroutine computes the concentration of aqueous variable activity species after mixing
!function compute_c_tilde(this,mixing_ratios,conc_old) result(c_tilde)
!    implicit none
!!> Argumnets
!    class(aqueous_chemistry_c), intent(in) :: this
!    real(kind=8), intent(in) :: mixing_ratios(:)  !> the first element corresponds to argument "this", and the rest correspond to argument "conc_old" (in the same order)
!    real(kind=8), intent(in) :: conc_old(:,:) !> concentration of aqueous variable activity species before mixing
!    real(kind=8), allocatable :: c_tilde(:) !> concentration of aqueous variable activity species after mixing
!!> Variables
!    integer(kind=4) :: i !> index aqueous variable activity species
!    integer(kind=4) :: j !> index mixing waters
!!> Pre-process
!    allocate(c_tilde(size(conc_old,1)))
!!> Process
!    do i=1,size(conc_old,1)
!        c_tilde(i)=mixing_ratios(1)*conc_old(i,1)
!        do j=1,size(conc_old,2)-1
!            c_tilde(i)=c_tilde(i)+mixing_ratios(j+1)*conc_old(i,j+1)
!        end do
!    end do
!end function 

!> This subroutine computes the concentration of aqueous variable activity species
!! we assume all mixing waters have the same aqueous species
subroutine compute_c_tilde(mix_waters,mixing_ratios_conc,c_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
    !real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> same order as "mix_waters"
    real(kind=8), intent(out) :: c_tilde(:) !> concentration of aqueous variable activity species after mixing (must be already allocated)
    !real(kind=8), intent(out) :: rk_tilde(:) !> aqueous kinetic reaction rates after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: c_tilde_aux(:) !> auxiliary variable to store the concentration of aqueous variable activity species after mixing
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    logical :: flag !> flag to check if an index is in "perm"
!> Pre-process
    !allocate(c_tilde(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
!> Process
    c_tilde=0d0 !> initialise concentration of aqueous variable activity species after mixing
    !rk_tilde=0d0 !> initialise aqueous kinetic reaction rates after mixing
    num_perm=0 !< counter in "perm"
    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
        !> Autentica chapuza
        if (mix_waters(1)%ind_var_act_species(i)/=mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i))) then
            if (num_perm>0) then
                call is_int_in_1D_array(i,perm(1:num_perm,2),flag)
                if (flag .eqv. .true.) then
                    continue
                else
                    num_perm=num_perm+1
                    perm(num_perm,1)=i
                    perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
                        mix_waters(1)%ind_var_act_species)
                end if
            else
                num_perm=num_perm+1
                perm(num_perm,1)=i
                perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
                    mix_waters(1)%ind_var_act_species)
            end if
        end if
        do j=1,num_mix_wat
            !print *, mix_waters(j)%indices_aq_species(mix_waters(j)%aq_phase%ind_diss_solids(i))
            c_tilde(i)=c_tilde(i)+mixing_ratios_conc(j)*mix_waters(j)%concentrations(mix_waters(j)%indices_aq_species(&
                mix_waters(j)%aq_phase%ind_diss_solids(i)))
        end do
    end do
    !do i=1,mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
    !    do j=1,num_mix_wat
    !        rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%rk(i)
    !    end do
    !end do
    c_tilde_aux=c_tilde !> we store the concentration of aqueous variable activity species after mixing (chapuza)
    !rk_tilde_aux=rk_tilde !> we store the aqueous kinetic reaction rates after mixing (chapuza)
    do i=1,num_perm
        !> we permute the indices in "perm"
        c_tilde(perm(i,1))=c_tilde_aux(perm(i,2))
        c_tilde(perm(i,2))=c_tilde_aux(perm(i,1))
        !rk_tilde(perm(i,1))=rk_tilde_aux(perm(i,2))
        !rk_tilde(perm(i,2))=rk_tilde_aux(perm(i,1))
    end do
end subroutine

!> This subroutine computes the concentration of aqueous variable activity species & aqueous kinetic reaction rates after water mixing
!! We assume all kinetic reactions are aqueous
!! we assume all mixing waters have the same kinetic reactions
subroutine compute_c_and_rk_tilde(mix_waters,mixing_ratios_conc,mixing_ratios_Rk,c_tilde,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
    real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
    real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> same order as "mix_waters"
    real(kind=8), intent(out) :: c_tilde(:) !> concentration of aqueous variable activity species after mixing
    real(kind=8), intent(out) :: rk_tilde(:) !> aqueous kinetic reaction rates after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: c_tilde_aux(:) !> auxiliary variable to store the concentration of aqueous variable activity species after mixing
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    logical :: flag !> flag to check if an index is in "perm"
!> Pre-process
    !allocate(c_tilde(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
!> Process
    c_tilde=0d0 !> initialise concentration of aqueous variable activity species after mixing
    rk_tilde=0d0 !> initialise aqueous kinetic reaction rates after mixing
    num_perm=0 !< counter in "perm"
    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
        !> Autentica chapuza
        if (mix_waters(1)%ind_var_act_species(i)/=mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i))) then
            if (num_perm>0) then
                call is_int_in_1D_array(i,perm(1:num_perm,2),flag)
                if (flag .eqv. .true.) then
                    continue
                else
                    num_perm=num_perm+1
                    perm(num_perm,1)=i
                    perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
                        mix_waters(1)%ind_var_act_species)
                end if
            else
                num_perm=num_perm+1
                perm(num_perm,1)=i
                perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
                    mix_waters(1)%ind_var_act_species)
            end if
        end if
        do j=1,num_mix_wat
            !print *, mix_waters(j)%indices_aq_species(mix_waters(j)%aq_phase%ind_diss_solids(i))
            c_tilde(i)=c_tilde(i)+mixing_ratios_conc(j)*mix_waters(j)%concentrations(mix_waters(j)%indices_aq_species(&
                mix_waters(j)%aq_phase%ind_diss_solids(i)))
            
        end do
    end do
    do i=1,mix_waters(1)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        do j=1,num_mix_wat
            rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%rk(i)
        end do
    end do
    !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
    !j=1
        !do !j=1,num_mix_wat
            j=1
            k=1
            i=1
            do
            !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
                if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
                    mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
                    rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%solid_chemistry%rk(i)
                    if (j<num_mix_wat) then
                        j=j+1
                        i=1
                        k=1
                    else
                        exit
                    end if
                else if (k<mix_waters(j)%solid_chemistry%mineral_zone%num_minerals_kin) then
                    k=k+1
                else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
                    i=i+1
                    k=1
                else if (j<num_mix_wat) then
                    j=j+1
                    i=1
                    k=1
                else
                    exit
                end if
            end do
            !rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%solid_chemistry%rk(i)
        !end do
    !end do
    end if
    c_tilde_aux=c_tilde !> we store the concentration of aqueous variable activity species after mixing (chapuza)
    !rk_tilde_aux=rk_tilde !> we store the aqueous kinetic reaction rates after mixing (chapuza)
    do i=1,num_perm
        !> we permute the indices in "perm"
        c_tilde(perm(i,1))=c_tilde_aux(perm(i,2))
        c_tilde(perm(i,2))=c_tilde_aux(perm(i,1))
        !rk_tilde(perm(i,1))=rk_tilde_aux(perm(i,2))
        !rk_tilde(perm(i,2))=rk_tilde_aux(perm(i,1))
    end do
end subroutine 

!> This subroutine computes the kinetic reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters, for the first time step
!> The estimation of the downstream waters reaction rates is done by taking the reaction rate of the closest upstream water
subroutine compute_rk_tilde_impl_opt1(mix_waters,mixing_ratios_Rk,num_up,num_down,theta,Delta_t,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(inout) :: mix_waters(:) !> waters that mix with each other
    !real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
    real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    !integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    !integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: Sk_nc_old(:,:) !> kinetic stoichiometric matrix associated to the previous target
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    !logical :: flag !> flag to check if an index is in "perm"
    real(kind=8), parameter :: eps=1d-16 !> machine epsilon
    
    !allocate(c_tilde(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    !allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_rk_tilde_impl_opt1"
    
!> Process
    !c_tilde=0d0 !> initialise concentration of aqueous variable activity species after mixing
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    !rk_tilde_down=0d0 !> initialise kinetic reaction contributions after mixing
    rk_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !num_perm=0 !< counter in "perm"
    !do i=1,mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
    !    !> Autentica chapuza
    !    if (mix_waters(1)%ind_var_act_species(i)/=mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i))) then
    !        if (num_perm>0) then
    !            call is_int_in_1D_array(i,perm(1:num_perm,2),flag)
    !            if (flag .eqv. .true.) then
    !                continue
    !            else
    !                num_perm=num_perm+1
    !                perm(num_perm,1)=i
    !                perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
    !                    mix_waters(1)%ind_var_act_species)
    !            end if
    !        else
    !            num_perm=num_perm+1
    !            perm(num_perm,1)=i
    !            perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
    !                mix_waters(1)%ind_var_act_species)
    !        end if
    !    end if
    !    do j=1,num_mix_wat
    !        !print *, mix_waters(j)%indices_aq_species(mix_waters(j)%aq_phase%ind_diss_solids(i))
    !        c_tilde(i)=c_tilde(i)+mixing_ratios_conc(j)*mix_waters(j)%concentrations(mix_waters(j)%indices_aq_species(&
    !            mix_waters(j)%aq_phase%ind_diss_solids(i)))
    !        
    !    end do
    !end do
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%rk+(1d0-theta)*&
                        mix_waters(j)%rk_old)
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%solid_chemistry%rk+&
                        (1d0-theta)*mix_waters(j)%solid_chemistry%rk_old)
                end if
                Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc)
            end if
        end do
        !Sk_nc_old=0d0 !> chapuza
        !> Downstream waters
        do j=1,num_down
            if (mix_waters(num_up+j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(num_up+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                !if (Sk_nc==Sk_nc_old) then
                !    rk_tilde=rk_tilde+
                !end if
                if (mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    if (norm_mat_inf(Sk_nc_old-Sk_nc)<eps) then !> we check if the kinetic stoichiometric matrix is the same
                        mix_waters(num_up+j)%Rk_est=Delta_t*(theta*mix_waters(num_up)%rk+&
                            (1d0-theta)*mix_waters(num_up)%rk_old)
                    else
                        mix_waters(num_up+j)%Rk_est=mix_waters(num_up+j)%rk_old*Delta_t
                    end if
                    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),mix_waters(num_up+j)%Rk_est/Delta_t)
                end if
                if (allocated(mix_waters(num_up+j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    if (norm_mat_inf(Sk_nc_old-Sk_nc)<eps) then
                        mix_waters(num_up+j)%solid_chemistry%Rk_est=Delta_t*(&
                            theta*mix_waters(num_up)%solid_chemistry%rk+(1d0-theta)*mix_waters(num_up)%solid_chemistry%rk_old)
                    else
                        mix_waters(num_up+j)%solid_chemistry%Rk_est=Delta_t*mix_waters(num_up+j)%solid_chemistry%rk_old
                    end if
                    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),&
                        mix_waters(num_up+j)%solid_chemistry%Rk_est/Delta_t)
                end if
                ! if (size(Sk_nc_old,1) /= size(Sk_nc,1)) then !> chapuza
                !     deallocate(Sk_nc_old)
                !     allocate(Sk_nc_old(size(Sk_nc,1),size(Sk_nc,2)))
                ! end if
                !Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc)
            end if
        end do
    !    do j=1,num_down
    !        Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
    !        rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),mix_waters(num_up+j)%rk_old)
    !        deallocate(Sk_nc)
    !    end do
    !end if
    !!do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    !if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
    !!j=1
    !    !do !j=1,num_mix_wat
    !    !> Upstream waters
    !    !j=1 !> index upstream waters
    !    !k=1 !> index mineral kinetic reactions
    !    !i=1 !> index mineral kinetic reactions
    !    do j=1,num_up !> loop upstream waters
    !        Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
    !        !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    !            !if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
    !            !    mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
    !                rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%solid_chemistry%rk)
    !                !if (j<num_up) then
    !                !    j=j+1
    !                !    i=1
    !                !    k=1
    !                !else
    !                !    exit
    !                !end if
    !            !else if (k<mix_waters(j)%solid_chemistry%mineral_zone%num_minerals_kin) then
    !            !    k=k+1
    !            !else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
    !            !    i=i+1
    !            !    k=1
    !            !else if (j<num_up) then
    !            !    j=j+1
    !            !    i=1
    !            !    k=1
    !            !else
    !            !    exit
    !            !end if
    !                deallocate(Sk_nc)
    !        end do
    !        !> Downstream waters
    !        !j=1 !> index downstream waters
    !        !k=1 !> index mineral kinetic reactions
    !        !i=1 !> index mineral kinetic reactions
    !        do j=1,num_down !> loop downstream waters
    !            Sk_nc=mix_waters(num_up+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th downstream water (chapuza)
    !        !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    !            !if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
    !            !    mix_waters(num_up+j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
    !                rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),&
    !                    mix_waters(num_up+j)%solid_chemistry%rk_old)
    !            !    if (j<num_down) then
    !            !        j=j+1
    !            !        i=1
    !            !        k=1
    !            !    else
    !            !        exit
    !            !    end if
    !            !else if (k<mix_waters(num_up+j)%solid_chemistry%mineral_zone%num_minerals_kin) then
    !            !    k=k+1
    !            !else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
    !            !    i=i+1
    !            !    k=1
    !            !else if (j<num_down) then
    !            !    j=j+1
    !            !    i=1
    !            !    k=1
    !            !else
    !            !    exit
    !            !end if
    !                deallocate(Sk_nc)
    !        end do
            !rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%solid_chemistry%rk(i)
        !end do
    !end do
    !end if
    !c_tilde_aux=c_tilde !> we store the concentration of aqueous variable activity species after mixing (chapuza)
    !!rk_tilde_aux=rk_tilde !> we store the aqueous kinetic reaction rates after mixing (chapuza)
    !do i=1,num_perm
    !    !> we permute the indices in "perm"
    !    c_tilde(perm(i,1))=c_tilde_aux(perm(i,2))
    !    c_tilde(perm(i,2))=c_tilde_aux(perm(i,1))
    !    !rk_tilde(perm(i,1))=rk_tilde_aux(perm(i,2))
    !    !rk_tilde(perm(i,2))=rk_tilde_aux(perm(i,1))
    !end do
    !rk_tilde=rk_tilde*Delta_t
end subroutine 

!> This subroutine computes the kinetic reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters
!> The estimation of the downstream waters reaction rates is done by extrapolation with the two previous reaction rates
subroutine compute_rk_tilde_impl_opt2(mix_waters,mixing_ratios_Rk,num_up,num_down,theta,Delta_t,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(inout) :: mix_waters(:) !> waters that mix with each other
    real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    !real(kind=8), allocatable :: Sk_nc_old(:,:) !> kinetic stoichiometric matrix associated to the previous target
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), parameter :: eps=1d-16 !> machine epsilon
    !real(kind=8) :: mu=1d0 !> parameter for extrapolation of downstream waters rk
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_rk_tilde_impl_opt2"
    
!> Process
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    !rk_tilde_down=0d0 !> initialise kinetic reaction contributions after mixing
    rk_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%rk+&
                        (1d0-theta)*mix_waters(j)%rk_old)
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%solid_chemistry%rk+&
                        (1d0-theta)*mix_waters(j)%solid_chemistry%rk_old)
                end if
                deallocate(Sk_nc)
            end if
        end do
    !> Downstream waters
        do j=1,num_down
            if (mix_waters(num_up+j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(num_up+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                if (mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    !if (allocated(Sk_nc_old) .and. norm_mat_inf(Sk_nc_old-Sk_nc)<eps) then
                    !    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),mix_waters(num_up)%rk)
                    !else
                    mix_waters(num_up+j)%Rk_est=Delta_t*((1d0+&
                            mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%mu*theta)*mix_waters(num_up+j)%rk_old-&
                            mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%mu*theta*mix_waters(num_up+j)%rk_old_old)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),mix_waters(num_up+j)%Rk_est/Delta_t)
                    !end if
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(num_up+j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    mix_waters(num_up+j)%solid_chemistry%Rk_est=Delta_t*(&
                        (1d0+mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%mu*theta)*&
                        mix_waters(num_up+j)%solid_chemistry%rk_old-&
                        mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%mu*theta*&
                        mix_waters(num_up+j)%solid_chemistry%rk_old_old)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),&
                        mix_waters(num_up+j)%solid_chemistry%Rk_est/Delta_t)
                end if
                !Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc)
            end if
        end do
end subroutine 

!> This subroutine computes the kinetic reaction rate contributions after water mixing, implicitly for upstream waters and explicitly for downstream waters
!> The estimation of the downstream waters reaction rates is done by assuming there is no mixing
subroutine compute_rk_tilde_impl_opt4(mix_waters,mixing_ratios_Rk,num_up,num_down,theta,Delta_t,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(inout) :: mix_waters(:) !> waters that mix with each other
    real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    real(kind=8), intent(out) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    real(kind=8), allocatable :: Rk_est(:) !> mean kinetic reaction amounts
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: conc_nc(:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: conc_nc_old(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    !real(kind=8) :: mu=0d0 !> parameter for extrapolation of downstream waters rk
    
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_rk_tilde_impl_opt3"
    
!> Process
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    !rk_tilde_down=0d0 !> initialise kinetic reaction contributions after mixing
    rk_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%rk+&
                        (1d0-theta)*mix_waters(j)%rk_old)
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%solid_chemistry%rk+&
                        (1d0-theta)*mix_waters(j)%solid_chemistry%rk_old)
                end if
                deallocate(Sk_nc)
            end if
        end do
    !> Downstream waters (we estimate them)
        do j=1,num_down
            if (mix_waters(num_up+j)%indices_rk%num_cols>0) then
                allocate(Rk_est(mix_waters(num_up+j)%indices_rk%num_cols))
                allocate(conc_nc(mix_waters(num_up+j)%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species))
                Sk_nc=mix_waters(num_up+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                call mix_waters(num_up+j)%initialise_conc_nc_iterative_method_ideal(&
                    mix_waters(num_up+j)%solid_chemistry%reactive_zone%CV_params%mu,conc_nc) !> we estimate the concentration of aqueous variable activity species in j-th downstream water (chapuza)
                !conc_nc=mix_waters(num_up+j)%get_conc_nc()
                conc_nc_old=mix_waters(num_up+j)%get_conc_nc_old() !> we get the old concentration of aqueous variable activity species in j-th downstream water (chapuza)
                call LU_lin_syst(matmul(Sk_nc,transpose(Sk_nc)),matmul(Sk_nc,conc_nc-conc_nc_old),tol,Rk_est) !> we compute the estimated mean reaction amount with ordianry least squares in j-th downstream water
                call mix_waters(num_up+j)%set_Rk_est(Rk_est) !> we set the estimated reaction rate in j-th downstream water
                rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),Rk_est/Delta_t) !> we update the reaction rate contributiuons
                !call mix_waters(num_up+j)%set_conc_var_act_species(conc_nc_old) !> we set the estimated reaction rate in j-th downstream water
                !     mix_waters(num_up+j)%rk(i)=(Rk_mean(i)-Delta_t*(1d0-theta)*mix_waters(num_up+j)%rk_old(i))/(Delta_t*theta)
                ! end do
                ! do i=1,mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
                !     mix_waters(num_up+j)%rk(mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                !         (Rk_mean(mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                !         mix_waters(num_up+j)%solid_chemistry%mineral_zone%num_minerals_kin+i)-Delta_t*(1d0-theta)*mix_waters(num_up+j)%rk_old(&
                !         mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i))/(Delta_t*theta)
                ! end do
                !if (allocated(mix_waters(num_up+j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    !rk_tilde_down=rk_tilde_down+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),Rk_mean/Delta_t)
                    ! do i=1,mix_waters(num_up+j)%solid_chemistry%mineral_zone%num_minerals_kin
                    !     mix_waters(num_up+j)%solid_chemistry%rk(i)=(Rk_mean(&
                    !         mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)-Delta_t*(1d0-theta)*&
                    !         mix_waters(num_up+j)%solid_chemistry%rk_old(i))/(Delta_t*theta)
                    ! end do
                !end if
                !Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc,conc_nc,conc_nc_old,Rk_est)
            end if
        end do
end subroutine

!> This subroutine computes the kinetic reaction rate contributions after water mixing implicitly for both upstream and downstream waters
!> The estimation of the downstream waters reaction amoutns is done by extrapolation with the second time derivative of the reaction rates
!> We assume uniform time step
subroutine compute_rk_tilde_impl_opt3(mix_waters,mixing_ratios_Rk,num_up,num_down,theta,Delta_t,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(inout) :: mix_waters(:) !> waters that mix with each other
    !real(kind=8), intent(in) :: c_tilde(:)  !> same order as "mix_waters"
    real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> ordered in upstream & downstream waters, respectively
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step
    real(kind=8), intent(out) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    !integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    !integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each target
    real(kind=8), allocatable :: conc_nc(:) !> kinetic stoichiometric matrix associated to the previous target
    real(kind=8), allocatable :: Rk_est(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    !logical :: flag !> flag to check if an index is in "perm"
    real(kind=8), parameter :: eps=1d-16 !> machine epsilon
    !real(kind=8) :: mu=5d-1 !> parameter for extrapolation of downstream waters rk
    real(kind=8), parameter :: tol=1d-16 !> tolerance for LU decomposition
    
    !allocate(c_tilde(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    !allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_rk_tilde_impl_opt4"
    
!> Process
    !c_tilde=0d0 !> initialise concentration of aqueous variable activity species after mixing
    !rk_tilde_up=0d0 !> initialise kinetic reaction contributions after mixing
    rk_tilde=0d0 !> initialise kinetic reaction contributions after mixing
    !> Upstream waters
        do j=1,num_up
            if (mix_waters(j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
                !deallocate(Sk_nc)
                if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%rk+&
                        (1d0-theta)*mix_waters(j)%rk_old)
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),theta*mix_waters(j)%solid_chemistry%rk+&
                        (1d0-theta)*mix_waters(j)%solid_chemistry%rk_old)
                end if
                deallocate(Sk_nc)
            end if
        end do
    !> Downstream waters (we estimate them)
        do j=1,num_down
            if (mix_waters(num_up+j)%indices_rk%num_cols>0) then
                Sk_nc=mix_waters(num_up+j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                if (mix_waters(num_up+j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                    !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                    !if (allocated(Sk_nc_old) .and. norm_mat_inf(Sk_nc_old-Sk_nc)<eps) then
                    !    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),mix_waters(num_up)%rk)
                    !else
                    mix_waters(num_up+j)%Rk_est=Delta_t*((1d0+&
                        2d0*theta)*mix_waters(num_up+j)%rk_old-3d0*theta*mix_waters(num_up+j)%rk_old_old+theta*&
                        mix_waters(num_up+j)%rk_old_old_old)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),mix_waters(num_up+j)%Rk_est/Delta_t)
                    !end if
                    !deallocate(Sk_nc)
                end if
                if (allocated(mix_waters(num_up+j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                    mix_waters(num_up+j)%solid_chemistry%Rk_est=Delta_t*(&
                        (1d0+2d0*theta)*&
                        mix_waters(num_up+j)%solid_chemistry%rk_old-3d0*theta*&
                        mix_waters(num_up+j)%solid_chemistry%rk_old_old+theta*&
                        mix_waters(num_up+j)%solid_chemistry%rk_old_old_old)
                    rk_tilde=rk_tilde+mixing_ratios_Rk(num_up+j)*matmul(transpose(Sk_nc),&
                        mix_waters(num_up+j)%solid_chemistry%Rk_est/Delta_t)
                end if
                !Sk_nc_old=Sk_nc !> chapuza
                deallocate(Sk_nc)
            end if
        end do
end subroutine 



!> This subroutine computes the kinetic reaction rate contributions after water mixing, explicitly for all target waters
subroutine compute_rk_tilde_expl(mix_waters,mixing_ratios_Rk,num_up,num_down,theta,Delta_t,rk_tilde)
    implicit none
!> Argumnets
    type(aqueous_chemistry_c), intent(inout) :: mix_waters(:) !> waters that mix with each other
    !real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
    real(kind=8), intent(in) :: mixing_ratios_Rk(:)
    integer(kind=4), intent(in) :: num_up  !> number of upstream waters ( we don't care, this is for polymorphism only)
    integer(kind=4), intent(in) :: num_down  !> number of downstream waters ( we don't care, this is for polymorphism only)
    real(kind=8), intent(in) :: theta !> time weighting factor (=0)
    real(kind=8), intent(in) :: Delta_t !> time step (not used)
    !real(kind=8), intent(out) :: rk_tilde_up(:) !> kinetic reaction rate contributions after mixing
    !real(kind=8), intent(out) :: rk_tilde_down(:) !> kinetic reaction rate contributions after mixing
    real(kind=8), intent(out) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
!> Variables
    integer(kind=4) :: i !> index in "ind_var_act_species"
    integer(kind=4) :: j !> index mixing waters
    integer(kind=4) :: k !> index mineral kinetic reactions
    !integer(kind=4) :: num_perm !> counter in "perm"
    integer(kind=4) :: num_mix_wat !> number of mixing waters
    !integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each mixing water
    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
    !logical :: flag !> flag to check if an index is in "perm"
!> Pre-process
    !allocate(c_tilde(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
    !allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
    num_mix_wat=size(mix_waters)
    !if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_rk_tilde_expl"
!> Process
    !c_tilde=0d0 !> initialise concentration of aqueous variable activity species after mixing
    !rk_tilde_down=0d0 !> initialise aqueous kinetic reaction rates after mixing
    !rk_tilde_up=0d0 !> initialise aqueous kinetic reaction rates after mixing
    rk_tilde=0d0 !> initialise kinetic reaction rates contributions after mixing
    !num_perm=0 !< counter in "perm"
    !do i=1,mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species
    !    !> Autentica chapuza
    !    if (mix_waters(1)%ind_var_act_species(i)/=mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i))) then
    !        if (num_perm>0) then
    !            call is_int_in_1D_array(i,perm(1:num_perm,2),flag)
    !            if (flag .eqv. .true.) then
    !                continue
    !            else
    !                num_perm=num_perm+1
    !                perm(num_perm,1)=i
    !                perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
    !                    mix_waters(1)%ind_var_act_species)
    !            end if
    !        else
    !            num_perm=num_perm+1
    !            perm(num_perm,1)=i
    !            perm(num_perm,2)=get_ind_int(mix_waters(1)%indices_aq_species(mix_waters(1)%aq_phase%ind_diss_solids(i)),&
    !                mix_waters(1)%ind_var_act_species)
    !        end if
    !    end if
    !    do j=1,num_mix_wat
    !        !print *, mix_waters(j)%indices_aq_species(mix_waters(j)%aq_phase%ind_diss_solids(i))
    !        c_tilde(i)=c_tilde(i)+mixing_ratios_conc(j)*mix_waters(j)%concentrations(mix_waters(j)%indices_aq_species(&
    !            mix_waters(j)%aq_phase%ind_diss_solids(i)))
    !        
    !    end do
    !end do
    do j=1,num_mix_wat
        if (mix_waters(j)%indices_rk%num_cols>0) then
            Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
            !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
            !deallocate(Sk_nc)
            if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
                !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
                rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk_old)
                !deallocate(Sk_nc)
            end if
            if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
                rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%solid_chemistry%rk_old)
            end if
            deallocate(Sk_nc)
        end if
    end do

        !do j=1,num_up
        !    rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%rk(i)
        !end do
        !do j=1,num_down
        !    rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(num_up+j)*mix_waters(num_up+j)%rk(i)
        !end do
    !end do
    !!do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    !if (allocated(mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
    !!> Mixing waters
    !    j=1 !> index mixing waters
    !    k=1 !> index mineral kinetic reactions
    !    i=1 !> index mineral kinetic reactions
    !        do
    !        !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
    !            if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
    !                mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
    !                rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%solid_chemistry%rk_old(i)
    !                if (j<num_mix_wat) then
    !                    j=j+1
    !                    i=1
    !                    k=1
    !                else
    !                    exit
    !                end if
    !            else if (k<mix_waters(j)%solid_chemistry%mineral_zone%num_minerals_kin) then
    !                k=k+1
    !            else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
    !                i=i+1
    !                k=1
    !            else if (j<num_mix_wat) then
    !                j=j+1
    !                i=1
    !                k=1
    !            else
    !                exit
    !            end if
    !        end do
    !    !!> Upstream waters
        !j=1 !> index upstream waters
        !k=1 !> index mineral kinetic reactions
        !i=1 !> index mineral kinetic reactions
        !    do
        !    !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
        !        if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
        !            mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
        !            rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%solid_chemistry%rk(i)
        !            if (j<num_up) then
        !                j=j+1
        !                i=1
        !                k=1
        !            else
        !                exit
        !            end if
        !        else if (k<mix_waters(j)%solid_chemistry%mineral_zone%num_minerals_kin) then
        !            k=k+1
        !        else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
        !            i=i+1
        !            k=1
        !        else if (j<num_up) then
        !            j=j+1
        !            i=1
        !            k=1
        !        else
        !            exit
        !        end if
        !    end do
        !    !> Downstream waters
        !    j=1 !> index downstream waters
        !    k=1 !> index mineral kinetic reactions
        !    i=1 !> index mineral kinetic reactions
        !    do
        !    !do i=1,mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin
        !        if (mix_waters(1)%solid_chemistry%mineral_zone%ind_min_chem_syst(i)==&
        !            mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst(k)) then !> j=1 no es necesario
        !            rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(num_up+j)*mix_waters(num_up+j)%solid_chemistry%rk_old(i)
        !            if (j<num_down) then
        !                j=j+1
        !                i=1
        !                k=1
        !            else
        !                exit
        !            end if
        !        else if (k<mix_waters(num_up+j)%solid_chemistry%mineral_zone%num_minerals_kin) then
        !            k=k+1
        !        else if (i<mix_waters(1)%solid_chemistry%mineral_zone%num_minerals_kin) then
        !            i=i+1
        !            k=1
        !        else if (j<num_down) then
        !            j=j+1
        !            i=1
        !            k=1
        !        else
        !            exit
        !        end if
        !    end do
            !rk_tilde(i)=rk_tilde(i)+mixing_ratios_Rk(j)*mix_waters(j)%solid_chemistry%rk(i)
        !end do
    !end do
    !end if
    !c_tilde_aux=c_tilde !> we store the concentration of aqueous variable activity species after mixing (chapuza)
    !!rk_tilde_aux=rk_tilde !> we store the aqueous kinetic reaction rates after mixing (chapuza)
    !do i=1,num_perm
    !    !> we permute the indices in "perm"
    !    c_tilde(perm(i,1))=c_tilde_aux(perm(i,2))
    !    c_tilde(perm(i,2))=c_tilde_aux(perm(i,1))
    !    !rk_tilde(perm(i,1))=rk_tilde_aux(perm(i,2))
    !    !rk_tilde(perm(i,2))=rk_tilde_aux(perm(i,1))
    !end do
    !rk_tilde=rk_tilde*Delta_t
end subroutine 

!!> This subroutine computes the kinetic reaction rate contributions after water mixing, explicitly for all target waters
!subroutine compute_rk_tilde_expl(mix_waters,mixing_ratios_Rk,num_up,num_down,rk_tilde)
!    implicit none
!!> Argumnets
!    type(aqueous_chemistry_c), intent(in) :: mix_waters(:) !> waters that mix with each other
!    !real(kind=8), intent(in) :: mixing_ratios_conc(:)  !> same order as "mix_waters"
!    real(kind=8), intent(in) :: mixing_ratios_Rk(:)  !> ordered in upstream & downstream waters, respectively
!    integer(kind=4), intent(in) :: num_up  !> number of upstream waters ( we don't care, this is for polymorphism only)
!    integer(kind=4), intent(in) :: num_down  !> number of downstream waters ( we don't care, this is for polymorphism only)
!    !real(kind=8), intent(in) :: Delta_t !> time step
!    real(kind=8), intent(out) :: rk_tilde(:) !> kinetic reaction rate contributions after mixing
!!> Variables
!    integer(kind=4) :: i !> index in "ind_var_act_species"
!    integer(kind=4) :: j !> index mixing waters
!    integer(kind=4) :: k !> index mineral kinetic reactions
!    !integer(kind=4) :: num_perm !> counter in "perm"
!    integer(kind=4) :: num_mix_wat !> number of mixing waters
!    !integer(kind=4), allocatable :: perm(:,:) !> indices that have to be permuted
!    real(kind=8), allocatable :: Sk_nc(:,:) !> kinetic stoichiometric matrix associated to each mixing water
!    !real(kind=8), allocatable :: rk_tilde_aux(:) !> auxiliary variable to store the aqueous kinetic reaction rates after mixing
!    !logical :: flag !> flag to check if an index is in "perm"
!!> Pre-process
!    !allocate(c_tilde(mix_waters(1)%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
!    !allocate(perm(mix_waters(1)%aq_phase%num_species-mix_waters(1)%aq_phase%wat_flag,2))
!    num_mix_wat=size(mix_waters)
!    if (num_mix_wat/=num_up+num_down) error stop "Dimension error in compute_rk_tilde_expl"
!!> Process
!    rk_tilde=0d0 !> initialise kinetic reaction rates contributions after mixing
!    do j=1,num_mix_wat
!        if (mix_waters(j)%indices_rk%num_cols>0) then
!            Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
!            !rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk)
!            !deallocate(Sk_nc)
!            if (mix_waters(j)%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts>0) then
!                !Sk_nc=mix_waters(j)%get_Sk_nc() !> we get the kinetic stoichiometric matrix in j-th upstream water (chapuza)
!                rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%rk_old)
!                !deallocate(Sk_nc)
!            end if
!            if (allocated(mix_waters(j)%solid_chemistry%mineral_zone%ind_min_chem_syst)) then
!                rk_tilde=rk_tilde+mixing_ratios_Rk(j)*matmul(transpose(Sk_nc),mix_waters(j)%solid_chemistry%rk_old)
!            end if
!            deallocate(Sk_nc)
!        end if
!    end do
!end subroutine

!> This subroutine initialises concentrations for any iterative method
!> It takes a linear combination of the concentrations in the two previous time steps
subroutine initialise_iterative_method(conc_old_old,conc_old,param,initial_guess)
    real(kind=8), intent(in) :: conc_old_old(:) !> concentrations at time step k-1
    real(kind=8), intent(in) :: conc_old(:) !> concentrations at time step k
    real(kind=8), intent(in) :: param !> parameter for linear combination
    real(kind=8), intent(out) :: initial_guess(:) !> initial guess concentrations at time step k+1
    
    integer(kind=4) :: i !> loop index
    
    if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
    do i=1,size(conc_old)
        if ((1d0+param)*conc_old(i)-param*conc_old_old(i)<0d0) then
            initial_guess(i)=conc_old(i)
        else
            initial_guess(i)=(1d0+param)*conc_old(i)-param*conc_old_old(i)
        end if
    end do
end subroutine

!> This subroutine initialises aqueous primary concentrations for any iterative method
!> It takes a linear combination of the aqueous primary concentrations in the two previous time steps
subroutine initialise_c1_aq_iterative_method(this,param)
    class(aqueous_chemistry_c) :: this
    !real(kind=8), intent(in) :: conc_old_old(:) !> concentrations at time step k-1
    !real(kind=8), intent(in) :: conc_old(:) !> concentrations at time step k
    real(kind=8), intent(in) :: param !> parameter for linear combination
    !real(kind=8), intent(out) :: initial_guess(:) !> initial guess concentrations at time step k+1
    
    integer(kind=4) :: i !> loop index
    
    if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
    do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
        if ((1d0+param)*this%conc_old(this%ind_var_act_species(i))-param*this%conc_old_old(this%ind_var_act_species(i))<0d0) then
            this%concentrations(this%ind_var_act_species(i))=this%conc_old(this%ind_var_act_species(i))
        else
            this%concentrations(this%ind_var_act_species(i))=(1d0+param)*this%conc_old(this%ind_var_act_species(i))-&
                param*this%conc_old_old(this%ind_var_act_species(i))
        end if
    end do
end subroutine

    subroutine set_ind_species(this)!,ind_prim_species,ind_sec_var_act_species)
        !> This subroutine sets the indices of primary and secondary species in state variables
        !! It also sets the indices of variable activity species
    class(aqueous_chemistry_c) :: this
    !integer(kind=4), intent(in), optional :: ind_prim_species(:), ind_sec_var_act_species(:) !> if one of them is present, the other must be present as well
    
    integer(kind=4) :: i,j,k,l, cntr_cst_act_sp
    integer(kind=4), allocatable :: aux_ind_sp(:)
    integer :: n_p_aq_old !> number of primary aqueous species before changing definition of primary species

    !if (present(ind_prim_species) .and. present(ind_sec_var_act_species)) then
    !    this%ind_prim_species=ind_prim_species
    !else !> we assume aq_species array in aq_phase pointer is already sorted in primary and secondary species
        ! if (.not. allocated(this%ind_prim_species)) then
        !     allocate(this%ind_prim_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
        ! else
        !     deallocate(this%ind_prim_species)
        !     allocate(this%ind_prim_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
        ! end if
        ! if (.not. allocated(this%ind_sec_var_act_species)) then
        !     allocate(this%ind_sec_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species))
        ! else
        !     deallocate(this%ind_sec_var_act_species)
        !     allocate(this%ind_sec_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species))
        ! end if
        ! if (.not. allocated(this%ind_sec_species)) then
        !     allocate(this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species))
        ! else
        !     deallocate(this%ind_sec_species)
        !     allocate(this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species))
        ! end if
        j=1 !> counter for primary species
        k=1 !> counter for secondary variable activity species
        l=1 !> counter for secondary species
        if (this%solid_chemistry%reactive_zone%speciation_alg%flag_comp .eqv. .false.) then
            allocate(this%ind_prim_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
            allocate(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
            allocate(this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species))
            !allocate(this%ind_cst_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species))
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
                this%ind_prim_species(i)=this%indices_aq_species(i)
                if (this%aq_phase%aq_species(i)%cst_act_flag .eqv. .false.) then
                    this%ind_var_act_species(j)=this%indices_aq_species(i)
                    j=j+1
                end if
            end do
            do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species
                this%ind_sec_species(i)=this%indices_aq_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                if (this%aq_phase%aq_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)%cst_act_flag &
                    .eqv. .false.) then
                    this%ind_var_act_species(j)=this%indices_aq_species(&
                    this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i)
                    j=j+1
                end if
            end do
            ! do i=1,this%aq_phase%num_species
            !     if (this%indices_aq_species(i)<=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
            !         this%ind_prim_species(this%indices_aq_species(i))=i
            !         this%ind_var_act_species(this%indices_aq_species(i))=i
            !         j=j+1
            !     else
            !         this%ind_sec_species(l)=i
            !         if (this%aq_phase%aq_species(this%indices_aq_species(i))%cst_act_flag .eqv. .false.) then
            !             this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+k)=i
            !             k=k+1
            !         end if
            !         l=l+1
            !     end if
            ! end do
        else
            cntr_cst_act_sp=0 !> counter for constant activity species in first num_aq_prim_species species
            ! allocate(aux_ind_sp(this%aq_phase%num_species))
            ! n_p_aq_old=size(this%ind_prim_species)
            ! aux_ind_sp(1:n_p_aq_old)=this%ind_prim_species
            ! aux_ind_sp(n_p_aq_old+1:this%aq_phase%num_species)=this%ind_sec_species
            if (allocated(this%ind_prim_species)) deallocate(this%ind_prim_species)
            if (allocated(this%ind_sec_species)) deallocate(this%ind_sec_species)
            if (allocated(this%ind_var_act_species)) deallocate(this%ind_var_act_species)
            allocate(this%ind_prim_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
            allocate(this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species))
            allocate(this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_sec_aq_species))
            do i=1,this%aq_phase%num_species
                if (this%aq_phase%aq_species(i)%cst_act_flag .eqv. .true.) then
                    this%ind_sec_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species+l)=&
                        this%indices_aq_species(i)
                    l=l+1
                    if (i<=this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species) then
                        cntr_cst_act_sp=cntr_cst_act_sp+1 !> counter for constant activity species in first num_aq_prim_species species
                    end if
                else if (i<=&
                        this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+cntr_cst_act_sp) then
                    this%ind_prim_species(j)=this%indices_aq_species(i)
                    this%ind_var_act_species(j)=this%indices_aq_species(i)
                    j=j+1
                else
                    this%ind_sec_species(k)=this%indices_aq_species(i)
                    this%ind_var_act_species(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+k)=&
                        this%indices_aq_species(i)
                    k=k+1
                end if
            end do
            !deallocate(aux_ind_sp)
        end if
    end subroutine

    subroutine set_conc_sec_aq_var_act_species(this,c2_aq_nc) !> 
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: c2_aq_nc(:)
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_sec_var_act_species
            this%concentrations(this%ind_var_act_species(&
                this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species+i))=&
            c2_aq_nc(i)
        end do
    end subroutine

    subroutine set_indices_rk(this) !> sets indices in chemical system for kinetic reactions
        !! OJO: FALTAN REACCIONES CINÉTICAS LINEALES
        class(aqueous_chemistry_c) :: this
        integer(kind=4) :: i,j,k,num_rk
        num_rk=0
        call this%indices_rk%allocate_array(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+&
            this%solid_chemistry%mineral_zone%num_minerals_kin+this%solid_chemistry%reactive_zone%chem_syst%num_gas_kin_reacts)
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        if (this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts>0) then !> chapuza
            do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
                call this%indices_rk%cols(num_rk+i)%allocate_vector(&
                    this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%num_aq_rk)
                j=1 !> counter for aqueous variable activity species
                k=1 !> counter for relevant aqueous species in kinetic reaction rates
                do
                    if (this%solid_chemistry%reactive_zone%chem_syst%species(j)%name == &
                        this%aq_phase%aq_species(&
                        this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%indices_aq_phase(k))%name) then
                        this%indices_rk%cols(num_rk+i)%col_1(k)=j
                        if (k<this%solid_chemistry%reactive_zone%chem_syst%redox_kin_reacts(i)%num_aq_rk) then
                            k=k+1
                            j=1
                        else
                            exit
                        end if
                    else if (j<this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species) then
                        j=j+1
                    else
                        error stop "Error in set_indices_rk"
                    end if
                end do
            end do
        end if
        num_rk=num_rk+this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        if (this%solid_chemistry%mineral_zone%num_minerals_kin>0) then !> chapuza
            do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
                call this%indices_rk%cols(num_rk+i)%allocate_vector(&
                    this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                    this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk)
                j=1 !> counter for aqueous variable activity species
                k=1 !> counter for relevant aqueous species in kinetic reaction rates
                do
                    if (this%solid_chemistry%reactive_zone%chem_syst%species(j)%name == &
                        this%aq_phase%aq_species(&
                        this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%indices_aq_phase(k))%name) then
                        this%indices_rk%cols(num_rk+i)%col_1(k)=j
                        if (k<this%solid_chemistry%reactive_zone%chem_syst%min_kin_reacts(&
                            this%solid_chemistry%mineral_zone%ind_min_chem_syst(i))%num_aq_rk) then
                            k=k+1
                            j=1
                        else
                            exit
                        end if
                    else if (j<this%solid_chemistry%reactive_zone%speciation_alg%num_aq_var_act_species) then
                        j=j+1
                    else
                        error stop "Error in set_indices_rk"
                    end if
                end do
            end do
        end if
    end subroutine

    function get_rk_new(this) result(rk_new)
    !!> This function returns the new kinetic reaction rates ordered in aqueous & mineral reactions (ESTA MAL)
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: rk_new(:)
    integer(kind=4) :: i
    allocate(rk_new(this%indices_rk%num_cols))
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
        rk_new(i)=this%rk_new(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        rk_new(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%rk_new(i)
    end do
    !> Faltan gases
    end function

    function get_rk_old(this) result(rk_old)
        !!> This function returns the old kinetic reaction rates (chapuza) ordered in aqueous & mineral reactions (ESTA MAL)
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: rk_old(:)
        integer(kind=4) :: i
        allocate(rk_old(this%indices_rk%num_cols))
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts
            rk_old(i)=this%rk_old(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_aq_kin_reacts+i)=this%solid_chemistry%rk_old(i)
        end do
        !> Faltan gases
    end function

    function get_Sk_nc(this) result(Sk_nc) !> returns variable activity stoichiometric matrix of kinetic reactions (chapuza)
        class(aqueous_chemistry_c) :: this
        real(kind=8), allocatable :: Sk_nc(:,:)
        integer(kind=4) :: i
        allocate(Sk_nc(this%indices_rk%num_cols,this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species))
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            Sk_nc(i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(i,&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            Sk_nc(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i,:)=&
                this%solid_chemistry%reactive_zone%chem_syst%Sk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%ind_min_chem_syst(i),&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            Sk_nc(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+&
                i,:)=this%solid_chemistry%reactive_zone%chem_syst%Sk(&
                this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%reactive_zone%chem_syst%num_minerals_kin+i,&
                1:this%solid_chemistry%reactive_zone%speciation_alg%num_var_Act_species)
        end do
        !> Faltan gases
    end function
    
    subroutine compute_Rk(this,theta,Delta_t) !> computes kinetic reaction amounts during a time step
    class(aqueous_chemistry_c) :: this !> aqueous chemistry object
    real(kind=8), intent(in) :: theta !> time weighting factor
    real(kind=8), intent(in) :: Delta_t !> time step
    integer(kind=4) :: i
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
        this%Rk(i)=theta*this%rk_new(i)+(1d0-theta)*this%rk_old(i)
    end do
    do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
        this%solid_chemistry%Rk(i)=theta*this%solid_chemistry%rk_new(i)+(1d0-theta)*this%solid_chemistry%rk_old(i)
    end do
    do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
        this%Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
            theta*this%rk_new(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)+&
            (1d0-theta)*this%rk_old(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
    end do
    this%Rk=this%Rk*Delta_t
    this%solid_chemistry%Rk=this%solid_chemistry%Rk*Delta_t
    end subroutine

    subroutine set_rk_mean(this,rk_mean)
        !!> This subroutine sets the mean kinetic reaction rates during a time step
        !! We assume that rk_mean attribute is already allocated
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: rk_mean(:) !> mean kinetic reaction rates during a time step
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%rk_mean(i)=rk_mean(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%rk_mean(i)=rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                rk_mean(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+i)
        end do
    end subroutine

    subroutine set_Rk(this,Rk)
        !!> This subroutine sets the kinetic reaction amounts during a time step
        !! We assume that Rk is already allocated
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: Rk(:) !> kinetic reaction amounts during a time step
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%Rk(i)=Rk(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%Rk(i)=Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                Rk(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+i)
        end do
        !this%Rk_mean=this%Rk_mean*Delta_t
        !this%solid_chemistry%Rk_mean=this%solid_chemistry%Rk_mean*Delta_t
    end subroutine

    subroutine initialise_conc_nc_iterative_method_ideal(this,param,conc_nc)
        !!> This subroutine initialises variable activity concentrations for any iterative method
        !! It takes a linear combination of the primary concentrations in the two previous time steps
        !! We assume all primary species are aqueous

        class(aqueous_chemistry_c) :: this !> aqueous chemistry object
        real(kind=8), intent(in) :: param !> parameter for linear combination
        real(kind=8), intent(out) :: conc_nc(:) !> concentrations at time step k+1 (must be already allocated)

        real(kind=8), allocatable :: c2nc(:) !> secondary variable activity concentrations (must be already allocated)
        integer(kind=4) :: i !> loop index

        
        
        if (param<0d0 .or. param>1d0) error stop "Initialisation parameter in iterative method must be in [0,1]"
        call this%initialise_c1_aq_iterative_method(param)
        !> chapuza
        do i=1,this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
            conc_nc(i)=this%concentrations(this%ind_var_act_species(i))
        end do
        !> autentica chapuza
        if (this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions>0) then
            allocate(c2nc(this%solid_chemistry%reactive_zone%speciation_alg%num_eq_reactions))
            call this%compute_c2nc_from_c1_aq_ideal(c2nc)
            conc_nc(this%solid_chemistry%reactive_zone%speciation_alg%num_prim_species+1:&
                this%solid_chemistry%reactive_zone%speciation_alg%num_var_act_species)=c2nc
            deallocate(c2nc)
        end if
    end subroutine
    
    subroutine set_Rk_est(this,Rk_est)
    !!> This subroutine sets the estimated kinetic reaction rates in the consistent WMA approach
        class(aqueous_chemistry_c) :: this
        real(kind=8), intent(in) :: Rk_est(:) !> estimated kinetic reaction rates
        integer(kind=4) :: i
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts
            this%Rk_est(i)=Rk_est(i)
        end do
        do i=1,this%solid_chemistry%mineral_zone%num_minerals_kin
            this%solid_chemistry%Rk_est(i)=Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)
        end do
        do i=1,this%solid_chemistry%reactive_zone%chem_syst%num_redox_kin_reacts
            this%Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+i)=&
                Rk_est(this%solid_chemistry%reactive_zone%chem_syst%num_lin_kin_reacts+&
                this%solid_chemistry%mineral_zone%num_minerals_kin+i)
        end do
    end subroutine
    
    function get_u_aq(this) result(u_aq)
    !> This function returns the aqueous component concentrations
    class(aqueous_chemistry_c) :: this
    real(kind=8), allocatable :: u_aq(:)
    ! print *, "Getting u_aq"
    ! print *, "Number of aqueous primary species: ", this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species
    ! print *, "Indices of variable activity species: ", this%ind_var_act_species
    ! print *, "Concentrations: ", this%concentrations(this%ind_var_act_species)
    ! print *, "Comp mat aq: ", this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq
    allocate(u_aq(this%solid_chemistry%reactive_zone%speciation_alg%num_aq_prim_species))
    u_aq=matmul(this%solid_chemistry%reactive_zone%speciation_alg%comp_mat_aq,&
        this%concentrations(this%ind_var_act_species))
    end function get_u_aq
    
    subroutine set_spec_alg_aq_chem(this,flag_comp,flag_surf,flag_Se,swap)
    !!> This subroutine sets all necessary variables in the speciation algebra object of the reactive zone pointer
        class(aqueous_chemistry_c) :: this
        logical, intent(in) :: flag_comp !> flag to indicate component matrix definition
        logical, intent(out) :: flag_surf !> flag to indicate if there are surface complexation reactions
        logical, intent(out) :: flag_Se !> TRUE if stoichiometric matrix has been modified, FALSE otherwise
        integer(kind=4), intent(out) :: swap(:) !> species or reactions to swap in stoichiometric matrix (already allocated)

        !logical :: flag_surf !> flag to indicate if there are surface complexation reactions

        call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(flag_comp)
        if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
            flag_surf=.true.
        else
            flag_surf=.false.
            call this%solid_chemistry%allocate_conc_solids()
            call this%solid_chemistry%allocate_activities()
            call this%solid_chemistry%allocate_log_act_coeffs_solid_chem()
        end if
        call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf)
        call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions()
        call this%set_ind_species()
        call this%solid_chemistry%reactive_zone%set_ind_eq_reacts()
        call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
        call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat()
        call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat()
        !allocate(swap(2))
        if (associated(this%gas_chemistry)) then
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap,this%gas_chemistry%activities)
        else
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
        end if
    end subroutine

    subroutine change_spec_alg_aq_chem(this)
    !!> This subroutine eliminates the constant activity species from the component matrix in the speciation algebra object of the &
        !! reactive zone pointer
        class(aqueous_chemistry_c) :: this !> aqueous chemistry object
        !logical, intent(in) :: flag_comp !> flag to indicate component matrix definition
        !logical, intent(out) :: flag_surf !> flag to indicate if there are surface complexation reactions
        !logical, intent(out) :: flag_Se !> TRUE if stoichiometric matrix has been modified, FALSE otherwise
        !integer(kind=4), intent(out) :: swap(:) !> species or reactions to swap in stoichiometric matrix (already allocated)

        logical :: flag_Se !> flag to indicate if stoichiometric matrix has been modified
        integer(kind=4), allocatable :: swap(:) !> species or reactions to swap in stoichiometric matrix

        !call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_comp(flag_comp)
        ! if (this%solid_chemistry%reactive_zone%cat_exch_zone%num_surf_compl>0) then
        !     flag_surf=.true.
        ! else
        !     flag_surf=.false.
        !     call this%solid_chemistry%allocate_conc_solids()
        !     call this%solid_chemistry%allocate_activities()
        !     call this%solid_chemistry%allocate_log_act_coeffs_solid_chem()
        ! end if
        ! call this%solid_chemistry%reactive_zone%speciation_alg%set_flag_cat_exch(flag_surf)
        call this%solid_chemistry%reactive_zone%set_speciation_alg_dimensions(.true.)
        call this%set_ind_species()
        call this%solid_chemistry%reactive_zone%set_ind_eq_reacts()
        call this%solid_chemistry%reactive_zone%set_stoich_mat_react_zone()
        call this%solid_chemistry%reactive_zone%set_ind_mins_stoich_mat()
        call this%solid_chemistry%reactive_zone%set_ind_gases_stoich_mat()
        allocate(swap(2))
        if (associated(this%gas_chemistry)) then
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap,this%gas_chemistry%activities)
        else
            call this%solid_chemistry%reactive_zone%compute_speciation_alg_arrays(flag_Se,swap)
        end if
    end subroutine

    subroutine set_pos(this,pos) !> set position in space
        !!> This subroutine sets the position of the aqueous chemistry object in the global state vector
        class(aqueous_chemistry_c) :: this !> aqueous chemistry object
        real(kind=8), intent(in) :: pos(:) !> position in space
        this%pos=pos
    end subroutine
    
end module
