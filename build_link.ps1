$ErrorActionPreference = 'Stop'
$gfortran = 'C:\Users\jordi\OneDrive\Documentos\fortran\mingw64\bin\gfortran.exe'
$root = 'C:\Users\jordi\source\repos\jordipg10\interfaz_remix'
$obj = @(
  'array_ops_m.o','strings_m.o','vectors_m.o','arrays_m.o','metodos_sist_lin_m.o','penta_cholesky_inv_m.o',
  'polynomials_m.o','special_fcts_m.o','time_discr_m.o','target_m.o','spatial_discr_m.o','spatial_discr_1D_m.o',
  'spatial_discr_2D_m.o','spatial_discr_rad_m.o','aq_phase_m.o','aq_species_m.o','aqueous_chemistry_m.o','biofilm_m.o',
  'chem_out_options_m.o','chem_system_m.o','chemistry_interface_m.o','compare_react_zones.o',
  'compute_c_nc_from_u_Newton_ideal.o','compute_c_nc_from_u_Newton.o','compute_c2_from_c1_ideal.o','compute_c2_from_c1_Picard.o',
  'compute_c2nc_from_c1_ideal.o','compute_c2nc_from_c1_Picard.o','compute_conc_surf_ideal_Newton.o','compute_conc_surf_ideal_Picard.o',
  'compute_conc_surf_ideal_bin.o','compute_d_log_gamma_d_I_aq_chem.o','compute_dc2_dc1_ideal.o','compute_dc2_dc1.o',
  'compute_dc2nc_dc1_ideal.o','compute_dc2nc_dc1.o','compute_dfk_dc1_EI_ideal.o','compute_dfk_dc_EI_ideal.o',
  'compute_drk_dc_mineral.o','compute_drk_dc_Monod.o','compute_Jacobian_rk_anal.o','compute_mass_bal_mins.o',
  'compute_Re.o','compute_Re_kin.o','compute_res_init.o','compute_res_Jac_res_anal_ideal.o','compute_res_Jac_res_anal_exch.o',
  'compute_res_Jac_res_anal.o','compute_res_Jac_res_incr_coeff.o','compute_rk_Jac_rk_anal.o','compute_rk_Jac_rk_incr_coeff.o',
  'compute_rk_mineral.o','compute_rk_Monod.o','compute_rk_new.o','compute_rk_drk_dc_Monod.o','convergence_params_m.o',
  'eq_reaction_m.o','exch_sites_conv_m.o','Gaines_Thomas_m.o','Gapon_m.o','gas_chemistry_m.o','gas_species_m.o','gas_phase_m.o',
  'get_indices_reaction.o','initialise_conc_anal_ideal.o','initialise_conc_anal_exch.o','initialise_conc_anal.o','initialise_conc_incr_coeff.o',
  'interfaz_comps_arch.o','interfaz_comps_arch_mix.o','interfaz_esp_arch.o','interfaz_comps_vars.o','kin_mineral_m.o',
  'kin_mineral_params_m.o','kin_params_m.o','kin_reaction_m.o','lin_kin_reaction_m.o','link_target_gases_reactive_zone.o',
  'link_target_solids_reactive_zone.o','link_target_solids_mineral_zone.o','link_target_waters_mineral_zone.o',
  'link_target_waters_reactive_zone.o','link_target_waters_target_gases.o','link_target_waters_target_solids.o',
  'local_chemistry_m.o','mineral_m.o','mineral_zone_m.o','mixing_iter_comp_ideal.o','mixing_iter_comp.o','mixing_iter_species.o',
  'Monod_params_m.o','Newton_EI_eq_ideal.o','Newton_EI_eq_kin_anal_ideal_opt2.o','Newton_EI_kin_anal_ideal_opt2.o',
  'params_act_coeff_m.o','params_aq_sol_m.o','params_spec_vol_m.o','phase_m.o','compute_react_term_EE_eq_kin.o','compute_react_term_EE_kin.o',
  'reaction_m.o','reactive_zone_m.o','read_association_react_PHREEQC.o','read_chem_system_CHEPROO.o','read_chemistry_interface_CHEPROO.o',
  'read_chemistry_interface.o','read_dissolution_react_PHREEQC.o','read_exchange_react_PHREEQC.o','read_init_bd_wat_types_CHEPROO.o',
  'read_init_gas_zones_CHEPROO.o','read_init_min_zones_CHEPROO.o','read_init_cat_exch_zones_CHEPROO.o','read_kinetics_DB.o',
  'read_master25.o','read_Monod_DB.o','read_PFLOTRAN_chemistry.o','read_PFLOTRAN_DB.o','read_PHREEQC_DB_opc1.o','read_PHREEQC_DB_opc2.o',
  'read_species.o','read_waters_init.o','read_tar_gas.o','read_tar_sol.o','read_wat_type_CHEPROO.o','redox_kin_reaction_m.o',
  'set_stoich_mat_gas.o','set_stoich_mat_react_zone.o','set_stoich_mat_sol.o','set_stoich_mat.o','solid_chemistry_m.o','solid_m.o',
  'solve_reactive_mixing_ideal_lump_iter.o','speciation_algebra_m.o','species_m.o','surf_compl_m.o','update_eq_reactions.o','Vanselow_m.o',
  'reactive_mixing_iter_EE_eq_kin_ideal.o','reactive_mixing_iter_EE_eq_kin.o','reactive_mixing_iter_EE_kin.o',
  'reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2.o','reactive_mixing_iter_EI_kin_anal_ideal_opt2.o','reactive_mixing_iter_EI_eq_anal_ideal.o',
  'write_chemistry.o','main_interfaz.o'
)
Push-Location (Join-Path $root 'obj')
& $gfortran -o ../bin/interfaz_remix.exe -static @obj (Join-Path $root 'lib/liblapack.a') (Join-Path $root 'lib/librefblas.a')
Write-Output "link exit=$LASTEXITCODE"
Pop-Location
