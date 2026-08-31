#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Portable build script for interfaz_remix on Linux and macOS.
# (On Windows, use the VS Code tasks or see BUILD_GUIDE.md.)
#
# IMPORTANT: the bundled lib/liblapack.a and lib/librefblas.a are Windows/MinGW
# archives and are NOT used here. This script links the platform-native BLAS/LAPACK:
#   - Linux : system reference/OpenBLAS  ->  -llapack -lblas
#             (install first, e.g.  sudo apt-get install -y liblapack-dev libblas-dev)
#   - macOS : the built-in Accelerate framework  ->  -framework Accelerate
#             (nothing to install; gfortran comes from  brew install gcc)
#
# The source list below mirrors the compile order in .vscode/tasks.json
# (compile-discr -> compile-chem -> compile-main). Keep the two in sync when
# you add or remove source files.
#
# Usage:
#   chmod +x build.sh
#   ./build.sh                                   # default (debug) build
#   FFLAGS="-O2 -fbacktrace" ./build.sh          # optimised/release build
#   FC=gfortran-14 ./build.sh                    # pick a specific gfortran
#
# Output: bin/interfaz_remix
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$ROOT/src"
OBJ="$ROOT/obj"
MOD="$ROOT/mod"
BIN="$ROOT/bin"
mkdir -p "$OBJ" "$MOD" "$BIN"

FC="${FC:-gfortran}"
# Default matches the flags the project is validated with. Override via env for release.
FFLAGS="${FFLAGS:--g -O0 -fcheck=all -fbacktrace -ffree-line-length-none -fno-range-check}"

# Ordered source list (base names, no .f90), mirroring .vscode/tasks.json.
SRCS=(
  # --- discretisation / utilities (compile-discr) ---
  array_ops_m
  strings_m
  vectors_m
  arrays_m
  metodos_sist_lin_m
  penta_cholesky_inv_m
  polynomials_m
  special_fcts_m
  target_m
  spatial_discr_m
  spatial_discr_1D_m
  spatial_discr_2D_m
  spatial_discr_rad_m
  time_discr_m
  # --- chemistry layer + routines (compile-chem) ---
  convergence_params_m
  params_aq_sol_m
  params_spec_vol_m
  params_act_coeff_m
  species_m
  solid_m
  gas_species_m
  aq_species_m
  phase_m
  biofilm_m
  reaction_m
  aq_phase_m
  gas_phase_m
  eq_reaction_m
  kin_params_m
  kin_reaction_m
  lin_kin_reaction_m
  kin_mineral_params_m
  mineral_m
  kin_mineral_m
  Monod_params_m
  redox_kin_reaction_m
  speciation_algebra_m
  exch_sites_conv_m
  Gaines_Thomas_m
  Gapon_m
  Vanselow_m
  surf_compl_m
  chem_system_m
  mineral_zone_m
  reactive_zone_m
  local_chemistry_m
  gas_chemistry_m
  solid_chemistry_m
  aqueous_chemistry_m
  chem_out_options_m
  chemistry_interface_m
  compare_react_zones
  compute_c_nc_from_u_Newton_ideal
  compute_c_nc_from_u_Newton
  compute_c2_from_c1_ideal
  compute_c2_from_c1_Picard
  compute_c2nc_from_c1_ideal
  compute_c2nc_from_c1_Picard
  compute_conc_surf_ideal_Newton
  compute_conc_surf_ideal_Picard
  compute_conc_surf_ideal_bin
  compute_d_log_gamma_d_I_aq_chem
  compute_dc2_dc1_ideal
  compute_dc2_dc1
  compute_dc2nc_dc1_ideal
  compute_dc2nc_dc1
  compute_dfk_dc1_EI_ideal
  compute_dfk_dc_EI_ideal
  compute_drk_dc_mineral
  compute_drk_dc_Monod
  compute_Jacobian_rk_anal
  compute_mass_bal_mins
  compute_Re
  compute_Re_kin
  compute_res_init
  compute_res_Jac_res_anal_ideal
  compute_res_Jac_res_anal_exch
  compute_res_Jac_res_anal
  compute_res_Jac_res_incr_coeff
  compute_rk_Jac_rk_anal
  compute_rk_Jac_rk_incr_coeff
  compute_rk_mineral
  compute_rk_Monod
  compute_rk_new
  compute_rk_drk_dc_Monod
  get_indices_reaction
  initialise_conc_anal_ideal
  initialise_conc_anal_exch
  initialise_conc_anal
  initialise_conc_incr_coeff
  interfaz_comps_arch
  interfaz_comps_arch_mix
  interfaz_esp_arch
  interfaz_comps_vars
  link_target_gases_reactive_zone
  link_target_solids_reactive_zone
  link_target_solids_mineral_zone
  link_target_waters_mineral_zone
  link_target_waters_reactive_zone
  link_target_waters_target_gases
  link_target_waters_target_solids
  mixing_iter_comp_ideal
  mixing_iter_comp
  mixing_iter_species
  Newton_EI_eq_ideal
  Newton_EI_kin_anal_ideal_opt2
  Newton_EI_eq_kin_anal_ideal_opt2
  compute_react_term_EE_eq_kin
  compute_react_term_EE_kin
  read_association_react_PHREEQC
  read_chem_system_CHEPROO
  read_chemistry_interface_CHEPROO
  read_chemistry_interface
  read_dissolution_react_PHREEQC
  read_exchange_react_PHREEQC
  read_init_bd_wat_types_CHEPROO
  read_init_gas_zones_CHEPROO
  read_init_min_zones_CHEPROO
  read_init_cat_exch_zones_CHEPROO
  read_kinetics_DB
  read_master25
  read_Monod_DB
  read_PFLOTRAN_chemistry
  read_PFLOTRAN_DB
  read_PHREEQC_DB_opc1
  read_PHREEQC_DB_opc2
  read_species
  read_waters_init
  read_tar_gas
  read_tar_sol
  read_wat_type_CHEPROO
  set_stoich_mat_gas
  set_stoich_mat_react_zone
  set_stoich_mat_sol
  set_stoich_mat
  solve_reactive_mixing_ideal_lump_iter
  update_eq_reactions
  reactive_mixing_iter_EE_eq_kin_ideal
  reactive_mixing_iter_EE_eq_kin
  reactive_mixing_iter_EE_kin
  reactive_mixing_iter_EI_eq_kin_anal_ideal_opt2
  reactive_mixing_iter_EI_kin_anal_ideal_opt2
  reactive_mixing_iter_EI_eq_anal_ideal
  write_chemistry
  # --- driver (compile-main) ---
  main_interfaz
)

echo ">> Compiling ${#SRCS[@]} files with: $FC $FFLAGS"
for f in "${SRCS[@]}"; do
  src="$SRC/$f.f90"
  [ -f "$src" ] || { echo "ERROR: source not found: $src"; exit 1; }
  echo "   [FC] $f.f90"
  "$FC" $FFLAGS -J "$MOD" -c "$src" -o "$OBJ/$f.o"
done

# Object list in the same order.
OBJS=()
for f in "${SRCS[@]}"; do OBJS+=("$OBJ/$f.o"); done

case "$(uname -s)" in
  Darwin)
    echo ">> Linking (macOS, Accelerate framework)"
    "$FC" -o "$BIN/interfaz_remix" "${OBJS[@]}" -framework Accelerate
    ;;
  *)
    echo ">> Linking (Linux, system LAPACK/BLAS)"
    # Embed the gfortran/GCC runtimes so the binary runs without gfortran
    # installed; LAPACK/BLAS stay dynamic (need liblapack.so.3 / libblas.so.3).
    "$FC" -static-libgfortran -static-libgcc -o "$BIN/interfaz_remix" "${OBJS[@]}" -llapack -lblas
    ;;
esac

chmod +x "$BIN/interfaz_remix"
echo ">> Done: $BIN/interfaz_remix"
