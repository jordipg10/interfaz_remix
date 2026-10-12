!> \file biofilm_m.f90
!> \brief Biofilm module for modeling biofilm processes in reactive transport
!> \details This module defines the biofilm class and related functionality for simulating 
!> biofilm growth, decay, and interaction with aqueous chemistry in porous media.
!>
!> The biofilm is modeled as an immobile zone (extending phase_c). This class holds ONLY the
!> intrinsic, medium-independent PARAMETERS of the biofilm (constitutive properties and
!> kinetic/model coefficients). The LOCAL STATE variables that depend on position and evolve in
!> time — biofilm thickness, internal composition fractions (water/cells/EPS), medium porosities,
!> morphology descriptors, the chemical microenvironment and the per-microorganism volumetric
!> fractions — live in solid_chemistry_c (see solid_chemistry_m).
!>
!> As a parameter set, the biofilm type:
!> - Hosts multiple microorganism species with individual growth/decay kinetics
!> - Provides mechanical, EPS, transport, detachment, temperature and bioclogging parameters
!> - Is pointed to by each solid_chemistry_c object that carries the biofilm's local state
!>
!> Parameter groups held by this class:
!> - Physical/mechanical parameters (max thickness, dry density, elastic/storage/loss moduli,
!>   yield stress, viscosity, cohesive energy, matrix permeability)
!> - EPS composition and production (polysaccharides, proteins, eDNA, lipids; Y_eps, alpha, beta)
!> - Transport parameters (effective diffusion ratios, tortuosity)
!> - Detachment/erosion model coefficients
!> - Temperature correction via Arrhenius (T_ref, theta, activation energy)
!> - Porosity-permeability coupling for porous media bioclogging (min porosity, KC exponent)
!> - Microorganism definitions (count and species array with kinetics)
!>
!> Property references grouped by category:
!>
!> Physical/mechanical:
!>   - Zhang & Bishop (1994) Water Research 28(11):2267 — density, porosity stratification
!>   - Billings et al. (2015) Rep. Prog. Phys. 78:036601 — mechanical review
!>   - Lieleg et al. (2011) Soft Matter 7:3307 — rheology P. aeruginosa
!>   - Boudarel et al. (2018) npj Biofilms Microbiomes 4:17 — standardized mechanics
!>   - Charlton et al. (2019) J. Bacteriol. 201:e00101-19 — viscoelasticity
!>
!> Chemical/EPS:
!>   - Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623 — EPS matrix review
!>   - Quan et al. (2022) Crit. Rev. Microbiol. 48:283 — water types in biofilm
!>
!> Transport:
!>   - Stewart (2003) J. Bacteriol. 185:1485 — diffusion in biofilms
!>   - Berg et al. (2020) Biotechnol. Bioeng. 117:2997 — diffusion measurement
!>
!> Kinetic:
!>   - Henze et al. (1987) ASM1, IAWPRC Report No. 1
!>   - Picioreanu et al. (2004) AEM 70:3024 — multispecies kinetics
!>   - Kreft et al. (2001) Microbiology 147:2897 — Luedeking-Piret EPS
!>
!> Structural:
!>   - Costerton et al. (2003) J. Clin. Invest. 112:1466
!>   - Hall-Stoodley et al. (2022) Nat. Rev. Microbiol. 20:608
!>
!> Detachment:
!>   - Wanner & Gujer (1986) Biotechnol. Bioeng. 28:314
!>   - Pitts et al. (2007) AEM 73:2870 — cohesive energy
!>   - Abbas et al. (2012) Math. Biosci. Eng. 9:215 — shear-dependent models
!>
!> Porosity-permeability:
!>   - Hommel et al. (2018) Transport in Porous Media — Kozeny-Carman review
!>   - Thullner et al. (2002) Water Resources Research — bioclogging models
!>
!> \see phase_m For base class
!> \see species_m For microorganism definitions
!> \see redox_kin_reaction_m For microbial reaction kinetics
!> \see Monod_params_m For substrate limitation parameters
!>
!> \author Jordi Petchame Guerrero
!> \date 2025-2026

module biofilm_m                                                                !< Module defining biofilm class for microbial biofilm processes
    use phase_m, only: phase_c
    use species_m, only: species_c, microorganism_c                             !< Import microorganism class for defining microbial species
    implicit none                                                               !< Enforce explicit variable declarations (Fortran best practice)
    save                                                                        !< Preserve module variable values between procedure calls
    private                                                                     !< Default accessibility is private; public entities must be explicitly declared
    
    !> \brief Biofilm class extending phase functionality
    !> \details Represents the medium-independent PARAMETER set of a biofilm (a special type of
    !> immobile zone containing microorganisms). Local state (thickness, composition fractions,
    !> porosities, morphology, chemical microenvironment, microorganism volume fractions) is stored
    !> per location in solid_chemistry_c, which holds a pointer to a biofilm_c parameter object.
    !>
    !> Physical interpretation (the state symbols below live in solid_chemistry_c):
    !> - Biofilm occupies a volumetric fraction of the porous medium (phi_imm from parent class)
    !> - Within the biofilm, water occupies fraction f_w, biomass occupies the remainder
    !> - Each microorganism type has its own volumetric fraction within the biofilm
    !>
    !> Mass balance considerations (using solid_chemistry_c state f_w, vol_fracts_microorgs):
    !> - Water in biofilm: V_water = phi_imm * f_w * V_total
    !> - Biomass volume: V_biomass = phi_imm * (1 - f_w) * V_total
    !> - Microorganism i volume: V_microorg_i = vol_fracts_microorgs(i) * V_biomass
    !>
    !> The parameters are organized into logical groups:
    !>   1. Physical/mechanical
    !>   2. EPS composition
    !>   3. Transport
    !>   4. Detachment/erosion
    !>   5. Temperature correction
    !>   6. Porosity-permeability coupling (bioclogging)
    !>   7. Biological/kinetic (microorganisms array)
    !>
    !> \note The immobile porosity (phi_imm) represents the total biofilm volumetric fraction
    type, public, extends(phase_c) :: biofilm_c                                 !< Biofilm class inheriting from phase base class
        
        !=====================================================================
        ! 1. PHYSICAL / MECHANICAL PROPERTIES
        !=====================================================================
        
        !> \brief Maximum allowed biofilm thickness [m]
        !> \details Upper limit for biofilm growth; used in detachment models.
        !> Sloughing occurs when thickness exceeds Lf_max:
        !>   r_det = k_sl * (Lf - Lf_max) for Lf > Lf_max
        !>
        !> Ref: Wanner & Gujer (1986) Biotechnol. Bioeng. 28:314
        real(kind=8) :: thickness_max                                           !< Maximum biofilm thickness Lf_max [m] (typ. 5e-4 to 2e-3)
        
        !> \brief Dry biomass density [g/m^3]
        !> \details Volumetric density of the dry biomass (VSS) within the biofilm.
        !> Stratified: base layers 5-10x denser than surface layers.
        !> Typical model values: 10,000-100,000 g VSS/m^3 (10-100 kg/m^3).
        !> Picioreanu et al. (1999) use CXm = 70 kg VSS/m^3.
        !> Wet biofilm density is near water: ~1000-1010 kg/m^3.
        !>
        !> Ref: Zhang & Bishop (1994) Water Research 28:2267
        !> Ref: Rittmann & McCarty (2001) Environmental Biotechnology, McGraw-Hill
        !> Ref: Picioreanu et al. (1999) Water Sci. Tech. 39(7):115
        real(kind=8) :: density_dry                                             !< Dry biomass density rho_X [g/m^3] (typ. 1e4-1e5)
        
        !> \brief Elastic modulus (Young's modulus) [Pa]
        !> \details Characterizes biofilm stiffness in response to deformation.
        !> Highly variable: from 0.1 Pa (very soft) to 9 MPa (basal layers).
        !> Typical range for soft biofilms: 100-10,000 Pa (0.1-10 kPa).
        !> Drinking water biofilms measured by AFM: 200-9000 kPa.
        !>   - Soft top layers: < 600 kPa
        !>   - Harder basal layers: up to 9000 kPa
        !>
        !> Ref: Billings et al. (2015) Rep. Prog. Phys. 78:036601
        !> Ref: Francius et al. (2011) Biofouling 27:729
        real(kind=8) :: elastic_modulus                                         !< Young's modulus E [Pa] (typ. 1e2-1e4)
        
        !> \brief Storage (shear) modulus G' [Pa]
        !> \details Elastic component of viscoelastic response.
        !> P. aeruginosa plateau modulus at 1 Hz: (2.0 +/- 0.9) kPa.
        !> General range across species: 0.1-1000 Pa (macrorheology).
        !> Biofilms are viscoelastic solids: G' > G'' at low frequencies (< 10 Hz).
        !>
        !> Ref: Lieleg et al. (2011) Soft Matter 7:3307
        !> Ref: Boudarel et al. (2018) npj Biofilms Microbiomes 4:17
        real(kind=8) :: storage_modulus                                         !< Storage modulus G' [Pa] (typ. 1e0-1e3)
        
        !> \brief Loss (viscous) modulus G'' [Pa]
        !> \details Viscous component of viscoelastic response.
        !> At low frequencies: G'' < G' (elastic-dominated, tan(delta) < 1).
        !> At high frequencies (> 30 Hz): G'' may exceed G'.
        !>
        !> Ref: Charlton et al. (2019) J. Bacteriol. 201:e00101-19
        real(kind=8) :: loss_modulus                                            !< Loss modulus G'' [Pa]
        
        !> \brief Yield stress [Pa]
        !> \details Stress threshold above which biofilm deforms plastically.
        !> Onset of non-linear behavior at a few percent strain.
        !> After yielding: modulus drops by > 3 orders of magnitude.
        !> Self-healing occurs within minutes.
        !> Range: 0.1-100 Pa.
        !>
        !> Ref: Lieleg et al. (2011) Soft Matter 7:3307
        !> Ref: Boudarel et al. (2018) npj Biofilms Microbiomes 4:17
        real(kind=8) :: yield_stress                                            !< Yield stress sigma_y [Pa] (typ. 0.1-100)
        
        !> \brief Dynamic viscosity [Pa*s]
        !> \details Frequency-dependent viscosity of the biofilm matrix.
        !> Range: 0.1-1000 Pa*s (100-1e6 mPa*s).
        !> S. aureus: ~10 kPa*s.
        !>
        !> Ref: Shaw et al. (2004); Rupp et al. (2005) — reviewed in Boudarel et al. (2018)
        real(kind=8) :: dynamic_viscosity                                       !< Dynamic viscosity eta [Pa*s] (typ. 0.1-1000)
        
        !> \brief Cohesive energy density [J/m^3 or nJ/um^3]
        !> \details Energy per unit volume required to detach biofilm via abrasion.
        !> Without Ca2+: 0.10 +/- 0.07 nJ/um^3 = 1e5 J/m^3
        !> With 10 mM Ca2+: 1.98 +/- 0.34 nJ/um^3 = ~2e6 J/m^3
        !>
        !> Ref: Pitts et al. (2007) AEM 73:2870
        real(kind=8) :: cohesive_energy                                         !< Cohesive energy density [J/m^3] (typ. 1e5-2e6)
        
        !> \brief Permeability of the biofilm matrix [m^2]
        !> \details Hydraulic permeability of the biofilm gel itself.
        !> Range: 1e-15 to 1e-7 m^2 depending on density and structure.
        !> Typical starting value for porous media: ~1e-12 m^2.
        !> High-density base layers are nearly impermeable;
        !> the loose surface layer may allow convective flow in channels.
        !>
        !> Ref: Thullner et al. (2002) Water Resources Research
        !> Ref: Doro & Stocker (2023) ETH Zurich — morphogenesis review
        real(kind=8) :: permeability                                            !< Biofilm matrix permeability kappa [m^2] (typ. 1e-15 to 1e-7)
        
        !=====================================================================
        ! 2. EPS COMPOSITION (mass fractions of dry EPS)
        !=====================================================================
        
        !> \brief Mass fraction of polysaccharides in dry EPS [-]
        !> \details Major structural component of biofilm matrix; species-specific.
        !> Range: 0.40-0.95 of dry EPS mass.
        !> Species-specific examples:
        !>   - P. aeruginosa (mucoid): alginate
        !>   - P. aeruginosa (non-mucoid): Psl, Pel
        !>   - S. aureus / S. epidermidis: PNAG
        !>   - B. subtilis: exopolysaccharide + TasA protein
        !>   - E. coli K-12: cellulose + curli amyloid fibers
        !>
        !> Ref: Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623
        !> Ref: Vu et al. (2009) Molecules 14:2535
        real(kind=8) :: eps_polysaccharide_frac                                 !< Polysaccharide mass fraction of dry EPS [-] (typ. 0.40-0.95)
        
        !> \brief Mass fraction of proteins in dry EPS [-]
        !> \details Functional enzymes and structural filaments.
        !> Range: 0.01-0.60 of dry EPS mass.
        !> Temperature effect: protein-to-polysaccharide ratio increases
        !> at elevated temperatures (lower sugar production).
        !>
        !> Ref: Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623
        real(kind=8) :: eps_protein_frac                                        !< Protein mass fraction of dry EPS [-] (typ. 0.01-0.60)
        
        !> \brief Mass fraction of extracellular DNA (eDNA) in dry EPS [-]
        !> \details Structural and informational role in the biofilm matrix.
        !> Range: 0.01-0.10.
        !>
        !> Ref: Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623
        real(kind=8) :: eps_edna_frac                                           !< eDNA mass fraction of dry EPS [-] (typ. 0.01-0.10)
        
        !> \brief Mass fraction of lipids in dry EPS [-]
        !> \details Membrane vesicles, hydrophobic fractions.
        !> Range: 0.01-0.10.
        !>
        !> Ref: Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623
        real(kind=8) :: eps_lipid_frac                                          !< Lipid mass fraction of dry EPS [-] (typ. 0.01-0.10)
        
        !> \brief EPS yield coefficient [-]
        !> \details Fraction of substrate consumed that is diverted to EPS production.
        !> Y_EPS = 0.12-0.22 gCOD_EPS/gCOD_substrate.
        !> Luedeking-Piret model: r_EPS = alpha * mu * X + beta * X
        !>   - alpha: growth-associated EPS production coefficient
        !>   - beta: non-growth-associated EPS production coefficient [1/T]
        !> Yield-based model: r_EPS = f * Y * r_S / (1 + f)
        !>   where f = fraction of substrate diverted to EPS.
        !>
        !> Ref: Kreft et al. (2001) Microbiology 147:2897
        !> Ref: Xavier & Foster (2007) PNAS 104:876
        !> Ref: Jayathilake et al. (2017) PLOS ONE 12:e0181965
        real(kind=8) :: Y_eps                                                   !< EPS yield coefficient [-] (typ. 0.12-0.22)
        
        !> \brief Growth-associated EPS production coefficient [-]
        !> \details alpha in Luedeking-Piret: r_EPS = alpha * mu * X + beta * X.
        !>
        !> Ref: Kreft et al. (2001) Microbiology 147:2897
        real(kind=8) :: eps_alpha                                               !< Luedeking-Piret alpha [-] (growth-associated)
        
        !> \brief Non-growth-associated EPS production rate [1/T]
        !> \details beta in Luedeking-Piret: r_EPS = alpha * mu * X + beta * X.
        !>
        !> Ref: Kreft et al. (2001) Microbiology 147:2897
        real(kind=8) :: eps_beta                                                !< Luedeking-Piret beta [1/s] (non-growth-associated)
        
        !=====================================================================
        ! 3. TRANSPORT PROPERTIES
        !=====================================================================
        
        !> \brief Effective diffusion ratio D_e/D_aq [-]
        !> \details Ratio of effective diffusion coefficient within biofilm
        !> to aqueous diffusion coefficient.
        !> General ranges (Stewart 2003):
        !>   - Light gases (O2, CO2, N2O, CH4): ~0.60
        !>   - Most organic solutes: ~0.25
        !>   - Overall range: 0.2-0.8, mean ~0.4
        !>   - O2 average across studies: 0.57
        !> If per-species ratios are needed, use D_ratio_array(:).
        !>
        !> Related by tortuosity: D_e = (porosity / tau^2) * D_aq
        !>
        !> Ref: Stewart (2003) J. Bacteriol. 185:1485
        !> Ref: Berg et al. (2020) Biotechnol. Bioeng. 117:2997
        real(kind=8) :: D_e_ratio                                               !< Effective diffusion ratio D_e/D_aq [-] (typ. 0.2-0.8, mean ~0.4)
        
        !> \brief Per-solute effective diffusion ratios D_e/D_aq [-]
        !> \details Array of ratios for each dissolved species.
        !> Examples at 25 degC (Stewart 2003):
        !>   Solute        | D_aq (m^2/s)   | D_e/D_aq | D_e (m^2/s)
        !>   O2            | 2.0e-9         | 0.57     | 1.14e-9
        !>   Glucose       | 6.7e-10        | 0.25     | 1.68e-10
        !>   CO2           | 1.92e-9        | ~0.60    | ~1.15e-9
        !>   NH4+          | 1.97e-9        | ~0.40    | ~7.9e-10
        !>   Sucrose       | 5.2e-10        | ~0.25    | ~1.3e-10
        !>   HSL (QS)      | ~5e-10         | ~0.20    | ~1.0e-10
        !>
        !> Ref: Stewart (2003) J. Bacteriol. 185:1485
        real(kind=8), allocatable :: D_ratio_array(:)                           !< Per-solute D_e/D_aq ratios [-] (size = num_solutes)
        
        !> \brief Tortuosity factor [-]
        !> \details Relates effective diffusivity to porosity:
        !>   D_e = (epsilon / tau^2) * D_aq
        !> where epsilon is the biofilm porosity.
        !> Range: ~1.2-2.5. Higher in dense/compact biofilms.
        !>
        !> Ref: Berg et al. (2020) Biotechnol. Bioeng. 117:2997
        real(kind=8) :: tortuosity                                              !< Tortuosity factor tau [-] (typ. 1.2-2.5)
        
        !=====================================================================
        ! 4. DETACHMENT / EROSION PARAMETERS
        !=====================================================================
        
        !> \brief Continuous erosion rate coefficient [1/s]
        !> \details First-order erosion: r_det = k_det * Lf.
        !> Proportional to wall shear stress in flow-cell experiments.
        !> Average detachment rate increases ~2-fold per 2x thickness increase.
        !>
        !> Ref: Wanner & Gujer (1986) Biotechnol. Bioeng. 28:314
        !> Ref: Liu et al. (2022) AEM 88:e01072-22
        real(kind=8) :: k_det                                                   !< Continuous erosion coefficient k_det [1/s]
        
        !> \brief Shear-dependent detachment coefficient [1/(Pa*s)]
        !> \details For shear-dependent detachment models:
        !>   r_det = delta * tau^gamma * Lf
        !> where tau = wall shear stress [Pa], gamma ~ 0.58.
        !> The shear stress tau is provided externally by the flow solver.
        !>
        !> Ref: Abbas et al. (2012) Math. Biosci. Eng. 9:215
        real(kind=8) :: k_det_shear                                             !< Shear-dependent detachment coefficient delta [varies]
        
        !> \brief Shear exponent for detachment [-]
        !> \details Exponent gamma in: r_det = delta * tau^gamma * Lf.
        !> Default: gamma = 0.58 (derived from experimental data).
        !>
        !> Ref: Abbas et al. (2012) Math. Biosci. Eng. 9:215
        real(kind=8) :: det_shear_exponent                                      !< Shear exponent gamma [-] (typ. ~0.58)
        
        !> \brief Sloughing rate coefficient [1/s]
        !> \details Rate coefficient for sloughing when thickness exceeds maximum:
        !>   r_det = k_sl * max(0, Lf - Lf_max)
        !>
        !> Ref: Wanner & Gujer (1986) Biotechnol. Bioeng. 28:314
        real(kind=8) :: k_sloughing                                             !< Sloughing coefficient k_sl [1/s]
        
        !=====================================================================
        ! 5. TEMPERATURE CORRECTION
        !=====================================================================
        
        !> \brief Reference temperature for kinetic parameters [K]
        !> \details Temperature at which all kinetic parameters (mu_max, K_s, etc.)
        !> are defined. Typically 293.15 K (20 degC) for ASM1 parameters.
        real(kind=8) :: T_ref                                                   !< Reference temperature T_ref [K] (typ. 293.15 = 20 degC)
        
        !> \brief Arrhenius temperature correction factor theta [-]
        !> \details Van't Hoff-Arrhenius temperature correction:
        !>   k(T) = k(T_ref) * theta^(T - T_ref)
        !> Typical values for biological systems: 1.01-1.10.
        !> For heterotrophic growth: theta ~ 1.07 (Q10 ~ 2).
        !> For nitrification: theta ~ 1.10 (more temperature-sensitive).
        !> Alternatively, Arrhenius activation energy E_a [J/mol]:
        !>   k(T) = k(T_ref) * exp(E_a/R * (1/T_ref - 1/T))
        !> Generic E_a for biodegradation: ~65.4 kJ/mol (ECHA recommendation).
        !>
        !> Ref: Henze et al. (1987) ASM1 — temperature correction
        !> Ref: ECHA guidance — E_a = 65.4 kJ/mol for biodegradation
        real(kind=8) :: theta                                                   !< Temperature correction factor theta [-] (typ. 1.01-1.10)
        
        !> \brief Arrhenius activation energy [J/mol]
        !> \details Alternative to theta for temperature correction:
        !>   k(T) = k(T_ref) * exp(E_a/R * (1/T_ref - 1/T))
        !> Generic value for biodegradation: ~65,400 J/mol (65.4 kJ/mol).
        !> Q10 values of 2-3 are typical for microbial processes.
        !>
        !> Ref: ECHA (2021) — recommended value
        real(kind=8) :: E_activation                                            !< Activation energy E_a [J/mol] (typ. ~65400)
        
        !=====================================================================
        ! 6. POROSITY-PERMEABILITY COUPLING (BIOCLOGGING)
        !=====================================================================
        
        !> \brief Minimum porosity (clogging limit) [-]
        !> \details Pores will not be completely blocked; permeability
        !> will not drop to zero. A residual porosity remains for flow.
        !> Typical: n_min ~ 0.01-0.10.
        real(kind=8) :: porosity_min                                            !< Minimum residual porosity n_min [-] (typ. 0.01-0.10)
        
        !> \brief Kozeny-Carman exponent for bioclogging [-]
        !> \details Exponent in generalized porosity-permeability relation:
        !>   k_b/k_0 = (n_b/n_0)^m
        !> Standard Kozeny-Carman: m = 3 (with correction for specific surface).
        !> Empirical range for bioclogging: m = 6-12 (stronger clogging effect).
        !>
        !> Ref: Hommel et al. (2018) Transport in Porous Media
        real(kind=8) :: KC_exponent                                             !< Kozeny-Carman exponent m [-] (typ. 3-12)
        
        !=====================================================================
        ! 7. BIOLOGICAL — MICROORGANISMS
        !=====================================================================
        
        !> \brief Number of microorganism species in the biofilm [-]
        integer(kind=4) :: num_microorgs                                        !< Count of distinct microorganism types [-]
        
        !> \brief Array of microorganism objects defining each species
        !> \details Each microorganism has properties such as:
        !> - Decay rate and decay flag
        !> - Growth kinetics parameters (mu_max, K_s, K_o2)
        !> - Stoichiometric coefficients (yield)
        type(microorganism_c), allocatable :: microorganisms(:)                 !< Microorganism species array (size = num_microorgs)
        
    contains
        !> Lifecycle and computation methods
        procedure :: compute_biomass_decay                                      !< First-order decay rates [M/(L^3*T)]
        procedure :: compute_D_eff                                              !< Effective diffusion from D_aq and ratio [m^2/s]
        procedure :: compute_eps_production                                     !< Luedeking-Piret EPS production rate [M/(L^3*T)]
        procedure :: correct_rate_for_temperature                               !< Arrhenius/theta correction [-]
    end type
    
contains
    
    !=========================================================================
    ! DECAY
    !=========================================================================
    
    !> \brief Compute first-order decay rates for all microorganism species
    !> \details Calculates biomass decay (death/lysis) rates:
    !>   rk_d(i) = k_decay(i) * X(i) if decay_flag = .true., else 0
    !>
    !> \param[in] this Biofilm object
    !> \param[in] conc Biomass concentrations [M/L^3] (size = num_microorgs)
    !> \return rk_d Decay rates [M/(L^3*T)] (size = num_microorgs)
    function compute_biomass_decay(this, conc) result(rk_d)
        class(biofilm_c) :: this                                                !< Biofilm object
        real(kind=8), intent(in) :: conc(:)                                     !< Biomass concentrations X(i) [M/L^3]
        real(kind=8), allocatable :: rk_d(:)                                    !< Decay rates [M/(L^3*T)]
        integer(kind=4) :: i
        
        allocate(rk_d(this%num_microorgs))
        do i = 1, this%num_microorgs
            if (this%microorganisms(i)%decay_flag) then
                rk_d(i) = this%microorganisms(i)%decay_rate * conc(i)
            else
                rk_d(i) = 0d0
            end if
        end do
    end function
    
    !=========================================================================
    ! EFFECTIVE DIFFUSION
    !=========================================================================
    
    !> \brief Compute effective diffusion coefficient for a solute in the biofilm
    !> \details D_e = D_aq * (D_e/D_aq)
    !> Uses per-solute ratio if available, otherwise the bulk D_e_ratio.
    !>
    !> Typical values at 25 degC (Stewart 2003):
    !>   O2:      D_e = 2.0e-9 * 0.57 = 1.14e-9 m^2/s
    !>   Glucose: D_e = 6.7e-10 * 0.25 = 1.68e-10 m^2/s
    !>
    !> \param[in] this Biofilm object
    !> \param[in] D_aq Aqueous diffusion coefficient [m^2/s]
    !> \param[in] i_solute Solute index (optional; uses D_ratio_array if allocated)
    !> \return D_eff Effective diffusion coefficient [m^2/s]
    !>
    !> Ref: Stewart (2003) J. Bacteriol. 185:1485
    pure real(kind=8) function compute_D_eff(this, D_aq, i_solute) result(D_eff)
        class(biofilm_c), intent(in) :: this                                    !< Biofilm object
        real(kind=8), intent(in) :: D_aq                                        !< Aqueous diffusion coeff [m^2/s]
        integer(kind=4), intent(in), optional :: i_solute                       !< Solute index [-]
        
        if (present(i_solute) .and. allocated(this%D_ratio_array)) then
            D_eff = D_aq * this%D_ratio_array(i_solute)
        else
            D_eff = D_aq * this%D_e_ratio
        end if
    end function
    
    !=========================================================================
    ! EPS PRODUCTION (Luedeking-Piret)
    !=========================================================================
    
    !> \brief Compute EPS production rate via Luedeking-Piret model
    !> \details r_EPS = alpha * mu * X + beta * X
    !>   - alpha: growth-associated term [-]
    !>   - beta: non-growth-associated maintenance term [1/T]
    !>   - mu: specific growth rate [1/T]
    !>   - X: biomass concentration [M/L^3]
    !>
    !> Alternative yield-based model:
    !>   r_EPS = f * Y * r_S / (1 + f)   with f ~ Y_eps
    !>
    !> \param[in] this Biofilm object
    !> \param[in] mu Specific growth rate [1/s]
    !> \param[in] X Biomass concentration [M/L^3]
    !> \return r_eps EPS volumetric production rate [M/(L^3*s)]
    !>
    !> Ref: Kreft et al. (2001) Microbiology 147:2897
    !> Ref: Xavier & Foster (2007) PNAS 104:876
    pure real(kind=8) function compute_eps_production(this, mu, X) result(r_eps)
        class(biofilm_c), intent(in) :: this                                    !< Biofilm object
        real(kind=8), intent(in) :: mu                                          !< Specific growth rate [1/s]
        real(kind=8), intent(in) :: X                                           !< Biomass concentration [M/L^3]
        
        r_eps = this%eps_alpha * mu * X + this%eps_beta * X
    end function
    
    !=========================================================================
    ! TEMPERATURE CORRECTION
    !=========================================================================
    
    !> \brief Correct a rate constant for temperature via theta factor
    !> \details Van't Hoff-Arrhenius correction:
    !>   k(T) = k(T_ref) * theta^(T - T_ref)
    !> where T and T_ref are in [K] or [degC] (same units for both).
    !>
    !> If theta is not set but E_activation is available, uses Arrhenius:
    !>   k(T) = k(T_ref) * exp(E_a/R * (1/T_ref - 1/T))
    !> where R = 8.314 J/(mol*K).
    !>
    !> Typical theta values:
    !>   - Heterotrophic growth: ~1.07
    !>   - Nitrification: ~1.10
    !>   - General biodegradation: ~1.04-1.08
    !>
    !> \param[in] this Biofilm object
    !> \param[in] k_ref Rate constant at T_ref [any units]
    !> \param[in] T Current temperature [K]
    !> \return k_T Temperature-corrected rate constant [same units as k_ref]
    !>
    !> Ref: Henze et al. (1987) ASM1
    pure real(kind=8) function correct_rate_for_temperature(this, k_ref, T) result(k_T)
        class(biofilm_c), intent(in) :: this                                    !< Biofilm object
        real(kind=8), intent(in) :: k_ref                                       !< Rate constant at T_ref
        real(kind=8), intent(in) :: T                                           !< Current temperature [K]
        real(kind=8), parameter :: R_gas = 8.314d0                              !< Universal gas constant [J/(mol*K)]
        
        if (this%theta > 0d0 .and. this%theta /= 1d0) then
            ! Theta-factor method
            k_T = k_ref * this%theta ** (T - this%T_ref)
        else if (this%E_activation > 0d0) then
            ! Arrhenius method
            k_T = k_ref * exp(this%E_activation / R_gas * (1d0/this%T_ref - 1d0/T))
        else
            ! No temperature correction
            k_T = k_ref
        end if
    end function
    
end module biofilm_m
