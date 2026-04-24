!> \file biofilm_m.f90
!> \brief Biofilm module for modeling biofilm processes in reactive transport
!> \details This module defines the biofilm class and related functionality for simulating 
!> biofilm growth, decay, and interaction with aqueous chemistry in porous media.
!>
!> The biofilm is modeled as an immobile zone (extending phase_c) that:
!> - Contains water and biomass with specific volumetric fractions
!> - Hosts multiple microorganism species with individual growth/decay kinetics
!> - Interacts with mobile aqueous phase through mass transfer
!> - Has mechanical, chemical, transport, structural, and detachment properties
!>
!> Key features:
!> - Biofilm volumetric fraction represented by immobile porosity (inherited from phase_c)
!> - Water content within biofilm specified by f_w parameter
!> - Multiple microorganism types with individual volumetric fractions
!> - Biomass decay kinetics with optional decay rates per microorganism
!> - EPS composition (polysaccharides, proteins, eDNA, lipids)
!> - Mechanical properties (elastic modulus, yield stress, viscosity)
!> - Transport properties (effective diffusion, tortuosity, boundary layer)
!> - Structural properties (roughness, fractal dimension, surface coverage)
!> - Detachment/erosion model parameters
!> - Temperature correction via Arrhenius
!> - Porosity-permeability coupling for porous media bioclogging
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
    !> \details Represents a biofilm as a special type of immobile zone containing microorganisms.
    !>
    !> Physical interpretation:
    !> - Biofilm occupies a volumetric fraction of the porous medium (phi_imm from parent class)
    !> - Within the biofilm, water occupies fraction f_w, biomass occupies the remainder
    !> - Each microorganism type has its own volumetric fraction within the biofilm
    !>
    !> Mass balance considerations:
    !> - Water in biofilm: V_water = phi_imm * f_w * V_total
    !> - Biomass volume: V_biomass = phi_imm * (1 - f_w) * V_total
    !> - Microorganism i volume: V_microorg_i = vol_fracs(i) * V_biomass
    !>
    !> The class collects ALL known measurable properties of biofilm into a single
    !> derived type, organized into logical groups:
    !>   1. Volumetric/porosity
    !>   2. Physical/mechanical
    !>   3. EPS composition
    !>   4. Transport
    !>   5. Structural/morphological
    !>   6. Detachment/erosion
    !>   7. Temperature correction
    !>   8. Porosity-permeability coupling (bioclogging)
    !>   9. Biological/kinetic (microorganisms array)
    !>
    !> \note The immobile porosity (phi_imm) represents the total biofilm volumetric fraction
    type, public, extends(phase_c) :: biofilm_c                                 !< Biofilm class inheriting from phase base class
        
        !=====================================================================
        ! 1. VOLUMETRIC / POROSITY PROPERTIES
        !=====================================================================
        
        !> \brief Volumetric fraction of water within the biofilm [-]
        !> \details Represents the porosity of the biofilm itself (not the bulk porous medium).
        !> Values typically range from 0.80 to 0.99 for hydrated biofilms.
        !> Water content by wet weight is typically > 70 wt% (Flemming & Wingender 2010),
        !> commonly cited as ~97 wt% (Costerton group).
        !> The remaining fraction (1 - f_w) is occupied by biomass (cells + EPS).
        !>
        !> Within the biofilm water, three types exist (Quan et al. 2022):
        !>   - Free water (highest diffusivity, in macropores/channels)
        !>   - Confined water (in channels, intermediate mobility)
        !>   - Bound water (adsorbed to EPS polymers, reduced mobility)
        !>
        !> Ref: Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623
        !> Ref: Quan et al. (2022) Crit. Rev. Microbiol. 48:283
        real(kind=8) :: f_w                                                     !< Water volumetric fraction in biofilm: f_w = V_water / V_biofilm [-] (typ. 0.80-0.99)
        
        !> \brief Volume fraction of cells within the biofilm [-]
        !> \details Fraction of biofilm volume occupied by bacterial cells.
        !> Typically 10-20% of total biofilm volume (0.10-0.20).
        !> Costerton et al. (2003) report ~15% cells, ~85% water channels.
        !> Volumetric bacterial density: 0.2-0.4 bacteria/um^3 (Quan et al. 2022).
        !>
        !> Ref: Costerton et al. (2003) J. Clin. Invest. 112:1466
        !> Ref: Quan et al. (2022) Crit. Rev. Microbiol. 48:283
        real(kind=8) :: f_cells                                                 !< Cell volume fraction in biofilm [-] (typ. 0.10-0.20)
        
        !> \brief Volume fraction of EPS within the biofilm [-]
        !> \details Fraction of biofilm volume occupied by extracellular polymeric substances.
        !> EPS can comprise 50-90% of total organic matter in biofilms.
        !> Constraint: f_w + f_cells + f_eps ~= 1.0
        !>
        !> Ref: Flemming & Wingender (2010) Nat. Rev. Microbiol. 8:623
        real(kind=8) :: f_eps                                                   !< EPS volume fraction in biofilm [-] (typ. 0.02-0.20)
        
        !=====================================================================
        ! 2. PHYSICAL / MECHANICAL PROPERTIES
        !=====================================================================
        
        !> \brief Biofilm thickness [m]
        !> \details Average thickness of the biofilm layer.
        !> Ranges from ~5 um (thin clinical) to >1 mm (industrial reactors).
        !> Classification (Flemming et al. 2023):
        !>   - Thin: < 10 um
        !>   - Medium: 10-200 um
        !>   - Thick: > 200 um
        !> Industrial MBBR/RBC: 200-1000 um.
        !>
        !> Ref: Flemming et al. (2023) npj Biofilms Microbiomes 9:36
        !> Ref: Hall-Stoodley et al. (2022) Nat. Rev. Microbiol. 20:608
        real(kind=8) :: thickness                                               !< Biofilm thickness Lf [m] (typ. 1e-5 to 1e-3)
        
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
        ! 3. EPS COMPOSITION (mass fractions of dry EPS)
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
        ! 4. TRANSPORT PROPERTIES
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
        
        !> \brief Oxygen penetration depth [m]
        !> \details Depth from biofilm surface at which O2 concentration reaches zero.
        !> P. aeruginosa flat slab, cell density 12 mg/cm^3, bulk O2 = 6 mg/L: ~77 um.
        !> Typical range for aerobic biofilms: 50-200 um.
        !> Biofilms > 100-200 um commonly have anoxic interior.
        !>
        !> Ref: Stewart (2003) J. Bacteriol. 185:1485
        !> Ref: de Beer et al. (1993) — O2 microelectrode studies
        real(kind=8) :: o2_penetration_depth                                    !< O2 penetration depth [m] (typ. 50e-6 to 200e-6)
        
        !> \brief External mass transfer boundary layer thickness [m]
        !> \details Thickness of the liquid film at the biofilm-bulk interface
        !> through which diffusion dominates substrate transport.
        !> Range: 10-300 um (depends on hydrodynamic conditions).
        !>   - Vigorous mixing/turbulent: ~10 um (minimal resistance)
        !>   - Laminar flow: up to 300 um (significant resistance)
        !> The boundary layer thickness determines the Sherwood number:
        !>   Sh = Lf / L_bl (local mass transfer coefficient)
        !>
        !> Ref: Picioreanu (1998) ISME Symposium — multidimensional modelling
        !> Ref: Stewart (2004) Antimicrob. Agents Chemother. 48:— (model)
        real(kind=8) :: bl_thickness                                            !< External boundary layer thickness L_bl [m] (typ. 10e-6 to 300e-6)
        
        !=====================================================================
        ! 5. STRUCTURAL / MORPHOLOGICAL PROPERTIES
        !=====================================================================
        
        !> \brief Surface roughness coefficient (Ra, dimensionless) [-]
        !> \details Beyenal-Lewandowski roughness:
        !>   Ra = sigma_L / L_mean
        !> where sigma_L = std. dev. of thickness, L_mean = mean thickness.
        !>   0 = perfectly smooth flat surface
        !>   1 = very rough surface
        !> Shewanella oneidensis thin biofilm: Ra = 0.26 +/- 0.05
        !> Rough finger-shaped protrusions: Ra = 0.41
        !>
        !> Ref: Beyenal et al. (2004) Biotechnol. Bioeng. 78:359
        !> Ref: Liu et al. (2022) AEM 88:e01072-22
        real(kind=8) :: roughness_Ra                                            !< Roughness coefficient Ra [-] (0 to 1)
        
        !> \brief Fractal dimension (2D cross-section) [-]
        !> \details Characterizes the geometric complexity of the biofilm surface.
        !> Range 2D: 1.5-2.5.
        !> Biofilm clusters in porous media (B. subtilis): Df ~ 1.84-1.93.
        !> High nutrient -> low Df (compact). Low nutrient -> high Df (fingering).
        !>
        !> Ref: Picioreanu et al. (1999) Water Sci. Tech. 39:115
        !> Ref: Doro et al. (2023) ETH Zurich — morphogenesis in porous media
        real(kind=8) :: fractal_dim_2d                                          !< 2D fractal dimension Df [-] (typ. 1.5-2.5)
        
        !> \brief Fractal dimension (3D) [-]
        !> \details 3D counterpart. Range: 2-3.
        !>
        !> Ref: Hermanowicz et al. (2001)
        real(kind=8) :: fractal_dim_3d                                          !< 3D fractal dimension Df [-] (typ. 2-3)
        
        !> \brief Surface coverage fraction [-]
        !> \details Fraction of substratum surface covered by biofilm.
        !> Range: 0.10-1.0 (10-100%). Function of attachment time and nutrient level.
        !>
        !> Ref: Costerton et al. (2003) J. Clin. Invest. 112:1466
        real(kind=8) :: surface_coverage                                        !< Surface coverage fraction [-] (0 to 1)
        
        !> \brief Specific biovolume per unit area [m^3/m^2 = m]
        !> \details Biofilm biovolume per substratum area. Measured by CLSM/COMSTAT.
        !> Range: 0.1-50 um^3/um^2 = 1e-7 to 5e-5 m.
        !>
        !> Ref: Renslow et al. (2011) Biotechnol. Bioeng. 108:1383
        real(kind=8) :: specific_biovolume                                      !< Specific biovolume per area [m^3/m^2]
        
        !=====================================================================
        ! 6. DETACHMENT / EROSION PARAMETERS
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
        ! 7. TEMPERATURE CORRECTION
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
        ! 8. POROSITY-PERMEABILITY COUPLING (BIOCLOGGING)
        !=====================================================================
        
        !> \brief Initial (clean) porosity of porous medium [-]
        !> \details Porosity before biofilm colonization.
        !> Used in Kozeny-Carman type relations to compute permeability reduction:
        !>   k_b/k_0 = (n_b/n_0)^3 * ((1-n_0)/(1-n_b))^2  (Kozeny-Carman)
        !> where n_b = n_0 - f_biofilm (biofilm-colonized porosity).
        !>
        !> Ref: Hommel et al. (2018) Transport in Porous Media — review
        real(kind=8) :: porosity_initial                                        !< Initial clean porosity of porous medium n_0 [-]
        
        !> \brief Biofilm-affected porosity of porous medium [-]
        !> \details Reduced porosity due to biofilm occupying pore space:
        !>   n_b = n_0 - phi_biofilm
        !> where phi_biofilm is the volumetric fraction of biofilm.
        !>
        !> Ref: Thullner et al. (2002) Water Resources Research
        real(kind=8) :: porosity_biofilm                                        !< Biofilm-affected porosity n_b [-]
        
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
        ! 9. CHEMICAL ENVIRONMENT (internal gradients)
        !=====================================================================
        
        !> \brief Interior pH range — minimum [-]
        !> \details Biofilm interior can have significant pH gradients.
        !> Examples:
        !>   - Dental plaque (S. mutans): pH 4.0-5.5 (lactic acid)
        !>   - Aerobic heterotrophic surface: near bulk 7.0-8.0
        !>   - Deep anoxic zones: can drop to 6.0-6.5
        !>   - Nitrifying biofilm: local pH decrease (H+ from nitrification)
        !>
        !> Ref: Stewart & Franklin (2008) Nat. Rev. Microbiol. 6:199
        real(kind=8) :: pH_min                                                  !< Minimum interior pH [-] (system-dependent)
        
        !> \brief Interior pH range — maximum [-]
        real(kind=8) :: pH_max                                                  !< Maximum interior pH [-] (system-dependent)
        
        !> \brief Redox potential range in aerobic zone [mV vs SHE]
        !> \details Aerobic zone near surface: +200 to +400 mV.
        !> Anaerobic core: -100 to -200 mV.
        !> Chemical heterogeneity leads to simultaneous aerobic respiration,
        !> denitrification, and sulfate reduction in thick (>200 um) biofilms.
        !>
        !> Ref: Stewart & Franklin (2008) Nat. Rev. Microbiol. 6:199
        real(kind=8) :: Eh_aerobic                                              !< Redox potential aerobic zone [mV vs SHE] (typ. +200 to +400)
        
        !> \brief Redox potential in anaerobic core [mV vs SHE]
        real(kind=8) :: Eh_anaerobic                                            !< Redox potential anaerobic zone [mV vs SHE] (typ. -100 to -200)
        
        !=====================================================================
        ! 10. BIOLOGICAL — MICROORGANISMS
        !=====================================================================
        
        !> \brief Number of microorganism species in the biofilm [-]
        integer(kind=4) :: num_microorgs                                        !< Count of distinct microorganism types [-]
        
        !> \brief Array of microorganism objects defining each species
        !> \details Each microorganism has properties such as:
        !> - Decay rate and decay flag
        !> - Growth kinetics parameters (mu_max, K_s, K_o2)
        !> - Stoichiometric coefficients (yield)
        type(microorganism_c), allocatable :: microorganisms(:)                 !< Microorganism species array (size = num_microorgs)
        
        !> \brief Volumetric fraction of each microorganism within the biomass [-]
        !> \details Specifies the spatial distribution of different microorganism types.
        !> Sum should equal 1.0: sum(vol_fracs) = 1.0
        !> These fractions are relative to the biomass volume, NOT total biofilm.
        !> Biomass volume = phi_biofilm * (1 - f_w) * V_total
        real(kind=8), allocatable :: vol_fracs(:)                               !< Microorganism volumetric fractions [-] (size = num_microorgs)
        
    contains
        !> Lifecycle and computation methods
        procedure :: compute_biomass_decay                                      !< First-order decay rates [M/(L^3*T)]
        procedure :: compute_D_eff                                              !< Effective diffusion from D_aq and ratio [m^2/s]
        procedure :: compute_permeability_reduction                             !< Kozeny-Carman permeability ratio [-]
        procedure :: compute_eps_production                                     !< Luedeking-Piret EPS production rate [M/(L^3*T)]
        procedure :: correct_rate_for_temperature                               !< Arrhenius/theta correction [-]
        procedure :: compute_detachment_rate                                    !< Detachment rate [m/s]
        procedure :: compute_thiele_modulus                                     !< Thiele modulus phi^2 [-]
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
    ! KOZENY-CARMAN PERMEABILITY REDUCTION
    !=========================================================================
    
    !> \brief Compute relative permeability reduction due to bioclogging
    !> \details Generalized Kozeny-Carman for bioclogging:
    !>   k_b/k_0 = (n_b / n_0)^m
    !> where n_b = max(porosity_biofilm, porosity_min),
    !> n_0 = porosity_initial, m = KC_exponent.
    !>
    !> Standard Kozeny-Carman (m=3, with specific surface correction):
    !>   k_b/k_0 = (n_b/n_0)^3 * ((1-n_0)/(1-n_b))^2
    !>
    !> Bioclogging models show m = 6-12 may better fit experiments.
    !>
    !> \param[in] this Biofilm object
    !> \return k_ratio Permeability ratio k_b/k_0 [-]
    !>
    !> Ref: Hommel et al. (2018) Transport in Porous Media
    !> Ref: Thullner et al. (2002) Water Resources Research
    pure real(kind=8) function compute_permeability_reduction(this) result(k_ratio)
        class(biofilm_c), intent(in) :: this                                    !< Biofilm object
        real(kind=8) :: n_b                                                     !< Effective biofilm porosity
        
        n_b = max(this%porosity_biofilm, this%porosity_min)
        
        if (this%porosity_initial > 0d0) then
            k_ratio = (n_b / this%porosity_initial) ** this%KC_exponent
        else
            k_ratio = 1d0
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
    
    !=========================================================================
    ! DETACHMENT RATE
    !=========================================================================
    
    !> \brief Compute total detachment rate [m/s]
    !> \details Combines continuous erosion and sloughing:
    !>   r_det = k_det * Lf                                    (erosion)
    !>         + k_sl * max(0, Lf - Lf_max)                    (sloughing)
    !>         + k_det_shear * tau^gamma * Lf                   (shear-dependent, optional)
    !>
    !> The result is the velocity at which the biofilm surface recedes [m/s].
    !> Thickness evolution: dLf/dt = u(Lf) - r_det
    !> where u(Lf) is the expansion velocity from growth.
    !>
    !> \param[in] this Biofilm object
    !> \param[in] tau Wall shear stress [Pa] (optional; for shear-dependent term)
    !> \return r_det Total detachment rate [m/s]
    !>
    !> Ref: Wanner & Gujer (1986) Biotechnol. Bioeng. 28:314
    !> Ref: Abbas et al. (2012) Math. Biosci. Eng. 9:215
    pure real(kind=8) function compute_detachment_rate(this, tau) result(r_det)
        class(biofilm_c), intent(in) :: this                                    !< Biofilm object
        real(kind=8), intent(in), optional :: tau                               !< Wall shear stress [Pa]
        
        ! Continuous erosion: r = k_det * Lf
        r_det = this%k_det * this%thickness
        
        ! Sloughing: r = k_sl * max(0, Lf - Lf_max)
        if (this%thickness > this%thickness_max .and. this%k_sloughing > 0d0) then
            r_det = r_det + this%k_sloughing * (this%thickness - this%thickness_max)
        end if
        
        ! Shear-dependent: r = delta * tau^gamma * Lf
        if (present(tau) .and. this%k_det_shear > 0d0) then
            r_det = r_det + this%k_det_shear * (tau ** this%det_shear_exponent) * this%thickness
        end if
    end function
    
    !=========================================================================
    ! THIELE MODULUS
    !=========================================================================
    
    !> \brief Compute the Thiele modulus (phi^2) for the biofilm
    !> \details Dimensionless number indicating diffusion vs. reaction limitation:
    !>   phi^2 = mu_max * X * Lf^2 / (D_eff * S_bulk * Y)
    !>
    !> Interpretation:
    !>   - Low phi^2 (metabolic limitation): compact, smooth, homogeneous biofilm
    !>   - High phi^2 (diffusion limitation): rough, porous, channeled biofilm
    !>
    !> This determines biofilm morphology (Picioreanu 1998).
    !>
    !> \param[in] this Biofilm object
    !> \param[in] mu_max Maximum specific growth rate [1/s]
    !> \param[in] X Biomass density [M/L^3]
    !> \param[in] D_eff Effective diffusion coefficient [m^2/s]
    !> \param[in] S_bulk Bulk substrate concentration [M/L^3]
    !> \param[in] Y Yield coefficient [-]
    !> \return phi2 Thiele modulus phi^2 [-]
    !>
    !> Ref: Picioreanu et al. (1998) Biotechnol. Bioeng. 57:718
    pure real(kind=8) function compute_thiele_modulus(this, mu_max, X, D_eff, S_bulk, Y) result(phi2)
        class(biofilm_c), intent(in) :: this                                    !< Biofilm object
        real(kind=8), intent(in) :: mu_max                                      !< Max specific growth rate [1/s]
        real(kind=8), intent(in) :: X                                           !< Biomass density [M/L^3]
        real(kind=8), intent(in) :: D_eff                                       !< Effective diffusion [m^2/s]
        real(kind=8), intent(in) :: S_bulk                                      !< Bulk substrate concentration [M/L^3]
        real(kind=8), intent(in) :: Y                                           !< Yield coefficient [-]
        
        if (D_eff > 0d0 .and. S_bulk > 0d0 .and. Y > 0d0) then
            phi2 = mu_max * X * this%thickness**2 / (D_eff * S_bulk * Y)
        else
            phi2 = 0d0
        end if
    end function
    
end module biofilm_m
