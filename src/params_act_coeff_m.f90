!> @file params_act_coeff_m.f90
!> @brief Activity coefficient parameters module for aqueous species
!> @details This module contains parameters that control the activity coefficient calculation
!> of aqueous species using various models (Debye-Hückel, Davies, Truesdell-Jones, b-dot).
!>
!> @par Activity Coefficient Theory:
!> The activity coefficient γᵢ relates activity aᵢ to molality mᵢ: aᵢ = γᵢ·mᵢ
!> For ionic species, log₁₀(γᵢ) depends on ionic strength I according to various models.
!>
!> @par General Form:
!> log₁₀(γᵢ) = -α·A·zᵢ²·√I / (1 + β·√I) + γ·I
!> where A, B are Debye-Hückel constants, zᵢ is valence, and α, β, γ are model-specific parameters.
!>
!> @par Supported Models:
!> 0. Ideal (γ = 1)
!> 1. Debye-Hückel restricted: α=1, β=0, γ=0
!> 2. Debye-Hückel extended: α=1, β=a₀·B, γ=0
!> 3. Davies: α=1, β=1, γ=0.3·A·z²
!> 4. Truesdell-Jones: α=1, β=B·aTJ, γ=bTJ
!> 5. b-dot (Helgeson): α=1, β=B·aTJ, γ=ḃ(T)
!>
!> @see params_aq_sol_m For Debye-Hückel constants A and B
!> @author Jordi Petchamé-Guerrero
!> @date November 2025

module params_act_coeff_m
    use params_aq_sol_m, only: params_aq_sol_s !< Import aqueous solution parameters (Debye-Hückel constants A, B)
    implicit none !< Enforce explicit variable declarations for type safety
    save !< Preserve module variables between procedure calls
    private !< Default visibility is private; expose only selected types/procedures
    !> @brief Parameters for activity coefficient models
    !> @details This type stores all parameters required to compute activity coefficients
    !> using various models. The general form is:
    !> log₁₀(γᵢ) = -α·A·zᵢ²·√I / (1 + β·√I) + γ·I
    !>
    !> @par Parameter Usage by Model:
    !> - **ion_size_param (a₀)**: Ion size parameter [Angstroms] used in extended Debye-Hückel
    !> - **a_TJ, b_TJ**: Truesdell-Jones fit parameters (species-specific)
    !> - **b_dot (ḃ)**: Temperature-dependent b-dot parameter (Helgeson model)
    !> - **alpha (α)**: Multiplier for Debye-Hückel term (typically 1)
    !> - **beta (β)**: Controls short-range interactions (varies by model)
    !> - **gamma (γ)**: Linear ionic strength term (varies by model)
    type, public :: params_act_coeff_c !< Activity coefficient parameters for various models
        real(kind=8) :: ion_size_param !< Effective hydrated ion diameter a₀ [Angstroms] (Debye-Hückel extended)
        real(kind=8) :: a_TJ !< Ion-specific fit parameter a_{TJ} [Angstroms] (Truesdell-Jones model)
        real(kind=8) :: b_TJ !< Ion-specific fit parameter b_{TJ} [kg/mol] (Truesdell-Jones model)
        real(kind=8) :: b_dot !< Temperature-dependent b-dot term ḃ(T) [kg/mol] (Helgeson 1969)
        real(kind=8) :: alpha=0d0 !< Multiplier for Debye-Hückel A term (α) [-] (default: 0 for ideal)
        real(kind=8) :: beta=0d0 !< Multiplier for short-range interaction term (β) [-] (default: 0)
        real(kind=8) :: gamma=0d0 !< Linear ionic strength coefficient (γ) [kg/mol] (default: 0)
    contains
    !> @name Setter Methods
    !> @{
        procedure :: set_ion_size_param !< Set ion size parameter a₀
    !> @}
    !> @name Computational Methods
    !> @{
        procedure :: compute_b_dot !< Compute temperature-dependent b-dot term (Helgeson 1969)
        procedure :: compute_csts !< Compute α, β, γ coefficients for specified model
    !> @}
    end type
    
    !> @brief Reserved interface block for future extensions
    !> @details Currently empty - may be used for operator overloading or generic interfaces
    interface
    end interface
    
    contains
        !> @brief Set ion size parameter for extended Debye-Hückel model
        !> @details Assigns the effective hydrated ion diameter a₀ used in extended Debye-Hückel:
        !> log₁₀(γᵢ) = -A·zᵢ²·√I / (1 + a₀·B·√I)
        !>
        !> @par Typical Values:
        !> - Small ions (Li⁺, H⁺): a₀ ≈ 3-6 Å
        !> - Medium ions (Na⁺, Cl⁻): a₀ ≈ 4-5 Å
        !> - Large ions (K⁺, I⁻): a₀ ≈ 5-9 Å
        !>
        !> @param[inout] this Activity coefficient parameters object
        !> @param[in] ion_size_param Effective hydrated ion diameter a₀ [Angstroms]
        subroutine set_ion_size_param(this,ion_size_param)
            implicit none !< Enforce explicit variable declarations
            class(params_act_coeff_c), intent(inout) :: this !< Activity coefficient parameters object [-]
            real(kind=8), intent(in) :: ion_size_param !< Ion size parameter a₀ [Angstroms]
            this%ion_size_param=ion_size_param !< Assign ion size parameter [Angstroms]
        end subroutine
        
        !> @brief Compute temperature-dependent b-dot parameter
        !> @details Calculates the b-dot term ḃ(T) for Helgeson (1969) activity coefficient model:
        !> ḃ(T) = 15698.4/T + 41.8808·ln(T) - 0.0367626·T - 974168/T² - 268.902
        !> where T is temperature [K].
        !>
        !> @par Helgeson Model:
        !> log₁₀(γᵢ) = -A·zᵢ²·√I / (1 + aTJ·B·√I) + ḃ(T)·I
        !>
        !> @par Reference:
        !> Helgeson, H.C. (1969). Thermodynamics of hydrothermal systems at elevated
        !> temperatures and pressures. Am. J. Sci., 267, 729-804.
        !>
        !> @param[inout] this Activity coefficient parameters object
        !> @param[in] temp Temperature [K]
        subroutine compute_b_dot(this,temp) !< Helgeson (1969) b-dot calculation
            implicit none !< Enforce explicit variable declarations
            class(params_act_coeff_c), intent(inout) :: this !< Activity coefficient parameters object [-]
            real(kind=8), intent(in) :: temp !< Temperature [K]
            !> Helgeson (1969) empirical formula for b-dot as function of temperature
            !> ḃ(T) = c₁/T + c₂·ln(T) + c₃·T + c₄/T² + c₅ [kg/mol]
            this%b_dot=15698.4/temp + 41.8808*log(temp) - 0.0367626*temp - 974168/(temp**2) - 268.902
        end subroutine
        
        !> @brief Compute α, β, γ coefficients for specified activity coefficient model
        !> @details Calculates the three parameters (alpha, beta, gamma) in the generalized form:
        !> log₁₀(γᵢ) = -α·A·zᵢ²·√I / (1 + β·√I) + γ·I
        !> based on the selected model and species valence.
        !>
        !> @par Model Implementations:
        !> - **0 (Ideal)**: γ = 1 always → α=0, β=0, γ=0
        !> - **1 (Debye-Hückel restricted)**: log₁₀(γ) = -A·z²·√I → α=1, β=0, γ=0
        !> - **2 (Debye-Hückel extended)**: β = a₀·B → α=1, β=a₀·B, γ=0
        !> - **3 (Davies)**: log₁₀(γ) = -A·z²·[√I/(1+√I) - 0.3·I] → α=1, β=1, γ=0.3·A·z²
        !> - **4 (Truesdell-Jones)**: Ion-specific parameters → α=1, β=B·a_TJ, γ=b_TJ
        !> - **5 (b-dot/Helgeson)**: Temperature-dependent γ → α=1, β=B·a_TJ, γ=ḃ(T)
        !>
        !> @param[inout] this Activity coefficient parameters object
        !> @param[in] z Ion valence (charge number) zᵢ [-]
        !> @param[in] params_aq_sol Aqueous solution parameters (Debye-Hückel A, B constants)
        !> @param[in] model Model selector: 0=ideal, 1=DH restricted, 2=DH extended, 3=Davies, 4=Truesdell-Jones, 5=b-dot
        subroutine compute_csts(this,z,params_aq_sol,model) !< Compute α, β, γ for activity model
            implicit none !< Enforce explicit variable declarations
            class(params_act_coeff_c), intent(inout) :: this !< Activity coefficient parameters object [-]
            integer(kind=4), intent(in) :: z !< Ion valence (charge number) [-]
            class(params_aq_sol_s), intent(in) :: params_aq_sol !< Aqueous solution parameters (A, B) [-]
            integer(kind=4), intent(in) :: model !< Model selector [0-5]
            
            !> Select model and compute corresponding α, β, γ parameters
            select case (model)
            case(0) !< Model 0: Ideal solution (γ = 1, no ionic strength effects)
                this%alpha=0d0 !< No Debye-Hückel term [-]
                this%beta=0d0 !< No denominator correction [-]
                this%gamma=0d0 !< No linear I term [kg/mol]
            case(1) !< Model 1: Debye-Hückel restricted (limiting law)
                this%alpha=1d0 !< Full Debye-Hückel term [-]
                this%beta=0d0 !< No denominator correction [-]
                this%gamma=0d0 !< No linear I term [kg/mol]
            case(2) !< Model 2: Debye-Hückel extended (ion size correction)
                this%alpha=1d0 !< Full Debye-Hückel term [-]
                this%beta=this%ion_size_param*params_aq_sol%B !< Ion size × B parameter [(Å)·(kg/mol)^(1/2)]
                this%gamma=0d0 !< No linear I term [kg/mol]
            case(3) !< Model 3: Davies equation
                this%alpha=1d0 !< Full Debye-Hückel term [-]
                this%beta=1d0 !< Davies uses β = 1 [-]
                this%gamma=0.3*params_aq_sol%A*(z**2) !< Davies correction: 0.3·A·z² [kg/mol]
            case(4) !< Model 4: Truesdell-Jones (ion-specific parameters)
                this%alpha=1d0 !< Full Debye-Hückel term [-]
                this%beta=params_aq_sol%B*this%a_TJ !< Species-specific β = B·aTJ [(kg/mol)^(1/2)]
                this%gamma=this%b_TJ !< Species-specific γ = bTJ [kg/mol]
            case(5) !< Model 5: b-dot (Helgeson) - temperature-dependent
                this%alpha=1d0 !< Full Debye-Hückel term [-]
                this%beta=params_aq_sol%B*this%a_TJ !< Species-specific β = B·aTJ [(kg/mol)^(1/2)]
                this%gamma=this%b_dot !< Temperature-dependent γ = ḃ(T) [kg/mol]
            case default !< Invalid model number
                error stop "Model not implemented yet" !< Terminate for unsupported model
            end select
        end subroutine
end module