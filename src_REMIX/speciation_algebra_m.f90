!> This structure:
!!>   contains algebraic information relevant for speciation: dimensions and arrays
module speciation_algebra_m
    use metodos_sist_lin_m, only: inv_matrix, compute_det, id_matrix
    implicit none
    save
    type, public :: speciation_algebra_c
        logical :: flag_comp !> TRUE if component matrix has no constant activity species (De Simoni et al, 2005)
        logical :: flag_cat_exch !> TRUE if there are cation exchange reactions in reactive zone, FALSE otherwise
        integer(kind=4) :: num_species !> number of species
        integer(kind=4) :: num_eq_reactions !> number of equilibrium reactions
        integer(kind=4) :: num_cst_act_species !> number of constant activity species
        integer(kind=4) :: num_var_act_species !> number of variable activity species
        integer(kind=4) :: num_sec_aq_species !> number of secondary aqueous species
        integer(kind=4) :: num_aq_sec_var_act_species !> number of aqueous secondary variable activity species
        integer(kind=4) :: num_aq_var_act_species !> number of aqueous variable activity species
        integer(kind=4) :: num_prim_species !> number of primary species
        integer(kind=4) :: num_aq_prim_species !> number of aqueous primary species
        integer(kind=4) :: num_sec_species !> number of secondary species
        real(kind=8), allocatable :: inv_Se_2(:,:) !> inverse of equilibrium secondary species stoichiometric submatrix
        real(kind=8), allocatable :: Se_1_star(:,:) !> =-inv_Se_2*S_e_1
        real(kind=8), allocatable :: inv_Se_nc_2(:,:) !> inverse of equilibrium variable activity secondary species stoichiometric submatrix
        real(kind=8), allocatable :: Se_nc_1_star(:,:) !> =-inv_Se_nc_2*S_e_1
        real(kind=8), allocatable :: logK_star(:) !> =inv_Se_nc_2*log(K)
        real(kind=8), allocatable :: logK_tilde(:) !> =inv_Se_2*log(K)
        real(kind=8), allocatable :: comp_mat_cst_act(:,:) !> component matrix with constant activity species (Saaltink et al, 1998)
        real(kind=8), allocatable :: comp_mat(:,:) !> component matrix without constant activity species (De Simoni et al, 2005)
        real(kind=8), allocatable :: comp_mat_aq(:,:) !> aqueous component matrix without constant activity species (De Simoni et al, 2005)
    contains
    !> Set
        procedure, public :: set_flag_comp
        procedure, public :: set_flag_cat_exch
        procedure, public :: set_dimensions
    !> Compute
        procedure, public :: compute_num_prim_species
        procedure, public :: compute_num_sec_species
        procedure, public :: compute_num_sec_aq_species
        procedure, public :: compute_num_aq_var_act_species
        procedure, public :: compute_num_aq_sec_var_act_species
        procedure, public :: compute_arrays
        procedure, private :: compute_Se_1_star
        procedure, private :: compute_Se_nc_1_star
        procedure, private :: compute_comp_mat
        procedure, private :: compute_comp_mat_aq
        procedure, private :: compute_comp_mat_cst_act
        procedure, public :: compute_logK_star
        procedure, public :: compute_logK_tilde
        procedure, public :: compute_inv_Se_nc_2
        procedure, public :: compute_inv_Se_2
    end type
    
    contains
        subroutine set_dimensions(this,n_sp,n_e,n_c,n_aq,n_nc_aq,num_min_kin,num_gas_kin)
            implicit none
            class(speciation_algebra_c) :: this
            integer(kind=4), intent(in) :: n_sp !> number of species
            integer(kind=4), intent(in) :: n_e !> number of equilibrium reactions
            integer(kind=4), intent(in) :: n_c !> number of constant activity species
            integer(kind=4), intent(in) :: n_aq !> number of aqueous species
            integer(kind=4), intent(in) :: n_nc_aq !> number of aqueous variable activity species
            integer(kind=4), intent(in), optional :: num_min_kin !> number of mineral kinetic reactions
            integer(kind=4), intent(in), optional :: num_gas_kin !> number of gases not in equilibrium
            this%num_species=n_sp
            this%num_eq_reactions=n_e
            this%num_cst_act_species=n_c
            this%num_var_act_species=this%num_species-this%num_cst_act_species
            this%num_aq_var_act_species=n_nc_aq
            call this%compute_num_prim_species(num_min_kin,num_gas_kin)
            call this%compute_num_sec_species()
            call this%compute_num_sec_aq_species(n_aq)            
            call this%compute_num_aq_sec_var_act_species()
        end subroutine

        subroutine compute_num_sec_aq_species(this,n_aq)
            implicit none
            class(speciation_algebra_c) :: this
            integer(kind=4), intent(in) :: n_aq !> number of aqueous species
            this%num_sec_aq_species=n_aq-this%num_aq_prim_species
        end subroutine

        subroutine compute_num_sec_species(this)
            implicit none
            class(speciation_algebra_c) :: this
            if (this%flag_comp.eqv..true.) then
                this%num_sec_species=this%num_eq_reactions+this%num_cst_act_species
            else
                this%num_sec_species=this%num_eq_reactions
            end if
        end subroutine
        
        subroutine compute_comp_mat(this) !> see De Simoni et al (2005) for defintion of component matrix
            implicit none
            class(speciation_algebra_c) :: this
            if (allocated(this%comp_mat)) then
                deallocate(this%comp_mat)
            end if
            allocate(this%comp_mat(this%num_prim_species,this%num_var_act_species))
            this%comp_mat(:,1:this%num_prim_species)=id_matrix(this%num_prim_species)
            this%comp_mat(:,this%num_prim_species+1:this%num_var_act_species)=transpose(this%Se_nc_1_star)
        end subroutine
        
        subroutine compute_comp_mat_aq(this) !> see De Simoni et al (2005) for defintion of component matrix
            implicit none
            class(speciation_algebra_c) :: this
            if (allocated(this%comp_mat_aq)) then
                deallocate(this%comp_mat_aq)
            end if
            allocate(this%comp_mat_aq(this%num_aq_prim_species,this%num_aq_var_act_species))
            this%comp_mat_aq(:,1:this%num_aq_prim_species)=id_matrix(this%num_aq_prim_species)
            this%comp_mat_aq(:,this%num_aq_prim_species+1:this%num_aq_var_act_species)=this%comp_mat(:,this%num_prim_species+1:&
                this%num_prim_species+this%num_aq_sec_var_act_species)
        end subroutine
        
        subroutine compute_comp_mat_cst_act(this) !> computes component matrix with constant activity species (Saaltink et al, 1998)
            implicit none
            class(speciation_algebra_c) :: this
            if (allocated(this%comp_mat_cst_act)) then
                deallocate(this%comp_mat_cst_act)
            end if
            allocate(this%comp_mat_cst_act(this%num_prim_species,this%num_species))
            this%comp_mat_cst_act(:,1:this%num_prim_species)=id_matrix(this%num_prim_species)
            this%comp_mat_cst_act(:,this%num_prim_species+1:this%num_species)=transpose(this%Se_1_star)
        end subroutine

        subroutine compute_inv_Se_2(this,Se_2,tol)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(in) :: Se_2(:,:) !> secondary equilibrium stoichiometic sbmatrix
            real(kind=8), intent(in) :: tol !> for computing inverse matrix
            if (size(Se_2,1)/=size(Se_2,2)) then
                error stop "Matrix must be square (compute_inv_Se_2)"
            end if
            if (allocated(this%inv_Se_2)) then
                deallocate(this%inv_Se_2)
            end if
            allocate(this%inv_Se_2(this%num_eq_reactions,this%num_eq_reactions))
            if (size(Se_2,1)==1) then
                this%inv_Se_2(1,1)=1d0/Se_2(1,1)
            else
                call inv_matrix(Se_2,tol,this%inv_Se_2)
            end if
        end subroutine
        
        subroutine compute_inv_Se_nc_2(this,Se_nc_2,tol)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(in) :: Se_nc_2(:,:) !> secondary variable activity equilibrium stoichiometic sbmatrix
            real(kind=8), intent(in) :: tol !> for computing inverse matrix
            if (allocated(this%inv_Se_nc_2)) then
                deallocate(this%inv_Se_nc_2)
            end if
            allocate(this%inv_Se_nc_2(this%num_eq_reactions,this%num_eq_reactions))
            call inv_matrix(Se_nc_2,tol,this%inv_Se_nc_2)
        end subroutine
        
        subroutine compute_arrays(this,Se,K,tol,flag,swap)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(inout) :: Se(:,:) !> equilibrium stoichiometric matrix
            real(kind=8), intent(inout) :: K(:) !> equilibrium constants
            real(kind=8), intent(in) :: tol !> for computing inverse matrix
            logical, intent(out) :: flag !> TRUE if stoichiometric matrix has been modified, FALSE otherwise
            integer(kind=4), intent(out) :: swap(:) !> indices of swapped species or reactions in reactive zone
            
            integer(kind=4) :: i
            real(kind=8) :: det
            logical :: error
            real(kind=8), allocatable :: Se_nc_2(:,:),Se_2(:,:),aux_Se(:,:),aux_K(:)
            
            flag=.false. !> by default
            !aux_Se=Se
            
            !allocate(swap(2)) !> chapuza
            !swap=0 !> no columns are swapped by default
            
            i=0
            !> esto de abajo deberia ir dentro del if
            aux_Se=Se
            aux_K=K
            if (this%flag_comp .eqv. .true.) then
                allocate(Se_nc_2(this%num_eq_reactions,this%num_eq_reactions))
                do
                    Se_nc_2=Se(:,this%num_prim_species+1:this%num_var_act_species) !> secondary variable activity equilibrium stoichiometric submatrix
                    call compute_det(Se_nc_2,tol,det,error)
                    if ((error .eqv. .true.) .or. (abs(det)<tol)) then
                        flag=.true.
                        i=i+1
                        Se=aux_Se
                        Se(:,this%num_prim_species+1)=aux_Se(:,this%num_prim_species-i+1)
                        Se(:,this%num_prim_species-i+1)=aux_Se(:,this%num_prim_species+1)
                    else
                        call this%compute_inv_Se_nc_2(Se_nc_2,tol)
                        swap(1)=this%num_prim_species-i+1 !> index of primary species to swap
                        swap(2)=this%num_prim_species+1 !> index of secondary variable activity species to swap
                        call this%compute_Se_nc_1_star(Se(:,1:this%num_prim_species))
                        call this%compute_logK_star(K)
                        call this%compute_comp_mat()
                        call this%compute_comp_mat_aq()
                        exit
                    end if
                end do
                deallocate(Se_nc_2)
            else
                allocate(Se_2(this%num_eq_reactions,this%num_eq_reactions))
                do
                    Se_2=Se(:,this%num_prim_species+1:this%num_species) !> secondary equilibrium stoichiometric submatrix
                    call compute_det(Se_2,tol,det,error)
                    if ((error .eqv. .true.) .or. (abs(det)<tol)) then
                        flag=.true.
                        i=i+1
                        Se=aux_Se
                        !> We swap equilibrium reactions
                        Se(this%num_eq_reactions,:)=aux_Se(this%num_eq_reactions-i,:) !> index of first equilibrium reaction to swap
                        Se(this%num_eq_reactions-i,:)=aux_Se(this%num_eq_reactions,:) !> index of second equilibrium reaction to swap
                        !Se(:,this%num_species)=aux_Se(:,this%num_species-i)
                        !Se(:,this%num_species-i)=aux_Se(:,this%num_species)
                    else
                        call this%compute_inv_Se_2(Se_2,tol)
                        swap(1)=this%num_eq_reactions-i
                        swap(2)=this%num_eq_reactions
                        !swap(1)=this%num_species-i
                        !swap(2)=this%num_species
                        call this%compute_Se_1_star(Se(:,1:this%num_prim_species))
                        K(swap(1))=aux_K(swap(2))
                        K(swap(2))=aux_K(swap(1))
                        call this%compute_logK_tilde(K)
                        call this%compute_comp_mat_cst_act()
                        exit
                    end if
                end do
                deallocate(Se_2)
            end if
        end subroutine 
        
        
        subroutine compute_Se_nc_1_star(this,Se_1)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(in) :: Se_1(:,:) !> primary equilibrium stoichiometric submatrix
            if (allocated(this%Se_nc_1_star)) then
                deallocate(this%Se_nc_1_star)
            end if
            allocate(this%Se_nc_1_star(this%num_eq_reactions,this%num_prim_species))
            this%Se_nc_1_star=-matmul(this%inv_Se_nc_2,Se_1)
        end subroutine
        
        subroutine compute_Se_1_star(this,Se_1)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(in) :: Se_1(:,:) !> primary equilibrium stoichiometric submatrix
            if (allocated(this%Se_1_star)) then
                deallocate(this%Se_1_star)
            end if
            allocate(this%Se_1_star(this%num_eq_reactions,this%num_prim_species))
            this%Se_1_star=-matmul(this%inv_Se_2,Se_1)
        end subroutine
        
        subroutine compute_logK_star(this,K)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(in) :: K(:) !> equilibrium constants
            this%logK_star=matmul(this%inv_Se_nc_2,log10(K))
        end subroutine
        
        subroutine compute_logK_tilde(this,K)
            implicit none
            class(speciation_algebra_c) :: this
            real(kind=8), intent(in) :: K(:) !> equilibrium constants
            this%logK_tilde=matmul(this%inv_Se_2,log10(K))
        end subroutine
        
        subroutine compute_num_prim_species(this,num_min_kin,num_gas_kin)
            implicit none
            class(speciation_algebra_c) :: this
            integer(kind=4), intent(in), optional :: num_min_kin
            integer(kind=4), intent(in), optional :: num_gas_kin
            if (this%flag_comp.eqv..true.) then
                this%num_prim_species=this%num_species-this%num_eq_reactions-this%num_cst_act_species
            else
                this%num_prim_species=this%num_species-this%num_eq_reactions
            end if
            this%num_aq_prim_species=this%num_prim_species
            if (this%flag_cat_exch.eqv..true.) then
                this%num_aq_prim_species=this%num_aq_prim_species-1
            end if
            if (present(num_min_kin) .and. present(num_gas_kin)) then
                if (num_min_kin>0 .and. this%flag_comp.eqv..false.) then
                    this%num_aq_prim_species=this%num_aq_prim_species-num_min_kin
                end if
                if (num_gas_kin>0 .and. this%flag_comp.eqv..false.) then
                    this%num_aq_prim_species=this%num_aq_prim_species-num_gas_kin
                end if
            end if
        end subroutine
        
        subroutine compute_num_aq_var_act_species(this)
            implicit none
            class(speciation_algebra_c) :: this
            if (this%flag_comp.eqv..true.) then
                this%num_aq_var_act_species=this%num_aq_prim_species+this%num_aq_sec_var_act_species
            end if
        end subroutine
        
        subroutine compute_num_aq_sec_var_act_species(this)
            implicit none
            class(speciation_algebra_c) :: this
            !integer(kind=4), intent(in) :: wat_flag !> cst act water flag
            if (this%flag_comp.eqv..true.) then !> chapuza
                this%num_aq_sec_var_act_species=this%num_aq_var_act_species-this%num_aq_prim_species
            else
                this%num_aq_sec_var_act_species=this%num_aq_var_act_species-this%num_aq_prim_species+1
            end if
        end subroutine
        
        subroutine set_flag_comp(this,flag_comp)
            implicit none
            class(speciation_algebra_c) :: this
            logical, intent(in) :: flag_comp
            this%flag_comp=flag_comp
        end subroutine
        
        subroutine set_flag_cat_exch(this,flag_cat_exch)
            implicit none
            class(speciation_algebra_c) :: this
            logical, intent(in) :: flag_cat_exch
            this%flag_cat_exch=flag_cat_exch
        end subroutine

    

end module 