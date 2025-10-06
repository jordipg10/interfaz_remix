!> Aqueous phase class (subclass of phase class):
!!      contains the aqueous species in a solution
!!      computes log_10 activity coefficients and its log_10-Jacobian
!!      rearranges aqueous species in variable & constant activity, respectively
module aq_phase_m
    use phase_m, only : phase_c
    use aq_species_m, only : aq_species_c, species_c
    use matrices_m, only : diag_matrix_c
    use params_aq_sol_m, only : params_aq_sol_t
    !use vectors_m
    implicit none
    save
    type, public, extends(phase_c) :: aq_phase_c
        integer(kind=4) :: num_aq_complexes=0 !> number of aqueous complexes
        type(aq_species_c), allocatable :: aq_species(:) !> aqueous species in solution (we assume they are sorted in primary, secondary variable activity & secondary constant activity
        integer(kind=4) :: wat_flag=0 !> 1 if water is present as a species, 0 otherwise
        integer(kind=4) :: ind_wat=0 !> index of water in 'aq_species' array
        integer(kind=4), allocatable :: ind_diss_solids(:) !> indices of dissolved solids
        integer(kind=4) :: ind_prot=0 !> index 'H+' in 'aq_species' array
        integer(kind=4) :: ind_oh=0 !> index 'OH-' in 'aq_species' array
        integer(kind=4) :: ind_bicarb=0 !> index 'HCO3-' in 'aq_species' array
        integer(kind=4) :: ind_carb=0 !> index 'CO3-2' in 'aq_species' array
    contains
    !> Set
        procedure, public :: set_ind_prot
        procedure, public :: set_ind_wat
        procedure, public :: set_ind_diss_solids
        procedure, public :: set_indices_aq_phase
        procedure, public :: set_num_aq_complexes
        procedure, public :: set_aq_species
        procedure, public :: set_single_aq_species
    !> Allocate
        procedure, public :: allocate_aq_species
        procedure, public :: allocate_ind_diss_solids
    !> Compute
        procedure, public :: compute_log_act_coeffs_aq_phase
        procedure, public :: compute_log_Jacobian_act_coeffs_aq_phase
    !> Rearrange
        procedure, public :: rearrange_aq_species
    !> Is
        procedure, public :: is_species_in_aq_phase
        procedure, public :: is_water_in_aq_phase
    !> Copy
        procedure, public :: copy_attributes
    !> Get
        procedure, public :: get_valences
    end type
    
    interface
        ! subroutine compute_log_act_coeffs_aq_phase(this,ionic_act,params_aq_sol,log_act_coeffs)
        !     import aq_phase_c
        !     import params_aq_sol_t
        !     implicit none
        !     class(aq_phase_c) :: this
        !     real(kind=8), intent(in) :: ionic_act
        !     class(params_aq_sol_t), intent(in) :: params_aq_sol
        !     real(kind=8), intent(out) :: log_act_coeffs(:) !> must be allocated
        ! end subroutine
        
        ! subroutine compute_log_Jacobian_act_coeffs_aq_phase(this,out_prod,conc,log_Jacobian_act_coeffs)
        !     import aq_phase_c
        !     implicit none
        !     class(aq_phase_c) :: this
        !     real(kind=8), intent(in) :: out_prod(:,:) !> outer product between d_log_gamma_d_I and z^2
        !     real(kind=8), intent(in) :: conc(:) !> concentration of species in a given target
        !     real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> must be allocated
        ! end subroutine        
        
    end interface
    
    contains
        
        subroutine set_num_aq_complexes(this,num_aq_complexes)
            implicit none
            class(aq_phase_c) :: this
            integer(kind=4), intent(in) :: num_aq_complexes
            this%num_aq_complexes=num_aq_complexes
        end subroutine
        
        subroutine set_ind_prot(this,ind_proton)
           implicit none
           class(aq_phase_c) :: this
           integer(kind=4), intent(in) :: ind_proton
           this%ind_prot=ind_proton
        end subroutine
        
        subroutine set_ind_wat(this,ind_wat)
           implicit none
           class(aq_phase_c) :: this
           integer(kind=4), intent(in) :: ind_wat
           this%ind_wat=ind_wat
        end subroutine
        
        subroutine allocate_aq_species(this,num_species)
            implicit none
            class(aq_phase_c) :: this
            integer(kind=4), intent(in), optional :: num_species
            if (present(num_species)) then
                this%num_species=num_species
            end if
            if (allocated(this%aq_species)) then
                deallocate(this%aq_species)
            end if
            allocate(this%aq_species(this%num_species))
        end subroutine
        
        subroutine allocate_ind_diss_solids(this)
           implicit none
           class(aq_phase_c) :: this
           allocate(this%ind_diss_solids(this%num_species-this%wat_flag))
        end subroutine

        
        subroutine set_aq_species(this,aq_species)
            implicit none
            class(aq_phase_c) :: this
            class(aq_species_c), intent(in) :: aq_species(:)
            this%aq_species=aq_species
            this%num_species=size(aq_species)
        end subroutine
        
        subroutine set_single_aq_species(this,aq_species,index)
            implicit none
            class(aq_phase_c) :: this
            class(aq_species_c), intent(in) :: aq_species
            integer(kind=4), intent(in) :: index
            if (.not. allocated(this%aq_species)) then
                error stop "Aqueous species must be allocated"
            else if (index>this%num_species) then
                error stop "Index is higher than number of aqueous species"
            else if (index<1) then
                error stop "Index must be positive"
            else
                this%aq_species(index)=aq_species
            end if
        end subroutine
        
       subroutine rearrange_aq_species(this)
    !> Rearranges aqueous species in variable and constant activity, respectively
           class(aq_phase_c) :: this
           
           type(aq_species_c), allocatable :: aux(:)
           integer(kind=4) :: i,j,k,l,m
           
           aux=this%aq_species !> chapuza
           j=1 !> counter aqueous species variable activity
           l=1 !> counter dissolved solids variable activity
           k=0 !> counter aqueous species constant activity
           m=0 !> counter dissolved solids constant activity
           do i=1,this%num_species
               if (aux(i)%cst_act_flag .eqv. .false. .and. aux(i)%name/='h2o') then
                    this%ind_diss_solids(l)=j
                    l=l+1
                    call this%aq_species(j)%assign_species(aux(i))
                    if (j<this%num_species) then
                        j=j+1
                    end if
               else if (aux(i)%name/='h2o(p)') then
                    this%ind_diss_solids(this%num_species-m)=this%num_species-k
                    m=m+1
                    call this%aq_species(this%num_species-k)%assign_species(aux(i))
                    k=k+1
                else if (aux(i)%cst_act_flag .eqv. .false.) then
                    call this%set_ind_wat(j)
                    call this%aq_species(j)%assign_species(aux(i))
                    if (j<this%num_species) then
                        j=j+1
                    end if
                else
                    call this%set_ind_wat(this%num_species-k)
                    call this%aq_species(this%num_species-k)%assign_species(aux(i))
                    k=k+1
               end if
           end do
       end subroutine
        
        ! subroutine append_aq_species(this,aq_species)
        !     implicit none
        !     class(aq_phase_c) :: this
        !     type(aq_species_c), intent(in) :: aq_species
            
        !     integer(kind=4) :: i
        !     type(aq_species_c), allocatable :: aux_array(:)
            
        !     aux_array=this%aq_species
        !     deallocate(this%aq_species)
        !     if (size(aux_array)>0) then
        !         call this%allocate_aq_species(size(aux_array)+1)
        !         do i=1,this%num_species-1
        !             this%aq_species(i)=aux_array(i)
        !         end do
        !         this%aq_species(this%num_species)=aq_species
        !     else
        !         call this%allocate_aq_species(1)
        !         call this%set_aq_species([aq_species]) !> chapuza
        !     end if
        ! end subroutine
        
        subroutine is_species_in_aq_phase(this,species,flag,species_ind)
            implicit none
            class(aq_phase_c) :: this
            class(species_c), intent(in) :: species
            logical, intent(out) :: flag
            integer(kind=4), intent(out), optional :: species_ind
            
            integer(kind=4) :: i
            
            flag=.false.
            if (present(species_ind)) then
                species_ind=0
            end if
            do i=1,this%num_species
                if (species%name==this%aq_species(i)%name) then
                    flag=.true.
                    !if (species%name=='h2o') then
                    !    this%ind_wat=i
                    !end if
                    if (present(species_ind)) then
                        species_ind=i
                    end if
                    exit
                end if
            end do
        end subroutine
        
        subroutine is_water_in_aq_phase(this,flag)
            implicit none
            class(aq_phase_c) :: this
            logical, intent(out) :: flag
            
            type(aq_species_c) :: water
            
            water%name='h2o'
            call this%is_species_in_aq_phase(water,flag)
        end subroutine
        
        subroutine set_indices_aq_phase(this)
           implicit none
           class(aq_phase_c) :: this
           
           integer(kind=4) :: i
           
           do i=1,this%num_species
               if (this%aq_species(i)%name=='h2o(p)') then
                  this%ind_wat=i
               else if (this%aq_species(i)%name=='h+') then
                   this%ind_prot=i
               else if (this%aq_species(i)%name=='oh-') then
                   this%ind_oh=i
               else if (this%aq_species(i)%name=='hco3-') then
                   this%ind_bicarb=i
               else if (this%aq_species(i)%name=='co3-2') then
                   this%ind_carb=i
               else
                   continue
               end if
           end do
        end subroutine
        
        subroutine set_ind_diss_solids(this)
           implicit none
           class(aq_phase_c) :: this
           
           integer(kind=4) :: i,j
           
           j=1
           do i=1,this%num_species
               if (this%aq_species(i)%name/='h2o') then
                   this%ind_diss_solids(j)=i
                   if (j<this%num_species-this%wat_flag) then
                       j=j+1
                   else
                       exit
                   end if
               end if
           end do
        end subroutine
                
        subroutine copy_attributes(this,aq_phase)
            implicit none
            class(aq_phase_c) :: this
            class(aq_phase_c), intent(in) :: aq_phase
            
            integer(KIND=4) :: i
            
            if (allocated(this%aq_species)) then
                continue
            else
                call this%allocate_aq_species(aq_phase%num_species)
            end if
            this%wat_flag=aq_phase%wat_flag
            this%num_aq_complexes=aq_phase%num_aq_complexes
            this%ind_wat=aq_phase%ind_wat
            !this%ind_proton=aq_phase%ind_proton
            this%ind_diss_solids=aq_phase%ind_diss_solids
            do i=1,this%num_species
                call this%aq_species(i)%assign_species(aq_phase%aq_species(i))
            end do
        end subroutine
        
        function get_valences(this) result(valences)
            implicit none
            class(aq_phase_c), intent(in) :: this
            integer(kind=4), allocatable :: valences(:)
            
            integer(kind=4) :: i
            
            allocate(valences(this%num_species))
            do i=1,this%num_species
                valences(i)=this%aq_species(i)%valence
            end do
        end function

         !> Computes log_10 activity coefficients aqueous species variable activity
subroutine compute_log_act_coeffs_aq_phase(this,ionic_act,params_aq_sol,log_act_coeffs)
    class(aq_phase_c) :: this
    real(kind=8), intent(in) :: ionic_act !> ionic activity
    class(params_aq_sol_t), intent(in) :: params_aq_sol
    real(kind=8), intent(out) :: log_act_coeffs(:) !> log_10 activity coefficients (must be already allocated)
    
    integer(kind=4) :: i
    
        
    !if (ionic_act<1d-2) then
    !    p_compute_log_act_coeff=>Debye_Huckel_restr
    !else if (ionic_act<=1d-1) then
    !    p_compute_log_act_coeff=>Debye_Huckel_ampl
    !else if (ionic_act<=7d-1) then
    !    p_compute_log_act_coeff=>Davies
    !else
    !    error stop
    !end if
    !p_compute_log_act_coeff=>Davies !> lo impongo (chapuza)
    do i=1,this%num_species-this%wat_flag
        !call p_compute_log_act_coeff(this%aq_species(i),ionic_act,log_act_coeffs(i))
        log_act_coeffs(this%ind_diss_solids(i))=-this%aq_species(this%ind_diss_solids(i))%params_act_coeff%alpha*(params_aq_sol%A*&
        this%aq_species(this%ind_diss_solids(i))%valence**2)*sqrt(ionic_act)/(1d0+this%aq_species(&
        this%ind_diss_solids(i))%params_act_coeff%beta*sqrt(ionic_act)) +&
         this%aq_species(this%ind_diss_solids(i))%params_act_coeff%gamma*ionic_act
    end do
    !log_act_coeffs(this%ind_wat)
end subroutine

!> Computes Jacobian log_10 activity coefficients aqueous species with respect to log_10 concentrations
!! We assume all primary species are aqueous
subroutine compute_log_Jacobian_act_coeffs_aq_phase(this,out_prod,conc,log_Jacobian_act_coeffs)
    implicit none
    class(aq_phase_c) :: this
    real(kind=8), intent(in) :: out_prod(:,:) !> subset outer product between d_log_gamma_d_I and z^2
    real(kind=8), intent(in) :: conc(:) !> subset concentrations aqueous species
    real(kind=8), intent(out) :: log_Jacobian_act_coeffs(:,:) !> (must be allocated)

    type(diag_matrix_c) :: conc_diag
    
    call conc_diag%set_diag_matrix(conc) !> chapuza
        
    log_Jacobian_act_coeffs=5d-1*log(1d1)*conc_diag%prod_mat_diag_mat(out_prod)
     
end subroutine
end module